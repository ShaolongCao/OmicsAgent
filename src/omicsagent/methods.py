"""Optional reference annotation and a transparent RCTD-free SPLIT core port."""
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .analysis import counts
from .config import RunConfig


def reference_data(adata: ad.AnnData, cfg: RunConfig) -> tuple[ad.AnnData, ad.AnnData]:
    ref = ad.read_h5ad(cfg.reference)
    key = cfg.reference_label_key
    if not ref.var_names.is_unique or not ref.obs_names.is_unique:
        raise ValueError("Reference gene and cell IDs must be unique")
    if key not in ref.obs or ref.obs[key].isna().any():
        raise ValueError("Reference labels missing")
    genes = adata.var_names.intersection(ref.var_names, sort=False)
    if len(genes) < cfg.min_shared_genes:
        raise ValueError("Insufficient shared genes")
    query = ad.AnnData(counts(adata, cfg.counts_layer), obs=adata.obs[[]].copy(),
                       var=adata.var[[]].copy())[:, genes].copy()
    reference = ad.AnnData(counts(ref, cfg.reference_counts_layer),
                          obs=ref.obs[[key]].copy(), var=ref.var[[]].copy())[:, genes].copy()
    reference.obs[key] = reference.obs[key].astype(str)
    if reference.obs[key].nunique() < 2 or (reference.obs[key].str.len() == 0).any():
        raise ValueError("Reference requires at least two nonempty cell types")
    if (np.asarray(query.X.sum(axis=1)).ravel() == 0).any():
        raise ValueError("Query has zero counts over shared genes")
    if (np.asarray(reference.X.sum(axis=1)).ravel() == 0).any():
        raise ValueError("Reference has zero counts over shared genes")
    return query, reference


def validate_weights(weights: pd.DataFrame, index: pd.Index) -> pd.DataFrame:
    if not weights.index.equals(index) or not weights.columns.is_unique:
        raise ValueError("Weight cell order/type names do not match")
    values = weights.to_numpy(dtype=float)
    if values.ndim != 2 or values.shape[1] < 2 or not np.isfinite(values).all():
        raise ValueError("Invalid weights")
    if (values < 0).any() or (values.sum(axis=1) <= 0).any():
        raise ValueError("Weights must be nonnegative with positive row sums")
    return pd.DataFrame(values / values.sum(axis=1, keepdims=True),
                        index=index, columns=weights.columns.astype(str))


def annotate_reference(adata: ad.AnnData, cfg: RunConfig) -> dict[str, Any]:
    query, reference = reference_data(adata, cfg)
    key = cfg.reference_label_key
    if cfg.annotation_method == "tacco":
        try:
            import tacco as tc
        except ImportError as exc:
            raise RuntimeError("Install OmicsAgent[tacco] to use TACCO") from exc
        tc.tl.annotate(query, reference, annotation_key=key, result_key="weights",
                       method="OT", platform_iterations=0, bisections=None, verbose=0)
        weights = query.obsm["weights"]
    else:
        try:
            import scvi
        except ImportError as exc:
            raise RuntimeError("Install OmicsAgent[scanvi] to use scANVI") from exc
        scvi.settings.seed = cfg.seed
        unknown = "__omicsagent_unlabeled__"
        if unknown in set(reference.obs[key]):
            raise ValueError("Reference label collides with unlabeled sentinel")
        query.obs["_label"] = unknown
        reference.obs["_label"] = reference.obs[key].astype(str)
        joint = ad.concat({"reference": reference, "query": query}, label="_origin",
                          index_unique="::", join="inner")
        scvi.model.SCVI.setup_anndata(joint, batch_key="_origin", labels_key="_label")
        model = scvi.model.SCVI(joint, n_latent=10, n_layers=1, gene_likelihood="nb")
        train = {"batch_size": cfg.scanvi_batch_size, "accelerator": "cpu", "devices": 1,
                 "enable_progress_bar": False, "early_stopping": False}
        model.train(max_epochs=cfg.scanvi_epochs, **train)
        classifier = scvi.model.SCANVI.from_scvi_model(
            model, unlabeled_category=unknown, labels_key="_label")
        classifier.train(max_epochs=cfg.scanvi_classifier_epochs, **train)
        indices = np.flatnonzero(joint.obs["_origin"].to_numpy() == "query")
        weights = classifier.predict(joint, indices=indices, soft=True)
        if not weights.index.equals(joint.obs_names[indices]):
            raise ValueError("scANVI prediction order mismatch")
        weights.index = query.obs_names
        classifier.save(str(cfg.out / "scanvi_model"), overwrite=False)
    weights = validate_weights(weights, adata.obs_names)
    adata.obsm["oa_weights"] = weights
    values = weights.to_numpy()
    ordered = np.sort(values, axis=1)
    adata.obs["oa_candidate"] = weights.idxmax(axis=1).astype(str)
    adata.obs["oa_weight_max"] = values.max(axis=1)
    adata.obs["oa_weight_margin"] = ordered[:, -1] - ordered[:, -2]
    adata.obs["oa_weight_entropy"] = -(values * np.log(np.maximum(values, 1e-30))).sum(axis=1)
    weights.to_csv(cfg.out / "annotation_weights.csv.gz", index_label="cell_id")
    result: dict[str, Any] = {"method": cfg.annotation_method, "shared_genes": query.n_vars,
                              "reference_cells": reference.n_obs,
                              "types": list(weights.columns),
                              "weights_are_calibrated_confidence": False}
    if cfg.split_enabled:
        result["split"] = apply_split(adata, reference, cfg, weights)
    return result


def split_core(raw: sparse.csr_matrix, weights: np.ndarray, profiles: np.ndarray,
               max_dense_elements: int = 20_000_000) -> sparse.csr_matrix:
    """Standard (no residual-removal) SPLIT v0.3.0 formula, cells x genes.

    Keep primary-type expected contribution / all-types expected contribution.
    Profiles are type means of library-normalized reference counts. No rounding or calibrated uncertainty claim.
    """
    raw = sparse.csr_matrix(raw, dtype=float)
    if weights.shape[0] != raw.shape[0] or profiles.shape != (weights.shape[1], raw.shape[1]):
        raise ValueError("SPLIT shapes disagree")
    if any(not np.isfinite(x).all() or (x < 0).any() for x in (raw.data, weights, profiles)):
        raise ValueError("SPLIT inputs must be finite and nonnegative")
    if (weights.sum(axis=1) <= 0).any():
        raise ValueError("Empty weights")
    if raw.shape[1] > max_dense_elements:
        raise ValueError("Gene dimension exceeds memory budget")
    weights = weights / weights.sum(axis=1, keepdims=True)
    primary = weights.argmax(axis=1)
    chunk = max(1, max_dense_elements // (3 * raw.shape[1]))
    blocks = []
    for start in range(0, raw.shape[0], chunk):
        stop = min(start + chunk, raw.shape[0])
        w = weights[start:stop]
        p = primary[start:stop]
        denominator = w @ profiles
        numerator = w[np.arange(len(p)), p, None] * profiles[p]
        ratio = np.divide(numerator, denominator, out=np.zeros_like(numerator),
                          where=denominator > 0)
        blocks.append(raw[start:stop].multiply(np.clip(ratio, 0, 1)).tocsr())
    return sparse.vstack(blocks, format="csr")


def apply_split(adata: ad.AnnData, reference: ad.AnnData, cfg: RunConfig,
                weights: pd.DataFrame) -> dict[str, Any]:
    if reference.n_vars * weights.shape[1] > cfg.max_dense_elements:
        raise ValueError("SPLIT reference profile exceeds memory budget")
    normalized = sparse.diags(1 / np.asarray(reference.X.sum(axis=1)).ravel()) @ reference.X
    profiles = np.vstack([
        np.asarray(normalized[reference.obs[cfg.reference_label_key].to_numpy() == label]
                   .mean(axis=0)).ravel() for label in weights.columns
    ])
    raw = counts(adata, cfg.counts_layer)
    locations = adata.var_names.get_indexer(reference.var_names)
    corrected = split_core(raw[:, locations], weights.to_numpy(), profiles,
                           cfg.max_dense_elements)
    # Preserve genes absent from the reference, explicitly flag their unchanged status.
    coo = corrected.tocoo()
    mapped = sparse.csr_matrix((coo.data, (coo.row, locations[coo.col])), shape=raw.shape)
    covered = np.zeros(adata.n_vars, dtype=bool)
    covered[locations] = True
    mapped += raw.multiply(~covered).tocsr()
    adata.layers["counts_split"] = mapped
    adata.var["oa_split_reference_covered"] = covered
    before = np.asarray(raw.sum(axis=1)).ravel()
    after = np.asarray(mapped.sum(axis=1)).ravel()
    adata.obs["oa_split_removed_fraction"] = 1 - after / np.maximum(before, 1)
    adata.obs[["oa_candidate", "oa_split_removed_fraction"]].to_csv(cfg.out / "split_qc.csv")
    return {"implementation": "python port of SPLIT v0.3.0 standard RCTD-free core",
            "official_R_validated": False, "residual_removal": False,
            "corrected_genes": int(covered.sum()), "unchanged_genes": int((~covered).sum()),
            "removed_fraction": float(1 - after.sum() / max(float(before.sum()), 1.0)),
            "layer": "counts_split", "downstream_uses_original_counts": True}
