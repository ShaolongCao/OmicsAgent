from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler

from .config import RunConfig


def counts(adata: ad.AnnData, layer: str) -> sparse.csr_matrix:
    if layer != "X" and layer not in adata.layers:
        raise ValueError(
            f"Missing raw counts layer: {layer!r}; specify X explicitly if appropriate"
        )
    matrix = sparse.csr_matrix(adata.X if layer == "X" else adata.layers[layer], dtype=float)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    values = matrix.data
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Counts must be finite and nonnegative")
    if not np.allclose(values, np.round(values), atol=1e-6, rtol=0):
        raise ValueError("Expected integer-like raw counts, not normalized/log expression")
    return matrix


def validate(adata: ad.AnnData, cfg: RunConfig) -> dict[str, Any]:
    if adata.n_obs < 4 or adata.n_vars < 3:
        raise ValueError("At least 4 observations and 3 genes are required")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Observation and gene IDs must be unique; resolve duplicates upstream")
    counts(adata, cfg.counts_layer)
    if cfg.annotation_key and cfg.annotation_key not in adata.obs:
        raise ValueError(f"Missing annotation_key: {cfg.annotation_key}")
    if cfg.normalization == "area":
        if cfg.area_key not in adata.obs:
            raise ValueError("Missing area column")
        area = np.asarray(adata.obs[cfg.area_key], dtype=float)
        if not np.isfinite(area).all() or (area <= 0).any():
            raise ValueError("Area values must be finite and positive")
    if cfg.modality == "spatial":
        if cfg.spatial_key not in adata.obsm:
            raise ValueError(f"Missing obsm[{cfg.spatial_key!r}]")
        xy = np.asarray(adata.obsm[cfg.spatial_key])
        if xy.shape != (adata.n_obs, 2) or not np.isfinite(xy).all():
            raise ValueError("Spatial coordinates must be a finite n_obs x 2 array")
        if cfg.sample_key not in adata.obs or adata.obs[cfg.sample_key].isna().any():
            raise ValueError("sample_key must exist without missing values")
        if cfg.center_group_key and cfg.center_group_key not in adata.obs:
            raise ValueError("Missing center_group_key")
    return {
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "obs_columns": list(adata.obs.columns),
        "layers": list(adata.layers),
    }


def preprocess(adata: ad.AnnData, cfg: RunConfig) -> tuple[ad.AnnData, dict[str, Any]]:
    """Keep full-gene log expression; PCA runs on a separate HVG working copy."""
    matrix = counts(adata, cfg.counts_layer)
    prefix = cfg.mt_prefix or ("MT-" if cfg.species == "human" else "mt-")
    mt = np.array(adata.var_names.str.startswith(prefix))
    total = np.asarray(matrix.sum(axis=1)).ravel()
    detected = np.asarray((matrix > 0).sum(axis=1)).ravel()
    pct_mt = 100 * np.asarray(matrix[:, mt].sum(axis=1)).ravel() / np.maximum(total, 1)
    obs = adata.obs.copy()
    obs["oa_total_counts"] = total
    obs["oa_n_genes"] = detected
    obs["oa_pct_mt"] = pct_mt
    cell_mask = (detected >= cfg.min_genes) & (pct_mt <= cfg.max_pct_mt) & (total > 0)
    gene_mask = np.asarray((matrix[cell_mask] > 0).sum(axis=0)).ravel() >= cfg.min_cells
    if cell_mask.sum() < 4 or gene_mask.sum() < 3:
        raise ValueError("QC left fewer than 4 observations or 3 genes; inspect thresholds")
    result = adata[cell_mask, gene_mask].copy()
    result.obs = obs.loc[result.obs_names]
    result.layers["counts"] = matrix[cell_mask][:, gene_mask].copy()
    if (np.asarray(result.layers["counts"].sum(axis=1)).ravel() == 0).any():
        raise ValueError("Gene filtering created zero-library observations; adjust QC thresholds")
    result.X = result.layers["counts"].copy()
    # A prior log1p tag in a source file does not describe the explicitly supplied counts.
    result.uns.pop("log1p", None)
    if cfg.normalization == "total":
        sc.pp.normalize_total(result, target_sum=1e4)
    elif cfg.normalization == "area":
        area = result.obs[cfg.area_key].to_numpy(dtype=float)
        result.X = result.X.multiply((np.median(area) / area)[:, None]).tocsr()
    sc.pp.log1p(result)
    result.layers["log1p"] = result.X.copy()
    work = result.copy()
    if cfg.platform != "imaging" and cfg.n_hvg < result.n_vars:
        sc.pp.highly_variable_genes(work, n_top_genes=cfg.n_hvg, flavor="seurat")
        mask = work.var["highly_variable"].to_numpy()
        if mask.sum() < 3:
            raise ValueError("Too few variable genes for PCA")
        result.var["highly_variable"] = mask
        work = work[:, mask].copy()
    else:
        result.var["highly_variable"] = True
    n_pcs = min(cfg.n_pcs, work.n_obs - 1, work.n_vars - 1)
    sc.pp.pca(work, n_comps=n_pcs, zero_center=True, svd_solver="arpack", random_state=cfg.seed)
    result.obsm["X_pca"] = work.obsm["X_pca"]
    sc.pp.neighbors(
        result,
        n_neighbors=min(cfg.n_neighbors, result.n_obs - 1),
        use_rep="X_pca",
        random_state=cfg.seed,
    )
    sc.tl.leiden(
        result,
        resolution=cfg.resolution,
        key_added="oa_cluster",
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=cfg.seed,
    )
    sc.tl.umap(result, random_state=cfg.seed)
    summary = {
        "before_cells": adata.n_obs,
        "after_cells": result.n_obs,
        "before_genes": adata.n_vars,
        "after_genes": result.n_vars,
        "low_gene_cells": int((detected < cfg.min_genes).sum()),
        "high_mt_cells": int((pct_mt > cfg.max_pct_mt).sum()),
        "zero_count_cells": int((total == 0).sum()),
        "mt_genes_found": int(mt.sum()),
        "n_pcs": n_pcs,
        "n_clusters": result.obs["oa_cluster"].nunique(),
        "normalization": cfg.normalization,
    }
    return result, summary


def markers(adata: ad.AnnData, cfg: RunConfig) -> pd.DataFrame:
    sizes = adata.obs["oa_cluster"].value_counts()
    eligible = list(sizes[sizes >= 2].index)
    columns = [
        "cluster",
        "gene",
        "score",
        "logfoldchanges",
        "pvals",
        "pvals_adj",
        "pct_in",
        "pct_out",
    ]
    if len(sizes) < 2 or not eligible:
        return pd.DataFrame(columns=columns)
    sc.tl.rank_genes_groups(
        adata,
        "oa_cluster",
        groups=eligible,
        method="wilcoxon",
        layer="log1p",
        use_raw=False,
        pts=True,
        n_genes=min(cfg.top_markers, adata.n_vars),
    )
    frame = sc.get.rank_genes_groups_df(adata, group=None)
    frame = frame.rename(
        columns={
            "group": "cluster",
            "names": "gene",
            "scores": "score",
            "pct_nz_group": "pct_in",
            "pct_nz_reference": "pct_out",
        }
    )
    return frame.reindex(columns=columns)


def reference_transfer(adata: ad.AnnData, cfg: RunConfig) -> dict[str, Any]:
    reference = ad.read_h5ad(cfg.reference)
    if not reference.var_names.is_unique:
        raise ValueError("Reference gene IDs must be unique")
    if (
        cfg.reference_label_key not in reference.obs
        or reference.obs[cfg.reference_label_key].isna().any()
    ):
        raise ValueError("Reference labels missing")
    shared = adata.var_names.intersection(reference.var_names)
    if len(shared) < cfg.min_shared_genes:
        raise ValueError(
            f"Only {len(shared)} shared genes; min_shared_genes={cfg.min_shared_genes}"
        )
    if (reference.n_obs + adata.n_obs) * len(shared) > cfg.max_dense_elements:
        raise ValueError("Reference PCA exceeds max_dense_elements; use a smaller reference")
    ref = ad.AnnData(
        counts(reference, cfg.reference_counts_layer)[:, reference.var_names.get_indexer(shared)]
    )
    query = ad.AnnData(adata.layers["counts"][:, adata.var_names.get_indexer(shared)])
    for obj in (ref, query):
        if (np.asarray(obj.X.sum(axis=1)).ravel() == 0).any():
            raise ValueError("Zero library in shared panel; reference mapping is undefined")
        sc.pp.normalize_total(obj, target_sum=1e4)
        sc.pp.log1p(obj)
    # Adapted from the user's reference-PCA template: fit scaling on reference ONLY.
    scaler = StandardScaler()
    ref_x = scaler.fit_transform(ref.X.toarray())
    query_x = scaler.transform(query.X.toarray())
    pca = PCA(n_components=min(cfg.n_pcs, ref.n_obs - 1, len(shared) - 1), random_state=cfg.seed)
    ref_pca = pca.fit_transform(np.clip(ref_x, -10, 10))
    query_pca = pca.transform(np.clip(query_x, -10, 10))
    knn = KNeighborsClassifier(n_neighbors=min(cfg.knn_k, ref.n_obs))
    knn.fit(ref_pca, reference.obs[cfg.reference_label_key].astype(str))
    votes = knn.predict_proba(query_pca).max(axis=1)
    candidate = knn.predict(query_pca).astype(object)
    adata.obs["oa_reference_candidate"] = candidate
    adata.obs["oa_reference_vote"] = votes
    candidate[votes < cfg.min_vote] = "Uncertain"
    adata.obs["oa_candidate"] = candidate
    return {
        "shared_genes": len(shared),
        "reference_cells": ref.n_obs,
        "low_vote_cells": int((votes < cfg.min_vote).sum()),
        "warning": "KNN vote is not calibrated confidence or an out-of-reference detector",
    }
