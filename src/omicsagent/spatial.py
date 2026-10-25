"""Sample-separated spatial graphs, permutation enrichment, Moran I and niches.

Niche core ported/adapted from the supplied spatial-niche-enrichment R reference.
No spatial edge or permutation crosses samples. No spot is called a single cell.
"""

from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

from .config import RunConfig


def bh(pvalues: np.ndarray) -> np.ndarray:
    result = np.full(len(pvalues), np.nan)
    valid = np.flatnonzero(np.isfinite(pvalues))
    if not len(valid):
        return result
    order = valid[np.argsort(pvalues[valid])]
    q = pvalues[order] * len(order) / np.arange(1, len(order) + 1)
    result[order] = np.minimum(1, np.minimum.accumulate(q[::-1])[::-1])
    return result


def radius_graph(
    xy: np.ndarray, samples: np.ndarray, radius: float, max_edges: int
) -> sparse.csr_matrix:
    sources: list[np.ndarray] = []
    targets: list[np.ndarray] = []
    total = 0
    for sample in np.unique(samples):
        indices = np.flatnonzero(samples == sample)
        tree = cKDTree(xy[indices])
        # Count BEFORE materializing a dense radius graph.
        n_directed = int(tree.count_neighbors(tree, radius)) - len(indices)
        total += n_directed
        if total > max_edges:
            raise ValueError("Spatial graph exceeds max_edges; reduce radius or subset data")
        pairs = tree.query_pairs(radius, output_type="ndarray")
        if len(pairs):
            left, right = indices[pairs[:, 0]], indices[pairs[:, 1]]
            sources.extend([left, right])
            targets.extend([right, left])
    if not sources:
        return sparse.csr_matrix((len(xy), len(xy)), dtype=float)
    row, col = np.concatenate(sources), np.concatenate(targets)
    return sparse.csr_matrix((np.ones(len(row)), (row, col)), shape=(len(xy), len(xy)))


def _null_stats(observed: np.ndarray, null: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean, sd = null.mean(axis=0), null.std(axis=0, ddof=1)
    z = np.divide(observed - mean, sd, out=np.full(observed.shape, np.nan), where=sd > 0)
    p = (1 + (np.abs(null - mean) >= np.abs(observed - mean) - 1e-12).sum(axis=0)) / (len(null) + 1)
    p = np.where(sd > 0, p, np.nan)  # degenerate null is not significant or testable
    return z, p


def spatial_analysis(adata: ad.AnnData, cfg: RunConfig) -> dict[str, Any]:
    assert cfg.radius is not None
    xy = np.asarray(adata.obsm[cfg.spatial_key], dtype=float)
    samples = adata.obs[cfg.sample_key].astype(str).to_numpy()
    graph = radius_graph(xy, samples, float(cfg.radius), cfg.max_edges)
    adata.obsp["oa_spatial_connectivities"] = graph
    # Ungated calls are kept in the audit but cannot silently drive cell-type enrichment.
    labels = adata.obs["oa_label"].astype(str).to_numpy(dtype=object)
    if cfg.platform == "imaging":
        labels[adata.obs["oa_gate_status"].to_numpy() != "marker_supported"] = "Uncertain"
    elif cfg.platform == "visium":
        labels = adata.obs["oa_cluster"].astype(str).to_numpy(dtype=object)
    if cfg.hierarchy:
        labels = adata.obs["oa_hierarchy_label"].astype(str).to_numpy(dtype=object)
    adata.obs["oa_spatial_label"] = pd.Categorical(labels)
    classes = sorted(set(labels))
    lookup = {label: i for i, label in enumerate(classes)}
    rng = np.random.default_rng(cfg.seed)
    enrichment, moran, niche, sample_rows = [], [], [], []
    for sample in sorted(set(samples)):
        idx = np.flatnonzero(samples == sample)
        local = graph[idx][:, idx]
        row, col = local.nonzero()
        code = np.array([lookup[x] for x in labels[idx]])
        sample_rows.append(
            {
                "sample": sample,
                "n_obs": len(idx),
                "directed_edges": local.nnz,
                "isolated_obs": int((np.asarray(local.sum(axis=1)).ravel() == 0).sum()),
            }
        )
        if not local.nnz:
            continue
        k = len(classes)

        def edge_counts(codes: np.ndarray) -> np.ndarray:
            return np.bincount(codes[row] * k + codes[col], minlength=k * k)

        observed = edge_counts(code)
        null = np.stack([edge_counts(rng.permutation(code)) for _ in range(cfg.n_permutations)])
        z, p = _null_stats(observed.astype(float), null)
        # Upper triangle avoids treating symmetric edge hypotheses as distinct.
        for i, first in enumerate(classes):
            for j in range(i, k):
                ix = i * k + j
                enrichment.append(
                    {
                        "sample": sample,
                        "label_a": first,
                        "label_b": classes[j],
                        "directed_count": int(observed[ix]),
                        "z": z[ix],
                        "p": p[ix],
                    }
                )
        x = sparse.csr_matrix(adata.layers["log1p"][idx])
        variance = (
            np.asarray(x.power(2).mean(axis=0)).ravel() - np.asarray(x.mean(axis=0)).ravel() ** 2
        )
        for gene_ix in np.argsort(variance)[::-1][: cfg.spatial_genes]:
            values = x[:, gene_ix].toarray().ravel()
            values -= values.mean()
            denom = values @ values
            if denom <= 1e-12:
                continue

            def stat(v: np.ndarray) -> float:
                return float(len(v) / local.sum() * (v @ local.dot(v)) / denom)

            obs_i = stat(values)
            null_i = np.array([stat(rng.permutation(values)) for _ in range(cfg.n_permutations)])
            # Positive autocorrelation is the declared one-sided alternative.
            pv = (1 + int((null_i >= obs_i).sum())) / (cfg.n_permutations + 1)
            moran.append(
                {"sample": sample, "gene": str(adata.var_names[gene_ix]), "moran_i": obs_i, "p": pv}
            )
        if cfg.center_group_key:
            assert cfg.center_groups is not None
            groups = (
                adata.obs[cfg.center_group_key]
                .iloc[idx]
                .astype(object)
                .fillna("")
                .astype(str)
                .to_numpy()
            )
            for group in cfg.center_groups:
                center_mask = groups == group
                if center_mask.sum() < cfg.min_centers:
                    sample_rows.append(
                        {
                            "sample": sample,
                            "skipped_center_group": group,
                            "reason": "fewer than min_centers",
                        }
                    )
                    continue
                # Counts are center-neighbor incidences: a neighbor of two centers counts twice.
                weights = np.asarray(local[center_mask].sum(axis=0)).ravel()
                if weights.sum() == 0:
                    sample_rows.append(
                        {"sample": sample, "skipped_center_group": group, "reason": "no neighbors"}
                    )
                    continue
                observed_n = np.bincount(code, weights=weights, minlength=k)
                null_n = np.stack(
                    [
                        np.bincount(rng.permutation(code), weights=weights, minlength=k)
                        for _ in range(cfg.n_permutations)
                    ]
                )
                zn, pn = _null_stats(observed_n, null_n)
                for j, label in enumerate(classes):
                    niche.append(
                        {
                            "sample": sample,
                            "group": group,
                            "neighbor_label": label,
                            "n_centers": int(center_mask.sum()),
                            "observed_incidences": int(observed_n[j]),
                            "z": zn[j],
                            "p": pn[j],
                        }
                    )
    frames = {}
    for name, rows, cols in [
        (
            "neighborhood_enrichment",
            enrichment,
            ["sample", "label_a", "label_b", "directed_count", "z", "p"],
        ),
        ("spatial_moran", moran, ["sample", "gene", "moran_i", "p"]),
        (
            "center_niches",
            niche,
            ["sample", "group", "neighbor_label", "n_centers", "observed_incidences", "z", "p"],
        ),
    ]:
        frame = pd.DataFrame(rows, columns=cols)
        frame["q"] = bh(frame["p"].to_numpy(dtype=float))
        frames[name] = frame
    return {
        "tables": frames,
        "samples": sample_rows,
        "directed_edges": graph.nnz,
        "radius": cfg.radius,
        "units": cfg.coordinate_units,
        "null": "labels shuffled within sample; fixed coordinates and graph",
        "warning": "Exploratory within-sample association, not replicated disease inference",
    }
