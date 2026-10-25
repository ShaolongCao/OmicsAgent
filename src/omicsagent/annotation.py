"""Conservative candidate generation and raw-count evidence audit.

Adapted from scRNAseq_cell_type_annotation_basic and xenium-annotation-qc-gate.
This is a partial, tested port, not a claim to reproduce the private worked examples.
"""

import json
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .config import RunConfig

# Illustration only: never presented as independently validated ground truth.
HUMAN_IMMUNE = {
    "T cell": ["CD3D", "CD3E", "TRAC"],
    "B cell": ["MS4A1", "CD79A", "CD79B"],
    "NK cell": ["NKG7", "GNLY", "KLRD1"],
    "Myeloid": ["LYZ", "LST1", "TYROBP"],
    "Platelet": ["PPBP", "PF4"],
}


def panel(cfg: RunConfig) -> dict[str, list[str]]:
    if cfg.marker_file:
        value = json.loads(cfg.marker_file.read_text())
    elif cfg.species == "human" and cfg.tissue.lower() in {"pbmc", "blood", "synthetic"}:
        value = HUMAN_IMMUNE
    else:
        return {}
    if not isinstance(value, dict) or not all(
        isinstance(k, str)
        and k
        and isinstance(v, list)
        and v
        and all(isinstance(g, str) and g for g in v)
        for k, v in value.items()
    ):
        raise ValueError("Marker JSON must map label -> nonempty list of gene symbols")
    return {k: list(dict.fromkeys(v)) for k, v in value.items()}


def marker_candidates(
    adata: ad.AnnData, marker_panel: dict[str, list[str]]
) -> list[dict[str, Any]]:
    """Heuristic broad candidates; tied/weak scores abstain. No probability claims."""
    labels = pd.Series("Uncertain", index=adata.obs_names)
    summaries = []
    for cluster in sorted(adata.obs["oa_cluster"].astype(str).unique()):
        mask = (adata.obs["oa_cluster"].astype(str) == cluster).to_numpy()
        scores = []
        for label, genes in marker_panel.items():
            present = [g for g in genes if g in adata.var_names]
            if len(present) < 2:
                continue
            idx = adata.var_names.get_indexer(present)
            detected = np.asarray((adata.layers["counts"][mask][:, idx] > 0).mean(axis=0)).ravel()
            supported = [g for g, d in zip(present, detected, strict=True) if d >= 0.2]
            if len(supported) >= 2:
                scores.append((float(detected.mean()), label, supported))
        scores.sort(reverse=True)
        candidate = "Uncertain"
        support: list[str] = []
        if scores and (len(scores) == 1 or scores[0][0] - scores[1][0] >= 0.1):
            candidate, support = scores[0][1], scores[0][2]
        labels.loc[mask] = candidate
        summaries.append(
            {
                "cluster": cluster,
                "candidate": candidate,
                "markers": support,
                "status": "provisional" if candidate != "Uncertain" else "unresolved",
                "evidence_ids": [],
                "method": "marker_detection_heuristic",
            }
        )
    adata.obs["oa_candidate"] = labels
    return summaries


def audit_labels(
    adata: ad.AnnData, marker_panel: dict[str, list[str]], threshold: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep original candidates, mark coverage gaps, and report threshold sensitivity."""
    labels = adata.obs["oa_candidate"].astype(str)
    support = np.full(adata.n_obs, np.nan)
    available = np.zeros(adata.n_obs, dtype=int)
    rows = []
    for label in sorted(labels.unique()):
        mask = (labels == label).to_numpy()
        requested = marker_panel.get(label, [])
        present = [g for g in requested if g in adata.var_names]
        available[mask] = len(present)
        if present:
            idx = adata.var_names.get_indexer(present)
            support[mask] = np.asarray(
                (adata.layers["counts"][mask][:, idx] >= 1).sum(axis=1)
            ).ravel()
        rows.append(
            {
                "label": label,
                "n_cells": int(mask.sum()),
                "markers_on_panel": len(present),
                "markers_requested": len(requested),
                "missing_markers": ";".join(g for g in requested if g not in present),
                "unsupported_at_1": int((support[mask] < 1).sum()) if present else None,
                "unsupported_at_2": int((support[mask] < 2).sum()) if len(present) >= 2 else None,
            }
        )
    gateable = (available >= threshold) & np.isfinite(support)
    gated = labels.to_numpy(dtype=object).copy()
    gated[gateable & (support < threshold)] = "Uncertain"
    # Preserve ungateable labels, but explicitly disallow interpreting them as passed.
    status = np.where(
        gateable, np.where(support >= threshold, "marker_supported", "unsupported"), "ungateable"
    )
    adata.obs["oa_marker_support"] = support
    adata.obs["oa_gate_status"] = status
    adata.obs["oa_label"] = pd.Categorical(gated)
    total = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
    bins = pd.qcut(pd.Series(total, index=adata.obs_names), 10, labels=False, duplicates="drop")
    table = pd.crosstab(bins, labels, normalize="index")
    diagnostics = []
    for label in labels.unique():
        rho = None
        if label in table and len(table) >= 3 and table[label].nunique() > 1:
            rho = float(spearmanr(table.index, table[label]).statistic)
        diagnostics.append(
            {
                "label": label,
                "rho_count_decile": rho,
                "sparsity_risk": rho is not None and rho < -0.6,
                "interpretation": "association only; biological abundance may also vary with depth",
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(diagnostics)
