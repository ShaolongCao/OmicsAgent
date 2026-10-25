"""Iterative, evidence-audited L1/L2 refinement; raw counts remain immutable."""
from pathlib import Path
from typing import Any, Literal

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from pydantic import BaseModel, ConfigDict, Field, model_validator
from scipy import sparse

from .analysis import counts, reference_transfer
from .config import RunConfig
from .trace import sha256, write_json


class DegeneratePool(ValueError):
    """A pool lacks usable variation; preserve it as unresolved."""


class HierarchyConfig(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    initial_l1_key: str
    compartments: dict[str, list[str]]
    l1_markers: dict[str, list[str]]
    l2_markers: dict[str, list[str]] = Field(default_factory=dict)
    l2_parent: dict[str, str] = Field(default_factory=dict)
    excluded_genes: list[str] = Field(default_factory=lambda: ["IRF7", "BST2"])
    batch_key: str | None = None
    reference: Path | None = None
    reference_l1_key: str = "l1"
    reference_l2_key: str = "l2"
    reference_counts_layer: str = "counts"
    min_shared_genes: int = Field(default=30, ge=3)
    max_dense_elements: int = Field(default=20_000_000, ge=1000)
    max_iterations: int = Field(default=3, ge=1, le=20)
    movement_tolerance: float = Field(default=.005, ge=0, le=1)
    min_pool_size: int = Field(default=12, ge=4)
    min_cluster_size: int = Field(default=5, ge=2)
    min_markers: int = Field(default=2, ge=1)
    detection_fraction: float = Field(default=.2, gt=0, le=1)
    margin: float = Field(default=.5, ge=0)
    subtype_margin: float = Field(default=.5, ge=0)
    atlas_purity: float = Field(default=.7, gt=.5, le=1)
    atlas_vote: float = Field(default=.6, gt=.5, le=1)
    n_pcs: int = Field(default=15, ge=2)
    n_neighbors: int = Field(default=15, ge=2)
    n_hvg: int = Field(default=2000, ge=3)
    resolution: float = Field(default=.5, gt=0, le=5)
    spatial_level: Literal["l1", "l2"] = "l1"
    mode: Literal["cell", "domain"] = "cell"

    @model_validator(mode="after")
    def ontology(self) -> "HierarchyConfig":
        labels = [label for pool in self.compartments.values() for label in pool]
        if not labels or len(labels) != len(set(labels)) or any(not x for x in labels):
            raise ValueError("Compartment L1 assignments must be nonempty and disjoint")
        if set(labels) != set(self.l1_markers) or len(labels) < 2:
            raise ValueError("Every L1 requires a marker set; at least two L1 types required")
        if set(self.l2_parent) != set(self.l2_markers) or not set(self.l2_parent.values()) <= set(labels):
            raise ValueError("L2 markers and parent map disagree")
        if any(not pool or not name for name, pool in self.compartments.items()):
            raise ValueError("Empty compartments are invalid")
        for markers in [self.l1_markers, self.l2_markers]:
            if any(not name or not genes or any(not gene for gene in genes)
                   for name, genes in markers.items()):
                raise ValueError("Marker names and lists must be nonempty")
        return self


def load_hierarchy(path: Path) -> HierarchyConfig:
    cfg = HierarchyConfig.model_validate_json(path.read_text())
    if cfg.reference:
        cfg.reference = (path.resolve().parent / cfg.reference).resolve()
    return cfg


def decide_cluster(current: str, best: str | None, strong: bool,
                   atlas: str | None, purity: float, vote: float,
                   min_purity: float = .7, min_vote: float = .6) -> tuple[str, str]:
    if best is None:
        return "Unknown", current
    if not strong:
        return "Ambiguous", current
    confident = atlas is not None and purity >= min_purity and vote >= min_vote
    if confident and atlas != best:
        return "Ambiguous", current
    if best == current:
        return "Keep", current
    if confident:
        return "Reassign", best
    return "Ambiguous", current


def stopping_reason(moves: list[float], tolerance: float, repeated: bool,
                    iteration: int, limit: int) -> str | None:
    if not moves:
        return "exact_stability"
    if repeated:
        return "cycle"
    if max(moves) <= tolerance:
        return "tolerance"
    if iteration >= limit:
        return "max_iterations"
    return None


def module_scores(obj: ad.AnnData, markers: dict[str, list[str]], excluded: set[str],
                  min_markers: int) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Fixed-query standardized marker minus expression-matched control means.

    Up to five nearest-mean nonmarker controls per marker, stable tie order.
    All marker sets are excluded from controls. At least two controls are needed.
    """
    x = sparse.csr_matrix(obj.X)
    means = np.asarray(x.mean(axis=0)).ravel()
    forbidden = set(g for genes in markers.values() for g in genes) | excluded
    controls = np.array([i for i, gene in enumerate(obj.var_names) if gene not in forbidden])
    scores = pd.DataFrame(index=obj.obs_names)
    coverage = {}
    for label, requested in markers.items():
        present = sorted(set(requested).intersection(obj.var_names) - excluded)
        coverage[label] = present
        values = np.full(obj.n_obs, np.nan)
        if len(present) >= min_markers and len(controls) >= 2:
            indices = obj.var_names.get_indexer(present)
            matched = np.unique(np.concatenate([
                controls[np.argsort(np.abs(means[controls] - means[i]), kind="stable")[:5]]
                for i in indices
            ]))
            raw = np.asarray(x[:, indices].mean(axis=1) - x[:, matched].mean(axis=1)).ravel()
            std = raw.std()
            if std > 1e-10:
                values = (raw - raw.mean()) / std
        scores[label] = values
    return scores, coverage


def pool_clusters(obj: ad.AnnData, cfg: HierarchyConfig, seed: int,
                  platform: str) -> tuple[np.ndarray, dict[str, Any]]:
    # Fresh AnnData excludes all prior PCA, graphs, normalization and cluster state.
    work = ad.AnnData(obj.layers["counts"].copy(), obs=obj.obs.copy(), var=obj.var[[]].copy())
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    if platform != "imaging" and work.n_vars > cfg.n_hvg:
        sc.pp.highly_variable_genes(work, n_top_genes=cfg.n_hvg, flavor="seurat")
        work = work[:, work.var.highly_variable].copy()
    variance = np.asarray(work.X.power(2).mean(axis=0)).ravel() - np.square(np.asarray(work.X.mean(axis=0)).ravel())
    if np.max(variance) < 1e-10:
        raise DegeneratePool("No variable genes in pool")
    n_pcs = min(cfg.n_pcs, work.n_obs - 1, work.n_vars - 1)
    sc.pp.pca(work, n_comps=n_pcs, svd_solver="arpack", random_state=seed)
    if not np.isfinite(work.obsm["X_pca"]).all() or np.std(work.obsm["X_pca"]) < 1e-10:
        raise DegeneratePool("Degenerate pool PCA")
    rep, harmony = "X_pca", "not_requested"
    if cfg.batch_key:
        if cfg.batch_key not in work.obs or work.obs[cfg.batch_key].isna().any():
            raise ValueError("Harmony batch metadata missing")
        if work.obs[cfg.batch_key].nunique() > 1:
            import harmonypy
            model = harmonypy.run_harmony(work.obsm["X_pca"], work.obs, cfg.batch_key,
                                          nclust=max(2, min(20, work.n_obs // 10)),
                                          random_state=seed, verbose=False)
            adjusted = np.asarray(model.Z_corr)
            expected = work.obsm["X_pca"].shape
            if adjusted.shape == expected[::-1]:
                adjusted = adjusted.T
            if adjusted.shape != expected or not np.isfinite(adjusted).all():
                raise ValueError("Harmony returned invalid coordinates")
            work.obsm["X_pca_harmony"] = adjusted
            rep, harmony = "X_pca_harmony", "applied"
        else:
            harmony = "single_batch_skipped"
    sc.pp.neighbors(work, use_rep=rep, n_neighbors=min(cfg.n_neighbors, work.n_obs-1),
                    random_state=seed)
    sc.tl.leiden(work, resolution=cfg.resolution, flavor="igraph", directed=False,
                 n_iterations=2, random_state=seed, key_added="cluster")
    return work.obs.cluster.astype(str).to_numpy(), {"harmony": harmony, "n_pcs": n_pcs,
                                                     "n_genes": work.n_vars}


def atlas_evidence(obj: ad.AnnData, cfg: HierarchyConfig, out: Path, seed: int) -> pd.DataFrame:
    evidence = pd.DataFrame({"l1": "Unknown", "l2": "Unknown", "vote": 0.}, index=obj.obs_names)
    if cfg.reference is None:
        return evidence
    reference = ad.read_h5ad(cfg.reference)
    if not reference.obs_names.is_unique or not reference.var_names.is_unique:
        raise ValueError("Reference IDs must be unique")
    for key in (cfg.reference_l1_key, cfg.reference_l2_key):
        if key not in reference.obs or reference.obs[key].isna().any():
            raise ValueError("Missing hierarchical atlas labels")
    labels = reference.obs[[cfg.reference_l1_key, cfg.reference_l2_key]].astype(str)
    if any(labels.groupby(cfg.reference_l2_key)[cfg.reference_l1_key].nunique() > 1):
        raise ValueError("A reference subtype has multiple L1 parents")
    for _, row in labels.drop_duplicates().iterrows():
        l1, l2 = row[cfg.reference_l1_key], row[cfg.reference_l2_key]
        if l2 in cfg.l2_parent and cfg.l2_parent[l2] != l1:
            raise ValueError("Reference subtype parent disagrees with configured ontology")
    query = obj.copy()
    transfer = RunConfig(input=cfg.reference, out=out, tissue="hierarchy",
                         reference=cfg.reference, reference_label_key=cfg.reference_l2_key,
                         reference_counts_layer=cfg.reference_counts_layer,
                         min_shared_genes=cfg.min_shared_genes, n_pcs=cfg.n_pcs, seed=seed,
                         max_dense_elements=cfg.max_dense_elements)
    reference_transfer(query, transfer)
    parent = labels.drop_duplicates().set_index(cfg.reference_l2_key)[cfg.reference_l1_key].to_dict()
    evidence["l2"] = query.obs.oa_reference_candidate.astype(str)
    evidence["l1"] = evidence.l2.map(parent).fillna("Unknown")
    evidence["vote"] = query.obs.oa_reference_vote
    evidence["ontology_known"] = evidence.l2.isin(cfg.l2_parent) & evidence.l1.isin(cfg.l1_markers)
    evidence.to_csv(out / "atlas_predictions.csv", index_label="cell_id")
    return evidence


def _best(scores: pd.DataFrame, indices: np.ndarray, coverage: dict[str, list[str]],
          obj: ad.AnnData, cfg: HierarchyConfig, margin: float) -> tuple[str | None, bool, float]:
    means = scores.iloc[indices].mean().dropna().sort_values(ascending=False, kind="stable")
    if means.empty:
        return None, False, 0.
    best = str(means.index[0])
    delta = float(means.iloc[0] - means.iloc[1]) if len(means) > 1 else 0.
    present = coverage[best]
    detected = np.asarray((obj.layers["counts"][indices][:, obj.var_names.get_indexer(present)] > 0)
                          .mean(axis=0)).ravel()
    supported = int((detected >= cfg.detection_fraction).sum()) >= cfg.min_markers
    # A single available module cannot establish a margin over unmeasurable alternatives.
    return best, bool(supported and len(means) >= 2 and delta > 1e-10 and delta >= margin and means.iloc[0] > 0), delta


def run_hierarchy(obj: ad.AnnData, cfg: HierarchyConfig, out: Path,
                  seed: int, platform: str) -> dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    if platform == "visium" and cfg.mode != "domain":
        raise ValueError("Visium hierarchy requires domain mode; spots are not cells")
    if cfg.initial_l1_key not in obj.obs or obj.obs[cfg.initial_l1_key].isna().any():
        raise ValueError("Initial L1 labels missing")
    if cfg.batch_key and (cfg.batch_key not in obj.obs or obj.obs[cfg.batch_key].isna().any()):
        raise ValueError("Harmony batch metadata missing")
    raw = counts(obj, "counts")
    fixed = ad.AnnData(raw.copy(), obs=obj.obs[[]].copy(), var=obj.var[[]].copy())
    fixed.layers["counts"] = raw.copy()
    sc.pp.normalize_total(fixed, target_sum=1e4)
    sc.pp.log1p(fixed)
    l1_scores, l1_coverage = module_scores(fixed, cfg.l1_markers, set(cfg.excluded_genes), cfg.min_markers)
    l2_scores, l2_coverage = module_scores(fixed, cfg.l2_markers, set(), cfg.min_markers)
    l1_scores.to_csv(out / "l1_scores.csv", index_label="cell_id")
    l2_scores.to_csv(out / "l2_scores.csv", index_label="cell_id")
    write_json(out / "marker_coverage.json", {"l1": l1_coverage, "l2": l2_coverage,
                                              "excluded": cfg.excluded_genes})
    atlas = atlas_evidence(obj, cfg, out, seed)
    current = obj.obs[cfg.initial_l1_key].astype(str).to_numpy().copy()
    obj.obs["oa_initial_l1"] = current.copy()
    l2 = np.full(obj.n_obs, "Unknown", dtype=object)
    status = np.full(obj.n_obs, "Unknown", dtype=object)
    cluster_ids = np.full(obj.n_obs, "unprocessed", dtype=object)
    pool_for = {label: pool for pool, labels in cfg.compartments.items() for label in labels}
    seen = {tuple(current)}
    history, all_decisions, dot_rows = [], [], []
    reason = "max_iterations"
    for iteration in range(1, cfg.max_iterations+1):
        pools = np.array([pool_for.get(x, "__unresolved__") for x in current])
        new = current.copy()
        moves = []
        round_rows, pool_rows = [], []
        l2[:] = "Unknown"
        for pool in sorted(set(pools)):
            members = np.flatnonzero(pools == pool)
            local = obj[members].copy()
            info: dict[str, Any] = {"pool": pool, "size": len(members)}
            if len(members) < cfg.min_pool_size:
                status[members] = "Unknown"
                cluster_ids[members] = f"r{iteration}:{pool}:tiny"
                info["reason"] = "tiny_pool"
                pool_rows.append(info)
                continue
            # Invalid batch/dependency errors must propagate, not silently skip Harmony.
            try:
                clusters, details = pool_clusters(local, cfg, seed+iteration, platform)
            except DegeneratePool:
                status[members] = "Unknown"
                cluster_ids[members] = f"r{iteration}:{pool}:degenerate"
                info["reason"] = "degenerate_pool"
                pool_rows.append(info)
                continue
            info.update(details)
            pool_rows.append(info)
            for cluster in sorted(set(clusters)):
                indices = members[clusters == cluster]
                cid = f"r{iteration}:{pool}:{cluster}"
                cluster_ids[indices] = cid
                best, strong, delta = _best(l1_scores, indices, l1_coverage, fixed, cfg, cfg.margin)
                if len(indices) < cfg.min_cluster_size:
                    best, strong = None, False
                predictions = atlas.iloc[indices]
                counts_l1 = predictions.l1.value_counts()
                atlas_l1: str | None = str(counts_l1.index[0])
                purity = float(counts_l1.iloc[0] / len(indices))
                votes = predictions.loc[predictions.l1 == atlas_l1, "vote"]
                vote = float(votes.median())
                if atlas_l1 not in cfg.l1_markers:
                    atlas_l1 = None
                counts_l2 = predictions.l2.value_counts()
                atlas_l2 = str(counts_l2.index[0])
                l2_purity = float(counts_l2.iloc[0] / len(indices))
                l2_vote = float(predictions.loc[predictions.l2 == atlas_l2, "vote"].median())
                changed = 0
                for previous in sorted(set(current[indices])):
                    subset = indices[current[indices] == previous]
                    decision, target = decide_cluster(previous, best, strong, atlas_l1, purity, vote,
                                                       cfg.atlas_purity, cfg.atlas_vote)
                    new[subset] = target
                    status[subset] = decision
                    changed += int(np.sum(current[subset] != target))
                    subtype = None
                    if decision == "Keep":
                        candidates = [name for name, parent in cfg.l2_parent.items() if parent == target]
                        if candidates:
                            subtype, supported, _ = _best(l2_scores[candidates], indices,
                                                          l2_coverage, fixed, cfg, cfg.subtype_margin)
                            subtype_conflict = (l2_purity >= cfg.atlas_purity and l2_vote >= cfg.atlas_vote
                                                and atlas_l2 != subtype)
                            if supported and not subtype_conflict:
                                l2[subset] = subtype
                    round_rows.append({"iteration": iteration, "cluster": cid, "pool": pool,
                                       "n_cells": len(subset), "source_l1": previous,
                                       "target_l1": target, "marker_best": best,
                                       "marker_margin": delta, "marker_strong": strong,
                                       "atlas_l1": atlas_l1, "atlas_purity": purity,
                                       "atlas_vote": vote, "atlas_l2": atlas_l2, "atlas_l2_purity": l2_purity,
                                       "atlas_l2_vote": l2_vote, "decision": decision})
                if changed:
                    moves.append(changed / len(members))
                for level, marker_sets in [("l1", l1_coverage), ("l2", l2_coverage)]:
                    for label, genes in marker_sets.items():
                        for gene in genes:
                            col = fixed.var_names.get_loc(gene)
                            dot_rows.append({"iteration": iteration, "cluster": cid,
                                             "level": level, "label": label, "gene": gene,
                                             "mean_log_expression": float(fixed.X[indices, col].mean()),
                                             "fraction_detected": float((raw[indices, col] > 0).mean())})
        moved = current != new
        # Moved cells cannot retain a subtype obtained in their old pool.
        l2[moved] = "Unknown"
        signature = tuple(new)
        reason_now = stopping_reason(moves, cfg.movement_tolerance, signature in seen,
                                     iteration, cfg.max_iterations)
        movement_by_pool = {pool: float(moved[pools == pool].mean()) for pool in set(pools)}
        history.append({"iteration": iteration, "moved_cells": int(moved.sum()),
                        "max_cluster_fraction": max(moves, default=0),
                        "total_fraction_by_pool": movement_by_pool,
                        "pools": pool_rows, "stop_reason": reason_now})
        pd.DataFrame({"cell_id": obj.obs_names, "pool": pools, "cluster": cluster_ids,
                      "before_l1": current, "after_l1": new, "l2": l2,
                      "decision": status}).to_csv(out / f"round_{iteration}_cells.csv", index=False)
        all_decisions.extend(round_rows)
        current = new
        seen.add(signature)
        if reason_now:
            reason = reason_now
            break
    obj.obs["oa_l1"] = current
    obj.obs["oa_l2"] = l2
    obj.obs["oa_hierarchy_status"] = status
    obj.obs["oa_hierarchy_cluster"] = cluster_ids
    # Only audited Keep calls can drive final enrichment, including at L1.
    selected = current if cfg.spatial_level == "l1" else l2
    obj.obs["oa_hierarchy_label"] = np.where((status == "Keep") & (selected != "Unknown"), selected, "Uncertain")
    decision_columns = ["iteration", "cluster", "pool", "n_cells", "source_l1", "target_l1",
                        "marker_best", "marker_margin", "marker_strong", "atlas_l1", "atlas_purity",
                        "atlas_vote", "atlas_l2", "atlas_l2_purity", "atlas_l2_vote", "decision"]
    pd.DataFrame(all_decisions, columns=decision_columns).to_csv(out / "decisions.csv", index=False)
    dot_columns = ["iteration", "cluster", "level", "label", "gene", "mean_log_expression", "fraction_detected"]
    pd.DataFrame(dot_rows, columns=dot_columns).to_csv(out / "dotplot_data.csv", index=False)
    write_json(out / "iterations.json", history)
    summary = {"iterations": len(history), "stop_reason": reason,
               "converged": reason in {"exact_stability", "tolerance"} and bool(np.any(status == "Keep")), "mode": cfg.mode,
               "status_counts": pd.Series(status).value_counts().to_dict(),
               "reference_sha256": sha256(cfg.reference) if cfg.reference else None,
               "reference_backend": "reference_fitted_pca_knn" if cfg.reference else "none",
               "real_atlas_validated": False, "spatial_level": cfg.spatial_level}
    write_json(out / "summary.json", summary)
    return summary
