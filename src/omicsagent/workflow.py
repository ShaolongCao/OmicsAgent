import importlib.metadata
import json
import platform
import shutil
import time
import uuid
from pathlib import Path
from typing import Any, cast

import anndata as ad
import pandas as pd

from . import __version__
from .analysis import markers, preprocess, reference_transfer, validate
from .annotation import audit_labels, marker_candidates, panel
from .config import RunConfig
from .evidence import EvidenceStore
from .hierarchy import load_hierarchy, run_hierarchy
from .methods import annotate_reference
from .report import make_report
from .spatial import spatial_analysis
from .trace import Trace, sha256, write_json


class Workflow:
    """Configured paths are fixed; the agent can only use this run's tools."""

    def __init__(self, config: RunConfig):
        self.cfg = config
        if config.out.exists() and any(config.out.iterdir()):
            raise ValueError("Output directory is not empty; choose a new directory")
        if config.out.resolve() == config.input.resolve():
            raise ValueError("Input and output must differ")
        config.out.mkdir(parents=True, exist_ok=True)
        self.run_id = uuid.uuid4().hex
        self.trace = Trace(config.out / "trace.jsonl", self.run_id)
        self.replay: Path | None = None
        self.stage = "new"
        self.adata: Any = None
        self.marker_table = pd.DataFrame()
        self.annotations: list[dict[str, Any]] = []
        self.summary: dict[str, Any] = {}
        self.warnings = [
            "Marker detection is not proof of identity; all annotations are provisional.",
            "No doublet removal or biological replicate DE is performed. Optional scANVI and SPLIT outputs require independent review.",
        ]
        if config.platform == "imaging" and config.normalization == "total":
            self.warnings.append(
                "Targeted-panel total normalization may distort expression; compare area normalization when cell area is available."
            )
        if config.platform == "visium":
            self.warnings.append(
                "Visium spots may contain multiple cells. Spatial enrichment uses expression clusters, not inferred single-cell labels."
            )
        self.evidence = EvidenceStore(config.out / "evidence", config.online_evidence)
        self.marker_panel = panel(config)
        write_json(config.out / "config.json", config.model_dump(mode="json"))
        self.trace.emit("created", version=__version__, config=config.model_dump(mode="json"))

    def require(self, stage: str) -> None:
        if self.stage != stage:
            raise ValueError(f"Expected stage={stage}; current={self.stage}")

    def inspect_data(self) -> dict[str, Any]:
        self.require("new")
        self.adata = ad.read_h5ad(self.cfg.input)
        self.summary["input"] = validate(self.adata, self.cfg)
        if 'oa_review_provenance_json' in self.adata.uns:
            write_json(self.cfg.out / 'provenance.json',
                       json.loads(self.adata.uns['oa_review_provenance_json']))
        self.stage = "loaded"
        return cast(dict[str, Any], self.summary["input"])

    def preprocess_data(self) -> dict[str, Any]:
        self.require("loaded")
        self.adata, self.summary["qc"] = preprocess(self.adata, self.cfg)
        if self.summary["qc"]["mt_genes_found"] == 0:
            self.warnings.append(
                "No mitochondrial-prefix genes found; mt QC is uninformative. Check gene identifiers/panel."
            )
        self.adata.write_h5ad(self.cfg.out / "checkpoint.h5ad", compression='gzip')
        self.stage = "processed"
        return cast(dict[str, Any], self.summary["qc"])

    def find_markers(self) -> dict[str, Any]:
        self.require("processed")
        self.marker_table = markers(self.adata, self.cfg)
        self.marker_table.to_csv(self.cfg.out / "markers.csv", index=False)
        if self.marker_table.empty:
            self.warnings.append(
                "Marker ranking unavailable: too few eligible clusters or observations."
            )
        self.stage = "markers"
        return {
            "rows": len(self.marker_table),
            "top_markers": self.marker_table.groupby("cluster", observed=True)
            .head(10)
            .fillna(0)
            .to_dict("records"),
        }

    def annotate(self) -> dict[str, Any]:
        self.require("markers")
        self.annotations = marker_candidates(self.adata, self.marker_panel)
        if self.cfg.reference:
            if self.replay and self.cfg.annotation_method != "knn":
                saved = ad.read_h5ad(self.replay / "result.h5ad")
                if not saved.obs_names.equals(self.adata.obs_names) or not saved.var_names.equals(self.adata.var_names):
                    raise ValueError("Replay dimensions/IDs differ")
                self.adata.obsm["oa_weights"] = saved.obsm["oa_weights"].copy()
                for column in saved.obs.columns:
                    if column.startswith("oa_weight") or column in {"oa_candidate", "oa_split_removed_fraction"}:
                        self.adata.obs[column] = saved.obs[column].copy()
                if self.cfg.split_enabled:
                    self.adata.layers["counts_split"] = saved.layers["counts_split"].copy()
                    self.adata.var["oa_split_reference_covered"] = saved.var["oa_split_reference_covered"].copy()
                    self.adata.obs[["oa_candidate", "oa_split_removed_fraction"]].to_csv(self.cfg.out / "split_qc.csv")
                self.summary["reference"] = json.loads((self.replay / "summary.json").read_text())["reference"]
                self.adata.obsm["oa_weights"].to_csv(self.cfg.out / "annotation_weights.csv.gz", index_label="cell_id")
            elif self.cfg.annotation_method == "knn":
                self.summary["reference"] = reference_transfer(self.adata, self.cfg)
            else:
                self.summary["reference"] = annotate_reference(self.adata, self.cfg)
        elif self.cfg.annotation_key:
            self.adata.obs["oa_candidate"] = (
                self.adata.obs[self.cfg.annotation_key]
                .astype(object)
                .fillna("Uncertain")
                .astype(str)
            )
        if self.cfg.reference or self.cfg.annotation_key:
            for record in self.annotations:
                mask = self.adata.obs["oa_cluster"].astype(str) == record["cluster"]
                candidates = self.adata.obs.loc[mask, "oa_candidate"].astype(str)
                record.update(
                    candidate=candidates.value_counts().index[0],
                    markers=[],
                    status="provisional",
                    method=f"reference_{self.cfg.annotation_method}" if self.cfg.reference else "provided_annotation",
                )
        self._audit()
        self.stage = "annotated"
        return {"annotations": self.annotations, "gate": self.summary["gate"]}

    def _audit(self) -> None:
        gate, sparsity = audit_labels(self.adata, self.marker_panel, self.cfg.min_marker_support)
        gate.to_csv(self.cfg.out / "annotation_gate.csv", index=False)
        sparsity.to_csv(self.cfg.out / "sparsity_diagnostics.csv", index=False)
        self.summary["gate"] = self.adata.obs["oa_gate_status"].value_counts().to_dict()

    def search_literature(self, query: str) -> dict[str, Any]:
        if self.stage not in {"markers", "annotated"}:
            raise ValueError(
                "Retrieve evidence after marker analysis and before spatial/report steps"
            )
        return self.evidence.search(query)

    def lookup_gene(self, symbol: str) -> dict[str, Any]:
        if self.stage not in {"markers", "annotated"}:
            raise ValueError("Gene lookup requires completed marker analysis")
        return self.evidence.lookup_gene(symbol, self.cfg.species)

    def propose_annotation(
        self, cluster: str, label: str, genes: list[str], evidence_ids: list[str], rationale: str
    ) -> dict[str, Any]:
        self.require("annotated")
        if self.cfg.reference or self.cfg.annotation_key or self.cfg.platform == "visium":
            raise ValueError(
                "v0.1 agent cannot overwrite provided/reference labels or label Visium spots as cells"
            )
        if cluster not in self.adata.obs["oa_cluster"].astype(str).unique():
            raise ValueError("Unknown cluster")
        if not isinstance(label, str) or not 1 <= len(label) <= 100 or len(rationale) > 2000:
            raise ValueError("Invalid label or rationale length")
        if label != "Uncertain":
            allowed = set(
                self.marker_table.loc[self.marker_table.cluster.astype(str) == cluster, "gene"]
            )
            if len(set(genes)) < 2 or not set(genes) <= allowed:
                raise ValueError("A proposal needs at least two distinct observed top markers")
            if not evidence_ids or any(e not in self.evidence.records for e in evidence_ids):
                raise ValueError("Citations must refer to retrieved evidence IDs")
            if any(not self.evidence.records[e].get("abstract") for e in evidence_ids):
                raise ValueError("Title-only sources cannot support a label proposal")
        self.adata.obs["oa_candidate"] = self.adata.obs["oa_candidate"].astype(str)
        mask = self.adata.obs["oa_cluster"].astype(str) == cluster
        self.adata.obs.loc[mask, "oa_candidate"] = label
        # Model-derived marker lists do NOT become independent gate rules.
        record = {
            "cluster": cluster,
            "candidate": label,
            "markers": genes,
            "status": "unresolved" if label == "Uncertain" else "provisional",
            "evidence_ids": evidence_ids,
            "rationale": rationale,
            "method": "claude_proposal_requires_expert_review",
        }
        self.annotations = [record if r["cluster"] == cluster else r for r in self.annotations]
        self._audit()
        return record

    def refine_hierarchy(self) -> dict[str, Any]:
        self.require("annotated")
        if self.cfg.hierarchy is None:
            raise ValueError("Hierarchy configuration is not enabled")
        if self.replay:
            saved = ad.read_h5ad(self.replay / "result.h5ad")
            if not saved.obs_names.equals(self.adata.obs_names) or not saved.var_names.equals(self.adata.var_names):
                raise ValueError("Hierarchy replay IDs differ")
            for column in ("oa_initial_l1", "oa_l1", "oa_l2", "oa_hierarchy_status", "oa_hierarchy_cluster", "oa_hierarchy_label"):
                self.adata.obs[column] = saved.obs[column].copy()
            shutil.copytree(self.replay / "hierarchy", self.cfg.out / "hierarchy")
            self.summary["hierarchy"] = json.loads((self.replay / "hierarchy" / "summary.json").read_text())
            self.stage = "hierarchy"
            return cast(dict[str, Any], self.summary["hierarchy"])
        hc = load_hierarchy(self.cfg.hierarchy)
        self.summary["hierarchy"] = run_hierarchy(
            self.adata, hc, self.cfg.out / "hierarchy", self.cfg.seed, self.cfg.platform)
        self.stage = "hierarchy"
        return cast(dict[str, Any], self.summary["hierarchy"])

    def analyze_spatial(self) -> dict[str, Any]:
        self.require("hierarchy" if self.cfg.hierarchy else "annotated")
        if self.cfg.modality != "spatial":
            raise ValueError("Spatial tool is unavailable for scrna modality")
        result = spatial_analysis(self.adata, self.cfg)
        for name, frame in result.pop("tables").items():
            frame.to_csv(self.cfg.out / f"{name}.csv", index=False)
        self.summary["spatial"] = result
        self.stage = "spatial"
        return result

    def generate_report(self) -> dict[str, Any]:
        self.require("spatial" if self.cfg.modality == "spatial" else ("hierarchy" if self.cfg.hierarchy else "annotated"))
        self.adata.write_h5ad(self.cfg.out / "result.h5ad", compression='gzip')
        self.adata.obs.to_csv(self.cfg.out / "cells.csv", index_label="cell_id")
        write_json(self.cfg.out / "annotations.json", self.annotations)
        write_json(self.cfg.out / "evidence_records.json", self.evidence.records)
        write_json(self.cfg.out / "summary.json", self.summary)
        make_report(
            self.adata,
            self.cfg,
            self.summary,
            self.annotations,
            self.evidence.records,
            self.warnings,
        )
        self.stage = "completed"
        return {"status": self.stage, "report": str(self.cfg.out / "report.md")}

    def call(self, name: str, args: dict[str, Any]) -> dict[str, Any]:
        allowed = {
            "inspect_data",
            "preprocess_data",
            "find_markers",
            "annotate",
            "search_literature",
            "lookup_gene",
            "propose_annotation",
            "analyze_spatial",
            "refine_hierarchy",
            "generate_report",
        }
        start = time.monotonic()
        try:
            if name not in allowed:
                raise ValueError("Unknown tool")
            value = getattr(self, name)(**args)
            result = {"ok": True, "data": value, "stage": self.stage}
        except Exception as exc:
            result = {
                "ok": False,
                "error": {"type": type(exc).__name__, "message": str(exc)},
                "stage": self.stage,
            }
        self.trace.emit(
            "tool", name=name, arguments=args, result=result, elapsed_s=time.monotonic() - start
        )
        self.manifest("completed" if self.stage == "completed" else "running")
        return result

    def manifest(self, status: str) -> None:
        dependencies = {
            name: importlib.metadata.version(name)
            for name in ["scanpy", "anndata", "numpy", "scipy", "pandas", "anthropic"]
        }
        for optional in ("scvi-tools", "tacco", "torch", "harmonypy"):
            try:
                dependencies[optional] = importlib.metadata.version(optional)
            except importlib.metadata.PackageNotFoundError:
                pass
        sources = {
            "input": {
                "path": str(self.cfg.input),
                "sha256": sha256(self.cfg.input) if self.cfg.input.is_file() else None,
            }
        }
        for name in ("reference", "marker_file", "hierarchy"):
            path = getattr(self.cfg, name)
            if path:
                sources[name] = {
                    "path": str(path),
                    "sha256": sha256(path) if path.is_file() else None,
                }
        if self.cfg.hierarchy:
            ref = load_hierarchy(self.cfg.hierarchy).reference
            if ref:
                sources["hierarchy_reference"] = {"path": str(ref), "sha256": sha256(ref)}
        artifacts = {
            str(p.relative_to(self.cfg.out)): sha256(p)
            for p in self.cfg.out.rglob("*")
            if p.is_file()
            and p.name != "manifest.json"
            and "replay-" not in str(p.relative_to(self.cfg.out))
        }
        code_dir = Path(__file__).parent
        write_json(
            self.cfg.out / "manifest.json",
            {
                "schema_version": 1,
                "run_id": self.run_id,
                "status": status,
                "stage": self.stage,
                "version": __version__,
                "python": platform.python_version(),
                "platform": platform.platform(),
                "dependencies": dependencies,
                "sources": sources,
                "artifacts": artifacts,
                "code_sha256": {p.name: sha256(p) for p in code_dir.glob("*.py")},
                "warnings": self.warnings,
            },
        )


def run_pipeline(cfg: RunConfig, replay: Path | None = None) -> dict[str, Any]:
    if replay:
        previous_config = json.loads((replay / "config.json").read_text())
        current_config = cfg.model_dump(mode="json")
        previous_config.pop("out", None)
        current_config.pop("out", None)
        if previous_config != current_config:
            raise ValueError("Replay config differs from saved analysis")
        saved_manifest = json.loads((replay / "manifest.json").read_text())
        for name, digest in saved_manifest["artifacts"].items():
            if sha256(replay / name) != digest:
                raise ValueError(f"Replay artifact changed: {name}")
        for record in saved_manifest["sources"].values():
            if record["sha256"] != sha256(Path(record["path"])):
                raise ValueError("Replay input/reference content changed")
    workflow = Workflow(cfg)
    workflow.replay = replay
    sequence = ["inspect_data", "preprocess_data", "find_markers", "annotate"]
    if cfg.hierarchy:
        sequence.append("refine_hierarchy")
    if cfg.modality == "spatial":
        sequence.append("analyze_spatial")
    sequence.append("generate_report")
    try:
        for name in sequence:
            result = workflow.call(name, {})
            if not result["ok"]:
                raise RuntimeError(f"{name}: {result['error']}")
            if name == "annotate" and replay:
                # Replay saved decisions, never rerun external retrieval/model calls.
                decisions = json.loads((replay / "annotations.json").read_text())
                previous = pd.read_csv(replay / "cells.csv", index_col="cell_id")
                if set(previous.index) != set(workflow.adata.obs_names):
                    raise ValueError("Replay cells differ; environment/input changed")
                if (
                    not previous.loc[workflow.adata.obs_names, "oa_cluster"]
                    .astype(str)
                    .equals(workflow.adata.obs["oa_cluster"].astype(str))
                ):
                    raise ValueError("Replay clusters differ; do not apply saved decisions")
                workflow.adata.obs["oa_candidate"] = previous.loc[
                    workflow.adata.obs_names, "oa_candidate"
                ].astype(str)
                workflow.annotations = decisions
                workflow.evidence.records = json.loads(
                    (replay / "evidence_records.json").read_text()
                )
                workflow._audit()
                workflow.trace.emit("replay", source=str(replay.resolve()))
        return {"status": workflow.stage, "out": str(cfg.out)}
    except Exception:
        workflow.manifest("failed")
        raise
