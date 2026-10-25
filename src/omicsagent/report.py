# mypy: disable-error-code="no-untyped-call"
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd

from .config import RunConfig


def scatter_by_label(
    xy: np.ndarray, labels: pd.Series, path: Path, title: str, spatial: bool = False
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    for label in sorted(labels.astype(str).unique()):
        keep = labels.astype(str).to_numpy() == label
        ax.scatter(xy[keep, 0], xy[keep, 1], s=5, alpha=0.7, label=label, rasterized=True)
    ax.set_title(title)
    ax.set_xlabel("x" if spatial else "UMAP 1")
    ax.set_ylabel("y" if spatial else "UMAP 2")
    if spatial:
        ax.set_aspect("equal")
        ax.invert_yaxis()
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=7, markerscale=2)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def make_report(
    adata: ad.AnnData,
    cfg: RunConfig,
    summary: dict[str, Any],
    annotations: list[dict[str, Any]],
    evidence: dict[str, Any],
    warnings: list[str],
) -> None:
    figures = cfg.out / "figures"
    figures.mkdir(exist_ok=True)
    scatter_by_label(
        adata.obsm["X_umap"],
        adata.obs["oa_cluster"],
        figures / "umap_clusters.png",
        "Expression clusters",
    )
    scatter_by_label(
        adata.obsm["X_umap"],
        adata.obs["oa_hierarchy_label"] if cfg.hierarchy else adata.obs["oa_label"],
        figures / "umap_labels.png",
        "Candidate labels after marker audit",
    )
    fig, axes = plt.subplots(1, 3, figsize=(12, 3))
    for ax, key in zip(axes, ["oa_total_counts", "oa_n_genes", "oa_pct_mt"], strict=True):
        ax.hist(adata.obs[key], bins=30)
        ax.set_title(key.removeprefix("oa_"))
    fig.suptitle("Post-QC distributions; filtering counts in summary.json")
    fig.tight_layout()
    fig.savefig(figures / "qc.png", dpi=130)
    plt.close(fig)
    image_links = [
        "![QC](figures/qc.png)",
        "![Clusters](figures/umap_clusters.png)",
        "![Labels](figures/umap_labels.png)",
    ]
    if cfg.modality == "spatial":
        # Separate panels/files prevent superimposing independent coordinate systems.
        for i, sample in enumerate(sorted(adata.obs[cfg.sample_key].astype(str).unique())):
            keep = adata.obs[cfg.sample_key].astype(str).to_numpy() == sample
            name = f"spatial_{i}.png"
            scatter_by_label(
                np.asarray(adata.obsm[cfg.spatial_key])[keep],
                adata.obs.loc[keep, "oa_spatial_label"],
                figures / name,
                f"Sample {sample} ({cfg.coordinate_units}); {cfg.platform}",
                True,
            )
            image_links.append(f"![Spatial sample {i}](figures/{name})")
    if cfg.hierarchy:
        table = pd.read_csv(cfg.out / "hierarchy" / "dotplot_data.csv")
        if not table.empty:
            last = table[table.iteration == table.iteration.max()]
            for level in ("l1", "l2"):
                frame = last[last.level == level].drop_duplicates(["cluster", "gene"])
                genes = sorted(frame.gene.unique())[:40]
                clusters = sorted(frame.cluster.unique())
                frame = frame[frame.gene.isin(genes)]
                if frame.empty:
                    continue
                fig, ax = plt.subplots(figsize=(max(7, len(genes)*.35), max(3, len(clusters)*.35)))
                dots = ax.scatter(frame.gene.map({g:i for i,g in enumerate(genes)}),
                                  frame.cluster.map({c:i for i,c in enumerate(clusters)}),
                                  s=frame.fraction_detected*100, c=frame.mean_log_expression,
                                  cmap="viridis")
                ax.set_xticks(range(len(genes)), genes, rotation=90)
                ax.set_yticks(range(len(clusters)), clusters, fontsize=7)
                ax.set_title(f"Final-round {level.upper()} markers; dot size = detection fraction")
                fig.colorbar(dots, ax=ax, label="Mean log expression")
                fig.tight_layout()
                name = f"hierarchy_{level}_dotplot.png"
                fig.savefig(figures / name, dpi=130)
                plt.close(fig)
                image_links.append(f"![Hierarchy {level} markers](figures/{name})")
    lines = [
        "# OmicsAgent analysis",
        "",
        f"Modality: **{cfg.modality} / {cfg.platform}**. Species: {cfg.species}. Tissue: {cfg.tissue}.",
        "",
        "This is an exploratory analysis. Candidate labels and retrieved evidence require expert review.",
        "",
        f"QC retained {adata.n_obs} observations and {adata.n_vars} genes; {summary['qc']['n_clusters']} expression clusters.",
        "",
        "## Interpretation limits",
        "",
    ]
    if cfg.hierarchy:
        hierarchy = summary["hierarchy"]
        lines += [f"Hierarchy mode: **{hierarchy['mode']}**; rounds: {hierarchy['iterations']}; "
                  f"stop reason: **{hierarchy['stop_reason']}**. See `hierarchy/decisions.csv`, "
                  "`hierarchy/iterations.json` and `hierarchy/dotplot_data.csv`. "
                  "Atlas transfer supports marker decisions; it is not a calibrated confidence measure. "
                  "Real ISD atlas validation is pending.", ""]
    if cfg.reference:
        lines += [f"Reference annotation: **{cfg.annotation_method}**. Soft weights are not calibrated confidence.", ""]
    if cfg.split_enabled:
        lines += ["SPLIT: Python port of the standard RCTD-free core (not validated against native R). "
                  "Corrected fractional counts are in `layers['counts_split']`; clustering, markers and "
                  "spatial statistics retain the original analysis. See `split_qc.csv` and `summary.json`.", ""]
    lines += [f"- {warning}" for warning in dict.fromkeys(warnings)]
    lines += [
        "",
        "## Cluster annotation decisions",
        "",
        "| Cluster | Candidate | Status | Markers | Evidence |",
        "|---|---|---|---|---|",
    ]
    for row in annotations:
        refs = ", ".join(
            f"[{i + 1}]({evidence[e]['url']})"
            for i, e in enumerate(row.get("evidence_ids", []))
            if e in evidence
        )

        def clean(value: Any) -> str:
            return str(value).replace("|", "/").replace("\n", " ")

        lines.append(
            f"| {clean(row['cluster'])} | {clean(row['candidate'])} | {row['status']} | {clean(', '.join(row.get('markers', [])))} | {refs or 'Not independently verified'} |"
        )
    lines += [
        "",
        "Per-cell reference/existing labels may differ from the dominant label in this summary; see cells.csv.",
        "",
        "## Figures",
        "",
        *image_links,
        "",
        "## Reproducibility",
        "",
        "See config.json, manifest.json, trace.jsonl and analysis.ipynb. Counts are preserved in result.h5ad layers['counts']; markers use full-gene log1p expression.",
        "Notebook reruns numerical tools offline from the original input, then replays saved annotation decisions. No API key is needed for replay.",
        "",
        "## Tables",
        "",
        "- markers.csv: exploratory post-clustering rank tests, not replicate-aware disease DE.",
        "- annotation_gate.csv and sparsity_diagnostics.csv: marker availability, gate sensitivity, depth association.",
    ]
    if cfg.modality == "spatial":
        lines += [
            "- neighborhood_enrichment.csv: sample-specific label permutations; BH across emitted hypotheses.",
            "- spatial_moran.csv: positive spatial autocorrelation for a variance-selected gene subset, not a genome-wide SVG screen.",
            "- center_niches.csv: optional fixed-center neighborhood incidence test; no paired disease DE is implemented.",
            "- Visium labels describe expression domains, not individual cells; no spot deconvolution is performed.",
        ]
    (cfg.out / "report.md").write_text("\n".join(lines) + "\n")
    notebook = nbformat.v4.new_notebook()
    notebook.metadata.kernelspec = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook.cells = [
        nbformat.v4.new_markdown_cell(
            "# OmicsAgent reproducible replay\nRun with OmicsAgent installed. Original input files must remain accessible. This reruns numerical analysis and replays recorded decisions; it does not call Claude."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\nimport json, uuid\nfrom omicsagent.config import RunConfig\nfrom omicsagent.workflow import run_pipeline\n"
            "base = Path.cwd()\nassert (base / 'config.json').exists(), 'Run notebook from its output directory'\n"
            "config = RunConfig.model_validate_json((base / 'config.json').read_text())\n"
            "config.out = base / ('replay-' + uuid.uuid4().hex[:8])\nconfig.online_evidence = False\n"
            "result = run_pipeline(config, replay=base)\nresult"
        ),
        nbformat.v4.new_code_cell(
            "import pandas as pd\npd.read_csv(config.out / 'markers.csv').head(20)"
        ),
        nbformat.v4.new_code_cell(
            "from IPython.display import Image, display\ndisplay(Image(filename=str(config.out / 'figures/umap_labels.png')))"
        ),
    ]
    nbformat.write(notebook, cfg.out / "analysis.ipynb")
