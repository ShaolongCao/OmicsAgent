# Iterative hierarchy: a complete worked example

## Run synthetic scRNA and imaging data

From the repository root:

```bash
uv sync --frozen --no-editable --extra dev --extra hierarchy
source .venv/bin/activate
omicsagent hierarchy-demo --out data/my-hierarchy
omicsagent run --config data/my-hierarchy/scrna.json
omicsagent run --config data/my-hierarchy/imaging.json
```

The generator creates 240 query cells and 240 independent reference cells with two L1 groups and four L2 types. Sixty query cells begin in the wrong L1 group. The input has 116 genes and two spatial samples. Truth labels exist only for evaluating the synthetic demonstration; hierarchy reads `initial_l1` and the separately generated reference.

Run `visium.json` for a synthetic spot-domain mechanics example. It uses the same simulated expression solely to exercise the domain-mode code path; it is not a realistic mixture/spot benchmark.

## Read the outputs

- `report.md`: global analysis, final audited hierarchy labels and marker dotplots.
- `cells.csv`: `oa_initial_l1`, `oa_l1`, `oa_l2`, `oa_hierarchy_status`, `oa_hierarchy_cluster`, `oa_hierarchy_label`.
- `hierarchy/round_N_cells.csv`: original round pool, cluster, before/after L1 and decision for every observation.
- `hierarchy/decisions.csv`: marker winner/margin, atlas L1/L2 majority/purity/vote and decision.
- `hierarchy/iterations.json`: pool sizes, Harmony status, movement fractions and stopping reasons.
- `hierarchy/l1_scores.csv`, `l2_scores.csv`, `marker_coverage.json`: fixed-query module evidence and unavailable markers.
- `hierarchy/atlas_predictions.csv`: per-cell independent reference predictions. Votes are not calibrated confidence.
- `hierarchy/dotplot_data.csv`: per-round mean log expression and detection fraction. Figures show the final round and at most 40 genes per level; complete values remain in the table.
- `result.h5ad`: preserved raw counts, all annotation columns and sample-separated spatial graph when applicable.
- `analysis.ipynb`: replay of validated saved decisions, without refitting hierarchy or optional models. Run it in the same package environment with original inputs available.

```python
import anndata as ad
x = ad.read_h5ad("data/my-hierarchy/imaging-output/result.h5ad")
print(x.obs[["oa_initial_l1", "oa_l1", "oa_l2", "oa_hierarchy_status"]].head())
wrong_initially = x.obs["initial_l1"] != x.obs["truth_l1"]
print("Initially misplaced:", int(wrong_initially.sum()))
print("Corrected:", int((x.obs.loc[wrong_initially, "oa_l1"] ==
                          x.obs.loc[wrong_initially, "truth_l1"]).sum()))
```

## Configure a real reference and marker ontology

A run JSON includes `"hierarchy": "hierarchy.json"`. Its initial `annotation_key` can coexist with the independent hierarchy reference. Hierarchy paths resolve against the hierarchy JSON location.

```json
{
  "initial_l1_key": "cell_type_level1",
  "compartments": {
    "myeloid": ["Mac", "DC", "pDC", "Mast", "Langerhans"],
    "lymphoid": ["T_NK", "B", "Plasma"]
  },
  "l1_markers": {},
  "l2_markers": {},
  "l2_parent": {},
  "excluded_genes": ["IRF7", "BST2"],
  "reference": "atlas.h5ad",
  "reference_l1_key": "atlas_l1",
  "reference_l2_key": "atlas_l2",
  "reference_counts_layer": "counts",
  "max_iterations": 3,
  "movement_tolerance": 0.005,
  "margin": 0.5,
  "spatial_level": "l1"
}
```

This is a schema illustration, not a runnable skin annotation configuration: replace the empty dictionaries with curated marker sets and an L2-to-L1 mapping. Every configured L1 needs markers. The software intentionally rejects incomplete ontologies. Use `hierarchy-demo` for a fully runnable JSON example.

Add `"batch_key": "donor"` to request Harmony within each pool. Supply real donor/batch metadata; do not invent batch labels or assume sample integration is always scientifically appropriate. Harmony is skipped with an explicit status when a pool contains one batch. Missing dependencies or metadata fail rather than silently claiming integration ran.

Inputs require unique gene/cell IDs, integer raw counts and compatible gene identifiers/species. The run-level `species` declaration does not infer or verify the atlas species: verify compatibility before transfer. The default shared-gene minimum is 30. Atlas subtype-parent inconsistencies raise an error. Foreign atlas classes remain explicit evidence, not automatically mapped to known types.

## How decisions work

Module scores are marker means minus deterministic expression-matched background means. Each marker uses up to five closest-mean background genes; marker genes and configured exclusions do not act as controls within that scoring level. Scores are standardized once over the full post-QC query. This choice keeps the scoring baseline fixed across iterations but makes it query-composition-dependent. A zero-variance or poorly covered module is unavailable.

At least two markers must be detected in 20% of cluster cells by default. The best standardized cluster-mean lineage score must be positive and exceed the runner-up by at least 0.5; ties never pass, even if the configured margin is zero. At least two measurable alternatives are needed for a margin-based call. Subtype margins are evaluated among siblings.

A different L1 can be assigned only if markers pass and atlas majority purity is at least 0.7 with median supporting vote at least 0.6. Strong contradictions are Ambiguous. Initial L1 is retained for unresolved bookkeeping. L2 remains Unknown on conflicts or inadequate evidence; a moved cell cannot retain its old-compartment subtype. Only audited Keep calls drive enrichment; other calls become Uncertain.

All pools are built from the round-start assignment. Changes apply together at round end. Stop reasons distinguish exact stability, per-cluster movement tolerance, cycle and maximum iterations. Inspect total pool movement as well: many small moves may satisfy the per-cluster tolerance. Unresolved cells may remain even at numerical stability.

## Interpret results responsibly

The first public release validates synthetic mechanics only for this new method. It does not establish novelty, accuracy on ISD, superiority over existing methods or suitability for clinical use. Evaluate real data with held-out donors, known marker coverage, annotation disagreement, rare populations, depth dependence and batch effects before drawing biological conclusions.
