---
name: tacco-annotation
description: Run TACCO optimal-transport reference annotation in OmicsAgent for single-cell expression or imaging spatial data and inspect composition weights.
---

Set `annotation_method: "tacco"`, `reference`, `reference_label_key` and raw-count layers. Install the `tacco` extra. Verify species, gene identifiers, tissue and reference class coverage. The adapter validates counts, intersects genes, uses TACCO OT with platform_iterations=0 and bisections=None, and preserves full cell-by-type weights. These settings are explicit baselines, not optimized parameters for every platform.

Run the normal pipeline or its `annotate` tool. Inspect `annotation_weights.csv.gz`, maximum weight, margin, entropy and raw-count marker audit. TACCO weights describe assignments/composition; do not report them as calibrated confidence. Apple Silicon may use the slower non-MKL fallback. Do not force x86 MKL installation.

For imaging, optionally set `split_enabled: true` after checking candidate assignments. Visium spots remain mixed domains in spatial statistics even if reference candidates are present; this implementation does not claim validated Visium cell deconvolution. Preserve uncertain labels and review count-depth bias.

Source adaptation: supplied xenium-tacco-annotation skill. Upstream: https://github.com/simonwm/tacco . Original dataset-specific results are not general validation.
