---
name: split-decontamination
description: Apply and audit OmicsAgent's Python implementation of standard RCTD-free SPLIT after scANVI or TACCO on imaging single-cell spatial counts.
---

This adapter is a Python port of the standard SPLIT v0.3.0 RCTD-free purification formula, not the official R package. Native R numerical parity has not been validated. It excludes aggressive residual removal and RCTD fitting.

Set `split_enabled: true`, `platform: "imaging"` and `annotation_method: "scanvi"` or `"tacco"`; supply a compatible labeled raw-count reference. The pipeline aligns weights, cells, genes and types. Reference profiles are per-type means after each reference cell is library-normalized, matching the supplied export skill.

For cell c, gene g, primary type p=argmax(weight), retain raw[c,g] times (weight[c,p]*profile[p,g]) / sum_t(weight[c,t]*profile[t,g]). Zero modeled support yields zero corrected signal. Genes absent from the reference are unchanged and flagged. Contributions use every type, not just the top two. Weights/posteriors are modeling assumptions; their uncertainty propagates into correction.

Inspect `split_qc.csv`, `summary.json` and `var.oa_split_reference_covered`. Results are fractional `layers['counts_split']`. Original counts, expression analysis and marker gate remain unchanged. Do not round or feed corrected values into count-likelihood scANVI automatically. Compare removed signal by cell type/depth and known markers before using corrections; removal is not proof that the signal was ambient.

Source adaptation: supplied xenium-split-decontamination skill. Algorithm: https://github.com/bdsc-tds/SPLIT/blob/v0.3.0/R/method_agnostic_purification.R . No claim that original private examples were reproduced.
