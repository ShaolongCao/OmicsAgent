# Supplied skills: integration map

Methods adapted from the author's previously developed single-cell and spatial analysis skills. Original environment-specific templates and private case studies are not distributed or claimed as independently validated.

| Skill | Current disposition |
|---|---|
| scRNAseq_cell_type_annotation_basic | Partial adaptation: shared-gene reference-fitted PCA/KNN; counts validation; marker panels. Query scaling corrected to use reference mean/SD. Harmony, multi-resolution Louvain and rare-type subclustering not ported. |
| xenium-annotation-qc-gate | Partial adaptation of A/B: raw-count own-label support, all-label panel coverage, preserved candidates, 1/2-marker sensitivity, count-decile diagnostic. C/D reference sparsity/centroid and E control probes/segmentation purity not implemented. |
| spatial-niche-enrichment | Python port of core fixed-center radius neighborhood/per-sample label permutations. Graph is sample-isolated, preserves incidence multiplicity, excludes self by identity rather than dropping all zero-distance observations, uses plus-one permutation p values. Stouffer, paired contrasts and DESeq2 not implemented. |
| xenium-tacco-annotation | Working optional TACCO OT adapter; full weights and marker audit. |
| xenium-scanvi-annotation | Working optional scVI → scANVI adapter; CPU training, weights and saved model. |
| xenium-split-decontamination | Python port of standard RCTD-free core; optional imaging correction layer. No native R parity or aggressive residual removal. |
| spata2-isg-gradient | Catalogued; no SPATA2 gradient runtime. |
| monocle3-h5ad-trajectory | Catalogued; no trajectory inference. |
| rnaseq-quickomics | Catalogued; bulk RNA-seq outside v0.1. |
| eln-writing | Catalogued; Markdown reports are generated, but no claim of ELN rubric compliance. |

## Scientific changes from the reference material

- A negative association between label frequency and count depth is a diagnostic risk, not proof of an artifact: real cell populations may differ in RNA content. The original categorical wording is softened in code and reports.
- A class with fewer available panel markers than the requested threshold is ungateable, not automatically failed. Ungateable candidates are preserved for audit but excluded from supported imaging labels used for enrichment.
- KNN and another embedding-based classifier do not become independent confirmation merely by agreeing. v0.1 does not implement a consensus vote or call KNN vote fraction a calibrated probability.
- The numeric results in private worked examples have not been independently reproduced here and are not included as validation claims.
- A zero-variance permutation null returns undefined z/p, rather than substituting an arbitrary standard deviation of one.
- A bounded graph-size check occurs before materializing radius pairs; sparse full-gene expression is preserved outside explicitly bounded reference PCA.

The original scripts are notebook-like templates with hardcoded environments and datasets. Their core algorithms were refactored into functions rather than automatically executing templates in the user's working environment. The new implementation is not a drop-in reproduction of all original pipelines.

The three executable method skills live in `src/omicsagent/skills/`. The Claude agent loads only the skills selected by the config. They are packaged with the application, not installed globally. See [worked example](HIERARCHY_TUTORIAL.md).

## Iterative hierarchy

The new `hierarchical-annotation` package skill drives the optional fixed-config refinement stage. It documents synchronous raw-count reclustering, marker/atlas evidence, conservative reassignment, stopping reasons and spatial domain semantics. Real ISD validation is pending.
