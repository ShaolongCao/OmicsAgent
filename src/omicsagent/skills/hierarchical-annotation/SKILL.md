---
name: hierarchical-annotation
description: Refine coarse L1 and subtype L2 annotations through bounded raw-count compartment reclustering and external atlas evidence in OmicsAgent for scRNA-seq or spatial data.
---

Use a fixed hierarchy JSON via RunConfig.hierarchy. Require initial L1 labels, disjoint compartment mapping, lineage/subtype marker sets and subtype parent mapping. A labeled reference is independent supporting evidence. Use the normal pipeline; the agent calls refine_hierarchy after annotate and before spatial analysis/reporting. The synthetic hierarchy-demo is a mechanics example, not ISD validation.

Each round snapshots all pools, reconstructs expression/PCA/optional Harmony/Leiden from raw counts, audits every L1 lineage and eligible L2 subtype, then applies changes simultaneously. Exclude specified ISG/overlap genes from lineage scoring. Scores use fixed-query standardized marker-minus-expression-matched-control means, with raw marker support; unavailable scores remain missing. Atlas PCA/KNN votes are not calibrated confidence and cannot override unsupported markers.

Reassign only when strong markers and adequate atlas purity/vote agree on a different lineage. Conflicts or insufficient margins remain Ambiguous. Missing evidence is Unknown. A moved cell loses its old L2 until its destination pool is refined. Confident same-lineage atlas subtype conflicts block L2 calls. Tiny and constant pools retain unresolved cells.

Default stop: no moving cluster exceeds 0.5 percent of its source pool, or at three rounds. Inspect stop_reason, total movement and unresolved counts. Cycle and iteration cap are not convergence; numerical stability does not establish biological correctness. Compare original L1, current L1 and L2; inspect per-round membership, decisions, score coverage, dotplots and atlas predictions.

Imaging graphs remain sample-separated; choose audited L1 or L2 labels for enrichment. Visium requires domain mode and expresses spot domains, never single-cell identity claims. Hierarchy does not consume fractional counts_split or imply a prior SPLIT correction was recomputed. Never treat agreement between methods as independent biological proof. Real atlas validation requires tissue-compatible markers and held-out donors.
