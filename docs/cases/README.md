# Three measured annotation decisions

**Current conflict call: retain Myeloid; withhold subtype.**

In cluster `r1:immune:5`, the marker subtype winner is Neutrophil, while reference transfer gives Monocyte. Both support the coarse Myeloid identity. The available evidence does not justify an automatic subtype assignment.

[Open the report source](evidence/report.html) or download the repository ZIP, extract it, and open `docs/cases/evidence/report.html` locally. The report embeds its evidence and works offline; adjacent tables preserve inspectable source links.

| Case | Affected cells | Evidence cluster cells | L1 | L2 |
|---|---:|---:|---|---|
| [agreement](#agreement) | 121 | 121 | B cell | Plasma |
| [conflict](#conflict) | 238 | 238 | Myeloid | Unknown |
| [insufficient](#insufficient) | 142 | 144 | Myeloid | Unknown |

These are selected examples from one fixed run, not a biological benchmark. Final-round eligible decision rows: agreement=7, conflict=5, insufficient=17. Two additional marker-resolved rows lack the required dual support and are excluded from these case categories. Rows sharing a cluster are not independent observations.

The query contains 2,161 cells after QC (2,286 prepared); reference contains 1,198 cells. Reference donor: TSP2. Query donors: TSP11 and TSP14. A stable-ID cap of 150 cells per donor/published fine label was fixed before inspecting outcomes. Unmapped classes were excluded; no labels or counts were artificially corrupted.

## conflict

Cluster `r1:immune:5`: **Myeloid / Unknown**. L1 supported: True.

- L1 marker winner Myeloid; margin 1.164; detection gate passed.
- Reference L1 Myeloid: purity 1.000, median vote 1.000.
- L2 marker winner Neutrophil passes the marker gate.
- Counterevidence: Reference L2 Monocyte conflicts with the marker winner; purity 1.000, median vote 1.000.

No subtype assigned.

**Next check:** Across 2 donors with at least 5 evaluable cells, marker/reference disagreement ranges 99.6% to 100.0%. Disagreement is observed in both available donors. Compare the competing subtype markers in this cluster and its mapped reference cells before accepting either subtype. This can expose panel specificity or reference annotation mismatch; the donor counts are descriptive, not a replication test.

| Donor | Evaluable cells | Disagreeing paired winners | Fraction |
|---|---:|---:|---:|
| TSP11 | 6 | 6 | 1.000 |
| TSP14 | 232 | 231 | 0.996 |

Per-cell winners compare competing measured marker modules with reference labels; they are not calibrated correctness calls. Unequal and small donor counts limit inference.

<details><summary>Inspect marker measurements</summary>

| Panel | Gene | Mean log expression | Detection fraction |
|---|---|---:|---:|
| Myeloid | LST1 | 1.801 | 0.899 |
| Myeloid | LYZ | 1.806 | 0.899 |
| Myeloid | TYROBP | 2.434 | 0.950 |
| Monocyte | CTSS | 2.881 | 0.987 |
| Monocyte | FCN1 | 1.801 | 0.899 |
| Monocyte | LYZ | 1.806 | 0.899 |
| Neutrophil | CSF3R | 3.183 | 1.000 |
| Neutrophil | FCGR3B | 2.953 | 1.000 |
| Neutrophil | S100A8 | 6.216 | 1.000 |

</details>

[Decisions](evidence/hierarchy/decisions.csv) · [Marker measurements](evidence/hierarchy/dotplot_data.csv) · [Parameters](evidence/hierarchy/config.json) · [Execution trace](evidence/trace.jsonl) · [Manifest](evidence/manifest.json)

## agreement

Cluster `r1:immune:18`: **B cell / Plasma**. L1 supported: True.

- L1 marker winner B cell; margin 0.536; detection gate passed.
- Reference L1 B cell: purity 0.975, median vote 1.000.
- L2 Plasma: marker_supported.
- Counterevidence: Per-cell marker/reference winners differ in 10/118 evaluable cells; cluster agreement does not imply uniform cells.

Provisional annotation; independent biological validation remains open.

| Donor | Evaluable cells | Disagreeing paired winners | Fraction |
|---|---:|---:|---:|
| TSP11 | 83 | 4 | 0.048 |
| TSP14 | 35 | 6 | 0.171 |

Per-cell winners compare competing measured marker modules with reference labels; they are not calibrated correctness calls. Unequal and small donor counts limit inference.

<details><summary>Inspect marker measurements</summary>

| Panel | Gene | Mean log expression | Detection fraction |
|---|---|---:|---:|
| B cell | CD79A | 0.689 | 0.719 |
| B cell | CD79B | 0.388 | 0.744 |
| B cell | MS4A1 | 0.022 | 0.025 |
| Plasma | JCHAIN | 1.880 | 0.669 |
| Plasma | MZB1 | 2.832 | 1.000 |
| Plasma | SDC1 | 1.031 | 0.893 |

</details>

[Decisions](evidence/hierarchy/decisions.csv) · [Marker measurements](evidence/hierarchy/dotplot_data.csv) · [Parameters](evidence/hierarchy/config.json) · [Execution trace](evidence/trace.jsonl) · [Manifest](evidence/manifest.json)

## insufficient

Cluster `r1:immune:2`: **Myeloid / Unknown**. L1 supported: True.

- L1 marker winner Myeloid; margin 1.188; detection gate passed.
- Reference L1 Myeloid: purity 0.986, median vote 1.000.
- Counterevidence: L2 unresolved: insufficient_marker_support.

No subtype assigned. Evidence is aggregated over 144 cluster cells; this review affects 142 cells with the recorded source label.

**Next check:** L2 marker margin 0.213 is below the configured minimum 0.500. Inspect per-cell co-expression of the competing subtype markers to distinguish a mixed cluster from overlapping panels before accepting the reference subtype.

| Donor | Evaluable cells | Disagreeing paired winners | Fraction |
|---|---:|---:|---:|
| TSP11 | 6 | 6 | 1.000 |
| TSP14 | 136 | 83 | 0.610 |

Per-cell winners compare competing measured marker modules with reference labels; they are not calibrated correctness calls. Unequal and small donor counts limit inference.

<details><summary>Inspect marker measurements</summary>

| Panel | Gene | Mean log expression | Detection fraction |
|---|---|---:|---:|
| Myeloid | LST1 | 1.363 | 0.806 |
| Myeloid | LYZ | 2.882 | 0.972 |
| Myeloid | TYROBP | 2.007 | 0.896 |
| Monocyte | CTSS | 2.676 | 0.986 |
| Monocyte | FCN1 | 2.706 | 0.986 |
| Monocyte | LYZ | 2.882 | 0.972 |
| Neutrophil | CSF3R | 2.025 | 0.944 |
| Neutrophil | FCGR3B | 1.248 | 0.750 |
| Neutrophil | S100A8 | 6.633 | 1.000 |

</details>

[Decisions](evidence/hierarchy/decisions.csv) · [Marker measurements](evidence/hierarchy/dotplot_data.csv) · [Parameters](evidence/hierarchy/config.json) · [Execution trace](evidence/trace.jsonl) · [Manifest](evidence/manifest.json)

## Provenance and limits

[Tabula Sapiens Figshare release v5](https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219/5), CC BY 4.0. Cite The Tabula Sapiens Consortium, Science 376, eabl4896 (2022). [Full data policy](DATA_POLICY.md). Raw integer counts were verified; the normalized X slot was not used. Published labels are reference metadata, not independent ground truth. TCL1A is absent from this source; no replacement was introduced to force an outcome.

Numerical CSVs in the evidence bundle are unchanged from the completed run. JSON/JSONL absolute path fields are replaced with basenames for portability; original artifact hashes and the source run hash remain in the manifest. Count matrices and original figures are retrieved/generated separately with the documented preparation/run commands.

[Design rationale and reproduction](../SCIENTIST_REVIEW.md) · [Validation record](../validation-review.md)
