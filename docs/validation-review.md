# Scientist review validation record

Date: 2026-10-28. Local platform: Windows, Python 3.12.14, Scanpy 1.11.5, AnnData 0.12.19. This record covers a local run, not public CI results.

## Real data

One pinned Tabula Sapiens bone marrow file, Figshare release v5, CC BY 4.0, DOI 10.6084/m9.figshare.14267219.v5. Publisher archive MD5 matched `1b665025142b0c74f7b63d7310b379d4`; extracted SHA256 matched `4e08b45d5668c7cc741ccddafd2626f65303ce764be4f3ebc1b7778d60c30c64`. All 30,482,041 stored raw values were checked: no fractional or negative counts. Raw X and the published raw_counts layer were identical; normalized X was excluded.

Preparation used 10X only, TSP2 reference, TSP11/TSP14 query, cap 150 per donor/published fine label in stable ID order, and the declared ontology/panels. Reference: 1,198 cells. Prepared query: 2,286 cells; retained after QC: 2,161. No biological thresholds or labels were adjusted to force case categories. A missing TCL1A marker remains recorded as missing.

The completed run stopped after one round with the existing tolerance criterion. This is numerical stopping, not biological convergence. Final status counts: Keep 1,873; Ambiguous 270; Reassign 18. Original labels remain saved. Source run hashes, actual thresholds and complete decision tables accompany the report.

Eligible final-round decision rows: agreement 7, conflict 5, insufficient 17. Two marker-resolved rows without dual support were excluded from these categories. Rows may share clusters and are not independent observations. The largest eligible row per category was selected with stable cluster/source-label tie breaking:

| Case | Affected cells | Evidence cluster cells | Result |
|---|---:|---:|---|
| Agreement | 121 | 121 | B cell / Plasma |
| Conflict | 238 | 238 | Myeloid; L2 withheld, Neutrophil markers versus Monocyte transfer |
| Insufficient | 142 | 144 | Myeloid; L2 withheld, marker margin 0.213 below 0.500 |

Conflict paired-winner disagreement: TSP11 6/6 evaluable cells; TSP14 231/232. The next check compares competing markers in query/reference because disagreement was observed in both donors. These unequal counts are descriptive, not a replication or donor-effect test.

## Mechanism and engineering checks

- Full pytest: **73 passed, 2 skipped**. Optional TACCO absent; Windows symlink creation lacks privilege. Environment-file rejection is still asserted; the symlink branch executes on platforms with permission.
- Ruff: passed. Mypy: passed for 20 source files.
- Public source export and release scanner: passed (74 selected files); UTF-8 decoding is explicit across Windows locales.
- Five explicit policy controls passed, including agreement resolution, coarse conflict deferral, atlas-only rejection, supported reassignment and weak-marker deferral. The all-Unknown control fails; its resolved coverage is zero. This is policy evaluation, not biological accuracy.
- Regression tests cover L1/L2 separation, original-label summaries after earlier reassignment, final skipped rounds avoiding stale judgments, donor paired-measurement denominators and ties, source-scoped portable bundles, old schemas, empty tables and HTML escaping.
- Node engine and Edge browser checks cover case switching, evidence expansion, marker plot rendering, accept/undo, JSON/CSV downloads, valid history import, atomic rejection of malformed/cross-run histories, original-first exports, script/CSV injection handling and mobile layout.
- The portable evidence bundle preserves numerical CSV bytes. Absolute JSON/JSONL path fields are replaced with basenames; original artifact hashes and parent run hash remain in the bundle manifest.

Two important independent-review findings were reproduced and fixed with regression tests: immutable original-label summaries and stale supported calls after a wholly skipped final round. No deferred minor findings were reported.

## Operational limits

The first large-file save exhausted local disk; the completed run was repeated with gzip H5AD output and identical analysis settings. Preparation reads selected raw rows instead of materializing unrelated corrected layers. These are storage/memory changes, not changes to scientific gates. Generated count matrices remain outside the public source release.

No paid model/API call, clinical validation, real ISD skin validation, held-out accuracy benchmark, calibrated correctness probability or superiority claim is made. Published annotations supply reference metadata, not independent truth. The cases and suggested checks require scientific review.
