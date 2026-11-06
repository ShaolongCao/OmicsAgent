# A measured annotation review

OmicsAgent should resolve what its evidence supports and make unresolved questions easy for a scientist to inspect. This prototype focuses on one decision: whether a cluster's label should be accepted when markers and reference transfer disagree.

## Three-minute walkthrough

1. Open the [case entry](cases/README.md), download its report bundle, and open `report.html` locally. No API key or server is needed.
2. Start with **Evidence conflict**. Read the separate L1 and L2 calls, supporting evidence, counterevidence, and the observed reason for the next check.
3. Expand marker measurements, donor summaries, parameters and source artifacts. A reference vote is not a correctness probability. Retaining an original label does not validate it.
4. Compare **Agreement** and **Insufficient evidence**. Read the eligible case counts: these are selected examples, not representative performance estimates.
5. Enter a review reason, accept the proposal, then undo it. Export the JSON history to preserve both actions. Export CSV to inspect reviewed labels. Original analysis files remain unchanged; a reload starts a new draft unless a saved history is imported.

## Reproduce

```bash
omicsagent prepare-review --out data/scientist-review
omicsagent run --config data/scientist-review/run.json
omicsagent review --out data/scientist-review/run --donor-key donor
omicsagent evaluate-review
```

Preparation downloads the pinned original Tabula Sapiens bone marrow processed-count file (about 403 MB compressed, 1 GB extracted), validates checksums, uses the integer raw slot, restricts to 10X, preserves published labels, and separates TSP2 reference from TSP11/TSP14 query. Install the development and hierarchy extras as in the main README. The stable cell-ID cap is fixed before viewing outcomes. Real runs need sufficient local disk and memory; generated count matrices stay outside the public source release.

## Why this problem

The useful output is a defensible decision and a way to continue working. Coarse identity may be supported while subtype evidence conflicts. Compressing both into one confidence label hides that distinction. Showing all evidence at once also obscures the most consequential question, so the first view is deliberately short and details remain inspectable.

## Tradeoffs

- Reuse deterministic marker/reference gates and recorded tables. This makes the explanation traceable, but panel quality and reference compatibility still constrain it.
- Use an offline report. It is easy to share and replay, but review history must be exported; this is not a collaborative database.
- Select cases by declared rules after a fixed run. This makes the portfolio legible, but does not estimate population accuracy or robustness.
- Diagnose donor patterns descriptively. Disagreement fractions can guide inspection; they cannot establish a causal donor effect.
- Preserve originals and export reviewed labels separately. The scientist can reverse a choice without erasing the evidence or earlier work.

## What would change this design

Repeated false resolutions on held-out donors would require revisiting marker specificity and decision gates before adding more automation. A high rate of missing-marker cases would favor better panel coverage diagnostics. Scientists repeatedly needing joint review would justify persistent multi-user review state. If the suggested checks do not help scientists decide what to do next, the trigger rules need revision, even when they are technically reproducible.

Behavior controls, selected real examples and biological validation are reported separately. No clinical use, calibrated probabilities, superiority benchmark or real LLM-call validation is claimed.
