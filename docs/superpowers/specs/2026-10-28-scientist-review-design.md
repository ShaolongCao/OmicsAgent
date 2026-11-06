# OmicsAgent: evidence-aware scientist review

Status: approved by the user on 2026-10-28. Scope: real public data, three cases, an interactive report, and behavioral evaluation.

## Intent

Help a scientist understand why a provisional annotation is defensible, what remains unresolved, and which observation would justify changing it. The portfolio entry should take about three minutes to explore from the README. Every case statement must be grounded in saved numerical artifacts.

## Approaches and decisions

The user selected a complete small prototype over a synthetic-only first iteration. Build on the existing Python hierarchy workflow and its offline artifacts. Add self-contained HTML with embedded structured case data and small JavaScript interactions. This implementation choice is the builder's recommendation: it preserves offline replay and avoids introducing a hosted service, account system, or additional model API calls. A separate hosted application would add deployment and persistence work without improving this initial evidence review task.

## Goals and boundaries

- Show three measured case types: agreement, marker/reference conflict, and insufficient evidence, with the conflict case selected initially.
- Separate coarse L1 identity, subtype L2 identity, and retained original labels. Retaining an original label is never presented as evidence that it is validated.
- Show both support and counterevidence, threshold outcomes, unresolved questions, and at most one prioritized next check.
- Preserve original labels; expose proposed changes, affected cells, before/after differences, accept and undo, and an exportable review history.
- Evaluate appropriate resolution as well as abstention. An all-Unknown policy must fail evaluation.
- Keep the existing numerical workflow and raw counts intact. Do not add a new annotation model, clinical claims, or a general chat interface.

## Data and scientific honesty

Use one pinned public processed-count dataset with explicit reusable terms and donor metadata. The initial candidate is a Tabula Sapiens immune-rich tissue (prefer blood); the exact release, tissue, file URL, reuse terms, checksum, metadata schema and citation must be verified before use. Use a smaller single tissue if needed for local memory limits. Do not silently replace licensed raw-count data with a minified trained-model object.

Create disjoint reference and query donor partitions when the metadata supports them; document donors, cell IDs, filtering, label ontology and random seed. Published annotations are reference metadata, not independent biological ground truth. Keep numerical data outside the published source tree, following the existing release restrictions; publish reproducible fetching/preparation code and compact, privacy-appropriate derived case artifacts with attribution.

Run the analysis before selecting cases. Record prespecified selection predicates and include eligible and selected cluster counts, stable case IDs, run hashes and selection reasons. Do not alter marker panels or thresholds solely to force the three desired outcomes. If the first candidate lacks a case type, inspect another tissue of the same dataset under the same documented policy. If three genuine case types still cannot be found, report the unmet criterion and retain a clearly labeled synthetic mechanism example separately; do not claim the requested real-case deliverable is complete.

## Decision evidence

Extend decisions with L1 and L2 marker winners, measured score margins, marker availability and detection support, reference winner/purity/vote, outcome reason codes and final retained labels. Include explicit reasons for tiny/degenerate pools and unavailable reference evidence. Preserve iteration and cluster identity in every link. Distinguish missing measurements from negative evidence, and transfer vote fractions from calibrated correctness probabilities.

Keep existing conservative reassignment semantics. A supported coarse identity can remain resolved while subtype evidence conflicts or is insufficient. A coarse-level conflict remains unresolved even if the old label is retained. Record proposed changes before they are applied numerically so the report can distinguish the pipeline proposal from the scientist's review choice.

## Triggered next check

Use deterministic, artifact-backed rules, not generated generic advice. For a subtype conflict with donor metadata and adequate representation across donors, calculate donor-stratified marker/reference summaries and recommend examining the observed concentration or reproducibility of the disagreement. Show denominators and the numerical trigger. If donor metadata or representation is inadequate, state that limitation and choose a relevant measured alternative, such as marker coverage when markers are unavailable. Do not call sample metadata donor metadata without verifying its meaning. A resolved agreement case need not receive a follow-up recommendation. These are diagnostic suggestions, not causal conclusions.

## Report and review behavior

The first view contains the current L1/L2 judgment, the principal supporting and opposing evidence, the unresolved question, and the next check when triggered. Three case selectors lead to expandable evidence tables, marker plots and parameters. Include links to decisions, round-cell membership, configuration, manifest and execution trace; provide a small self-contained evidence subset so the shared portfolio report remains inspectable offline.

Accept/undo changes the report's review state, never the original analysis outputs. Show affected cell counts and original/proposed/reviewed labels. A later accepted action can be undone with a new history event. State that browser review state is a local draft; allow JSON export/import scoped to the dataset and run hash, and CSV export of reviewed labels. Validate imported reviews against the run and allowed proposals. Exported records include event sequence, action, timestamp, proposal ID, reason/evidence references and prior/new state. If a rerun is desired, reviewed labels are a separate explicit input, not a mutation of the original h5ad. Escape data in HTML and avoid external scripts or network-dependent report assets.

## Integration and reproducibility

Add report generation to the existing pipeline and preserve Markdown/notebook outputs. Add CLI entry points for reproducible data preparation, case report generation and behavior evaluation as needed. Include generated evidence and report files in manifest hashing without a self-referential hash. Ensure cached replay reuses recorded decisions and regenerates equivalent evidence instead of refitting. Update release allowlists deliberately for the new documentation and static demo; do not relax private-data or credential screening broadly.

README links to a committed, directly viewable case entry and explains how to open the HTML locally; use a hosted link only if an existing authorized hosting route is available. Provide a three-minute walkthrough and a one-page design rationale with tradeoffs and observations that would change the design. Keep repository source and documentation in English, consistent with the existing release policy.

## Evaluation and acceptance

Use meaningful behavior tests plus real-data descriptive results, with these reported separately. Cover: supported agreement resolves; supported L1 plus conflicting L2 retains L1 and defers L2; strong reference cannot override missing/failed markers; inadequate evidence abstains; missing metadata cannot trigger a donor claim; supported reassignment preserves originals; accept/undo/export/import retain run-scoped traceability; and all-Unknown fails the resolution requirement. Evaluate coverage and inappropriate assignment alongside abstention, never one scalar score that rewards blanket abstention.

Verify report statements against source tables, local evidence links, thresholds and cell counts. Exercise browser case switching, evidence expansion, accept/undo and exports. Run relevant existing hierarchy/replay tests, full feasible project checks, lint/type checks and release scanning. Record actual results and skipped optional dependencies; do not claim a biological benchmark or a real LLM/API validation.

## Builder discretion and delivery

The builder chooses the exact tissue, deterministic sampling caps, visual styling, implementation modules and test organization after inspecting data and runtime constraints. Record these decisions and their reasons. No API key or paid model call is required. Deliver the code changes, reproducible commands, actual three-case report if supported, evaluation record and design rationale. Prepare a reviewable branch/patch; do not merge into the default branch as part of this scope.
