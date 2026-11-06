# Evidence rubric and loop control implementation plan

**Goal:** Implement the user's approved evaluation + runtime control approach and push verified changes to GitHub.

**Architecture:** A deterministic, versioned rubric consumes recorded cluster decisions. Runtime control uses cell-label transitions, unresolved evidence, repeated states and a fixed iteration budget. Reports display measured scores and unassessed dimensions without claiming calibrated confidence or biological validation.

**Tech stack:** Existing Python, pandas, offline HTML/JavaScript; no new dependencies or LLM calls.

**Approved design:** Eight dimensions from the conversation: applicability, support, counterevidence, granularity, selective resolution, next check, progress, audit. Scores 0/1/2 or null; no aggregate correctness score. Invalid assignments are hard failures. Accept at supported granularity, continue only for changed assignments requiring re-pooling, request scientific review for unresolved conflicts, stop unresolved when no executable progress remains. Existing numerical stop reasons remain visible separately. Historical runs may be assessed retrospectively but must never be represented as rubric-controlled runs.

## Tasks

- [x] Add failing policy tests: supported resolution, coarse-only resolution, conflict deferral, unsafe assignment, absent evidence, reassignment below movement tolerance, cycle, budget and all-Unknown controls.
- [x] Implement rubric.py, integrate per-round records and runtime control in hierarchy.py, retain original labels and old numerical stop diagnostics.
- [x] Display rubric dimensions and loop trajectory in the report; distinguish historical retrospective evaluation. Add report and workflow regression tests.
- [x] Document scoring anchors/limits, run full tests, lint/types, JavaScript report tests and public export checks; obtain independent review and fix findings.
- [ ] Push the verified change and update the existing PR, without merging.

## Review focus

Missing measurements are not passes; no truth labels are inferred from atlas votes. A stable unresolved run must not be called accepted. A small reassignment still needs validation in its new pool. Cycles and budget exhaustion cannot silently accept changed labels. Skipped pools count toward unresolved cells. Historical files remain immutable.
