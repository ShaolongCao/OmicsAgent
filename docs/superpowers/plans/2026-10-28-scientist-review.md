# Scientist Review Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Build an evidence-grounded three-case scientist review prototype.

**Architecture:** Extend existing hierarchy artifacts, derive structured case evidence from saved tables, and render a standalone HTML review report. Prepare one licensed real dataset reproducibly and evaluate policy behavior separately from biological accuracy.

**Tech Stack:** Python 3.12, Scanpy/AnnData, pandas, existing pytest/ruff/mypy, standalone JavaScript/HTML.

**Spec:** ../specs/2026-10-28-scientist-review-design.md

## Global Constraints

- No paid API calls; original counts and labels remain immutable.
- Real examples require measured evidence; synthetic fixtures remain labeled.
- Donor claims require verified donor metadata; reference/query donor IDs must be disjoint.
- Source/docs remain English; existing release screening remains enforced.
- Review actions create separate draft exports, never mutate analysis artifacts.

## Review Focus

- Missing/degenerate evidence must remain represented, with no invented confidence.
- Donor denominators and insufficient representation must not imply reproducibility.
- Imported events from another run or malformed histories must be rejected atomically.
- Embedded labels containing markup must render as text and not execute scripts.
- Cached replay and older decision tables must remain readable.

### Task 1: Record decision evidence

**Files:** modify src/omicsagent/hierarchy.py; test tests/test_hierarchy_workflow.py.
**Interfaces:** produces decisions.csv with l2 marker strength/margin, selected subtype, and explicit l1/l2 reasons; existing columns retained.
- [x] Add a failing integration assertion for supported L1/conflicting L2 and explicit reason fields.
- [x] Run pytest on the test and observe the missing-field failure.
- [x] Implement evidence fields and reason codes without weakening existing gates.
- [x] Run hierarchy tests; expected all applicable tests pass.

### Task 2: Build structured cases, behavioral evaluation and review report

**Files:** create src/omicsagent/review.py, src/omicsagent/review_html.py, tests/test_review.py; modify report.py and cli.py.
**Interfaces:** build_review(out: Path, donor_key: str | None = None) -> dict; write_review(out: Path, donor_key: str | None = None) -> dict. Consumes existing run artifacts; produces review.json/report.html with source links and a deterministic run hash. evaluate_behavior() -> dict produces nontrivial policy evaluation.
- [x] Test agreement/conflict/insufficient classification, missing metadata, all-Unknown rejection and HTML injection; observe failure.
- [x] Implement selection predicates, source hashes, support/counterevidence, triggered next-check and standalone report.
- [x] Add run-scoped accept/undo and validated JSON import/export plus reviewed-cell CSV export.
- [x] Exercise browser actions and malformed import handling; run pytest tests/test_review.py.

### Task 3: Prepare real data and publish measured cases

**Files:** create examples/prepare_review.py, docs/cases/*, docs/SCIENTIST_REVIEW.md; modify README.md and scripts/export_release.py.
**Interfaces:** reproducible preparation CLI with pinned dataset provenance; real run artifacts consumed by Task 2.
- [x] Verify reuse terms, count matrix and donor/schema before selecting a pinned tissue.
- [x] Add tests for disjoint donor splitting and refusal of normalized values as counts; observe failure, then implement preparation.
- [x] Execute fixed-policy real run; select measured cases, record eligible counts and provenance; never force outcomes.
- [x] Copy compact self-contained case artifacts and report into docs/cases; update explicit release allowlist and walkthrough/design rationale.

### Task 4: Integrate and verify

**Files:** tests and documentation as required; update this checklist and validation record.
- [x] Run full feasible pytest suite, ruff, mypy and release scanner; record failures/skips precisely.
- [x] Obtain independent whole-branch review; fix important findings with regression tests.
- [x] Export reviewable patch/source and user report; prepare branch without merging default branch.

## Execution decisions

The user explicitly requested implementation of the approved design; proceed inline without another scope/plan approval. This clone is already an isolated task-local checkout on scientist-review-prototype. Use that checkout rather than creating another worktree. Git author identity is currently unconfigured: preserve changes and provide a patch if a commit cannot be made without inventing user identity.
