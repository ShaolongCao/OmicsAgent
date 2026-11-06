# Cell annotation evidence rubric v1.0

This rubric evaluates whether a decision is defensible from recorded evidence and controls when the hierarchy loop continues. It is deterministic; an LLM does not grade its own output. It is not a calibrated probability, a biological accuracy benchmark, or a license to accept a label because its total score is high. There is deliberately no total score.

## Scoring anchors

Scores are 0 (criterion not met), 1 (partial), 2 (criterion met), or null (unassessed). Scores apply to the stated criterion only. Each output includes the measured fields or the reason assessment is unavailable.

| Dimension | Implemented anchor | Limits |
| --- | --- | --- |
| Applicability | null | Tissue/technology compatibility, out-of-reference status and reference quality need additional measurements. Metadata alone is not validation. |
| Support | 2 when marker gates support the proposed parent without confident reference opposition; 0 when they fail; null when absent | Gates include marker count, detection, positive module score and separation. Multiple markers need not be independent or specific. |
| Counterevidence | 2 for a detected conflict that is safely deferred; 1 when no measured conflict is detected; 0 for unsafe assignment | Absence of conflict does not rule out doublets, ambient RNA, batch effects or a shared error. |
| Granularity | 2 for supported L1/L2 with safe subtype deferral; 1 for unresolved identity; 0 for an unsafe assigned level | Retaining an original label never validates it. L1 can remain supported while L2 requires review. |
| Selective resolution | null per case | Evaluate coverage and inappropriate resolution across a balanced control set or independently adjudicated data. |
| Next check | 2 when an observed issue triggers a question with an explicit decision consequence; otherwise null | Only re-pooling changed assignments is executable automatically. Other checks require a scientist. This score does not prove a check was performed. |
| Progress | null per case; round record gives 1 for changed assignments, 0 for unchanged assignments | Movement creates a validation obligation, not biological evidence. There is no claim of new information from rerunning fixed markers/reference. |
| Audit | 1 when a cluster decision is recorded; otherwise null | Artifact hashes and UI accept/undo are verified separately; the policy does not self-certify integrity. |

## Hard failures and actions

Unsupported parent assignment, reference-unsupported reassignment, subtype assignment without marker support, assigned subtype contrary to confident reference evidence, or subtype outside the configured parent are hard failures. Runtime aborts before publishing those labels and leaves a rubric diagnostic. No other score can compensate.

Four actions are recorded independently of the legacy numerical stop reason:

1. `accept`: all assessed subsets support the reported granularity, with no skipped pools. This is policy acceptance, not external validation.
2. `continue`: some cell identities changed and require validation in newly constructed pools. This overrides movement tolerance, but cannot override a repeated state or iteration budget.
3. `review`: a stable unresolved marker/reference conflict requires scientific judgment. Supported parent identities are retained; disputed subtypes remain Unknown.
4. `stop_unresolved`: no executable progress remains, a pool was skipped, evidence is insufficient, a state cycles, or the budget is exhausted. Counts of unresolved and supported cells are retained.

Stable output is not necessarily resolved output. `summary.converged` now means policy acceptance; `numerically_stable` and the legacy `stop_reason` retain numerical diagnostics. Read `rubric_action` and `rubric_reason` for scientific workflow status. A terminal review can include both supported and unresolved populations. Changed cells at the budget limit remain unvalidated, with their subtypes cleared by the existing hierarchy logic.

## Artifacts and report

Each hierarchy run writes `hierarchy/rubric.json`, with rubric version, every round's decisions, evidence paths, skipped/moved cell counts, remaining budget and control action. `hierarchy/iterations.json` and `summary.json` also expose the action. The standard manifest hashes the rubric artifact. The offline report has per-case rubric details and a loop trajectory, including when no cases are eligible.

`omicsagent review --out EXISTING_RUN` displays retrospective case scoring if the run predates the rubric. It does not rewrite its decisions or claim that the historical loop followed the new policy. Changes to report evidence can change its run hash; old review histories are intentionally not imported across different evidence snapshots.

## Evaluation

Run `omicsagent evaluate-rubric` for synthetic decision and loop controls, including unsupported subtype, blind reference acceptance, blanket Unknown, small reassignments, repeated states and budget exhaustion. Run `omicsagent evaluate-review` for the original five decision controls. Unit/integration tests additionally exercise conflict deferral, skipped pools, absence of evidence and the offline UI.

These tests validate implementation behavior, not biological accuracy. The three real Tabula Sapiens examples are illustrative, not a representative benchmark. A future biological benchmark should use independently adjudicated labels and held-out donors/references; report L1/L2 accuracy, macro-averaged performance, accuracy versus coverage, unknown precision/recall, false confident assignment, per-donor behavior and computational cost. Freeze thresholds before evaluation. Published atlas labels used for transfer cannot simultaneously serve as independent truth.

The real Tabula Sapiens rerun on 2026-11-06 exercised the runtime controller: 18 changed cells in round 1 overrode the old tolerance stop, 15 changed in round 2, and 20 changed in round 3 exhausted the configured budget. The final action is `stop_unresolved`, not convergence. At the final round 1,047 cells meet the per-decision acceptance policy and 1,114 remain unresolved at the requested granularity; this is not an accuracy estimate. See the [recorded validation summary](cases/rubric-validation.json). Re-clustering uses a round-specific seed, so small changes can include algorithmic variation: further iterations are not necessarily better and the budget must remain bounded. Existing one-round case examples retain their original measurements and are labeled retrospective.

Validation on 2026-11-06: 87 Python tests passed, 2 skipped (optional TACCO unavailable and Windows symlink privilege absent). Ruff and mypy passed. The actual offline runtime report was checked in Edge for three recorded rounds, budget exhaustion and unassessed dimensions; the review interaction suite also passed. Independent review identified a rejected-round evidence persistence gap; its regression was reproduced and fixed, and the reviewer confirmed the fix. Public export scanning passed for 79 files.

## Design choice

An evaluation-only layer would describe mistakes after a run but could not prevent premature termination. The user selected evaluation plus control: deterministic evidence gates govern actions while the report exposes the remaining uncertainty. New biological checks should only become executable loop actions after their inputs, interpretation rules, resource budget and failure behavior are tested. No autonomous external analysis or LLM judge was added.

Hierarchical rejection motivates preserving a supported parent rather than forcing a fine label or rejecting every level: [Uncertainty-aware single-cell annotation with a hierarchical reject option, Bioinformatics (2024)](https://academic.oup.com/bioinformatics/article/40/3/btae128/7619283). The scoring anchors above are OmicsAgent engineering choices, not a published validated scale.
