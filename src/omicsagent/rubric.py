"""Versioned evidence rubric. Scores describe policy evidence, never probabilities."""
import math
from typing import Any

VERSION = '1.0'


def _true(value: Any) -> bool:
    return value is True or str(value).lower() == 'true'


def _label(value: Any) -> str | None:
    return None if value is None or str(value) in {'', 'Unknown', 'nan', 'None'} else str(value)


def _number(value: Any) -> float:
    try:
        number = float(value)
        return number if math.isfinite(number) else 0.
    except (TypeError, ValueError):
        return 0.


def assess_decision(row: dict[str, Any], thresholds: dict[str, Any] | None = None) -> dict[str, Any]:
    """Assess a recorded cluster/source-label subset, without manufacturing evidence.

    Applicability and external accuracy require independent measurements and remain
    unassessed. Existing marker gates include coverage, detection and score margin.
    """
    cfg = thresholds or {}
    def confident(level: str) -> bool:
        prefix = 'atlas' if level == 'l1' else 'atlas_l2'
        return (_label(row.get('atlas_' + level)) is not None
                and _number(row.get(prefix + '_purity')) >= cfg.get('atlas_purity', .7)
                and _number(row.get(prefix + '_vote')) >= cfg.get('atlas_vote', .6))

    best, target = _label(row.get('marker_best')), _label(row.get('target_l1'))
    subtype, subtype_best = _label(row.get('selected_l2')), _label(row.get('l2_marker_best'))
    strong, subtype_strong = _true(row.get('marker_strong')), _true(row.get('l2_marker_strong'))
    conflict = strong and confident('l1') and _label(row.get('atlas_l1')) != best
    subtype_conflict = subtype_strong and confident('l2') and _label(row.get('atlas_l2')) != subtype_best
    decision = row.get('decision')
    supported = bool(strong and best and best == target and not conflict)
    failures = []
    if decision in {'Keep', 'Reassign'} and not supported:
        failures.append('unsupported_l1_assignment')
    if decision == 'Reassign' and not (confident('l1') and _label(row.get('atlas_l1')) == best):
        failures.append('reassignment_without_reference_agreement')
    parents = cfg.get('l2_parent')
    if subtype and (decision != 'Keep' or not supported or not subtype_strong
                    or subtype != subtype_best or subtype_conflict
                    or (parents is not None and parents.get(subtype) != target)):
        failures.append('unsupported_l2_assignment')
    accepted = ('L2' if subtype else 'L1') if supported and decision == 'Keep' and not failures else None
    needs_subtype = row.get('l2_reason') not in {'no_subtype_panel', None} and not subtype
    if failures:
        action, reason = 'stop_unresolved', 'invalid_assignment'
    elif conflict or subtype_conflict:
        action, reason = 'review', 'marker_reference_conflict'
    elif decision == 'Reassign':
        action, reason = 'continue', 'validate_changed_assignments'
    elif accepted and not needs_subtype:
        action, reason = 'accept', 'supported_at_reported_granularity'
    else:
        action, reason = 'stop_unresolved', 'insufficient_evidence'

    check = None
    if action == 'review':
        check = {'kind': 'competing_markers', 'executor': 'scientist',
                 'trigger': reason, 'question': 'Do competing markers support a mixed population or a reference mismatch?',
                 'decision_effect': 'Resolve the contradiction before assigning the disputed level; otherwise retain deferral.'}
    elif action == 'continue':
        check = {'kind': 'repool', 'executor': 'pipeline', 'trigger': reason,
                 'question': 'Does the changed identity remain supported in its new pool?',
                 'decision_effect': 'Accept only if the new round independently passes the existing evidence gates.'}
    elif action == 'stop_unresolved' and not failures:
        check = {'kind': 'marker_coverage_and_separation', 'executor': 'scientist', 'trigger': reason,
                 'question': 'Are markers missing, poorly detected, or insufficiently separated?',
                 'decision_effect': 'Obtain discriminating evidence before assigning a finer identity; do not rerun unchanged inputs.'}

    def dim(score: int | None, evidence: str) -> dict[str, Any]:
        return {'score': score, 'evidence': evidence}

    dimensions = {
        'applicability': dim(None, 'Tissue/technology compatibility and out-of-reference assessment are not measured by this decision row.'),
        'support': dim(None if 'marker_strong' not in row else (2 if supported else 0),
                       f'Marker winner={best}; strong={strong}; proposed parent={target}.'),
        'counterevidence': dim(0 if failures else (2 if conflict or subtype_conflict else 1),
                              f'L1 conflict={conflict}; L2 conflict={subtype_conflict}. Other confounders are not ruled out.'),
        'granularity': dim(0 if failures else (2 if accepted else 1),
                           f'Accepted evidence level={accepted}; proposed subtype={subtype}.'),
        'selective_resolution': dim(None, 'Measured across controls or independently adjudicated cases, not inferred from one row.'),
        'next_check': dim(2 if check else None, str(check['trigger']) if check else 'No additional check required by measured gates.'),
        'progress': dim(None, 'Requires round transitions; see loop control record.'),
        'audit': dim(1 if row.get('cluster') else None, 'Decision fields recorded; artifact integrity and reversibility are separate checks.'),
    }
    return {'version': VERSION, 'cluster': row.get('cluster'), 'source_l1': row.get('source_l1'),
            'n_cells': int(_number(row.get('n_cells'))), 'accepted_level': accepted,
            'action': action, 'reason': reason, 'dimensions': dimensions,
            'hard_failures': failures, 'next_check': check}


def control_loop(assessments: list[dict[str, Any]], *, moved_cells: int, repeated: bool,
                 iteration: int, limit: int, skipped_cells: int,
                 numerical_stop: str | None) -> dict[str, Any]:
    """Only re-pooling changed assignments is currently executable by this pipeline."""
    if any(item['hard_failures'] for item in assessments):
        action, reason = 'stop_unresolved', 'invalid_assignment'
    elif moved_cells and repeated:
        action, reason = 'stop_unresolved', 'cycle'
    elif moved_cells and iteration >= limit:
        action, reason = 'stop_unresolved', 'budget_exhausted'
    elif moved_cells:
        action, reason = 'continue', 'validate_changed_assignments'
    elif any(item['action'] == 'review' for item in assessments):
        action, reason = 'review', 'unresolved_conflict_requires_scientist'
    elif skipped_cells or not assessments or any(item['action'] != 'accept' for item in assessments):
        action, reason = 'stop_unresolved', 'no_executable_progress'
    else:
        action, reason = 'accept', 'supported_at_reported_granularity'
    return {'version': VERSION, 'action': action, 'reason': reason,
            'iteration': iteration, 'budget_remaining': max(0, limit-iteration),
            'moved_cells': moved_cells, 'skipped_cells': skipped_cells,
            'numerical_stop_reason': numerical_stop,
            'progress': {'score': 1 if moved_cells else 0,
                         'meaning': 'Changed assignments require validation; movement is not new biological evidence.'},
            'accepted_cells': sum(item['n_cells'] for item in assessments if item['action'] == 'accept'),
            'unresolved_cells': skipped_cells + sum(item['n_cells'] for item in assessments if item['action'] != 'accept')}


def evaluate_rubric() -> dict[str, Any]:
    """Policy controls, not biological truth or a substitute for held-out validation."""
    from .review import evaluate_predictions
    base = dict(cluster='control', n_cells=10, source_l1='T', target_l1='T', decision='Keep',
                marker_best='T', marker_strong=True, atlas_l1='T', atlas_purity=.9, atlas_vote=.9,
                selected_l2='Unknown', l2_reason='no_subtype_panel')
    controls: list[tuple[str, dict[str, Any], str]] = [
        ('supported_parent', {}, 'accept'),
        ('coarse_conflict', {'atlas_l1': 'B', 'decision': 'Ambiguous'}, 'review'),
        ('no_evidence', {'marker_best': None, 'marker_strong': False, 'decision': 'Unknown'}, 'stop_unresolved'),
        ('blind_reference', {'marker_strong': False}, 'stop_unresolved'),
        ('unverified_subtype', {'selected_l2': 'CD4'}, 'stop_unresolved'),
        ('repool_required', {'decision': 'Reassign', 'source_l1': 'B'}, 'continue'),
    ]
    observed = [assess_decision(base | change)['action'] for _, change, _ in controls]
    vocabulary = {'accept': 'resolved', 'review': 'deferred', 'stop_unresolved': 'unknown', 'continue': 'reassigned'}
    expected = [vocabulary[expected] for _, _, expected in controls]
    result = evaluate_predictions(expected, [vocabulary[value] for value in observed])
    loop_controls = [
        ('small_move', 1, False, 1, 'tolerance', 'continue'),
        ('cycle', 1, True, 2, 'cycle', 'stop_unresolved'),
        ('budget', 1, False, 3, 'max_iterations', 'stop_unresolved'),
        ('stable_supported', 0, True, 1, 'exact_stability', 'accept'),
    ]
    loop_results = []
    for name, moved, repeated, iteration, numerical, expected_action in loop_controls:
        observed_action = control_loop([assess_decision(base)], moved_cells=moved, repeated=repeated,
                                       iteration=iteration, limit=3, skipped_cells=0,
                                       numerical_stop=numerical)['action']
        loop_results.append({'id': name, 'expected': expected_action, 'observed': observed_action})
    result['passed'] = result['passed'] and all(r['expected'] == r['observed'] for r in loop_results)
    return {'version': VERSION, 'scope': 'Synthetic policy controls, not biological accuracy', **result,
            'loop_controls': loop_results,
            'controls': [{'id': name, 'expected': expected, 'observed': actual}
                         for (name, _, expected), actual in zip(controls, observed, strict=True)],
            'all_unknown_control': evaluate_predictions(expected, ['unknown'] * len(expected))}
