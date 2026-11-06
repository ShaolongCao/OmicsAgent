import pytest

from omicsagent.rubric import assess_decision, control_loop, evaluate_rubric


def row(**changes):
    value = dict(cluster='c', source_l1='T', target_l1='T', n_cells=10,
                 decision='Keep', marker_best='T', marker_strong=True,
                 atlas_l1='T', atlas_purity=.95, atlas_vote=.9,
                 selected_l2='CD4', l2_marker_best='CD4', l2_marker_strong=True,
                 atlas_l2='CD4', atlas_l2_purity=.95, atlas_l2_vote=.9,
                 l2_reason='marker_supported')
    return value | changes


def test_supported_and_coarse_only_are_resolved():
    assert assess_decision(row())['action'] == 'accept'
    result = assess_decision(row(selected_l2='Unknown', l2_reason='no_subtype_panel'))
    assert result['action'] == 'accept'
    assert result['accepted_level'] == 'L1'
    assert result['dimensions']['applicability']['score'] is None


def test_conflict_retains_supported_parent_and_requests_review():
    result = assess_decision(row(selected_l2='Unknown', atlas_l2='CD8',
                                 l2_reason='marker_reference_conflict'))
    assert result['action'] == 'review'
    assert result['accepted_level'] == 'L1'
    assert not result['hard_failures']
    assert result['next_check']['kind'] == 'competing_markers'


@pytest.mark.parametrize('changes', [dict(atlas_l2='CD8'), dict(marker_strong=False),
                                     dict(l2_marker_strong=False), dict(target_l1='B')])
def test_unsupported_assignments_are_hard_failures(changes):
    assert assess_decision(row(**changes))['hard_failures']


def test_missing_evidence_is_not_a_pass():
    result = assess_decision({'decision': 'Unknown'})
    assert result['action'] == 'stop_unresolved'
    assert result['dimensions']['support']['score'] is None


def test_small_moves_continue_even_if_numerically_converged():
    assessment = assess_decision(row(decision='Reassign', source_l1='B', selected_l2='Unknown',
                                     l2_reason='awaiting_new_pool'))
    result = control_loop([assessment], moved_cells=1, repeated=False, iteration=1, limit=3,
                          skipped_cells=0, numerical_stop='tolerance')
    assert result['action'] == 'continue'
    assert result['reason'] == 'validate_changed_assignments'


@pytest.mark.parametrize('repeated,iteration,reason', [(True,2,'cycle'), (False,3,'budget_exhausted')])
def test_changed_labels_cannot_bypass_cycle_or_budget(repeated, iteration, reason):
    result = control_loop([assess_decision(row())], moved_cells=1, repeated=repeated,
                          iteration=iteration, limit=3, skipped_cells=0, numerical_stop=None)
    assert result['action'] == 'stop_unresolved'
    assert result['reason'] == reason


def test_stability_does_not_resolve_conflict_or_skipped_pools():
    conflict = assess_decision(row(selected_l2='Unknown', atlas_l2='CD8'))
    args = dict(moved_cells=0, repeated=True, iteration=1, limit=3, numerical_stop='exact_stability')
    assert control_loop([conflict], skipped_cells=0, **args)['action'] == 'review'
    assert control_loop([], skipped_cells=20, **args)['action'] == 'stop_unresolved'


def test_behavior_controls_reject_blanket_abstention():
    result = evaluate_rubric()
    assert result['passed']
    assert not result['all_unknown_control']['passed']
    assert result['all_unknown_control']['resolved_coverage'] == 0
