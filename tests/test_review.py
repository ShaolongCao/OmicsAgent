import json

import pandas as pd
import pytest


def fixture_run(tmp_path):
    h = tmp_path / 'hierarchy'
    h.mkdir()
    rows, cells = [], []
    for cluster, subtype, reason in [('a', 'CD4', 'marker_supported'),
                                     ('b', 'Unknown', 'marker_reference_conflict'),
                                     ('c', 'Unknown', 'insufficient_marker_support')]:
        rows.append(dict(iteration=1, cluster=cluster, pool='immune', n_cells=20,
                         source_l1='T cell', target_l1='T cell', marker_best='T cell',
                         marker_strong=cluster != 'c', marker_margin=1.2,
                         atlas_l1='T cell', atlas_purity=.95, atlas_vote=.9,
                         atlas_l2='CD8' if cluster == 'b' else 'CD4',
                         atlas_l2_purity=.9, atlas_l2_vote=.8, decision='Keep',
                         l1_reason='marker_supported_current', l2_marker_best='CD4',
                         l2_marker_strong=cluster != 'c', l2_marker_margin=.8,
                         l2_reason=reason, selected_l2=subtype))
        for i in range(20):
            cells.append(dict(cell_id=f'{cluster}{i}', cluster=cluster,
                              before_l1='T cell', after_l1='T cell', l2=subtype,
                              donor='d1' if i < 10 else 'd2'))
    pd.DataFrame(rows).to_csv(h/'decisions.csv', index=False)
    pd.DataFrame(cells).to_csv(h/'round_1_cells.csv', index=False)
    pd.DataFrame(cells).set_index('cell_id').to_csv(tmp_path/'cells.csv')
    pd.DataFrame([dict(iteration=1, cluster='b', level='l1', label='T cell', gene='CD3D',
                       mean_log_expression=2., fraction_detected=.8)]).to_csv(h/'dotplot_data.csv', index=False)
    (h/'marker_coverage.json').write_text(json.dumps({'l1': {'T cell': ['CD3D', 'CD3E']}, 'l2': {}}))
    (h/'config.json').write_text(json.dumps({'atlas_purity': .7, 'atlas_vote': .6, 'margin': .5}))
    (tmp_path/'config.json').write_text('{}')
    return tmp_path


def test_measured_cases_distinguish_coarse_support_and_subtype_conflict(tmp_path):
    from omicsagent.review import build_review
    report = build_review(fixture_run(tmp_path))
    by_kind = {c['kind']: c for c in report['cases']}
    assert set(by_kind) == {'agreement', 'conflict', 'insufficient'}
    assert by_kind['conflict']['l1'] == 'T cell'
    assert by_kind['conflict']['l2'] == 'Unknown'
    assert by_kind['conflict']['l1_supported'] is True
    assert by_kind['conflict']['counterevidence']
    assert by_kind['conflict']['next_check']['kind'] != 'donor'
    assert report['selection']['eligible']['conflict'] == 1


def test_report_requires_real_measurements_for_donor_suggestion(tmp_path):
    from omicsagent.review import build_review
    report = build_review(fixture_run(tmp_path), donor_key='donor')
    conflict = next(c for c in report['cases'] if c['kind'] == 'conflict')
    assert conflict['next_check']['kind'] != 'donor'


def test_marker_winner_opposing_retained_label_is_counterevidence(tmp_path):
    from omicsagent.review import build_review
    root = fixture_run(tmp_path)
    frame = pd.read_csv(root/'hierarchy'/'decisions.csv')
    frame.loc[frame.cluster == 'b', ['marker_best', 'decision']] = ['B cell', 'Ambiguous']
    frame.to_csv(root/'hierarchy'/'decisions.csv', index=False)
    conflict = next(c for c in build_review(root)['cases'] if c['kind'] == 'conflict')
    assert conflict['l1_supported'] is False
    assert any('L1 marker winner B cell' in text for text in conflict['counterevidence'])
    assert not any('L1 marker winner B cell' in text for text in conflict['support'])


@pytest.mark.parametrize('tied', [False, True])
def test_donor_trigger_uses_evaluable_pairs_and_excludes_tied_winners(tmp_path, tied):
    from omicsagent.review import build_review
    root = fixture_run(tmp_path)
    ids = [f'b{i}' for i in range(20)]
    pd.DataFrame({'CD4': 1., 'CD8': 1. if tied else 0., 'Other': -1.}, index=ids).to_csv(
        root/'hierarchy'/'l2_scores.csv')
    pd.DataFrame({'l2': ['CD8']*10 + ['CD4']*10}, index=ids).to_csv(
        root/'hierarchy'/'atlas_predictions.csv')
    params = json.loads((root/'hierarchy'/'config.json').read_text())
    params['l2_parent'] = {'CD4': 'T cell', 'CD8': 'T cell', 'Other': 'T cell'}
    (root/'hierarchy'/'config.json').write_text(json.dumps(params))
    conflict = next(c for c in build_review(root, 'donor')['cases'] if c['kind'] == 'conflict')
    if tied:
        assert conflict['next_check']['kind'] != 'donor'
    else:
        assert conflict['next_check']['kind'] == 'donor'
        assert [d['n_evaluable'] for d in conflict['donor_summary']] == [10, 10]
        assert [d['disagreement_fraction'] for d in conflict['donor_summary']] == [1., 0.]


def test_unknown_policy_cannot_pass_behavior_evaluation():
    from omicsagent.review import evaluate_predictions
    expected = ['resolved', 'deferred', 'unknown', 'reassigned']
    result = evaluate_predictions(expected, ['unknown'] * 4)
    assert result['passed'] is False
    assert result['resolved_coverage'] == 0
    assert result['inappropriate_assignment'] == 0
    with pytest.raises(ValueError):
        evaluate_predictions(expected, [])


def test_marker_only_resolved_label_is_not_an_insufficient_evidence_case(tmp_path):
    from omicsagent.review import build_review
    root = fixture_run(tmp_path)
    frame = pd.read_csv(root/'hierarchy'/'decisions.csv')
    frame.loc[frame.cluster == 'c', ['selected_l2', 'l2_reason', 'marker_strong',
                                     'l2_marker_strong', 'atlas_l2_purity']] = ['CD4', 'marker_supported', True, True, .4]
    frame.to_csv(root/'hierarchy'/'decisions.csv', index=False)
    report = build_review(root)
    assert report['selection']['eligible']['insufficient'] == 0
    assert 'insufficient' in report['selection']['missing']


def test_embedded_report_escapes_untrusted_labels(tmp_path):
    from omicsagent.review import write_review
    root = fixture_run(tmp_path)
    frame = pd.read_csv(root/'hierarchy'/'decisions.csv')
    frame.loc[0, 'target_l1'] = '</script><script>window.pwned=1</script>'
    frame.to_csv(root/'hierarchy'/'decisions.csv', index=False)
    write_review(root)
    html = (root/'report.html').read_text(encoding='utf-8')
    assert '</script><script>window.pwned' not in html
    assert '\\u003c/script\\u003e' in html


def test_old_run_without_new_evidence_is_not_claimed_resolved(tmp_path):
    from omicsagent.review import build_review
    root = fixture_run(tmp_path)
    frame = pd.read_csv(root/'hierarchy'/'decisions.csv').drop(
        columns=['l2_reason', 'l2_marker_strong', 'selected_l2'])
    frame.to_csv(root/'hierarchy'/'decisions.csv', index=False)
    assert all(c['kind'] != 'agreement' for c in build_review(root)['cases'])


def test_review_summary_uses_initial_labels_and_discloses_cluster_evidence_scope(tmp_path):
    from omicsagent.review import build_review
    root = fixture_run(tmp_path)
    cells = pd.read_csv(root/'cells.csv', index_col=0)
    cells['oa_initial_l1'] = 'B cell'
    cells.to_csv(root/'cells.csv')
    path = root/'hierarchy'/'round_1_cells.csv'
    members = pd.read_csv(path)
    extra = members[members.cluster == 'b'].iloc[:2].copy()
    extra['cell_id'] = ['extra1', 'extra2']
    extra['before_l1'] = 'B cell'
    pd.concat([members, extra]).to_csv(path, index=False)
    conflict = next(c for c in build_review(root)['cases'] if c['kind'] == 'conflict')
    assert conflict['original_l1'] == 'B cell'
    assert conflict['evidence_n_cells'] == 22
    assert conflict['n_cells'] == 20


def test_empty_decision_table_discloses_missing_case_types(tmp_path):
    from omicsagent.review import build_review
    root = fixture_run(tmp_path)
    frame = pd.read_csv(root/'hierarchy'/'decisions.csv').iloc[:0]
    frame.to_csv(root/'hierarchy'/'decisions.csv', index=False)
    report = build_review(root)
    assert report['selection']['complete'] is False
    assert report['selection']['missing'] == ['agreement', 'conflict', 'insufficient']


def test_skipped_final_round_cannot_resurrect_previous_supported_calls(tmp_path):
    from omicsagent.review import build_review
    root = fixture_run(tmp_path)
    members = pd.read_csv(root/'hierarchy'/'round_1_cells.csv')
    members['l2'] = 'Unknown'
    members['decision'] = 'Unknown'
    members.to_csv(root/'hierarchy'/'round_2_cells.csv', index=False)
    (root/'hierarchy'/'iterations.json').write_text(json.dumps([
        {'iteration': 1, 'pools': []},
        {'iteration': 2, 'pools': [{'pool': 'immune', 'reason': 'degenerate_pool'}]}]))
    report = build_review(root)
    assert report['cases'] == []
    assert report['selection']['final_iteration'] == 2
    assert report['selection']['skipped_pools'][0]['reason'] == 'degenerate_pool'
