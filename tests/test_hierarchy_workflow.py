import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from omicsagent.config import load_config
from omicsagent.demo import make_hierarchy_demo
from omicsagent.workflow import run_pipeline


@pytest.fixture
def hierarchy_demo(tmp_path):
    root = tmp_path / 'hierarchy'
    make_hierarchy_demo(root)
    return root


def test_full_hierarchy_preserves_counts_and_corrects_misplaced_cells(hierarchy_demo):
    cfg = load_config(hierarchy_demo / 'scrna.json')
    assert run_pipeline(cfg)['status'] == 'completed'
    result = ad.read_h5ad(cfg.out / 'result.h5ad')
    original = ad.read_h5ad(cfg.input)
    np.testing.assert_array_equal(result.layers['counts'].toarray(), original.layers['counts'].toarray())
    wrong = original.obs.initial_l1 != original.obs.truth_l1
    assert wrong.sum() > 0
    assert (result.obs.loc[wrong, 'oa_l1'] == original.obs.loc[wrong, 'truth_l1']).all()
    assert (cfg.out / 'hierarchy' / 'decisions.csv').exists()
    assert (cfg.out / 'report.html').exists()
    assert (cfg.out / 'review.json').exists()
    summary = json.loads((cfg.out / 'summary.json').read_text())['hierarchy']
    assert summary['stop_reason'] == 'exact_stability'
    assert summary['iterations'] >= 2
    assert (result.obs.oa_l2 != 'Unknown').any()


def test_hierarchy_spatial_replay_rejects_tampering(hierarchy_demo):
    cfg = load_config(hierarchy_demo / 'imaging.json')
    run_pipeline(cfg)
    first = cfg.out
    result = ad.read_h5ad(first / 'result.h5ad')
    row, col = result.obsp['oa_spatial_connectivities'].nonzero()
    assert np.all(result.obs['sample'].to_numpy()[row] == result.obs['sample'].to_numpy()[col])
    assert np.array_equal(result.obs.oa_spatial_label.astype(str), result.obs.oa_hierarchy_label.astype(str))
    cfg.out = hierarchy_demo / 'replay'
    run_pipeline(cfg, replay=first)
    replay = ad.read_h5ad(cfg.out / 'result.h5ad')
    assert np.array_equal(result.obs.oa_l1, replay.obs.oa_l1)
    with (first / 'cells.csv').open('a') as handle:
        handle.write('\ntampered\n')
    cfg.out = hierarchy_demo / 'tampered-replay'
    with pytest.raises(ValueError, match='artifact'):
        run_pipeline(cfg, replay=first)


def test_one_round_invalidates_moved_subtype(hierarchy_demo):
    hp = hierarchy_demo / 'hierarchy.json'
    config = json.loads(hp.read_text())
    config['max_iterations'] = 1
    hp.write_text(json.dumps(config))
    cfg = load_config(hierarchy_demo / 'scrna.json')
    run_pipeline(cfg)
    result = ad.read_h5ad(cfg.out / 'result.h5ad')
    moved = result.obs.oa_initial_l1 != result.obs.oa_l1
    assert moved.any()
    assert (result.obs.loc[moved, 'oa_l2'] == 'Unknown').all()
    assert json.loads((cfg.out / 'summary.json').read_text())['hierarchy']['stop_reason'] == 'max_iterations'


def test_replay_does_not_refit_hierarchy(hierarchy_demo, monkeypatch):
    cfg = load_config(hierarchy_demo/'scrna.json')
    run_pipeline(cfg)
    first = cfg.out
    def no_refit(*args, **kwargs):
        raise AssertionError('Hierarchy was refit during replay')
    monkeypatch.setattr('omicsagent.workflow.run_hierarchy',no_refit)
    cfg.out = hierarchy_demo/'cached-replay'
    assert run_pipeline(cfg,replay=first)['status'] == 'completed'


def test_conflicting_atlas_subtype_does_not_validate_l2(hierarchy_demo):
    reference_path = hierarchy_demo/'reference.h5ad'
    ref = ad.read_h5ad(reference_path)
    ref.obs['truth_l2'] = ref.obs.truth_l2.astype(str).replace({'Mono':'DC','DC':'Mono'})
    ref.write_h5ad(reference_path)
    cfg = load_config(hierarchy_demo/'scrna.json')
    run_pipeline(cfg)
    result = ad.read_h5ad(cfg.out/'result.h5ad')
    myeloid = result.obs.truth_l1 == 'Myeloid'
    assert (result.obs.loc[myeloid,'oa_l2'] == 'Unknown').all()
    decisions = pd.read_csv(cfg.out/'hierarchy'/'decisions.csv')
    kept = decisions[(decisions.source_l1 == 'Myeloid') & (decisions.decision == 'Keep')]
    assert not kept.empty
    assert (kept.l1_reason == 'marker_supported_current').all()
    assert (kept.l2_reason == 'marker_reference_conflict').all()
    assert (kept.l2_marker_strong).all()
    assert (kept.selected_l2 == 'Unknown').all()


def test_bounded_agent_calls_hierarchy_and_generates_report(hierarchy_demo):
    from test_agent_evidence import FakeClient

    from omicsagent.agent import run_agent
    cfg = load_config(hierarchy_demo/'scrna.json')
    cfg.model = 'test-model'
    cfg.max_total_tokens = 1000000
    client = FakeClient(['inspect_data','preprocess_data','find_markers','annotate',
                         'refine_hierarchy','generate_report'])
    assert run_agent(cfg,client=client)['status'] == 'completed'
    result = ad.read_h5ad(cfg.out/'result.h5ad')
    assert 'oa_l1' in result.obs


def test_harmony_records_applied_batches(hierarchy_demo):
    pytest.importorskip('harmonypy')
    path = hierarchy_demo/'hierarchy.json'
    hc = json.loads(path.read_text())
    hc['batch_key'] = 'sample'
    path.write_text(json.dumps(hc))
    cfg = load_config(hierarchy_demo/'scrna.json')
    run_pipeline(cfg)
    rounds = json.loads((cfg.out/'hierarchy'/'iterations.json').read_text())
    assert all(pool['harmony'] == 'applied' for row in rounds for pool in row['pools'])


def test_all_tiny_pools_produce_readable_empty_tables_and_report(hierarchy_demo):
    import pandas as pd
    path = hierarchy_demo/'hierarchy.json'
    hc = json.loads(path.read_text())
    hc['min_pool_size'] = 1000
    path.write_text(json.dumps(hc))
    cfg = load_config(hierarchy_demo/'scrna.json')
    run_pipeline(cfg)
    assert (cfg.out/'report.md').exists()
    assert pd.read_csv(cfg.out/'hierarchy'/'dotplot_data.csv').empty
    assert pd.read_csv(cfg.out/'hierarchy'/'decisions.csv').empty
    result = ad.read_h5ad(cfg.out/'result.h5ad')
    assert len(result) == 240
    assert (result.obs.oa_hierarchy_status == 'Unknown').all()
