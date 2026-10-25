import numpy as np
import pytest

from omicsagent.hierarchy import HierarchyConfig, decide_cluster, stopping_reason


def test_atlas_cannot_override_markers():
    result = decide_cluster('A', 'A', True, 'B', 0.99, 0.99)
    assert result == ('Ambiguous', 'A')
    assert decide_cluster('A', None, False, 'B', 1, 1) == ('Unknown', 'A')
    assert decide_cluster('A', 'B', True, 'B', .8, .9) == ('Reassign', 'B')
    assert decide_cluster('A', 'B', True, None, 0, 0) == ('Ambiguous', 'A')
    assert decide_cluster('A', 'A', False, 'A', .9, .9) == ('Ambiguous', 'A')


def test_stop_boundary_cycles_and_cap():
    assert stopping_reason([.005], .005, False, 1, 3) == 'tolerance'
    assert stopping_reason([.006], .005, False, 1, 3) is None
    assert stopping_reason([.02], .005, True, 2, 3) == 'cycle'
    assert stopping_reason([.02], .005, False, 3, 3) == 'max_iterations'
    assert stopping_reason([], .005, False, 1, 3) == 'exact_stability'


def test_invalid_ontology_rejected():
    with pytest.raises(ValueError):
        HierarchyConfig.model_validate(dict(initial_l1_key='l1', compartments={'p':['A'],'q':['A']},
            l1_markers={'A':['x','y']}, l2_markers={'a':['z']}, l2_parent={'a':'B'}))


def test_nonfinite_threshold_rejected():
    with pytest.raises(ValueError):
        HierarchyConfig.model_validate(dict(initial_l1_key='l1', compartments={'p':['A','B']},
            l1_markers={'A':['x','y'],'B':['z','w']}, l2_markers={}, l2_parent={}, margin=np.nan))


def test_unmeasurable_modules_are_not_zero_evidence():
    import anndata as ad
    import pandas as pd
    from scipy import sparse

    from omicsagent.hierarchy import module_scores
    obj = ad.AnnData(sparse.csr_matrix(np.ones((6, 5))),
                     var=pd.DataFrame(index=['a','b','c','d','e']))
    scores, coverage = module_scores(obj, {'A':['a','b'], 'B':['missing','other']}, set(), 2)
    assert scores.isna().all().all()
    assert coverage['B'] == []


def test_degenerate_pool_is_explicitly_unresolved(tmp_path):
    import anndata as ad
    import pandas as pd
    from scipy import sparse

    from omicsagent.hierarchy import run_hierarchy
    obj = ad.AnnData(sparse.csr_matrix(np.ones((30, 8))),
                     obs=pd.DataFrame({'l1':['A']*30},index=[f'c{i}' for i in range(30)]),
                     var=pd.DataFrame(index=list('abcdefgh')))
    obj.layers['counts'] = obj.X.copy()
    cfg = HierarchyConfig(initial_l1_key='l1',compartments={'p':['A'],'q':['B']},
                          l1_markers={'A':['a','b'],'B':['c','d']})
    result = run_hierarchy(obj, cfg, tmp_path, 42, 'scrna')
    assert (obj.obs.oa_hierarchy_status == 'Unknown').all()
    assert result['stop_reason'] == 'exact_stability'
    assert result['converged'] is False


def test_tied_scores_remain_ambiguous_at_zero_threshold():
    import anndata as ad
    import pandas as pd
    from scipy import sparse

    from omicsagent.hierarchy import _best
    obj = ad.AnnData(sparse.csr_matrix(np.ones((4,4))),var=pd.DataFrame(index=list('abcd')))
    obj.layers['counts'] = obj.X.copy()
    cfg = HierarchyConfig(initial_l1_key='l1',compartments={'p':['A','B']},
                          l1_markers={'A':['a','b'],'B':['c','d']},margin=0)
    scores = pd.DataFrame({'A':[1.]*4,'B':[1.]*4},index=obj.obs_names)
    _, strong, _ = _best(scores,np.arange(4),cfg.l1_markers,obj,cfg,0)
    assert not strong
