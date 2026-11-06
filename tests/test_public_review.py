import numpy as np
import pandas as pd
import pytest


def test_donor_split_is_disjoint_and_preserves_original_labels():
    from omicsagent.public_review import partition_cells
    obs = pd.DataFrame({'donor': ['TSP2', 'TSP2', 'TSP11', 'TSP14'],
                        'method': ['10X']*4,
                        'cell_ontology_class': ['B cell', 'CD4-positive, alpha-beta T cell',
                                                'B cell', 'CD8-positive, alpha-beta T cell'],
                        'free_annotation': ['naive B cell', 'CD4 T cell', 'memory B cell', 'CD8 T cell']},
                       index=['a', 'b', 'c', 'd'])
    ref, query = partition_cells(obs)
    assert list(ref.index) == ['a', 'b']
    assert list(query.index) == ['c', 'd']
    assert set(ref.donor).isdisjoint(query.donor)
    assert query.loc['c', 'initial_l1'] == 'B cell'
    assert query.loc['d', 'published_label'] == 'CD8-positive, alpha-beta T cell'
    assert obs.loc['d', 'cell_ontology_class'] == 'CD8-positive, alpha-beta T cell'


def test_normalized_matrix_is_rejected_as_counts():
    import anndata as ad
    from scipy import sparse

    from omicsagent.public_review import raw_source
    obj = ad.AnnData(sparse.csr_matrix(np.array([[0., 1.5], [1., 2.]])))
    with pytest.raises(ValueError, match='integer'):
        raw_source(obj)


def test_missing_donor_or_one_query_donor_is_explicit():
    from omicsagent.public_review import partition_cells
    with pytest.raises(ValueError, match='metadata'):
        partition_cells(pd.DataFrame({'method': ['10X']}))
