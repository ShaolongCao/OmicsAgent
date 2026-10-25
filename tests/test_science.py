import anndata as ad
import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError
from scipy import sparse

from omicsagent.analysis import counts, reference_transfer, validate
from omicsagent.annotation import audit_labels
from omicsagent.config import RunConfig, load_config
from omicsagent.spatial import _null_stats, bh, radius_graph


def test_spatial_config_requires_units_and_sample():
    with pytest.raises(ValidationError, match="sample_key"):
        RunConfig(
            input="a.h5ad", out="out", tissue="kidney", modality="spatial", platform="imaging"
        )


def test_paths_resolve_relative_to_config(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    assert cfg.input == demo_dir / "synthetic.h5ad"


@pytest.mark.parametrize("data", [[0.1, 1], [-1, 1], [float("nan"), 1]])
def test_counts_reject_invalid_values(data):
    obj = ad.AnnData(sparse.csr_matrix([data]))
    with pytest.raises(ValueError):
        counts(obj, "X")


def test_missing_counts_not_silently_using_x():
    obj = ad.AnnData(sparse.csr_matrix([[1, 2]]))
    with pytest.raises(ValueError, match="Missing raw counts"):
        counts(obj, "counts")


def test_duplicate_ids_rejected(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    obj = ad.read_h5ad(cfg.input)
    obj.var_names = ["duplicate"] * obj.n_vars
    with pytest.raises(ValueError, match="unique"):
        validate(obj, cfg)


def test_gate_preserves_original_and_marks_coverage_gap():
    obj = ad.AnnData(
        sparse.csr_matrix([[1, 0], [0, 0], [1, 0]]), var=pd.DataFrame(index=["A", "B"])
    )
    obj.layers["counts"] = obj.X.copy()
    obj.obs["oa_candidate"] = ["type1", "type1", "type2"]
    gate, _ = audit_labels(obj, {"type1": ["A"], "type2": ["absent"]}, 1)
    assert list(obj.obs.oa_candidate) == ["type1", "type1", "type2"]
    assert list(obj.obs.oa_label) == ["type1", "Uncertain", "type2"]
    assert list(obj.obs.oa_gate_status) == ["marker_supported", "unsupported", "ungateable"]
    assert gate.loc[gate.label == "type2", "markers_on_panel"].iloc[0] == 0


def test_gate_threshold_cannot_exceed_available_panel():
    obj = ad.AnnData(sparse.csr_matrix([[1, 0]]), var=pd.DataFrame(index=["A", "B"]))
    obj.layers["counts"] = obj.X.copy()
    obj.obs["oa_candidate"] = ["type1"]
    audit_labels(obj, {"type1": ["A", "absent"]}, 2)
    assert obj.obs.oa_gate_status.iloc[0] == "ungateable"


def test_graph_never_crosses_samples_and_preserves_colocated_other_cells():
    xy = np.array([[0, 0], [0, 0], [1, 0], [0, 0]])
    graph = radius_graph(xy, np.array(["a", "a", "a", "b"]), 2, 100)
    assert graph[0, 1] == 1
    assert graph[0, 3] == 0
    assert graph.diagonal().sum() == 0
    assert graph.nnz == 6


def test_graph_budget_checked_before_materializing():
    with pytest.raises(ValueError, match="max_edges"):
        radius_graph(np.zeros((100, 2)), np.repeat("a", 100), 10, 10)


def test_degenerate_null_is_not_reported_significant():
    z, p = _null_stats(np.array([5.0]), np.ones((19, 1)))
    assert np.isnan(z[0]) and np.isnan(p[0])


def test_permutation_pvalue_uses_plus_one():
    _, p = _null_stats(np.array([100.0]), np.arange(19).reshape(19, 1))
    assert p[0] == 1 / 20


def test_bh_preserves_nan_and_monotonic_adjustment():
    assert np.allclose(bh(np.array([0.01, 0.04, 0.03, np.nan]))[:3], [0.03, 0.04, 0.04])


def test_reference_projection_matches_self_and_rejects_low_overlap(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    cfg.reference = demo_dir / "reference.h5ad"
    cfg.reference_label_key = "known_label"
    obj = ad.read_h5ad(cfg.input)
    result = reference_transfer(obj, cfg)
    assert result["shared_genes"] == obj.n_vars
    assert (obj.obs.oa_candidate == obj.obs.known_label).mean() > 0.95
    # This is only a mechanics test, not independent biological validation.
    cfg.min_shared_genes = 1000
    with pytest.raises(ValueError, match="shared genes"):
        reference_transfer(obj, cfg)
