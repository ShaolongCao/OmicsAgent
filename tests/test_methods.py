import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from omicsagent.config import load_config
from omicsagent.methods import apply_split, reference_data, split_core, validate_weights
from omicsagent.workflow import run_pipeline


def test_split_hand_calculated_all_types_and_zero_support():
    raw = sparse.csr_matrix([[12., 10., 6.]])
    weights = np.array([[0.5, 0.3, 0.2]])
    profiles = np.array([[2., 0., 0.], [1., 1., 0.], [3., 1., 0.]])
    got = split_core(raw, weights, profiles).toarray()
    np.testing.assert_allclose(got, [[12 / 1.9, 0, 0]])
    np.testing.assert_allclose(raw.toarray(), [[12, 10, 6]])


def test_split_chunk_invariance_and_invalid_weights():
    raw = sparse.csr_matrix(np.ones((10, 4)))
    weights = np.tile([0.8, 0.2], (10, 1))
    profiles = np.array([[1, 2, 3, 4], [4, 3, 2, 1]])
    np.testing.assert_allclose(split_core(raw, weights, profiles, 12).toarray(),
                               split_core(raw, weights, profiles).toarray())
    with pytest.raises(ValueError):
        split_core(raw, -weights, profiles)
    with pytest.raises(ValueError):
        validate_weights(pd.DataFrame([[1, 0]], index=["wrong"]), pd.Index(["cell"]))


def test_split_preserves_uncovered_genes_and_input(demo_dir):
    cfg = load_config(demo_dir / "tacco-split.json")
    cfg.out.mkdir()
    obj = ad.read_h5ad(cfg.input)
    _, ref = reference_data(obj, cfg)
    ref = ref[:, :-1].copy()
    labels = sorted(ref.obs.known_label.unique())
    weights = pd.DataFrame(np.ones((obj.n_obs, 3)) / 3, index=obj.obs_names, columns=labels)
    obj.obs["oa_candidate"] = weights.idxmax(axis=1)
    before = obj.layers["counts"].copy()
    result = apply_split(obj, ref, cfg, weights)
    np.testing.assert_array_equal(before.toarray(), obj.layers["counts"].toarray())
    np.testing.assert_array_equal(before[:, -1].toarray(), obj.layers["counts_split"][:, -1].toarray())
    assert result["unchanged_genes"] == 1
    assert np.all(obj.layers["counts_split"].toarray() <= before.toarray() + 1e-9)


def test_config_requires_valid_split_scope(demo_dir):
    cfg = json.loads((demo_dir / "tacco-split.json").read_text())
    cfg["annotation_method"] = "knn"
    from omicsagent.config import RunConfig
    with pytest.raises(ValueError, match="SPLIT"):
        RunConfig.model_validate(cfg)


def test_tacco_split_pipeline_and_cached_replay(demo_dir, monkeypatch):
    pytest.importorskip("tacco")
    cfg = load_config(demo_dir / "tacco-split.json")
    assert run_pipeline(cfg)["status"] == "completed"
    original = cfg.out
    first = ad.read_h5ad(original / "result.h5ad")
    assert first.obsm["oa_weights"].shape == (180, 3)
    assert "counts_split" in first.layers
    def must_not_train(*args):
        raise AssertionError("Replay must not rerun reference annotation")
    monkeypatch.setattr("omicsagent.workflow.annotate_reference", must_not_train)
    cfg.out = demo_dir / "replay-advanced"
    run_pipeline(cfg, replay=original)
    second = ad.read_h5ad(cfg.out / "result.h5ad")
    np.testing.assert_allclose(first.layers["counts_split"].toarray(),
                               second.layers["counts_split"].toarray())
