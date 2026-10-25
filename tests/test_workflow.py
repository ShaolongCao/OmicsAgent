import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from omicsagent.config import load_config
from omicsagent.trace import sha256
from omicsagent.workflow import Workflow, run_pipeline


@pytest.mark.parametrize("config_name", ["scrna.json", "spatial.json", "visium.json"])
def test_complete_workflows(demo_dir, config_name):
    cfg = load_config(demo_dir / config_name)
    assert run_pipeline(cfg)["status"] == "completed"
    result = ad.read_h5ad(cfg.out / "result.h5ad")
    original = ad.read_h5ad(cfg.input)
    assert np.array_equal(result.layers["counts"].toarray(), original.layers["counts"].toarray())
    assert result.n_vars == original.n_vars  # HVG PCA must not discard marker genes.
    assert np.all(result.layers["log1p"].data >= 0)
    manifest = json.loads((cfg.out / "manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert manifest["artifacts"]["report.md"] == sha256(cfg.out / "report.md")
    assert (cfg.out / "analysis.ipynb").exists()
    if cfg.modality == "spatial":
        row, col = result.obsp["oa_spatial_connectivities"].nonzero()
        samples = result.obs[cfg.sample_key].to_numpy()
        assert np.all(samples[row] == samples[col])
        assert (cfg.out / "neighborhood_enrichment.csv").exists()
        if cfg.platform == "imaging":
            assert not pd.read_csv(cfg.out / "center_niches.csv").empty
        else:
            assert list(result.obs.oa_spatial_label) == list(result.obs.oa_cluster)


def test_replay_reproduces_labels(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    run_pipeline(cfg)
    previous = cfg.out
    cfg.out = demo_dir / "replay"
    run_pipeline(cfg, replay=previous)
    first = pd.read_csv(previous / "cells.csv")
    second = pd.read_csv(cfg.out / "cells.csv")
    assert first.oa_label.equals(second.oa_label)


def test_out_of_order_and_unknown_tools_are_structured_errors(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    wf = Workflow(cfg)
    assert not wf.call("preprocess_data", {})["ok"]
    assert not wf.call("__getattribute__", {})["ok"]
    assert wf.stage == "new"


def test_existing_output_not_overwritten(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    Workflow(cfg)
    with pytest.raises(ValueError, match="not empty"):
        Workflow(cfg)


def test_empty_qc_has_failed_manifest(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    cfg.min_genes = 10000
    with pytest.raises(RuntimeError, match="QC left"):
        run_pipeline(cfg)
    assert json.loads((cfg.out / "manifest.json").read_text())["status"] == "failed"


def test_missing_input_still_records_failure(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    cfg.input = demo_dir / "missing.h5ad"
    with pytest.raises(RuntimeError, match="inspect_data"):
        run_pipeline(cfg)
    result = json.loads((cfg.out / "manifest.json").read_text())
    assert result["status"] == "failed"
    assert result["sources"]["input"]["sha256"] is None


def test_fabricated_citations_rejected(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    wf = Workflow(cfg)
    for name in ["inspect_data", "preprocess_data", "find_markers", "annotate"]:
        assert wf.call(name, {})["ok"]
    cluster = str(wf.marker_table.cluster.iloc[0])
    genes = (
        wf.marker_table.loc[wf.marker_table.cluster.astype(str) == cluster, "gene"].head(2).tolist()
    )
    result = wf.call(
        "propose_annotation",
        {
            "cluster": cluster,
            "label": "T cell",
            "genes": genes,
            "evidence_ids": ["made-up"],
            "rationale": "test",
        },
    )
    assert not result["ok"] and "retrieved evidence" in result["error"]["message"]
