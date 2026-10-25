import anndata as ad
import pandas as pd
import pytest
from scipy import sparse

from omicsagent.ingest import import_10x


def test_metadata_aligned_by_barcode_not_row(tmp_path, monkeypatch):
    obj = ad.AnnData(
        sparse.csr_matrix([[1, 2], [3, 4]]),
        obs=pd.DataFrame(index=["a", "b"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    monkeypatch.setattr("omicsagent.ingest.sc.read_10x_h5", lambda *a, **k: obj.copy())
    meta = tmp_path / "cells.csv"
    pd.DataFrame({"cell_id": ["b", "a"], "x": [20, 10], "y": [2, 1]}).to_csv(meta, index=False)
    out = tmp_path / "imported.h5ad"
    import_10x(tmp_path / "fake.h5", meta, out, "cell_id", "x", "y")
    loaded = ad.read_h5ad(out)
    assert loaded.obsm["spatial"].tolist() == [[10, 1], [20, 2]]
    with pytest.raises(ValueError, match="overwrite"):
        import_10x(tmp_path / "fake.h5", meta, out, "cell_id", "x", "y")


def test_missing_metadata_barcode_rejected(tmp_path, monkeypatch):
    obj = ad.AnnData(sparse.csr_matrix([[1, 2]]), obs=pd.DataFrame(index=["a"]))
    monkeypatch.setattr("omicsagent.ingest.sc.read_10x_h5", lambda *a, **k: obj.copy())
    meta = tmp_path / "cells.csv"
    pd.DataFrame({"cell_id": ["b"]}).to_csv(meta, index=False)
    with pytest.raises(ValueError, match="missing 1"):
        import_10x(tmp_path / "fake.h5", meta, tmp_path / "out.h5ad", "cell_id", None, None)
