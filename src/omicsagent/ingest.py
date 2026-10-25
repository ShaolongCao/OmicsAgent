from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def import_10x(
    matrix: Path,
    metadata: Path | None,
    out: Path,
    id_column: str,
    x_column: str | None,
    y_column: str | None,
) -> dict[str, object]:
    """Join 10x gene-expression H5 to an explicitly specified metadata CSV by ID."""
    if out.exists():
        raise ValueError("Refusing to overwrite existing h5ad")
    if bool(x_column) != bool(y_column):
        raise ValueError("Specify both x and y columns")
    obj = sc.read_10x_h5(matrix, gex_only=True)
    obj.var["original_symbol"] = obj.var_names.astype(str)
    duplicates = int(obj.var_names.duplicated().sum())
    # Retain mapping: repeated symbols get suffixes, never silently sum ambiguous genes.
    obj.var_names_make_unique()
    if metadata:
        frame = pd.read_csv(metadata, dtype={id_column: str})
        if (
            id_column not in frame
            or frame[id_column].isna().any()
            or frame[id_column].duplicated().any()
        ):
            raise ValueError("Metadata requires unique, nonmissing IDs")
        frame = frame.set_index(id_column)
        missing = obj.obs_names.difference(frame.index)
        if len(missing):
            raise ValueError(f"Metadata missing {len(missing)} matrix barcodes")
        obj.obs = frame.loc[obj.obs_names].copy()
        if x_column:
            xy = obj.obs[[x_column, y_column]].to_numpy(dtype=float)
            if not np.isfinite(xy).all():
                raise ValueError("Coordinates must be finite")
            obj.obsm["spatial"] = xy
    elif x_column:
        raise ValueError("Coordinates require a metadata CSV")
    obj.layers["counts"] = obj.X.copy()
    out.parent.mkdir(parents=True, exist_ok=True)
    obj.write_h5ad(out)
    return {
        "status": "imported",
        "shape": list(obj.shape),
        "duplicate_symbols_suffixed": duplicates,
        "out": str(out),
    }
