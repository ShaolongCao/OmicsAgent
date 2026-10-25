"""Download a public 10x Visium example. No histology image is extracted or required.

Source: https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-1-0
10x Genomics, Human Lymph Node (Space Ranger 1.1.0), CC BY 4.0.
Run from repository root: python examples/fetch_visium.py
"""

import tarfile
from pathlib import Path

import httpx
import pandas as pd
from scipy.spatial import cKDTree

from omicsagent.ingest import import_10x
from omicsagent.trace import sha256, write_json

DIRECTORY = Path("data/visium-lymph-node")
BASE = "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/"
PREFIX = "V1_Human_Lymph_Node"


def main() -> None:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    sources = []
    for suffix in ("filtered_feature_bc_matrix.h5", "spatial.tar.gz"):
        name = f"{PREFIX}_{suffix}"
        path = DIRECTORY / name
        if not path.exists():
            partial = path.with_suffix(path.suffix + ".partial")
            with httpx.stream("GET", BASE + name, timeout=60, follow_redirects=True) as response:
                response.raise_for_status()
                with partial.open("wb") as handle:
                    for data in response.iter_bytes():
                        handle.write(data)
            partial.replace(path)
        sources.append({"url": BASE + name, "sha256": sha256(path)})
    with tarfile.open(DIRECTORY / f"{PREFIX}_spatial.tar.gz") as archive:
        members = [
            m
            for m in archive.getmembers()
            if m.name.endswith("tissue_positions_list.csv") and m.isfile()
        ]
        if len(members) != 1:
            raise ValueError("Expected exactly one tissue positions table")
        stream = archive.extractfile(members[0])
        frame = pd.read_csv(
            stream,
            header=None,
            names=["cell_id", "in_tissue", "array_row", "array_col", "y_pixel", "x_pixel"],
        )
    frame["sample_id"] = PREFIX
    frame.to_csv(DIRECTORY / "positions.csv", index=False)
    target = DIRECTORY / "visium.h5ad"
    if not target.exists():
        print(
            import_10x(
                DIRECTORY / f"{PREFIX}_filtered_feature_bc_matrix.h5",
                DIRECTORY / "positions.csv",
                target,
                "cell_id",
                "x_pixel",
                "y_pixel",
            )
        )
    # One ring: use measured nearest-neighbor spacing in full-resolution pixel coordinates.
    tissue = frame.loc[frame.in_tissue == 1, ["x_pixel", "y_pixel"]].to_numpy()
    distances, _ = cKDTree(tissue).query(tissue, k=2)
    radius = float(pd.Series(distances[:, 1]).median() * 1.15)
    write_json(
        DIRECTORY / "provenance.json",
        {
            "source": "10x Genomics Human Lymph Node, Space Ranger 1.1.0",
            "license": "CC BY 4.0",
            "sources": sources,
            "coordinate_units": "pixel",
        },
    )
    write_json(
        DIRECTORY / "run.json",
        {
            "input": "visium.h5ad",
            "out": "../../outputs/visium-lymph-node",
            "modality": "spatial",
            "platform": "visium",
            "species": "human",
            "tissue": "lymph node",
            "min_genes": 200,
            "max_pct_mt": 20,
            "sample_key": "sample_id",
            "coordinate_units": "pixel",
            "radius": radius,
            "n_permutations": 99,
            "spatial_genes": 10,
        },
    )
    print(f"Prepared {target}; radius={radius:.2f} pixels")


if __name__ == "__main__":
    main()
