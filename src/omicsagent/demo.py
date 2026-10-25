from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .annotation import HUMAN_IMMUNE
from .trace import write_json


def make_demo(directory: Path) -> None:
    """Synthetic, explicitly labeled integration fixtures, NOT biological benchmarks."""
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("Demo directory must be empty")
    directory.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    panel = {k: HUMAN_IMMUNE[k] for k in ["T cell", "B cell", "Myeloid"]}
    genes = (
        [g for values in panel.values() for g in values]
        + [f"GENE{i}" for i in range(70)]
        + ["MT-CO1"]
    )
    labels = np.repeat(list(panel), 60)
    x = rng.poisson(0.3, size=(len(labels), len(genes)))
    for i, label in enumerate(labels):
        idx = [genes.index(g) for g in panel[label]]
        x[i, idx] += rng.poisson(8, size=len(idx))
    obs = pd.DataFrame(index=[f"cell-{i}" for i in range(len(labels))])
    obs["known_label"] = labels
    obs["sample"] = np.tile(np.repeat(["sample_a", "sample_b", "sample_c"], 20), 3)
    obs["cell_area"] = rng.uniform(50, 150, len(labels))
    obs["center_group"] = np.tile(np.array(["A", "B", "", ""]), len(labels) // 4)
    matrix = sparse.csr_matrix(x)
    obj = ad.AnnData(matrix, obs=obs, var=pd.DataFrame(index=genes))
    obj.layers["counts"] = matrix.copy()
    offsets = {"T cell": (0, 0), "B cell": (35, 0), "Myeloid": (15, 35)}
    obj.obsm["spatial"] = np.array([rng.normal(offsets[label], 15) for label in labels])
    obj.uns["omicsagent_demo"] = "synthetic test data, not biological evidence"
    obj.write_h5ad(directory / "synthetic.h5ad")
    ref = obj.copy()
    reference_x = rng.poisson(0.3, size=x.shape)
    for i, label in enumerate(labels):
        idx = [genes.index(g) for g in panel[label]]
        reference_x[i, idx] += rng.poisson(8, size=len(idx))
    ref.X = sparse.csr_matrix(reference_x)
    ref.layers["counts"] = ref.X.copy()
    ref.obs_names = [f"reference-{i}" for i in range(len(labels))]
    ref.write_h5ad(directory / "reference.h5ad")
    write_json(directory / "markers.json", panel)
    common = {
        "input": "synthetic.h5ad",
        "tissue": "synthetic",
        "min_genes": 5,
        "min_cells": 2,
        "max_pct_mt": 100,
        "n_pcs": 10,
        "n_hvg": 50,
        "marker_file": "markers.json",
        "n_permutations": 19,
        "spatial_genes": 3,
    }
    write_json(directory / "scrna.json", {**common, "out": "scrna-output"})
    write_json(
        directory / "spatial.json",
        {
            **common,
            "out": "spatial-output",
            "modality": "spatial",
            "platform": "imaging",
            "sample_key": "sample",
            "coordinate_units": "micron",
            "radius": 20,
            "annotation_key": "known_label",
            "center_group_key": "center_group",
            "center_groups": ["A", "B"],
        },
    )
    write_json(
        directory / "visium.json",
        {
            **common,
            "out": "visium-output",
            "modality": "spatial",
            "platform": "visium",
            "sample_key": "sample",
            "coordinate_units": "pixel",
            "radius": 20,
        },
    )

    for method in ("tacco", "scanvi"):
        write_json(directory / f"{method}-split.json", {
            **common, "out": f"{method}-split-output", "modality": "spatial",
            "platform": "imaging", "sample_key": "sample", "coordinate_units": "micron",
            "radius": 20, "reference": "reference.h5ad", "reference_label_key": "known_label",
            "annotation_method": method, "split_enabled": True,
            "scanvi_epochs": 2, "scanvi_classifier_epochs": 2,
        })
    write_json(directory / "scanvi-scrna.json", {
        **common, "out": "scanvi-scrna-output", "reference": "reference.h5ad",
        "reference_label_key": "known_label", "annotation_method": "scanvi",
        "scanvi_epochs": 2, "scanvi_classifier_epochs": 2,
    })


def make_hierarchy_demo(directory: Path) -> None:
    """Independent synthetic atlas with deliberately misplaced initial L1 cells."""
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("Demo directory must be empty")
    directory.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(123)
    l1 = {'Myeloid': ['LYZ', 'LST1', 'TYROBP'], 'Lymphoid': ['CD3D', 'CD3E', 'TRAC']}
    l2 = {'Mono': ['S100A8', 'S100A9'], 'DC': ['FCER1A', 'CD1C'],
          'T_naive': ['CCR7', 'LEF1'], 'T_effector': ['GZMK', 'CCL5']}
    parent = {'Mono': 'Myeloid', 'DC': 'Myeloid', 'T_naive': 'Lymphoid', 'T_effector': 'Lymphoid'}
    genes = [g for values in [*l1.values(), *l2.values()] for g in values]
    genes += [f'BACKGROUND_{i}' for i in range(100)] + ['IRF7', 'BST2']
    types = np.repeat(list(l2), 60)
    def generate(prefix: str) -> ad.AnnData:
        x = rng.poisson(.5, (len(types), len(genes)))
        for i, subtype in enumerate(types):
            for gene in l1[parent[subtype]] + l2[subtype]:
                x[i, genes.index(gene)] += rng.poisson(12)
        obs = pd.DataFrame({'truth_l2': types, 'truth_l1': [parent[t] for t in types],
                            'sample': np.tile(['sample_a', 'sample_b'], len(types)//2)},
                           index=[f'{prefix}-{i}' for i in range(len(types))])
        obj = ad.AnnData(sparse.csr_matrix(x), obs=obs, var=pd.DataFrame(index=genes))
        obj.layers['counts'] = obj.X.copy()
        obj.obsm['spatial'] = rng.uniform(0, 100, (len(types), 2))
        return obj
    query, reference = generate('query'), generate('reference')
    query.obs['initial_l1'] = query.obs.truth_l1.copy()
    query.obs.loc[query.obs_names[:30], 'initial_l1'] = 'Lymphoid'
    query.obs.loc[query.obs_names[120:150], 'initial_l1'] = 'Myeloid'
    query.write_h5ad(directory / 'query.h5ad')
    reference.write_h5ad(directory / 'reference.h5ad')
    write_json(directory / 'markers.json', l1)
    hierarchy = {'initial_l1_key': 'initial_l1',
                 'compartments': {'myeloid_pool': ['Myeloid'], 'lymphoid_pool': ['Lymphoid']},
                 'l1_markers': l1, 'l2_markers': l2, 'l2_parent': parent,
                 'reference': 'reference.h5ad', 'reference_l1_key': 'truth_l1',
                 'reference_l2_key': 'truth_l2', 'resolution': .6}
    write_json(directory / 'hierarchy.json', hierarchy)
    base = {'input': 'query.h5ad', 'tissue': 'synthetic', 'min_genes': 5, 'min_cells': 2,
            'max_pct_mt': 100, 'n_hvg': 80, 'n_pcs': 10, 'annotation_key': 'initial_l1',
            'marker_file': 'markers.json', 'hierarchy': 'hierarchy.json',
            'n_permutations': 19, 'spatial_genes': 3}
    write_json(directory / 'scrna.json', {**base, 'out': 'scrna-output'})
    write_json(directory / 'imaging.json', {**base, 'out': 'imaging-output', 'modality': 'spatial',
               'platform': 'imaging', 'sample_key': 'sample', 'coordinate_units': 'micron', 'radius': 15})
    write_json(directory / 'hierarchy-domain.json', {**hierarchy, 'mode': 'domain'})
    write_json(directory / 'visium.json', {**base, 'out': 'visium-output', 'modality': 'spatial',
               'platform': 'visium', 'sample_key': 'sample', 'coordinate_units': 'pixel', 'radius': 15,
               'hierarchy': 'hierarchy-domain.json'})
