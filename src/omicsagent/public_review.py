"""Reproducible, donor-disjoint Tabula Sapiens bone marrow case preparation."""
import hashlib
import json
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .analysis import counts
from .trace import sha256, write_json

SOURCE = 'https://ndownloader.figshare.com/files/34701967'
SOURCE_MD5 = '1b665025142b0c74f7b63d7310b379d4'
SOURCE_SHA256 = '4e08b45d5668c7cc741ccddafd2626f65303ce764be4f3ebc1b7778d60c30c64'
L1_MARKERS = {'T cell': ['CD3D', 'CD3E', 'TRAC'],
              'B cell': ['CD79A', 'CD79B', 'MS4A1'],
              'Myeloid': ['LYZ', 'LST1', 'TYROBP'],
              'NK cell': ['NKG7', 'GNLY', 'KLRD1']}
L2_MARKERS = {'CD4 T': ['IL7R', 'CCR7', 'LTB'],
              'CD8 T': ['CD8A', 'CD8B', 'CCL5'],
              'Naive B': ['TCL1A', 'IGHD', 'IL4R'],
              'Memory B': ['CD27', 'AIM2', 'TNFRSF13B'],
              'Plasma': ['MZB1', 'JCHAIN', 'SDC1'],
              'Monocyte': ['FCN1', 'CTSS', 'LYZ'],
              'Macrophage': ['C1QA', 'C1QB', 'APOE'],
              'Neutrophil': ['FCGR3B', 'CSF3R', 'S100A8'],
              'NK': ['NKG7', 'GNLY', 'KLRD1']}
L2_PARENT = {name: ('T cell' if name in {'CD4 T', 'CD8 T'} else
                    'B cell' if name in {'Naive B', 'Memory B', 'Plasma'} else
                    'NK cell' if name == 'NK' else 'Myeloid') for name in L2_MARKERS}


def _map_label(ontology: str, free: str) -> tuple[str, str] | None:
    label = ontology.lower()
    detail = free.lower()
    if 'cd4' in label:
        return 'T cell', 'CD4 T'
    if 'cd8' in label:
        return 'T cell', 'CD8 T'
    if 'plasma' in label:
        return 'B cell', 'Plasma'
    if 'b cell' in label:
        if 'naive' in label + detail:
            return 'B cell', 'Naive B'
        if 'memory' in label + detail:
            return 'B cell', 'Memory B'
    if 'natural killer' in label or label == 'nk cell':
        return 'NK cell', 'NK'
    if 'macrophage' in label:
        return 'Myeloid', 'Macrophage'
    if 'monocyte' in label:
        return 'Myeloid', 'Monocyte'
    if 'neutrophil' in label or 'granulocyte' in label:
        return 'Myeloid', 'Neutrophil'
    return None


def partition_cells(obs: pd.DataFrame, cap: int = 150) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fixed donors; stable-ID capped per published fine label, not inferred outcome."""
    required = {'donor', 'method', 'cell_ontology_class', 'free_annotation'}
    if not required <= set(obs) or obs[list(required)].isna().any().any():
        raise ValueError('Required donor/method/annotation metadata is missing')
    frame = obs.loc[obs.method.astype(str).eq('10X')].copy().sort_index()
    pairs = [_map_label(str(a), str(b)) for a, b in
             zip(frame.cell_ontology_class, frame.free_annotation, strict=True)]
    keep = np.array([pair is not None for pair in pairs])
    frame = frame.loc[keep].copy()
    mapped = [pair for pair in pairs if pair is not None]
    frame['initial_l1'] = [pair[0] for pair in mapped]
    frame['reference_l1'] = frame.initial_l1.astype(str)
    frame['reference_l2'] = [pair[1] for pair in mapped]
    frame['published_label'] = frame.cell_ontology_class.astype(str)
    frame = frame.groupby(['donor', 'cell_ontology_class'], observed=True, sort=True).head(cap)
    ref = frame[frame.donor.astype(str) == 'TSP2'].copy()
    query = frame[frame.donor.astype(str).isin(['TSP11', 'TSP14'])].copy()
    if ref.empty or query.empty or set(ref.donor.astype(str)) & set(query.donor.astype(str)):
        raise ValueError('Required fixed donor partition unavailable')
    return ref, query


def raw_source(obj: ad.AnnData) -> ad.AnnData:
    """Prefer the published raw slot, validating integer-like counts explicitly."""
    source = obj.raw.to_adata() if obj.raw is not None else obj.copy()
    matrix = counts(source, 'X')
    source.X = matrix
    source.layers['counts'] = matrix.copy()
    return source


def fetch_source(directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    archive = directory/'TS_Bone_Marrow.h5ad.zip'
    if not archive.exists():
        request = urllib.request.Request(SOURCE, headers={'User-Agent': 'OmicsAgent/0.3'})
        with urllib.request.urlopen(request, timeout=120) as response, archive.open('wb') as handle:
            while chunk := response.read(1024 * 1024):
                handle.write(chunk)
    digest = hashlib.md5(usedforsecurity=False)
    with archive.open('rb') as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    if digest.hexdigest() != SOURCE_MD5:
        raise ValueError('Published archive checksum mismatch; remove failed download before retry')
    destination = directory/'TS_Bone_Marrow.h5ad'
    if not destination.exists():
        with zipfile.ZipFile(archive) as source:
            matches = [info for info in source.infolist() if info.filename.endswith('TS_Bone_Marrow.h5ad')
                       and not info.filename.startswith('__MACOSX/')]
            if len(matches) != 1:
                raise ValueError('Expected exactly one bone marrow AnnData file')
            with source.open(matches[0]) as reader, destination.open('wb') as writer:
                while chunk := reader.read(1024 * 1024):
                    writer.write(chunk)
    return destination


def prepare_public_review(source: Path, out: Path, cap: int = 150) -> dict[str, Any]:
    if out.exists() and any(out.iterdir()):
        raise ValueError('Preparation output must be empty')
    if cap < 5:
        raise ValueError('Cell cap must be at least 5')
    if sha256(source) != SOURCE_SHA256:
        raise ValueError('Input does not match the pinned public dataset SHA256')
    data = ad.read_h5ad(source, backed='r')
    try:
        ref_obs, query_obs = partition_cells(data.obs, cap)
        ids = ref_obs.index.union(query_obs.index)
        if data.raw is None:
            raise ValueError('Pinned source is missing its verified raw count slot')
        # Read only selected raw rows; to_memory() would load unrelated corrected layers.
        matrix = data.raw.X[data.obs_names.get_indexer(ids)].tocsr()
        selected = ad.AnnData(matrix, obs=data.obs.loc[ids].copy(), var=data.raw.var.copy())
    finally:
        data.file.close()
    raw = raw_source(selected)
    raw.obs = pd.concat([ref_obs, query_obs]).loc[raw.obs_names].copy()
    # This pinned file's verified raw.var index already contains unique symbols.
    # Ancillary gene_symbol contains duplicated aliases and must not replace it.
    if not raw.var_names.is_unique:
        raise ValueError('Duplicate gene symbols require explicit aggregation, not suffixing')
    genes = np.asarray((sparse.csr_matrix(raw.X) > 0).sum(axis=0)).ravel() >= 3
    raw = raw[:, genes].copy()
    provenance = {'dataset': 'Tabula Sapiens v1 bone marrow, Figshare release v5',
                  'source_url': 'https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219/5',
                  'doi': '10.6084/m9.figshare.14267219.v5', 'download_url': SOURCE,
                  'archive_md5': SOURCE_MD5, 'source_h5ad_sha256': sha256(source),
                  'license': 'CC BY 4.0', 'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                  'citation': 'The Tabula Sapiens Consortium, Science 376, eabl4896 (2022). '
                              'doi:10.1126/science.abl4896',
                  'data_kind': 'real', 'donor_key': 'donor', 'reference_donors': ['TSP2'],
                  'query_donors': sorted(query_obs.donor.astype(str).unique()),
                  'sampling': f'10X only; stable cell-ID order; cap {cap} per donor/published fine label',
                  'excluded': 'Unmapped published classes excluded; no artificial label corruption',
                  'reference_cells': len(ref_obs), 'query_cells': len(query_obs),
                  'n_genes': raw.n_vars,
                  'limits': 'Published labels are reference metadata, not independent ground truth. '
                            'Panels and thresholds are illustrative; no biological accuracy claim. '
                            'One NK module cannot establish a subtype margin.'}
    raw.uns['oa_review_provenance_json'] = json.dumps(provenance)
    out.mkdir(parents=True, exist_ok=True)
    raw[ref_obs.index].copy().write_h5ad(out/'reference.h5ad', compression='gzip')
    raw[query_obs.index].copy().write_h5ad(out/'query.h5ad', compression='gzip')
    write_json(out/'provenance.json', provenance)
    hierarchy = {'initial_l1_key': 'initial_l1', 'reference': 'reference.h5ad',
                 'reference_l1_key': 'reference_l1', 'reference_l2_key': 'reference_l2',
                 'compartments': {'immune': list(L1_MARKERS)},
                 'l1_markers': L1_MARKERS, 'l2_markers': L2_MARKERS, 'l2_parent': L2_PARENT,
                 'max_iterations': 3, 'n_hvg': 2000, 'n_pcs': 15, 'resolution': .5,
                 'max_dense_elements': 120_000_000}
    write_json(out/'hierarchy.json', hierarchy)
    write_json(out/'run.json', {'input': 'query.h5ad', 'out': 'run', 'tissue': 'bone marrow',
                              'counts_layer': 'counts', 'hierarchy': 'hierarchy.json',
                              'annotation_key': 'initial_l1', 'n_hvg': 2000,
                              'max_dense_elements': 120_000_000})
    return provenance
