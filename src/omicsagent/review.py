"""Derive inspectable review cases from recorded measurements, without refitting."""
import hashlib
import json
from pathlib import Path
from typing import Any, cast

import pandas as pd

from .trace import sha256, write_json

KINDS = ('agreement', 'conflict', 'insufficient')


def evaluate_predictions(expected: list[str], predicted: list[str]) -> dict[str, Any]:
    """Behavioral scoring: blanket abstention cannot satisfy resolution controls."""
    allowed = {'resolved', 'deferred', 'unknown', 'reassigned'}
    if not expected or len(expected) != len(predicted) or not set(expected + predicted) <= allowed:
        raise ValueError('Nonempty, aligned outcomes from the behavior vocabulary are required')
    resolve = [i for i, value in enumerate(expected) if value in {'resolved', 'reassigned'}]
    defer = [i for i, value in enumerate(expected) if value in {'deferred', 'unknown'}]
    coverage = sum(predicted[i] == expected[i] for i in resolve) / max(len(resolve), 1)
    inappropriate = sum(predicted[i] in {'resolved', 'reassigned'} for i in defer) / max(len(defer), 1)
    return {'passed': expected == predicted and bool(resolve) and bool(defer),
            'n_cases': len(expected), 'resolved_coverage': coverage,
            'inappropriate_assignment': inappropriate,
            'correct': sum(a == b for a, b in zip(expected, predicted, strict=True))}


def evaluate_behavior() -> dict[str, Any]:
    """Small explicit policy controls, separate from real-data biological validation."""
    from .hierarchy import decide_cluster
    controls = [
        ('supported_agreement', ('T', 'T', True, 'T', .95, .9), 'resolved'),
        ('coarse_conflict', ('T', 'T', True, 'B', .95, .9), 'deferred'),
        ('atlas_only', ('T', None, False, 'B', .99, .99), 'unknown'),
        ('supported_reassignment', ('T', 'B', True, 'B', .95, .9), 'reassigned'),
        ('weak_markers', ('T', 'B', False, 'B', .99, .99), 'deferred'),
    ]
    vocabulary = {'Keep': 'resolved', 'Ambiguous': 'deferred',
                  'Unknown': 'unknown', 'Reassign': 'reassigned'}
    expected = [c[2] for c in controls]
    predicted = [vocabulary[decide_cluster(*c[1])[0]] for c in controls]
    return {'scope': 'Deterministic mechanism controls; not biological accuracy',
            **evaluate_predictions(expected, predicted),
            'controls': [{'id': c[0], 'expected': c[2], 'observed': p}
                         for c, p in zip(controls, predicted, strict=True)],
            'all_unknown_control': evaluate_predictions(expected, ['unknown'] * len(controls))}


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {}


def _bool(value: Any) -> bool:
    return value is True or str(value).lower() == 'true'


def _text(value: Any, default: str = 'Unknown') -> str:
    return default if value is None or pd.isna(value) else str(value)


def _donor_evidence(out: Path, row: dict[str, Any], members: pd.DataFrame,
                    donor_key: str | None, metadata: pd.DataFrame,
                    thresholds: dict[str, Any]) -> list[dict[str, Any]]:
    """Per-cell marker/reference disagreement by verified donor, with denominators."""
    if not donor_key or donor_key not in metadata or not metadata.index.is_unique:
        return []
    scores_path = out/'hierarchy'/'l2_scores.csv'
    atlas_path = out/'hierarchy'/'atlas_predictions.csv'
    if not scores_path.is_file() or not atlas_path.is_file():
        return []
    scores = pd.read_csv(scores_path, index_col=0)
    atlas = pd.read_csv(atlas_path, index_col=0)
    parent = thresholds.get('l2_parent', {})
    labels = [k for k in scores if parent.get(k) == _text(row.get('target_l1'))]
    if len(labels) < 2 or 'l2' not in atlas:
        return []
    ids = members.cell_id.astype(str)
    values = scores.reindex(ids)[labels]
    # No winner is invented when competing modules are missing or tied.
    valid = values.notna().all(axis=1) & values.eq(values.max(axis=1), axis=0).sum(axis=1).eq(1)
    donor = metadata.reindex(ids)[donor_key]
    predicted = atlas.reindex(ids).l2
    valid &= donor.notna() & predicted.isin(labels)
    if not valid.any():
        return []
    best = values.loc[valid].idxmax(axis=1)
    data = pd.DataFrame({'donor': donor.loc[valid].astype(str),
                         'disagree': best != predicted.loc[valid]})
    return [{'donor': str(name), 'n_evaluable': len(group),
             'n_cluster_cells': int((donor.astype(str) == name).sum()),
             'n_disagree': int(group.disagree.sum()),
             'disagreement_fraction': float(group.disagree.mean())}
            for name, group in data.groupby('donor', sort=True)]


def build_review(out: Path, donor_key: str | None = None) -> dict[str, Any]:
    """Select one representative per observed category from the final recorded round.

    Selection policy is fixed: eligible final-round rows, largest cell count first,
    then cluster/source label. Every eligible count is disclosed, including absent types.
    """
    out = Path(out)
    hierarchy = out/'hierarchy'
    decisions = pd.read_csv(hierarchy/'decisions.csv')
    thresholds = _read(hierarchy/'config.json')
    provenance = _read(out/'provenance.json')
    if donor_key is None:
        donor_key = provenance.get('donor_key')
    paths = [out/'config.json', out/'cells.csv', out/'provenance.json',
             *sorted(hierarchy.glob('*.csv')), *sorted(hierarchy.glob('*.json'))]
    sources = {p.relative_to(out).as_posix(): sha256(p) for p in paths if p.is_file()}
    run_hash = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
    metadata = pd.read_csv(out/'cells.csv', index_col=0) if (out/'cells.csv').is_file() else pd.DataFrame()
    metadata.index = metadata.index.astype(str)
    dot_path = hierarchy/'dotplot_data.csv'
    dots = pd.read_csv(dot_path) if dot_path.is_file() else pd.DataFrame()
    coverage = _read(hierarchy/'marker_coverage.json')
    history_path = hierarchy/'iterations.json'
    history = json.loads(history_path.read_text(encoding='utf-8')) if history_path.is_file() else []
    recorded_rounds = [int(p.stem.split('_')[1]) for p in hierarchy.glob('round_*_cells.csv')]
    recorded_rounds += [int(r['iteration']) for r in history]
    final_iteration = max(recorded_rounds, default=int(decisions.iteration.max()) if not decisions.empty else 0)
    skipped_pools = [pool for entry in history if entry['iteration'] == final_iteration
                     for pool in entry.get('pools', []) if pool.get('reason')]
    final = decisions[decisions.iteration == final_iteration].copy()
    eligible: dict[str, list[dict[str, Any]]] = {kind: [] for kind in KINDS}
    excluded_resolved = 0
    for row in final.to_dict('records'):
        l1_conflict = (_bool(row.get('marker_strong')) and
                       _text(row.get('atlas_l1')) not in {'Unknown', _text(row.get('marker_best'))} and
                       float(row.get('atlas_purity', 0)) >= thresholds.get('atlas_purity', .7) and
                       float(row.get('atlas_vote', 0)) >= thresholds.get('atlas_vote', .6))
        l2_conflict = row.get('l2_reason') == 'marker_reference_conflict'
        if l1_conflict or l2_conflict:
            kind = 'conflict'
        elif (row.get('decision') == 'Keep' and _bool(row.get('marker_strong')) and
              row.get('l2_reason') == 'marker_supported' and
              _text(row.get('selected_l2')) != 'Unknown' and
              _text(row.get('atlas_l2')) == _text(row.get('selected_l2')) and
              float(row.get('atlas_l2_purity', 0)) >= thresholds.get('atlas_purity', .7) and
              float(row.get('atlas_l2_vote', 0)) >= thresholds.get('atlas_vote', .6)):
            kind = 'agreement'
        elif (_text(row.get('selected_l2')) == 'Unknown' or row.get('decision') != 'Keep'
              or not _bool(row.get('marker_strong'))):
            kind = 'insufficient'
        else:
            excluded_resolved += 1
            continue
        eligible[kind].append(row)
    cases = []
    round_cache: dict[int, pd.DataFrame] = {}
    for kind in KINDS:
        if not eligible[kind]:
            continue
        row = sorted(eligible[kind], key=lambda r: (-int(r['n_cells']), str(r['cluster']),
                                                    str(r['source_l1'])))[0]
        iteration = int(row['iteration'])
        if iteration not in round_cache:
            round_cache[iteration] = pd.read_csv(hierarchy/f'round_{iteration}_cells.csv')
        round_cells = round_cache[iteration]
        members = round_cells[(round_cells.cluster == row['cluster']) &
                              (round_cells.before_l1 == row['source_l1'])].copy()
        if len(members) != int(row['n_cells']) or not members.cell_id.is_unique:
            raise ValueError('Decision membership/count mismatch; refusing misleading review')
        l1 = _text(row['target_l1'])
        l2 = _text(row.get('selected_l2'))
        l1_supported = (row['decision'] == 'Keep' and _bool(row.get('marker_strong')))
        support: list[str] = []
        counter: list[str] = []
        if _bool(row.get('marker_strong')):
            destination = support if _text(row.get('marker_best')) == l1 else counter
            destination.append(f"L1 marker winner {_text(row.get('marker_best'))}; "
                               f"margin {float(row['marker_margin']):.3f}; detection gate passed.")
        else:
            counter.append('L1 markers do not pass the configured support and separation gates.')
        reference = (f"Reference L1 {_text(row.get('atlas_l1'))}: "
                     f"purity {float(row.get('atlas_purity', 0)):.3f}, "
                     f"median vote {float(row.get('atlas_vote', 0)):.3f}.")
        (support if _text(row.get('atlas_l1')) == l1 else counter).append(reference)
        if row.get('l2_reason') == 'marker_reference_conflict':
            support.append(f"L2 marker winner {_text(row.get('l2_marker_best'))} passes the marker gate.")
            counter.append(f"Reference L2 {_text(row.get('atlas_l2'))} conflicts with the marker winner; "
                           f"purity {float(row.get('atlas_l2_purity', 0)):.3f}, "
                           f"median vote {float(row.get('atlas_l2_vote', 0)):.3f}.")
        elif l2 != 'Unknown':
            support.append(f"L2 {l2}: {_text(row.get('l2_reason'))}.")
        else:
            counter.append(f"L2 unresolved: {_text(row.get('l2_reason'), 'not recorded in this run')}.")
        subset = dots[(dots.iteration == iteration) & (dots.cluster == row['cluster'])] if not dots.empty else dots
        donor_summary = _donor_evidence(out, row, members, donor_key, metadata, thresholds)
        adequate = [d for d in donor_summary if d['n_evaluable'] >= 5]
        if kind == 'agreement' and donor_summary:
            n_pairs = sum(d['n_evaluable'] for d in donor_summary)
            n_disagree = sum(d['n_disagree'] for d in donor_summary)
            if n_disagree:
                counter.append(f'Per-cell marker/reference winners differ in {n_disagree}/{n_pairs} '
                               'evaluable cells; cluster agreement does not imply uniform cells.')
        next_check = None
        if kind == 'conflict' and len(adequate) >= 2:
            rates = [d['disagreement_fraction'] for d in adequate]
            next_check = {'kind': 'donor',
                          'observation': f"Across {len(adequate)} donors with at least 5 evaluable cells, "
                                         f"marker/reference disagreement ranges {min(rates):.1%} to {max(rates):.1%}.",
                          'reason': 'Inspect donor-stratified competing markers to assess whether the '
                                    'observed disagreement is concentrated or recurs across donors; '
                                    'this is descriptive evidence, not a donor-effect test.'}
            if min(rates) >= .5 and max(rates) - min(rates) < .3:
                next_check['kind'] = 'cross_donor_conflict'
                next_check['reason'] = ('Disagreement is observed in both available donors. Compare '
                                        'the competing subtype markers in this cluster and its mapped '
                                        'reference cells before accepting either subtype. This can '
                                        'expose panel specificity or reference annotation mismatch; '
                                        'the donor counts are descriptive, not a replication test.')
        elif kind == 'conflict':
            next_check = {'kind': 'competing_markers',
                          'observation': 'Marker and reference winners conflict; adequate donor-level '
                                         'paired measurements are unavailable.',
                          'reason': 'Inspect competing marker detection and expression in the saved '
                                    'table before accepting a subtype. A donor explanation is not established.'}
        elif kind == 'insufficient':
            next_check = {'kind': 'coverage',
                          'observation': f"Recorded unresolved reason: {_text(row.get('l2_reason'))}; "
                                         f"L1 reason: {_text(row.get('l1_reason'))}.",
                          'reason': 'Inspect marker availability, detection fractions and pool size '
                                    'against the configured gates before requesting a new label.'}
            subtype_margin = float(row.get('l2_marker_margin', 0))
            min_margin = float(thresholds.get('subtype_margin', .5))
            if row.get('l2_reason') == 'insufficient_marker_support' and subtype_margin < min_margin:
                next_check = {'kind': 'subtype_separation',
                              'observation': f'L2 marker margin {subtype_margin:.3f} is below the '
                                             f'configured minimum {min_margin:.3f}.',
                              'reason': 'Inspect per-cell co-expression of the competing subtype '
                                        'markers to distinguish a mixed cluster from overlapping '
                                        'panels before accepting the reference subtype.'}
        cells = []
        for member in members.to_dict('records'):
            cell_id = str(member['cell_id'])
            original = (_text(metadata.loc[cell_id, 'oa_initial_l1'])
                        if cell_id in metadata.index and 'oa_initial_l1' in metadata else
                        _text(member['before_l1']))
            cells.append({'cell_id': cell_id, 'original_l1': original,
                          'proposed_l1': l1, 'proposed_l2': l2})
        cid = f"{iteration}:{row['cluster']}:{row['source_l1']}"
        originals = pd.Series([cell['original_l1'] for cell in cells]).value_counts()
        original_summary = (str(originals.index[0]) if len(originals) == 1 else
                            'Mixed: ' + ', '.join(f'{label} ({number})' for label, number in originals.items()))
        evidence_n = int((round_cells.cluster == row['cluster']).sum())
        scope_note = (f' Evidence is aggregated over {evidence_n} cluster cells; this review affects '
                      f'{len(members)} cells with the recorded source label.' if evidence_n != len(members) else '')
        links = {'Decision table': 'hierarchy/decisions.csv',
                 'Marker measurements': 'hierarchy/dotplot_data.csv',
                 'Marker coverage': 'hierarchy/marker_coverage.json',
                 'Iteration and stopping record': 'hierarchy/iterations.json',
                 'Parameters': 'hierarchy/config.json',
                 'Cell membership': f'hierarchy/round_{iteration}_cells.csv',
                 'Configuration': 'config.json', 'Manifest': 'manifest.json',
                 'Execution trace': 'trace.jsonl'}
        cases.append({'id': cid, 'kind': kind, 'iteration': iteration,
                      'cluster': str(row['cluster']), 'n_cells': len(members),
                      'evidence_n_cells': evidence_n,
                      'l1': l1, 'l2': l2, 'l1_supported': l1_supported,
                      'original_l1': original_summary, 'proposed_l1': l1,
                      'support': support, 'counterevidence': counter,
                      'unresolved': ('No subtype assigned.' if l2 == 'Unknown' else
                                     'Provisional annotation; independent biological validation remains open.')
                                    + ('' if l1_supported else ' Retained L1 is not validated.') + scope_note,
                      'next_check': next_check, 'decision': row,
                      'markers': subset.to_dict('records'), 'coverage': coverage,
                      'donor_summary': donor_summary, 'cells': cells,
                      'links': {label: path for label, path in links.items()
                                if (out/path).is_file() or path in {'manifest.json', 'trace.jsonl'}},
                      'evidence_ids': [f'{cid}:decision', f'{cid}:markers', f'{cid}:membership']})
    missing = [kind for kind in KINDS if not eligible[kind]]
    # Round-trip replaces NaN from older CSV schemas; no invalid JSON is embedded.
    result = {'schema_version': 1, 'run_hash': run_hash, 'provenance': provenance,
              'selection': {'policy': 'Final-round rows; largest n_cells, then cluster/source label',
                            'eligible': {kind: len(rows) for kind, rows in eligible.items()},
                            'excluded_resolved_without_dual_support': excluded_resolved,
                            'final_iteration': final_iteration, 'skipped_pools': skipped_pools,
                            'missing': missing, 'complete': not missing},
              'cases': cases, 'sources': sources, 'thresholds': thresholds}
    return cast(dict[str, Any], json.loads(json.dumps(result, default=str), parse_constant=lambda _: None))


def write_review(out: Path, donor_key: str | None = None) -> dict[str, Any]:
    from .review_html import render_review
    report = build_review(out, donor_key)
    write_json(Path(out)/'review.json', report)
    (Path(out)/'report.html').write_text(render_review(report), encoding='utf-8')
    manifest_path = Path(out)/'manifest.json'
    if manifest_path.is_file():
        manifest = _read(manifest_path)
        for name in ('report.html', 'review.json', 'provenance.json'):
            if (Path(out)/name).is_file():
                manifest.setdefault('artifacts', {})[name] = sha256(Path(out)/name)
        write_json(manifest_path, manifest)
    return {'report': str(Path(out)/'report.html'), 'selection': report['selection'],
            'run_hash': report['run_hash']}
