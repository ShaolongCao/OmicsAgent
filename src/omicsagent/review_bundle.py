"""Export a portable evidence snapshot without shipping count matrices or local paths."""
import json
import shutil
from pathlib import Path, PureWindowsPath
from typing import Any

from .review import build_review, write_review
from .trace import sha256, write_json


def _portable(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _portable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_portable(v) for v in value]
    if isinstance(value, str):
        if PureWindowsPath(value).is_absolute():
            return PureWindowsPath(value).name
        if Path(value).is_absolute():
            return Path(value).name
    return value


def export_review_bundle(source: Path, destination: Path) -> dict[str, Any]:
    """Keep numerical CSVs byte-for-byte, redact path-only JSON fields, record lineage."""
    if destination.exists():
        raise ValueError('Bundle destination must not exist')
    original = build_review(source)
    destination.mkdir(parents=True)
    files = [source/'config.json', source/'cells.csv', source/'trace.jsonl',
             *sorted((source/'hierarchy').glob('*.csv')),
             *sorted((source/'hierarchy').glob('*.json'))]
    redacted = []
    for file in files:
        if not file.is_file():
            continue
        name = file.relative_to(source).as_posix()
        target = destination/name
        target.parent.mkdir(parents=True, exist_ok=True)
        if file.suffix in {'.json', '.jsonl'}:
            lines = file.read_text(encoding='utf-8').splitlines()
            if file.suffix == '.json':
                value = json.loads('\n'.join(lines))
                cleaned = _portable(value)
                write_json(target, cleaned)
            else:
                value = [json.loads(line) for line in lines if line.strip()]
                cleaned = _portable(value)
                target.write_text(''.join(json.dumps(row, allow_nan=False)+'\n' for row in cleaned), encoding='utf-8')
            if value != cleaned:
                redacted.append(name)
        else:
            shutil.copyfile(file, target)
    provenance = dict(original['provenance'])
    provenance.update({'source_run_hash': original['run_hash'],
                       'bundle_note': 'Numerical CSVs are unchanged. Absolute path fields in JSON/JSONL '
                                      'are replaced by basenames; original hashes are retained in manifest. '
                                      'Count matrices are retrieved separately with prepare-review.'})
    write_json(destination/'provenance.json', provenance)
    original_manifest = (json.loads((source/'manifest.json').read_text(encoding='utf-8'))
                         if (source/'manifest.json').is_file() else {})
    write_json(destination/'manifest.json', {
        'status': 'derived_review_bundle', 'source_run_hash': original['run_hash'],
        'source_run_id': original_manifest.get('run_id'),
        'original_artifact_hashes': original_manifest.get('artifacts', {}),
        'redacted_files': sorted(redacted),
        'artifacts': {p.relative_to(destination).as_posix(): sha256(p)
                      for p in destination.rglob('*') if p.is_file() and p.name != 'manifest.json'}})
    return write_review(destination)
