import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


class Trace:
    def __init__(self, path: Path, run_id: str) -> None:
        self.path = path
        self.run_id = run_id
        self.sequence = 0

    def emit(self, event: str, **details: Any) -> None:
        self.sequence += 1
        record = dict(
            schema_version=1,
            run_id=self.run_id,
            event_id=self.sequence,
            timestamp=datetime.now(timezone.utc).isoformat(),
            event=event,
            **details,
        )
        with self.path.open("a") as handle:
            handle.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
