import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast

import httpx

from .trace import write_json


class EvidenceStore:
    """Fixed endpoint retrieval, cached records, no model-authored URLs accepted."""

    def __init__(self, directory: Path, online: bool, client: Any = None):
        self.directory = directory
        directory.mkdir(exist_ok=True)
        self.online = online
        self.client = client or httpx.Client(timeout=15, follow_redirects=False)
        self.records: dict[str, dict[str, Any]] = {}

    def lookup_gene(self, symbol: str, species: str) -> dict[str, Any]:
        if not isinstance(symbol, str) or not re.fullmatch(r"[A-Za-z0-9_.-]{1,80}", symbol):
            raise ValueError("Use a single gene symbol, not a query or URL")
        organism = {"human": "homo_sapiens", "mouse": "mus_musculus"}[species]
        url = f"https://rest.ensembl.org/lookup/symbol/{organism}/{symbol}"
        key = hashlib.sha256(url.encode()).hexdigest()[:20]
        path = self.directory / f"gene-{key}.json"
        if path.exists():
            return cast(dict[str, Any], json.loads(path.read_text()))
        if not self.online:
            return {"ok": False, "error": "offline"}
        for attempt in range(3):
            try:
                response = self.client.get(url, headers={"Accept": "application/json"})
                if response.status_code in (400, 404):
                    result = {
                        "ok": True,
                        "status": "not_found",
                        "symbol": symbol,
                        "species": species,
                    }
                    break
                response.raise_for_status()
                body = response.json()
                result = {
                    "ok": True,
                    "status": "found",
                    "symbol": symbol,
                    "species": species,
                    "ensembl_id": body.get("id"),
                    "description": body.get("description"),
                    "biotype": body.get("biotype"),
                    "url": url,
                    "retrieved_at": datetime.now(timezone.utc).isoformat(),
                    "scope": "gene identity/function only; not cell-type evidence",
                }
                break
            except httpx.HTTPError:
                if attempt == 2:
                    return {"ok": False, "error": "retrieval_failed"}
                time.sleep(0.25 * 2**attempt)
        write_json(path, result)
        return result

    def search(self, query: str) -> dict[str, Any]:
        if not isinstance(query, str) or not 3 <= len(query) <= 500:
            raise ValueError("query must contain 3–500 characters")
        key = hashlib.sha256(query.encode()).hexdigest()[:20]
        path = self.directory / f"{key}.json"
        if path.exists():
            result = json.loads(path.read_text())
        elif not self.online:
            return {"ok": False, "error": "offline", "records": []}
        else:
            response = None
            for attempt in range(3):
                try:
                    response = self.client.get(
                        "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                        params={
                            "query": query,
                            "format": "json",
                            "resultType": "core",
                            "pageSize": 3,
                        },
                    )
                    if response.status_code == 429 or response.status_code >= 500:
                        response.raise_for_status()
                    response.raise_for_status()
                    break
                except httpx.HTTPError:
                    if attempt == 2:
                        return {"ok": False, "error": "retrieval_failed", "records": []}
                    time.sleep(0.25 * 2**attempt)
            assert response is not None
            entries = []
            for item in response.json().get("resultList", {}).get("result", []):
                source, article_id = str(item.get("source", "")), str(item.get("id", ""))
                if not source.isalnum() or not article_id.replace("-", "").isalnum():
                    continue
                record_id = f"{key}:{source}:{article_id}"
                entries.append(
                    {
                        "id": record_id,
                        "source": "Europe PMC",
                        "title": item.get("title", ""),
                        "abstract": item.get("abstractText", "")[:6000],
                        "year": item.get("pubYear"),
                        "doi": item.get("doi"),
                        "url": f"https://europepmc.org/article/{source}/{article_id}",
                        "query": query,
                        "retrieved_at": datetime.now(timezone.utc).isoformat(),
                    }
                )
            result = {"ok": True, "records": entries}
            write_json(path, result)
        for record in result["records"]:
            self.records[record["id"]] = record
        return cast(dict[str, Any], result)
