from types import SimpleNamespace

import httpx
import pytest

from omicsagent.agent import run_agent
from omicsagent.config import load_config
from omicsagent.evidence import EvidenceStore


def test_evidence_cache_and_source_ids(tmp_path):
    calls = []

    def handler(request):
        calls.append(request)
        return httpx.Response(
            200,
            json={
                "resultList": {
                    "result": [
                        {
                            "source": "MED",
                            "id": "123",
                            "title": "Example",
                            "abstractText": "A supplied test abstract.",
                        }
                    ]
                }
            },
        )

    store = EvidenceStore(tmp_path, True, httpx.Client(transport=httpx.MockTransport(handler)))
    first = store.search("T cell CD3D")
    second = store.search("T cell CD3D")
    assert first == second and len(calls) == 1
    assert first["records"][0]["id"] in store.records
    assert first["records"][0]["url"] == "https://europepmc.org/article/MED/123"


def test_offline_evidence_never_fabricated(tmp_path):
    result = EvidenceStore(tmp_path, False).search("example query")
    assert not result["ok"] and result["records"] == []


def test_ensembl_identity_is_not_cell_type_evidence(tmp_path):
    client = httpx.Client(
        transport=httpx.MockTransport(
            lambda req: httpx.Response(200, json={"id": "ENSG_TEST", "description": "Test gene"})
        )
    )
    store = EvidenceStore(tmp_path, True, client)
    result = store.lookup_gene("CD3D", "human")
    assert result["ensembl_id"] == "ENSG_TEST"
    assert "not cell-type evidence" in result["scope"]
    assert store.records == {}


def test_ensembl_not_found_is_explicit(tmp_path):
    client = httpx.Client(transport=httpx.MockTransport(lambda req: httpx.Response(400)))
    assert (
        EvidenceStore(tmp_path, True, client).lookup_gene("FAKE", "human")["status"] == "not_found"
    )


def test_ensembl_rejects_url_input(tmp_path):
    with pytest.raises(ValueError, match="single gene"):
        EvidenceStore(tmp_path, False).lookup_gene("https://example.com", "human")


def test_model_and_credentials_required_before_output(demo_dir, monkeypatch):
    cfg = load_config(demo_dir / "scrna.json")
    with pytest.raises(ValueError, match="model"):
        run_agent(cfg)
    cfg.model = "test-model"
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    with pytest.raises(ValueError, match="ANTHROPIC_API_KEY"):
        run_agent(cfg)
    assert not cfg.out.exists()


class Block:
    def __init__(self, name):
        self.name = name

    def model_dump(self, **kwargs):
        return {"type": "tool_use", "id": self.name, "name": self.name, "input": {}}


class FakeClient:
    def __init__(self, tool_names):
        self.messages = self
        self.names = iter(tool_names)
        self.requests = []

    def create(self, **kwargs):
        self.requests.append(kwargs)
        return SimpleNamespace(
            content=[Block(next(self.names))],
            stop_reason="tool_use",
            usage=SimpleNamespace(model_dump=lambda: {"input_tokens": 100, "output_tokens": 20}),
        )


def test_fake_claude_completes_real_pipeline(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    cfg.model = "fake-test-model"
    client = FakeClient(
        ["inspect_data", "preprocess_data", "find_markers", "annotate", "generate_report"]
    )
    result = run_agent(cfg, client)
    assert result["status"] == "completed"
    assert result["tool_calls"] == 5
    assert (cfg.out / "report.md").exists()


def test_agent_stops_on_repeated_invalid_tools(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    cfg.model = "fake-test-model"
    result = run_agent(cfg, FakeClient(["invalid"] * 3))
    assert result["status"] == "incomplete" and result["reason"] == "consecutive_tool_errors"


def test_agent_obeys_tool_budget(demo_dir):
    cfg = load_config(demo_dir / "scrna.json")
    cfg.model = "fake-test-model"
    cfg.max_tool_calls = 1
    result = run_agent(cfg, FakeClient(["inspect_data", "preprocess_data"]))
    assert result["status"] == "incomplete"
    assert result["reason"] == "tool_or_token_budget_exhausted"
