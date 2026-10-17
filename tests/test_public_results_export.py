"""Publication checks use synthetic files, never experimental outcomes."""

import importlib.util
import json
from pathlib import Path

import pytest


@pytest.fixture
def exporter(tmp_path, monkeypatch):
    source = Path(__file__).resolve().parents[1] / "scripts/export_public_results.py"
    spec = importlib.util.spec_from_file_location("public_results_export", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "ROOT", tmp_path)
    monkeypatch.setattr(module, "DESTINATION", tmp_path / "docs/archive")
    return module


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def test_repeated_publication_exports_aggregates_and_leaves_private_sidecars(exporter):
    root = exporter.ROOT
    run = root / "reports/repeated-research/20261017-experiments-1-2"
    results = {"run_status": "completed", "candidate_count": 8,
               "methods": {"control": {"mean": {"macro_ap": .1}}},
               "usage": {"control": {"requests": 3}},
               "experiments": {"research": {"status": "inconclusive"}},
               "aligned_groups": ["private-group"],
               "records": [{"candidate_id": "private-variant", "label": True}]}
    write_json(run / "results.json", results)
    write_json(run / "private-scored-records.json", results["records"])
    (run / "requests.jsonl").write_text('{"messages":["private prompt"]}\n')
    audit = run.parent / "20261017-experiments-1-2-audit.json"
    write_json(audit, {"run_status": "completed", "valid_proposals": 3,
                       "verification": {"all_hashes_valid": True},
                       "independent_request_audit": {"requests": 3},
                       "request_payload": {"messages": ["private prompt"]}})
    report = run.parent / "20261017-experiments-1-2-review.md"
    report.write_text(f"[Results](<{run / 'results.json'}>)\n"
                      f"[Audit](<{audit}>)\n"
                      f"[Private](<{run / 'private-scored-records.json'}>)\n")
    originals = {p: p.read_bytes() for p in (run / "results.json", audit, report)}

    exporter.main()

    published = exporter.DESTINATION / "reports/repeated-research"
    public_result = json.loads((published / run.name / "results.json").read_text())
    public_audit = json.loads((published / audit.name).read_text())
    assert public_result["methods"] == results["methods"]
    assert public_result["candidate_count"] == 8
    assert "records" not in public_result
    assert "aligned_groups" not in public_result
    assert public_audit["valid_proposals"] == 3
    assert "request_payload" not in public_audit
    assert not (published / run.name / "private-scored-records.json").exists()
    assert not (published / run.name / "requests.jsonl").exists()
    public_text = (published / report.name).read_text()
    assert str(root) not in public_text
    assert f"[Results]({run.name}/results.json)" in public_text
    assert f"[Audit]({audit.name})" in public_text
    assert "[Private]" not in public_text
    assert "local artifact" in public_text
    assert all(p.read_bytes() == raw for p, raw in originals.items())
    manifest = json.loads((exporter.DESTINATION / "manifest.json").read_text())
    assert len(manifest["artifacts"]) == 3
    for entry in manifest["artifacts"]:
        assert entry["source_sha256"] == exporter.digest((root / entry["source"]).read_bytes())
        assert entry["public_sha256"] == exporter.digest(
            (root / entry["public_path"]).read_bytes())


@pytest.mark.parametrize("private", [
    {"records": [{"candidate_id": "variant", "label": True}]},
    {"requests": [{"messages": ["private prompt"]}]},
    {"observations": [{"variant_id": "variant", "resistant": True}]},
    {"mean": {"macro_ap": [0, 1, 0, 1]}},
])
def test_repeated_publication_rejects_private_content_inside_allowed_fields(exporter, private):
    run = exporter.ROOT / "reports/repeated-research/20261017-experiments-1-2"
    write_json(run / "results.json", {"methods": {"control": private}})
    with pytest.raises(ValueError, match="private|aggregate"):
        exporter.main()
    assert not (exporter.DESTINATION / run.relative_to(exporter.ROOT) / "results.json").exists()


@pytest.mark.parametrize("private", [
    {"api_key": "synthetic-secret"},
    {"response": {"body": "synthetic-agent-output"}},
    {"requests": ["synthetic-agent-input"]},
])
def test_repeated_usage_exports_counts_only(exporter, private):
    run = exporter.ROOT / "reports/repeated-research/20261017-experiments-1-2"
    write_json(run / "results.json", {"usage": {"control": private}})
    with pytest.raises(ValueError, match="private|aggregate"):
        exporter.main()
