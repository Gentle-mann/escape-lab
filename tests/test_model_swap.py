import json

import httpx
import pytest

from escape_lab.data import sha256
from escape_lab.improvement import write_json
from escape_lab.model_swap import (
    AuditedClient,
    released_measurements,
    score_released,
    verify_commitment,
)


def test_adapter_keeps_original_payload_and_records_only_safe_metadata(tmp_path):
    captured = []

    def respond(request):
        payload = json.loads(request.content)
        captured.append(payload)
        return httpx.Response(200, json={
            "id": "synthetic", "model": "glm", "usage": {
                "prompt_tokens": 4, "completion_tokens": 3, "total_tokens": 7},
            "choices": [{"finish_reason": "tool_calls", "message": {
                "content": None, "reasoning_content": "PRIVATE_NOT_FOR_LOGGING",
                "tool_calls": []}}]})

    from escape_lab.crusoe import ModelInfo
    with AuditedClient("fake", model_id="glm", transport=httpx.MockTransport(respond),
                       audit_path=tmp_path / "requests.jsonl", max_tokens=8192,
                       max_retries=0) as client:
        client.models = {"glm": ModelInfo(id="glm", supported_parameters=["tools"])}
        client.chat([{"role": "user", "content": "synthetic"}], [])
    assert set(captured[0]) == {"model", "messages", "tools", "tool_choice", "max_tokens", "stream"}
    assert captured[0]["max_tokens"] == 8192
    log = (tmp_path / "requests.jsonl").read_text()
    assert "PRIVATE_NOT_FOR_LOGGING" not in log
    assert json.loads(log)["response"]["usage"]["total_tokens"] == 7


def test_predictions_must_match_frozen_commitment_before_outcomes(tmp_path):
    p = tmp_path / "predictions.json"
    p.write_text("{}")
    freeze = tmp_path / "freeze.json"
    write_json(freeze, {"retrospective_reuse_authorized": True,
                        "hashes": {str(p): sha256(p)}})
    verify_commitment(freeze)
    p.write_text('{"changed": true}')
    with pytest.raises(ValueError, match="changed"):
        verify_commitment(freeze)


def test_retrospective_requires_explicit_label(tmp_path):
    p = tmp_path / "freeze.json"
    write_json(p, {"hashes": {}})
    with pytest.raises(ValueError, match="retrospective"):
        verify_commitment(p)


def test_released_measurements_require_complete_matching_original_methods():
    manifest = {"variants": [{"variant_id": "v1", "group_id": "g1"}]}
    methods = ("v1_fixed", "biology_fixed", "random", "grid", "agent")
    rows = [{"method": m, "variant_id": "v1", "drug": d, "group_id": "g1",
             "residual_activity": 10.0, "measurement_status": "observed", "source_cell": "B2"}
            for m in methods for d in ("pem", "fut")]
    assert len(released_measurements(rows, manifest)) == 2
    rows[-1]["residual_activity"] = 20.0
    with pytest.raises(ValueError, match="disagree"):
        released_measurements(rows, manifest)
    with pytest.raises(ValueError, match="coverage"):
        released_measurements(rows[:-1], manifest)


def test_retrospective_scores_released_records_without_workbook_or_receipt_write(
    tmp_path, monkeypatch
):
    import openpyxl

    def forbid_workbook(*args, **kwargs):
        pytest.fail("Retrospective evaluation must never open a workbook")

    monkeypatch.setattr(openpyxl, "load_workbook", forbid_workbook)
    reference, output = tmp_path / "reference", tmp_path / "new"
    receipt = tmp_path / "data/revision-v2/external-receipt.json"
    variants = [{"variant_id": f"v{i}", "group_id": f"g{i // 2}"} for i in range(6)]
    manifest = {"panel_id": "synthetic", "variants": variants}
    methods = ("v1_fixed", "biology_fixed", "random", "grid", "agent")
    rows = [{"method": m, "variant_id": v["variant_id"], "group_id": v["group_id"],
             "drug": d, "residual_activity": float(i), "measurement_status": "observed",
             "source_cell": f"B{i+1}"}
            for m in methods for i, v in enumerate(variants) for d in ("pem", "fut")]
    records = reference / "external/external-records.json"
    write_json(records, {"records": rows})
    report = reference / "external/external-report.json"
    write_json(report, {"records_path": str(records), "records_sha256": sha256(records)})
    write_json(receipt, {"status": "completed", "report_path": str(report),
                         "report_sha256": sha256(report), "records_sha256": sha256(records)})
    before = receipt.read_bytes()
    predictions = output / "external-predictions.json"
    write_json(predictions, {"predictions": [
        {"method": r["method"], "variant_id": r["variant_id"], "drug": r["drug"],
         "score": r["residual_activity"] / 6} for r in rows]})
    freeze = output / "external-freeze.json"
    write_json(freeze, {"retrospective_reuse_authorized": True,
                        "hashes": {str(predictions): sha256(predictions)},
                        "released_artifact_sha256": {str(p): sha256(p)
                                                     for p in (receipt, report, records)}})
    result = score_released(tmp_path, reference, output, manifest, freeze)
    assert result["retrospective_reuse"] is True
    assert result["methods"]["agent"]["macro_spearman"] == pytest.approx(1)
    assert receipt.read_bytes() == before
    records.write_text("{}")
    with pytest.raises(ValueError, match="changed"):
        score_released(tmp_path, reference, output, manifest, freeze)
