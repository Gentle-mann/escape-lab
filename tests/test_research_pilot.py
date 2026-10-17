import hashlib
import json
from copy import deepcopy
from pathlib import Path

import httpx
import pytest

from escape_lab.crusoe import CrusoeError, ModelInfo
from escape_lab.data import sha256
from escape_lab.experiment_planner import SYSTEM_PROMPT as ORIGINAL_PROMPT
from escape_lab.improvement import write_json
from escape_lab.research_pilot import (
    ARMS,
    SYSTEM_PROMPT,
    EvidenceClient,
    score_joint_released,
    verify_joint_commitment,
)


def test_both_arms_share_instructions_and_only_evidence_differs(tmp_path):
    captured = []

    def respond(request):
        captured.append(json.loads(request.content))
        return httpx.Response(200, json={
            "id": "synthetic", "model": "glm", "usage": {
                "prompt_tokens": 4, "completion_tokens": 3, "total_tokens": 7},
            "choices": [{"finish_reason": "tool_calls", "message": {
                "content": None, "tool_calls": []}}]})

    card = {"id": "geometry", "claim": "Distances are geometric descriptors.",
            "source_url": "https://example.org/paper", "publication_year": 2020}
    incoming = [{"role": "system", "content": ORIGINAL_PROMPT},
                {"role": "user", "content": '{"diagnostics":{"macro_ap":0.2}}'}]
    before = deepcopy(incoming)
    for arm, cards in (("control", []), ("research", [card])):
        with EvidenceClient("fake", model_id="glm", evidence_cards=cards,
                            audit_path=tmp_path / f"{arm}.jsonl", max_tokens=8192,
                            max_retries=0, transport=httpx.MockTransport(respond)) as client:
            client.models = {"glm": ModelInfo(id="glm", supported_parameters=["tools"])}
            client.chat(incoming, [{"type": "function"}], prompt_version="original")
    assert incoming == before
    assert ORIGINAL_PROMPT != SYSTEM_PROMPT
    assert captured[0]["messages"][0] == captured[1]["messages"][0]
    assert captured[0]["messages"][0]["content"] == SYSTEM_PROMPT
    assert "Only aggregate inner-cross-validation feedback is supplied." not in SYSTEM_PROMPT
    assert "Use only supplied catalog definitions and observed diagnostics;" not in SYSTEM_PROMPT
    assert json.loads(captured[0]["messages"][1]["content"])["biological_evidence"] == []
    assert json.loads(captured[1]["messages"][1]["content"])["biological_evidence"] == [card]
    left, right = deepcopy(captured[0]), deepcopy(captured[1])
    left["messages"].pop(1)
    right["messages"].pop(1)
    assert left == right
    audit = json.loads((tmp_path / "research.adapter.jsonl").read_text())
    assert audit["effective_system_prompt_sha256"] == hashlib.sha256(
        SYSTEM_PROMPT.encode()).hexdigest()
    assert audit["base_planner_prompt_version"] == "original"
    assert audit["effective_prompt_version"] != "original"


def test_evidence_context_limit_is_checked_before_transport(tmp_path):
    def no_transport(request):
        pytest.fail("Oversized evidence must fail before transport")

    cards = [{"id": "large", "claim": "x" * 40_001,
              "source_url": "https://example.org", "publication_year": 2020}]
    with (EvidenceClient("fake", model_id="glm", evidence_cards=cards,
                         audit_path=tmp_path / "calls.jsonl",
                         transport=httpx.MockTransport(no_transport)) as client,
          pytest.raises(CrusoeError, match="context")):
        client.chat([{"role": "system", "content": ORIGINAL_PROMPT}], [],
                    prompt_version="original")


def test_adapter_rejects_unknown_card_fields_and_unexpected_prompt(tmp_path):
    with pytest.raises(ValueError, match="card"):
        EvidenceClient("fake", evidence_cards=[{"id": "x", "outcomes": [1]}],
                       audit_path=tmp_path / "bad.jsonl")
    with (EvidenceClient("fake", evidence_cards=[], audit_path=tmp_path / "good.jsonl") as client,
          pytest.raises(CrusoeError, match="planner prompt")):
        client.chat([{"role": "system", "content": "Outcome-bearing replacement"}], [],
                    prompt_version="unknown")


def joint_fixture(tmp_path):
    root = tmp_path / "root"
    output = root / "run"
    protocol, evidence = output / "protocol.json", output / "evidence.json"
    write_json(protocol, {"synthetic": True})
    write_json(evidence, {"cards": []})
    hashes = {str(p): sha256(p) for p in (protocol, evidence)}
    arms = {}
    for arm in ARMS:
        arm_output = output / arm
        paths = {"prediction_path": arm_output / "external-predictions.json",
                 "search_path": arm_output / "full-search/search.json",
                 "fitted_models_path": arm_output / "fitted-models.json",
                 "external_freeze_path": arm_output / "external-freeze.json"}
        model = arm_output / "models/agent/pem.pkl"
        model.parent.mkdir(parents=True, exist_ok=True)
        model.write_bytes(b"synthetic model")
        model_hashes = {str(model): sha256(model)}
        for name, path in paths.items():
            write_json(path, {"synthetic": name})
            hashes[str(path)] = sha256(path)
        hashes.update(model_hashes)
        arms[arm] = {**{k: str(v) for k, v in paths.items()}, "model_hashes": model_hashes,
                     "successful_proposals": 1, "selected_configurations": {"agent": "fixed"}}
    freeze = output / "joint-freeze.json"
    write_json(freeze, {"schema_version": "research-pilot-joint-freeze-v1",
                        "retrospective_reuse_authorized": True,
                        "protocol_path": str(protocol), "evidence_path": str(evidence),
                        "arms": arms, "hashes": hashes})
    return root, output, freeze


@pytest.mark.parametrize("mutation", ["missing_arm", "changed_prediction", "missing_model", "zero_valid"])
def test_joint_gate_fails_before_any_released_scoring(tmp_path, monkeypatch, mutation):
    import escape_lab.research_pilot as pilot

    root, output, freeze = joint_fixture(tmp_path)
    payload = json.loads(freeze.read_text())
    if mutation == "missing_arm":
        payload["arms"].pop("research")
    elif mutation == "changed_prediction":
        Path(payload["arms"]["research"]["prediction_path"]).write_text("changed")
    elif mutation == "missing_model":
        payload["hashes"].pop(next(iter(payload["arms"]["research"]["model_hashes"])))
    else:
        payload["arms"]["research"]["successful_proposals"] = 0
    write_json(freeze, payload)

    def forbid_score(*args):
        pytest.fail("Scorer called before complete two-arm commitment")

    monkeypatch.setattr(pilot, "score_released", forbid_score)
    with pytest.raises(ValueError):
        score_joint_released(root, root / "reference", output, {}, freeze)


def test_joint_gate_revalidates_before_second_arm(tmp_path, monkeypatch):
    import escape_lab.research_pilot as pilot

    root, output, freeze = joint_fixture(tmp_path)
    assert verify_joint_commitment(freeze)["retrospective_reuse_authorized"] is True
    calls = []

    def score(*args):
        calls.append(args[2].name)
        payload = json.loads(freeze.read_text())
        Path(payload["arms"]["research"]["prediction_path"]).write_text("changed")
        return {"status": "completed"}

    monkeypatch.setattr(pilot, "score_released", score)
    with pytest.raises(ValueError, match="changed"):
        score_joint_released(root, root / "reference", output, {}, freeze)
    assert calls == ["control"]


def test_arm_runner_uses_fresh_clients_and_does_not_transfer_search_results(tmp_path, monkeypatch):
    import escape_lab.research_pilot as pilot

    clients, sent = [], []
    original_client = EvidenceClient

    def respond(request):
        if request.method == "GET":
            return httpx.Response(200, json={"data": [
                {"id": pilot.MODEL, "supported_parameters": ["tools"]}]})
        sent.append(json.loads(request.content))
        return httpx.Response(200, json={
            "id": "synthetic", "model": pilot.MODEL,
            "usage": {"prompt_tokens": 4, "completion_tokens": 3, "total_tokens": 7},
            "choices": [{"finish_reason": "tool_calls", "message": {"tool_calls": []}}]})

    def client_factory(*args, **kwargs):
        client = original_client(*args, **kwargs, transport=httpx.MockTransport(respond))
        clients.append(client)
        return client

    def search(training, labels, store, planner, seed, deadline, output):
        assert planner.proposals_made == 0
        assert planner.client.total_tokens == 0
        marker = "CONTROL_ONLY_HISTORY" if output.parent.name == "control" else "RESEARCH_ONLY_HISTORY"
        planner.client.chat([{"role": "system", "content": ORIGINAL_PROMPT},
                             {"role": "user", "content": marker}], [], prompt_version="base")
        return {"selections": {"agent": "fixed"}, "agent_status": "completed",
                "successful_proposals": 3, "failed_opportunities": 0,
                "agent_proposals": [{"calls": 1}],
                "inner_metrics": {"OTHER_ARM_GRID_RESULT": 123}}

    def predict(*args, **kwargs):
        return {"model_hashes": {}, "predictions": [], "diagnostics": []}

    monkeypatch.setattr(pilot, "credentials", lambda root: ("fake", None))
    monkeypatch.setattr(pilot, "EvidenceClient", client_factory)
    monkeypatch.setattr(pilot, "search_candidates", search)
    monkeypatch.setattr(pilot, "predict_selected", predict)
    card = {"id": "geometry", "claim": "Distances are geometric descriptors.",
            "source_url": "https://example.org/paper", "publication_year": 2020}
    for arm in ARMS:
        pilot._run_arm(tmp_path, tmp_path / arm, [] if arm == "control" else [card],
                       [], {"RAW_LABEL_SENTINEL": True}, None, [])
    assert len(clients) == 2 and clients[0] is not clients[1]
    assert "CONTROL_ONLY_HISTORY" not in json.dumps(sent[1])
    assert "RESEARCH_ONLY_HISTORY" not in json.dumps(sent[0])
    assert all("OTHER_ARM_GRID_RESULT" not in json.dumps(p) and
               "RAW_LABEL_SENTINEL" not in json.dumps(p) for p in sent)
    assert clients[0].total_tokens == clients[1].total_tokens == 7
