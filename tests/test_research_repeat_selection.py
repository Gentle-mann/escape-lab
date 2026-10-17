import hashlib
import json
from copy import deepcopy

import httpx
import pytest

from escape_lab.crusoe import CrusoeError, ModelInfo
from escape_lab.experiment_planner import SYSTEM_PROMPT as BASE_SYSTEM_PROMPT
from escape_lab.experiment_planner import ExperimentProposal
from escape_lab.research_repeat_selection import (
    ARMS,
    NUMERICAL_POLICY,
    SYSTEM_PROMPT,
    RepeatEvidenceClient,
    numerical_selection,
    select_from_cache,
)
from escape_lab.revision_evaluation import BIOLOGY_FIXED, REFERENCE
from escape_lab.revision_models import catalog


def metric(ap):
    return {"macro_ap": ap, "per_drug": {
        drug: {"average_precision": ap, "brier": .1} for drug in ("pem", "fut")
    }, "by_fold": {f"inner-{i}": {"macro_ap": ap} for i in range(3)}}


def cache():
    return {candidate: metric(.2) for candidate in catalog()}


def proposal(candidate, parent=REFERENCE):
    return ExperimentProposal(candidate_id=candidate, parent_id=parent,
                              hypothesis="A specified representation improves AP.",
                              expected_effect="Higher pooled AP.",
                              rejection_condition="Acceptance guards fail.",
                              status="proposed", error=None, calls=1, model_id="synthetic")


class Planner:
    def __init__(self, responses):
        self.responses = iter(responses)
        self.feedback = []
        self.attempts = []

    def propose(self, configurations, feedback):
        self.feedback.append(deepcopy(feedback))
        response = next(self.responses)
        if isinstance(response, Exception):
            raise response
        return response


def test_only_reference_and_own_evaluated_metrics_enter_feedback():
    metrics = cache()
    first, second, third = list(catalog())[2:5]
    metrics[first] = metric(.3)
    metrics[second] = metric(.4)
    metrics[third] = metric(.5)
    for value in metrics.values():
        value.update(outer_outcomes=["OUTER_SECRET"], source_path="OUTER_FILE")
        value["per_drug"]["pem"]["raw_labels"] = ["RAW_SECRET"]
    planner = Planner([proposal(first), proposal(second, first), proposal(third, second)])
    result = select_from_cache(metrics, planner, seed=12)
    assert len(planner.feedback) == 3
    for index, feedback in enumerate(planner.feedback):
        assert {r["candidate_id"] for r in feedback["history"]} == {
            REFERENCE, *[first, second, third][:index]}
        serialized = json.dumps(feedback)
        assert all(secret not in serialized for secret in ("OUTER_SECRET", "OUTER_FILE", "RAW_SECRET"))
    assert result["selected"] == third
    assert result["successful_proposals"] == 3
    assert result["failed_opportunities"] == 0
    assert result["status"] == "completed"
    assert "OUTER_SECRET" not in json.dumps(result)
    assert "OUTER_FILE" not in json.dumps(result)


def test_failures_consume_opportunities_and_keep_incumbent():
    metrics = cache()
    metrics[BIOLOGY_FIXED] = metric(.3)
    planner = Planner([RuntimeError("secret raw exception"), proposal(BIOLOGY_FIXED),
                       proposal(BIOLOGY_FIXED)])
    result = select_from_cache(metrics, planner, seed=12)
    assert len(planner.feedback) == len(result["proposals"]) == 3
    assert result["selected"] == BIOLOGY_FIXED
    assert result["successful_proposals"] == 1
    assert result["failed_opportunities"] == 2
    assert result["status"] == "degraded"
    assert result["proposals"][0]["error_category"] == "RuntimeError"
    assert result["proposals"][2]["status"] == "duplicate"
    assert "secret raw exception" not in json.dumps(result)
    assert len(result["attempts"]) >= 2


def test_deadline_records_all_remaining_opportunities_without_calls():
    planner = Planner([])
    result = select_from_cache(cache(), planner, seed=12, seconds=0)
    assert result["selected"] == REFERENCE
    assert result["status"] == "unavailable"
    assert len(result["proposals"]) == result["failed_opportunities"] == 3
    assert all(p["status"] == "exhausted" for p in result["proposals"])
    assert planner.feedback == []


def test_late_result_is_preserved_but_not_evaluated(monkeypatch):
    import escape_lab.research_repeat_selection as selection

    clock = {"now": 0.}
    monkeypatch.setattr(selection.time, "monotonic", lambda: clock["now"])

    class LatePlanner(Planner):
        def propose(self, configurations, feedback):
            clock["now"] = 20.
            return super().propose(configurations, feedback)

    result = select_from_cache(cache(), LatePlanner([proposal(BIOLOGY_FIXED)]),
                               seed=12, seconds=10)
    assert result["selected"] == REFERENCE
    assert result["successful_proposals"] == 0
    assert result["proposals"][0]["candidate_id"] == BIOLOGY_FIXED
    assert result["proposals"][0]["status"] == "exhausted"
    assert "metrics" not in result["proposals"][0]


def test_unrequested_parent_and_unknown_candidate_are_failed_opportunities():
    planner = Planner([proposal(BIOLOGY_FIXED, "interactions-logistic-c1"),
                       proposal("outer-result"), proposal(BIOLOGY_FIXED)])
    result = select_from_cache(cache(), planner, seed=12)
    assert result["successful_proposals"] == 1
    assert result["failed_opportunities"] == 2
    assert result["proposals"][0]["status"] == "invalid_response"
    assert result["proposals"][1]["status"] == "invalid_response"


@pytest.mark.parametrize("changed", ["missing", "extra"])
def test_cache_requires_exact_controller_catalog(changed):
    metrics = cache()
    if changed == "missing":
        metrics.pop(BIOLOGY_FIXED)
    else:
        metrics["outer_scores"] = metric(.9)
    with pytest.raises(ValueError, match="catalog"):
        select_from_cache(metrics, Planner([]), seed=12)


@pytest.mark.parametrize("accepted_representation", ["chemistry", "structure", "interactions"])
def test_numerical_policy_uses_three_queries_with_frozen_coarse_to_fine_rule(accepted_representation):
    metrics = cache()
    if accepted_representation in ("structure", "interactions"):
        metrics[BIOLOGY_FIXED] = metric(.3)
    if accepted_representation == "interactions":
        metrics["interactions-logistic-c0.1"] = metric(.4)
    third = f"{accepted_representation}-logistic-c0.01"
    metrics[third] = metric(.5)
    # A superior unvisited metric cannot steer the local policy.
    metrics["interactions-hgb-depth2"] = metric(.99)
    result = numerical_selection(metrics, seed=35)
    assert result["tried"]["numerical"] == [BIOLOGY_FIXED, "interactions-logistic-c0.1", third]
    assert result["selections"]["numerical"] == third
    assert result["selections"]["grid"] == "interactions-hgb-depth2"
    assert result["selections"]["fixed"] == BIOLOGY_FIXED
    assert result["budgets"] == {"fixed": 0, "random": 3, "numerical": 3, "grid": 11}
    assert result["grid_budget_matched"] is False
    assert result["tried"]["fixed"] == []
    assert len(set(result["tried"]["random"])) == 3
    assert REFERENCE not in result["tried"]["random"]
    assert len(set(result["tried"]["grid"])) == 11
    assert result == numerical_selection(metrics, seed=35)
    assert NUMERICAL_POLICY["proposal_budget"] == 3


def test_all_adapters_share_prompt_and_bound_card_payloads(tmp_path):
    captured = []

    def respond(request):
        captured.append(json.loads(request.content))
        return httpx.Response(200, json={"id": "synthetic", "model": "glm", "usage": {
            "prompt_tokens": 4, "completion_tokens": 3, "total_tokens": 7},
            "choices": [{"finish_reason": "tool_calls", "message": {
                "content": None, "tool_calls": []}}]})

    relevant = {"id": "geometry", "claim": "Distances describe ligand geometry.",
                "source_url": "https://pubmed.ncbi.nlm.nih.gov/123/", "publication_year": 2020}
    placebo = {**relevant, "id": "placebo", "claim": "Plants grow through cellular division."}
    cards_by_arm = {"control": [], "relevant": [relevant], "placebo": [placebo]}
    incoming = [{"role": "system", "content": BASE_SYSTEM_PROMPT},
                {"role": "user", "content": '{"diagnostics":{"macro_ap":0.2}}'}]
    before = deepcopy(incoming)
    for arm in ARMS:
        cards = cards_by_arm[arm]
        with RepeatEvidenceClient("fake", model_id="glm", evidence_cards=cards,
                                  audit_path=tmp_path / f"{arm}.jsonl", max_tokens=16384,
                                  max_retries=0, transport=httpx.MockTransport(respond)) as client:
            client.models = {"glm": ModelInfo(id="glm", supported_parameters=["tools"])}
            client.chat(incoming, [{"type": "function"}], prompt_version="original")
    assert incoming == before
    for arm, payload in zip(ARMS, captured, strict=True):
        assert payload["messages"][0]["content"] == SYSTEM_PROMPT
        assert payload["max_tokens"] == 16384
        assert json.loads(payload["messages"][1]["content"])["biological_evidence"] == cards_by_arm[arm]
        audit = json.loads((tmp_path / f"{arm}.adapter.jsonl").read_text())
        assert audit["effective_system_prompt_sha256"] == hashlib.sha256(SYSTEM_PROMPT.encode()).hexdigest()
        assert audit["base_planner_prompt_version"] == "original"
    for payload in captured:
        payload["messages"].pop(1)
    assert captured[0] == captured[1] == captured[2]
    assert "240" in SYSTEM_PROMPT


def test_evidence_adapter_rejects_extra_fields_and_oversized_context(tmp_path):
    with pytest.raises(ValueError, match="card"):
        RepeatEvidenceClient("fake", evidence_cards=[{"id": "x", "outcomes": [1]}],
                             audit_path=tmp_path / "bad.jsonl")

    def forbidden(request):
        pytest.fail("Invalid payload reached transport")

    cards = [{"id": "large", "claim": "x" * 40001,
              "source_url": "https://pubmed.ncbi.nlm.nih.gov/123/", "publication_year": 2020}]
    with (RepeatEvidenceClient("fake", model_id="glm", evidence_cards=cards,
                               audit_path=tmp_path / "calls.jsonl",
                               transport=httpx.MockTransport(forbidden)) as client,
          pytest.raises(CrusoeError, match="context")):
        client.chat([{"role": "system", "content": BASE_SYSTEM_PROMPT}], [])


def test_evidence_adapter_rejects_replacement_system_prompt(tmp_path):
    with (RepeatEvidenceClient("fake", evidence_cards=[], audit_path=tmp_path / "calls.jsonl") as client,
          pytest.raises(CrusoeError, match="planner prompt")):
        client.chat([{"role": "system", "content": "Outer results: secret"}], [])
