import json

import httpx
from test_crusoe import completion

from escape_lab.agent import NemotronInvestigator
from escape_lab.baseline import NumericalSearch
from escape_lab.crusoe import CrusoeClient

SHORTLISTS = {
    "exploit": [{"candidate_id": "v1", "numerical_probability": 0.7,
                 "gene": "FGFR2", "mutation": "N549K"}],
    "explore": [{"candidate_id": "v2", "numerical_probability": 0.4}],
}


def make_agent(responses):
    bodies = []

    def handler(request):
        if request.method == "GET":
            return httpx.Response(200, json={"data": [{
                "id": "nvidia/test-nemotron", "supported_parameters": ["tools"]
            }]})
        bodies.append(json.loads(request.content))
        return httpx.Response(200, json=responses.pop(0))

    client = CrusoeClient("unused-test-key", model_id="nvidia/test-nemotron",
                          transport=httpx.MockTransport(handler))
    client.discover_models()
    return NemotronInvestigator(client), bodies


def selection(candidate_id, **kwargs):
    return completion(arguments={"candidate_id": candidate_id,
                                 "reason": "Test measured resistance coverage.",
                                 "evidence_ids": [], **kwargs})


def test_selected_probability_comes_from_numerical_model():
    agent, bodies = make_agent([selection("v2")])
    decision = agent.choose_action({"budget_remaining": 20}, SHORTLISTS)
    assert decision.candidate_id == "v2"
    assert decision.numerical_probability == 0.4
    assert decision.fallback_used is False
    assert decision.prompt_version
    assert decision.model_id == "nvidia/test-nemotron"
    assert len(bodies) == 1


def test_invalid_candidate_gets_one_correction_with_matching_tool_response():
    agent, bodies = make_agent([selection("invented"), selection("v1")])
    decision = agent.choose_action({"budget_remaining": 20}, SHORTLISTS)
    assert decision.candidate_id == "v1"
    assert decision.fallback_used is False
    correction_messages = bodies[1]["messages"]
    assert correction_messages[-2]["role"] == "assistant"
    assert correction_messages[-1]["role"] == "tool"
    assert correction_messages[-1]["tool_call_id"] == "call-1"


def test_repeated_invalid_calls_fall_back_explicitly_and_do_not_fake_agent_success():
    agent, bodies = make_agent([selection("invented"), selection("invented")])
    decision = agent.choose_action({"budget_remaining": 20}, SHORTLISTS)
    assert decision.candidate_id == "v1"
    assert decision.fallback_used is True
    assert "invalid" in decision.fallback_reason.lower()
    assert len(bodies) == 2


def test_llm_cannot_supply_its_own_calibrated_probability():
    agent, _ = make_agent([selection("v1", probability=0.999),
                          selection("v1", probability=0.999)])
    decision = agent.choose_action({"budget_remaining": 20}, SHORTLISTS)
    assert decision.fallback_used
    assert decision.numerical_probability == 0.7


def test_revealed_contradictory_evidence_is_visible_but_hidden_state_is_not():
    agent, bodies = make_agent([selection("v2", evidence_ids=["ev1"])])
    state = {"budget_remaining": 19, "revealed_ids": ["v1"],
             "observations": [{"evidence_id": "ev1", "candidate_id": "v1",
                               "dual_resistant": False, "prediction_before": 0.9}],
             "hidden_outcomes": {"v2": "must never serialize"},
             "api_key": "must never serialize"}
    decision = agent.choose_action(state, SHORTLISTS)
    assert decision.candidate_id == "v2"
    payload = json.dumps(bodies[0]["messages"])
    assert "dual_resistant" in payload
    assert "false" in payload
    assert "must never serialize" not in payload
    assert "v1" not in bodies[0]["messages"][-1]["content"].split('"shortlists":')[1]


def test_premature_finish_is_rejected_for_coverage_audit():
    finish = completion("finish", {"status": "unresolved", "reason": "Enough.",
                                    "evidence_ids": []})
    agent, _ = make_agent([finish, selection("v2")])
    decision = agent.choose_action({"budget_remaining": 20,
                                    "objective": "coverage_audit"}, SHORTLISTS)
    assert decision.action == "select_candidate"
    assert decision.candidate_id == "v2"


def test_zero_budget_stops_without_inference():
    agent, bodies = make_agent([])
    decision = agent.choose_action({"budget_remaining": 0}, SHORTLISTS)
    assert decision.action == "finish"
    assert decision.status == "budget_exhausted"
    assert decision.fallback_used is False
    assert not bodies


def test_hallucinated_evidence_and_duplicate_candidates_are_rejected():
    agent, bodies = make_agent([selection("v1", evidence_ids=["invented-evidence"]),
                               selection("v1")])
    decision = agent.choose_action({"budget_remaining": 19, "revealed_ids": ["v1"]},
                                   SHORTLISTS)
    assert decision.fallback_used
    assert decision.candidate_id == "v2"
    assert len(bodies) == 2


def test_missing_tool_calls_require_correction_then_explicit_fallback():
    no_tool = {"choices": [{"message": {"role": "assistant", "content": "Trust me."}}]}
    agent, bodies = make_agent([no_tool, no_tool])
    decision = agent.choose_action({"budget_remaining": 20}, SHORTLISTS)
    assert decision.fallback_used
    assert len(bodies) == 2
    assert bodies[1]["messages"][-1]["role"] == "user"


def test_verified_falsification_can_finish_with_supporting_evidence():
    finish = completion("finish", {"status": "counterexample_verified",
                                    "reason": "The revealed assay refutes universal coverage.",
                                    "evidence_ids": ["ev1"]})
    agent, _ = make_agent([finish])
    decision = agent.choose_action({"budget_remaining": 19, "objective": "falsification",
                                    "observations": [{"evidence_id": "ev1", "candidate_id": "v1",
                                                      "dual_resistant": True}]}, SHORTLISTS)
    assert decision.action == "finish"
    assert decision.status == "counterexample_verified"
    assert not decision.fallback_used


def test_actual_baseline_features_survive_but_unrevealed_outcomes_do_not():
    candidate = {
        "candidate_id": "v2", "gene": "FGFR2", "protein_variant": "p.N549K",
        "wt": "N", "mut": "K", "context": "ABCDE", "aligned_position": 549,
        "group_id": "block-54", "dual_resistant": True,
        "pem_median": 999, "fut_median": 888, "hidden_outcomes": "do-not-serialize",
        "score": 999, "median": 888, "classification": "hidden-class",
        "replicates": [999, 888, 777],
    }
    shortlists = NumericalSearch().shortlists([candidate], revealed_ids=set())
    agent, bodies = make_agent([selection("v2", evidence_ids=["ev1"])])
    decision = agent.choose_action({
        "budget_remaining": 19,
        "observations": [{"evidence_id": "ev1", "candidate_id": "v1",
                          "dual_resistant": False, "pem_median": 1.2, "fut_median": 0.8,
                          "score": 1.2, "median": 0.8, "classification": "neutral",
                          "replicates": [1.1, 1.2, 1.3],
                          "hidden_outcomes": "do-not-serialize"}],
    }, shortlists)
    assert decision.candidate_id == "v2"
    sent = json.loads(bodies[0]["messages"][-1]["content"])
    for policy_rows in sent["shortlists"].values():
        row = policy_rows[0]
        for field in ("protein_variant", "wt", "mut", "context", "group_id",
                      "model_disagreement", "exploration_score", "gene_prior"):
            assert row[field] == shortlists["exploit"][0][field]
        assert not {"dual_resistant", "pem_median", "fut_median", "hidden_outcomes",
                    "score", "median", "classification", "replicates"} & row.keys()
    observation = sent["state"]["observations"][0]
    assert observation["pem_median"] == 1.2
    assert observation["fut_median"] == 0.8
    assert observation["score"] == 1.2
    assert observation["median"] == 0.8
    assert observation["classification"] == "neutral"
    assert observation["replicates"] == [1.1, 1.2, 1.3]
    assert "do-not-serialize" not in json.dumps(sent)


def test_decision_limit_bounds_inference_and_flags_fallback_continuation():
    agent, bodies = make_agent([selection("v1")])
    agent.max_decisions = 1
    first = agent.choose_action({"budget_remaining": 20}, SHORTLISTS)
    second = agent.choose_action({"budget_remaining": 19, "revealed_ids": ["v1"]}, SHORTLISTS)
    assert not first.fallback_used
    assert second.fallback_used
    assert second.candidate_id == "v2"
    assert "budget" in second.fallback_reason
    assert len(bodies) == 1


def test_token_stop_threshold_prevents_a_correction_call():
    agent, bodies = make_agent([selection("invented")])
    agent.client.max_total_tokens = 30
    decision = agent.choose_action({"budget_remaining": 20}, SHORTLISTS)
    assert decision.fallback_used
    assert "token budget" in decision.fallback_reason
    assert len(bodies) == 1


def test_oversized_correction_cannot_bypass_context_size_limit():
    oversized = completion()
    oversized["choices"][0]["message"]["content"] = "x" * 50_000
    agent, bodies = make_agent([oversized, selection("v2")])
    decision = agent.choose_action({"budget_remaining": 20}, SHORTLISTS)
    assert decision.fallback_used
    assert "context" in decision.fallback_reason.lower()
    assert len(bodies) == 1
