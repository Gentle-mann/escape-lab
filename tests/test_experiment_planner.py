import hashlib
import json
from copy import deepcopy
from dataclasses import asdict
from types import SimpleNamespace

import pytest

from escape_lab.crusoe import CrusoeError
from escape_lab.experiment_planner import (
    PROMPT_HASH,
    PROMPT_VERSION,
    ExperimentPlanner,
    sanitize_feedback,
)

BASE = "chemistry-logistic-c0.1"
NEXT = "structure-logistic-c0.1"
THIRD = "structure-hgb-depth2"
CATALOG = {
    BASE: {"representation": "chemistry", "model": "logistic", "C": 0.1, "max_depth": 2},
    NEXT: {"representation": "structure", "model": "logistic", "C": 0.1, "max_depth": 2},
    THIRD: {
        "representation": "structure",
        "model": "hist_gradient_boosting",
        "C": 0.1,
        "max_depth": 2,
    },
}
METRICS = {
    "macro_ap": 0.63,
    "macro_brier": 0.14,
    "mean_fold_macro_ap": 0.62,
    "fold_ap_std": 0.02,
    "coefficient_norm_mean": 2.3,
    "per_drug": {
        "pem": {
            "average_precision": 0.64,
            "brier": 0.13,
            "n": 100,
            "positive_count": 20,
            "prevalence": 0.2,
        },
        "fut": {
            "average_precision": 0.62,
            "brier": 0.15,
            "n": 100,
            "positive_count": 25,
            "prevalence": 0.25,
        },
    },
}


def feedback():
    return deepcopy(
        {
            "incumbent": {"candidate_id": BASE, "metrics": METRICS},
            "history": [{"candidate_id": BASE, "status": "baseline", "metrics": METRICS}],
            "diagnostics": {"by_fold": {"inner-0": METRICS}},
            "structural_coverage": {"FGFR1": 0.9, "FGFR2": 0.8},
            "remaining_seconds": 200,
            "seed": 13,
        }
    )


def response(candidate_id=NEXT, parent_id=BASE, **extra):
    arguments = {
        "candidate_id": candidate_id,
        "parent_id": parent_id,
        "hypothesis": "Structural features reduce errors in unfamiliar aligned positions.",
        "expected_effect": "Inner macro AP rises by at least 0.01 without higher Brier.",
        "rejection_condition": "Reject if paired inner macro AP is not higher or Brier exceeds guard.",
        **extra,
    }
    return {
        "role": "assistant",
        "content": None,
        "tool_calls": [
            {
                "id": "experiment-1",
                "type": "function",
                "function": {"name": "propose_experiment", "arguments": json.dumps(arguments)},
            }
        ],
    }


class FakeClient:
    model_id = "nvidia/test-nemotron"

    def __init__(self, responses):
        self.responses = responses
        self.requests = []
        self.total_tokens = 0
        self.usage_records = []

    def chat(self, messages, tools, *, prompt_version):
        self.requests.append(
            {
                "messages": json.loads(json.dumps(messages)),
                "tools": tools,
                "prompt_version": prompt_version,
            }
        )
        self.total_tokens += 31
        self.usage_records.append({"budget_charge_tokens": 31, "usage_estimated": False})
        value = self.responses.pop(0)
        if isinstance(value, Exception):
            raise value
        if isinstance(value, tuple):
            value, finish_reason = value
        else:
            finish_reason = "tool_calls"
        return SimpleNamespace(
            message=value,
            finish_reason=finish_reason,
            usage={"prompt_tokens": 20, "completion_tokens": 11, "total_tokens": 31},
        )


def test_valid_selection_preserves_hypothesis_and_usage_fingerprints():
    client = FakeClient([response()])
    planner = ExperimentPlanner(client)
    result = planner.propose(CATALOG, feedback())
    assert result.status == "proposed"
    assert result.candidate_id == NEXT
    assert result.parent_id == BASE
    assert result.hypothesis and result.expected_effect and result.rejection_condition
    assert result.calls == 1
    assert result.error is None
    assert result.model_id == client.model_id
    assert result.prompt_version == PROMPT_VERSION
    assert result.usage["budget_charge_tokens"] == 31
    assert result.fingerprints["prompt_hash"] == PROMPT_HASH
    assert len(result.fingerprints["catalog_hash"]) == 64
    assert planner.attempts[0]["candidate_id"] == NEXT
    system = client.requests[0]["messages"][0]["content"]
    assert hashlib.sha256(system.encode()).hexdigest() == PROMPT_HASH


def test_schema_is_flat_and_enumerates_only_remaining_candidates_and_evaluated_parents():
    client = FakeClient([response()])
    ExperimentPlanner(client).propose(CATALOG, feedback())
    schema = client.requests[0]["tools"][0]["function"]["parameters"]
    assert schema["additionalProperties"] is False
    assert schema["properties"]["candidate_id"]["enum"] == sorted([NEXT, THIRD])
    assert schema["properties"]["parent_id"]["enum"] == [BASE]
    assert set(schema["required"]) == {
        "candidate_id",
        "parent_id",
        "hypothesis",
        "expected_effect",
        "rejection_condition",
    }
    assert "$ref" not in json.dumps(schema)
    assert "$defs" not in schema


def test_fixed_semantics_and_acceptance_card_match_controller_protocol():
    client = FakeClient([response()])
    ExperimentPlanner(client).propose(CATALOG, feedback())
    system = client.requests[0]["messages"][0]["content"]
    assert "larger C means weaker shrinkage" in system
    assert "Representations are cumulative" in system
    assert "four explicit gene indicators" in " ".join(system.split())
    assert "four-paralog" not in system
    assert "training median" in system
    assert "reference Brier * 1.05 + 0.001" in system
    assert "at least two of three inner folds" in system
    assert "strictly above the\ncurrent incumbent" in system
    assert "Outer and external results are never proposal feedback" in system


@pytest.mark.parametrize(
    "invalid",
    [
        response("outside-catalog"),
        response(BASE),
        response(parent_id=THIRD),
        response(parent_id=None),
        response(hypothesis=" "),
        response(expected_effect=""),
        response(rejection_condition=""),
        response(code="arbitrary execution"),
        response(candidate_id=1),
    ],
)
def test_invalid_selection_repairs_once_with_matched_tool_messages(invalid):
    client = FakeClient([invalid, response()])
    planner = ExperimentPlanner(client)
    result = planner.propose(CATALOG, feedback())
    assert result.status == "proposed"
    assert result.calls == 2
    messages = client.requests[1]["messages"]
    assert messages[-2]["role"] == "assistant"
    assert messages[-1]["role"] == "tool"
    assert messages[-1]["tool_call_id"] == "experiment-1"
    assert [entry["status"] for entry in planner.attempts] == ["invalid_response", "proposed"]


def test_two_invalid_responses_are_failure_without_baseline_substitution():
    client = FakeClient([response("unknown"), response("unknown")])
    result = ExperimentPlanner(client).propose(CATALOG, feedback())
    assert result.status == "invalid_response"
    assert result.candidate_id is None and result.parent_id is None
    assert result.calls == 2
    assert result.usage["budget_charge_tokens"] == 62


def test_prior_proposals_are_not_repeated_and_are_not_implicitly_evaluated_parents():
    client = FakeClient([response(), response(THIRD, NEXT), response(THIRD, BASE)])
    planner = ExperimentPlanner(client)
    assert planner.propose(CATALOG, feedback()).candidate_id == NEXT
    result = planner.propose(CATALOG, feedback())
    assert result.candidate_id == THIRD
    assert result.parent_id == BASE
    assert result.calls == 2
    schema = client.requests[1]["tools"][0]["function"]["parameters"]
    assert schema["properties"]["candidate_id"]["enum"] == [THIRD]
    sent = json.loads(client.requests[1]["messages"][1]["content"])
    assert sent["previous_proposals"][0]["candidate_id"] == NEXT


def test_evaluated_previous_candidate_can_become_parent_and_remaining_catalog_exhausts():
    client = FakeClient([response(THIRD, NEXT)])
    state = feedback()
    state["history"].append({"candidate_id": NEXT, "status": "rejected", "metrics": METRICS})
    planner = ExperimentPlanner(client)
    assert planner.propose(CATALOG, state).parent_id == NEXT
    result = planner.propose(CATALOG, state)
    assert result.status == "exhausted"
    assert result.calls == 0


def test_aggregate_sanitization_removes_nested_data_and_unapproved_labels():
    state = feedback()
    secret = "SECRET-RAW-OUTCOME"
    state.update(raw_rows=secret, schema=secret, objective=secret, selection_rule=secret)
    state["incumbent"].update(hypothesis=secret, config={"raw": secret})
    state["diagnostics"]["by_fold"]["outer-0"] = {"macro_ap": 0.99}
    state["diagnostics"]["by_fold"]["inner-0"]["replicates"] = secret
    state["diagnostics"]["by_fold"]["inner-0"]["per_drug"]["pem"]["predictions"] = secret
    state["diagnostics"]["by_fold"]["inner-0"]["per_drug"][secret] = {"brier": 0.0}
    state["history"].append({"candidate_id": secret, "metrics": METRICS})
    state["structural_coverage"].update({secret: 1, "FGFR3": {"rows": secret}})
    state["diagnostics"]["macro_ap"] = float("nan")
    safe = sanitize_feedback(state, CATALOG)
    assert secret not in json.dumps(safe)
    assert "outer-0" not in safe["diagnostics"]["by_fold"]
    assert "macro_ap" not in safe["diagnostics"]
    assert safe["diagnostics"]["by_fold"]["inner-0"]["per_drug"]["pem"]["n"] == 100
    assert safe["structural_coverage"] == {"FGFR1": 0.9, "FGFR2": 0.8}


def test_catalog_unknown_fields_are_not_sent_and_invalid_known_fields_fail_closed():
    catalog = {key: dict(value) for key, value in CATALOG.items()}
    catalog[NEXT]["raw_outcomes"] = "SECRET-CATALOG"
    client = FakeClient([response()])
    assert ExperimentPlanner(client).propose(catalog, feedback()).status == "proposed"
    assert "SECRET-CATALOG" not in json.dumps(client.requests)
    catalog[NEXT]["model"] = "exec-python"
    result = ExperimentPlanner(FakeClient([])).propose(catalog, feedback())
    assert result.status == "invalid_catalog"
    assert result.calls == 0


def test_transport_failure_preserves_only_allowlisted_category_and_consumes_opportunity():
    client = FakeClient([CrusoeError("Crusoe transport failed (ReadTimeout).")])
    planner = ExperimentPlanner(client, max_proposals=1)
    first = planner.propose(CATALOG, feedback())
    assert first.status == "failed" and first.candidate_id is None
    assert "ReadTimeout" in first.error
    assert first.usage["budget_charge_tokens"] == 31
    assert planner.propose(CATALOG, feedback()).status == "exhausted"
    assert len(client.requests) == 1


def test_arbitrary_provider_errors_and_private_fields_do_not_reach_logs_or_results():
    client = FakeClient([CrusoeError("SECRET-PROVIDER-BODY")])
    planner = ExperimentPlanner(client)
    result = planner.propose(CATALOG, feedback())
    assert "SECRET" not in json.dumps({"result": asdict(result), "attempts": planner.attempts})
    assert "CrusoeError" in result.error


def test_truncation_immediately_fails_without_unchanged_cap_retry():
    value = response()
    value.update(content="SECRET-CONTENT", reasoning_content="SECRET-REASONING")
    client = FakeClient([(value, "length")])
    planner = ExperimentPlanner(client)
    result = planner.propose(CATALOG, feedback())
    assert result.status == "truncated_response"
    assert result.candidate_id is None and result.calls == 1
    assert planner.attempts[0]["finish_reason"] == "length"
    assert "SECRET" not in json.dumps({"result": asdict(result), "attempts": planner.attempts})


def test_missing_tool_calls_and_wrong_envelopes_fail_safely():
    client = FakeClient([{"content": "No tool"}, {"tool_calls": "malformed"}])
    result = ExperimentPlanner(client).propose(CATALOG, feedback())
    assert result.status == "invalid_response"
    assert result.calls == 2 and result.candidate_id is None


def test_zero_or_elapsed_time_prevents_extra_calls(monkeypatch):
    state = feedback()
    state["remaining_seconds"] = 0
    assert ExperimentPlanner(FakeClient([])).propose(CATALOG, state).status == "exhausted"
    clock = iter([100, 100, 102])

    def next_time():
        return next(clock)

    monkeypatch.setattr("escape_lab.experiment_planner.monotonic", next_time)
    state["remaining_seconds"] = 1
    client = FakeClient([response()])
    result = ExperimentPlanner(client).propose(CATALOG, state)
    assert result.status == "exhausted"
    assert result.calls == 1 and result.candidate_id is None


def test_unknown_finish_reason_is_not_logged_verbatim():
    client = FakeClient([(response(), "SECRET-FINISH")])
    planner = ExperimentPlanner(client)
    assert planner.propose(CATALOG, feedback()).status == "proposed"
    assert planner.attempts[0]["finish_reason"] == "unknown"
    assert "SECRET-FINISH" not in json.dumps(planner.attempts)


@pytest.mark.parametrize("limit", [True, 0, 4, -1, 1.5])
def test_maximum_opportunities_is_strictly_bounded(limit):
    with pytest.raises(ValueError):
        ExperimentPlanner(FakeClient([]), max_proposals=limit)


def test_new_planner_shares_client_token_accounting_but_not_experiment_memory():
    client = FakeClient([response(), response()])
    first = ExperimentPlanner(client).propose(CATALOG, feedback())
    second = ExperimentPlanner(client).propose(CATALOG, feedback())
    assert first.candidate_id == second.candidate_id == NEXT
    assert second.usage == {"budget_charge_tokens": 31, "client_total_tokens": 62}
    sent = json.loads(client.requests[1]["messages"][1]["content"])
    assert sent["previous_proposals"] == []


def test_oversized_bad_tool_arguments_are_not_echoed_into_repair_context():
    invalid = response()
    invalid["tool_calls"][0]["function"]["arguments"] = "UNBOUNDED" * 2000
    client = FakeClient([invalid, response()])
    result = ExperimentPlanner(client).propose(CATALOG, feedback())
    assert result.status == "proposed"
    assert result.calls == 2
    assert "UNBOUNDED" not in json.dumps(client.requests[1])
