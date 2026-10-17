import hashlib
import json
from dataclasses import asdict
from types import SimpleNamespace

import pytest

from escape_lab.crusoe import CrusoeError
from escape_lab.improvement_agent import (
    PROMPT_HASH,
    PROMPT_VERSION,
    TOOL_SCHEMA,
    PipelineProposer,
    ProposeArguments,
)
from escape_lab.prediction import PredictionConfig


def proposal_response(config=None, rationale="Improve cross-gene generalization."):
    if config is None:
        config = PredictionConfig().model_dump()
    return {"role": "assistant", "content": None, "tool_calls": [{
        "id": "proposal-call", "type": "function", "function": {
            "name": "propose_pipeline",
            "arguments": json.dumps({"config": config, "rationale": rationale}),
        },
    }]}


class FakeClient:
    model_id = "nvidia/test-nemotron"

    def __init__(self, responses):
        self.responses = responses
        self.requests = []

    def chat(self, messages, tools, *, prompt_version):
        self.requests.append({"messages": json.loads(json.dumps(messages)),
                              "tools": tools, "prompt_version": prompt_version})
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        if isinstance(response, tuple):
            response, finish_reason = response
        else:
            finish_reason = "tool_calls" if response.get("tool_calls") else "stop"
        return SimpleNamespace(message=response, finish_reason=finish_reason)


def test_valid_proposal_is_typed_and_logs_prompt_and_feedback_fingerprints():
    client = FakeClient([proposal_response()])
    proposer = PipelineProposer(client)
    result = proposer.propose({"objective": "macro_average_precision", "revision_index": 1,
                               "incumbent": {"config": PredictionConfig().model_dump(),
                                             "metrics": {"macro_average_precision": 0.5}}})
    assert isinstance(result.config, PredictionConfig)
    assert result.status == "proposed"
    assert result.error is None
    assert result.calls == 1
    assert result.model_id == client.model_id
    assert result.prompt_version == PROMPT_VERSION
    assert asdict(result)["status"] == "proposed"
    assert proposer.attempts[0]["prompt_hash"] == PROMPT_HASH
    assert len(proposer.attempts[0]["feedback_hash"]) == 64
    system = client.requests[0]["messages"][0]["content"]
    assert hashlib.sha256(system.encode()).hexdigest() == PROMPT_HASH


def test_bad_config_receives_one_repair_with_matching_tool_message():
    invalid = {**PredictionConfig().model_dump(), "model": "arbitrary-python"}
    client = FakeClient([proposal_response(invalid), proposal_response()])
    proposer = PipelineProposer(client)
    result = proposer.propose({"remaining_seconds": 200})
    assert result.status == "proposed"
    assert result.calls == 2
    messages = client.requests[1]["messages"]
    assert messages[-2]["role"] == "assistant"
    assert messages[-1]["role"] == "tool"
    assert messages[-1]["tool_call_id"] == "proposal-call"
    assert [attempt["status"] for attempt in proposer.attempts] == [
        "invalid_response", "proposed"]


def test_invalid_config_never_becomes_a_numerical_fallback():
    invalid = {**PredictionConfig().model_dump(), "split": "test"}
    client = FakeClient([proposal_response(invalid), proposal_response(invalid)])
    result = PipelineProposer(client).propose({})
    assert result.config is None
    assert result.status == "invalid_response"
    assert result.calls == 2
    assert result.error


def test_config_must_be_complete_even_if_prediction_model_has_defaults():
    client = FakeClient([proposal_response({"model": "logistic"}),
                         proposal_response({"model": "logistic"})])
    proposer = PipelineProposer(client)
    result = proposer.propose({})
    assert result.config is None
    parameters = client.requests[0]["tools"][0]["function"]["parameters"]
    assert set(parameters["properties"]["config"]["required"]) == set(
        PredictionConfig.model_fields
    )


def test_inlined_config_schema_preserves_exact_allowed_space():
    expected = ProposeArguments.model_json_schema()
    original_definition = expected.pop("$defs")["PredictionConfig"]
    original_definition["required"] = list(PredictionConfig.model_fields)
    expected["properties"]["config"] = original_definition
    assert TOOL_SCHEMA == expected
    assert TOOL_SCHEMA["properties"]["config"]["type"] == "object"
    assert TOOL_SCHEMA["properties"]["config"]["additionalProperties"] is False
    assert '"$ref"' not in json.dumps(TOOL_SCHEMA)
    assert "$defs" not in TOOL_SCHEMA


def test_stringified_config_is_still_rejected_without_silent_coercion():
    stringified = json.dumps(PredictionConfig().model_dump())
    client = FakeClient([proposal_response(stringified), proposal_response(stringified)])
    proposer = PipelineProposer(client)
    result = proposer.propose({})
    assert result.status == "invalid_response"
    assert result.config is None
    assert result.calls == 2
    assert [attempt["status"] for attempt in proposer.attempts] == [
        "invalid_response", "invalid_response"]


def test_transport_error_and_exhaustion_return_explicit_failures_without_substitution():
    client = FakeClient([CrusoeError("Crusoe returned HTTP 503.")])
    proposer = PipelineProposer(client, max_proposals=1)
    first = proposer.propose({})
    second = proposer.propose({})
    assert first.status == "failed"
    assert first.config is None
    assert first.calls == 1
    assert second.status == "exhausted"
    assert second.config is None
    assert second.calls == 0
    assert len(client.requests) == 1


def test_nested_feedback_only_exposes_approved_aggregate_fields():
    client = FakeClient([proposal_response()])
    feedback = {
        "objective": "macro_average_precision", "revision_index": 1, "seed": 23,
        "remaining_seconds": 180,
        "hidden_outcomes": "SECRET-CANDIDATE-DATA",
        "incumbent": {
            "config": PredictionConfig().model_dump(), "candidate_id": "SECRET-CANDIDATE-DATA",
            "metrics": {"macro_average_precision": 0.61,
                        "per_drug": {"pemigatinib": {"average_precision": 0.65,
                                                     "candidate_ids": ["SECRET-CANDIDATE-DATA"]},
                                     "futibatinib": {"average_precision": 0.57}},
                        "predictions": [0.9, 0.1]},
        },
        "history": [{"revision_index": 0, "accepted": False,
                     "raw_rows": "SECRET-CANDIDATE-DATA",
                     "metrics": {"macro_average_precision": 0.5}}],
        "error_analysis": {
            "by_gene": {"FGFR2": {"average_precision": 0.6, "n": 200,
                                    "replicates": "SECRET-CANDIDATE-DATA"},
                        "SECRET-CANDIDATE-DATA": {"average_precision": 0.99}},
            "by_fold": {"0": {"macro_average_precision": 0.6,
                                 "raw_rows": "SECRET-CANDIDATE-DATA"}},
            "candidate_id": "SECRET-CANDIDATE-DATA",
        },
        "allowed_config": {"exec": "SECRET-CANDIDATE-DATA", "split": "test"},
    }
    proposer = PipelineProposer(client)
    result = proposer.propose(feedback)
    assert result.config is not None
    serialized = json.dumps(client.requests)
    assert "SECRET-CANDIDATE-DATA" not in serialized
    sent = json.loads(client.requests[0]["messages"][1]["content"])
    assert sent["incumbent"]["metrics"]["per_drug"]["pemigatinib"]["average_precision"] == 0.65
    assert sent["error_analysis"]["by_gene"]["FGFR2"]["n"] == 200
    assert "exec" not in sent["allowed_config"]
    assert sent["seed"] == 23


def test_no_tool_call_is_repaired_once_then_returns_failure():
    client = FakeClient([{"role": "assistant", "content": "I refuse to select a config."},
                         {"role": "assistant", "content": "Still no structured proposal."}])
    result = PipelineProposer(client).propose({})
    assert result.config is None
    assert result.status == "invalid_response"
    assert result.calls == 2
    assert len(client.requests) == 2


def test_exact_controller_aggregate_layout_and_invalid_nested_types():
    client = FakeClient([proposal_response()])
    metrics = {"macro_ap": 0.62, "per_drug": {
        "pem": {"average_precision": 0.64, "brier": 0.12, "n": 300, "positive_count": 60},
        "fut": {"average_precision": 0.60, "brier": 0.14, "n": 300, "positive_count": 50},
    }}
    proposer = PipelineProposer(client)
    result = proposer.propose({
        "incumbent": {"config": PredictionConfig().model_dump(), "metrics": metrics,
                      "status": {"candidate_id": "DO-NOT-SEND"}},
        "error_analysis": {
            "by_gene": {"FGFR2": metrics["per_drug"]},
            "by_fold": {"development-fold-0": metrics,
                        "final-fold-0": {"macro_ap": 0.99}},
        },
    })
    assert result.config is not None
    sent = json.loads(client.requests[0]["messages"][1]["content"])
    assert sent["incumbent"]["metrics"] == metrics
    assert "status" not in sent["incumbent"]
    assert sent["error_analysis"]["by_gene"]["FGFR2"]["pem"]["positive_count"] == 60
    assert sent["error_analysis"]["by_fold"]["development-fold-0"] == metrics
    assert "final-fold-0" not in sent["error_analysis"]["by_fold"]
    assert "DO-NOT-SEND" not in json.dumps(client.requests)


def test_oversized_invalid_tool_arguments_do_not_trigger_an_unbounded_repair():
    message = proposal_response()
    message["tool_calls"][0]["function"]["arguments"] = "invalid-json" * 10_000
    client = FakeClient([message])
    result = PipelineProposer(client).propose({})
    assert result.config is None
    assert result.status == "invalid_response"
    assert result.calls == 1
    assert "context" in result.error


def test_elapsed_feedback_deadline_prevents_a_second_completion(monkeypatch):
    clock = iter([100, 100, 102])

    def next_time():
        return next(clock)

    monkeypatch.setattr("escape_lab.improvement_agent.monotonic", next_time)
    client = FakeClient([proposal_response({"model": "logistic"})])
    result = PipelineProposer(client).propose({"remaining_seconds": 1})
    assert result.config is None
    assert result.status == "exhausted"
    assert result.calls == 1
    assert len(client.requests) == 1


@pytest.mark.parametrize(("message", "category"), [
    ("Crusoe transport failed (ReadTimeout).", "ReadTimeout"),
    ("Crusoe transport failed (ConnectTimeout).", "ConnectTimeout"),
    ("Crusoe returned HTTP 503.", "HTTP503"),
    ("Crusoe returned HTTP 401.", "HTTP401"),
    ("Crusoe total token budget exhausted.", "TokenBudgetExhausted"),
])
def test_operational_failures_preserve_only_safe_categories(message, category):
    client = FakeClient([CrusoeError(message)])
    proposer = PipelineProposer(client)
    result = proposer.propose({})
    assert result.status == "failed"
    assert result.config is None
    assert result.calls == 1
    assert result.error == f"Crusoe request failed ({category}); no pipeline was proposed."
    assert proposer.attempts[0]["error_category"] == category


@pytest.mark.parametrize("message", [
    "provider body SECRET-CREDENTIAL-123 and arbitrary private details",
    "Crusoe transport failed (SECRET-CREDENTIAL-123).",
    "Crusoe returned HTTP 503. SECRET-CREDENTIAL-123",
    "Crusoe transport failed (ReadTimeout).\nSECRET-CREDENTIAL-123",
])
def test_unknown_error_messages_never_escape_into_result_or_attempt_log(message):
    client = FakeClient([CrusoeError(message)])
    proposer = PipelineProposer(client)
    result = proposer.propose({})
    assert result.error == "Crusoe request failed (CrusoeError); no pipeline was proposed."
    assert proposer.attempts[0]["error_category"] == "CrusoeError"
    serialized = json.dumps({"result": asdict(result), "attempts": proposer.attempts})
    assert "SECRET-CREDENTIAL-123" not in serialized


def test_token_truncation_is_explicit_and_never_retried_as_schema_repair():
    message = proposal_response()
    message["content"] = "PRIVATE-CONTENT-MUST-NOT-BE-LOGGED"
    message["reasoning_content"] = "PRIVATE-REASONING-MUST-NOT-BE-LOGGED"
    client = FakeClient([(message, "length")])
    proposer = PipelineProposer(client)
    result = proposer.propose({})
    assert result.status == "truncated_response"
    assert result.config is None
    assert result.calls == 1
    assert "completion token limit" in result.error
    assert len(client.requests) == 1
    assert proposer.attempts[0]["status"] == "truncated_response"
    assert proposer.attempts[0]["finish_reason"] == "length"
    assert "config" not in proposer.attempts[0]
    serialized = json.dumps({"result": asdict(result), "attempts": proposer.attempts})
    assert "PRIVATE-CONTENT" not in serialized
    assert "PRIVATE-REASONING" not in serialized


@pytest.mark.parametrize("finish_reason", ["stop", "tool_calls", "function_call", "content_filter"])
def test_valid_response_logs_allowlisted_finish_reason(finish_reason):
    client = FakeClient([(proposal_response(), finish_reason)])
    proposer = PipelineProposer(client)
    result = proposer.propose({})
    assert result.status == "proposed"
    assert proposer.attempts[0]["finish_reason"] == finish_reason


def test_unknown_finish_reason_is_not_echoed_in_logs_or_results():
    client = FakeClient([(proposal_response(), "PRIVATE-FINISH-REASON")])
    proposer = PipelineProposer(client)
    result = proposer.propose({})
    assert result.status == "proposed"
    assert proposer.attempts[0]["finish_reason"] == "unknown"
    assert "PRIVATE-FINISH-REASON" not in json.dumps({
        "result": {**asdict(result), "config": result.config.model_dump()},
        "attempts": proposer.attempts})
