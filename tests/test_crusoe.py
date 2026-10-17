import json

import httpx
import pytest

from escape_lab.crusoe import CrusoeClient, CrusoeError


def completion(name="select_candidate", arguments=None):
    return {
        "model": "nvidia/test-nemotron",
        "choices": [{"message": {
            "role": "assistant",
            "content": None,
            "reasoning_content": "private reasoning must not be retained",
            "tool_calls": [{"id": "call-1", "type": "function", "function": {
                "name": name, "arguments": json.dumps(arguments or {})
            }}],
        }, "finish_reason": "tool_calls"}],
        "usage": {"prompt_tokens": 20, "completion_tokens": 10, "total_tokens": 30},
    }


def client_for(handler, **kwargs):
    return CrusoeClient("test-secret-do-not-log", transport=httpx.MockTransport(handler),
                        sleep=lambda _: None, **kwargs)


def test_discovery_preserves_exact_ids_and_capability_metadata():
    def handler(request):
        assert request.url.path == "/v1/models"
        return httpx.Response(200, json={"data": [
            {"id": "nvidia/Nemotron-Exact", "supported_parameters": ["tools", "max_tokens"]},
            {"id": "other/model"}, {"id": "Nemotron-Unknown"},
        ]})

    client = client_for(handler)
    models = client.discover_models()
    assert [model.id for model in models] == ["nvidia/Nemotron-Exact", "Nemotron-Unknown"]
    assert models[0].tools_supported is True
    assert models[1].tools_supported is None


def test_unknown_capabilities_require_explicit_smoke_and_retain_no_reasoning():
    requests = []

    def handler(request):
        requests.append(request)
        if request.method == "GET":
            return httpx.Response(200, json={"data": [{"id": "nvidia/test-nemotron"}]})
        body = json.loads(request.content)
        assert body["model"] == "nvidia/test-nemotron"
        return httpx.Response(200, json=completion("report_ready", {"ready": True}))

    client = client_for(handler, model_id="nvidia/test-nemotron")
    client.discover_models()
    with pytest.raises(CrusoeError, match="smoke"):
        client.chat([], [])
    assert len(requests) == 1
    assert client.run_smoke().passed is True
    assert client.total_tokens == 30
    serialized = json.dumps(client.usage_records)
    assert "private reasoning" not in serialized
    assert "test-secret" not in serialized


def test_transient_status_retries_are_bounded_and_auth_does_not_retry():
    attempts = []

    def handler(request):
        attempts.append(request)
        return httpx.Response(429, headers={"retry-after": "9999"}, json={"error": "x"})

    client = client_for(handler, max_retries=2)
    with pytest.raises(CrusoeError, match="429"):
        client.discover_models()
    assert len(attempts) == 3

    attempts.clear()

    def unauthorized(request):
        attempts.append(request)
        return httpx.Response(401, text="do not echo test-secret-do-not-log")

    client = client_for(unauthorized)
    with pytest.raises(CrusoeError) as error:
        client.discover_models()
    assert len(attempts) == 1
    assert "test-secret" not in str(error.value)


def test_explicit_unsupported_tools_fails_before_paid_request():
    calls = []

    def handler(request):
        calls.append(request)
        return httpx.Response(200, json={"data": [
            {"id": "nvidia/test-nemotron", "supported_parameters": ["max_tokens"]}
        ]})

    client = client_for(handler, model_id="nvidia/test-nemotron")
    client.discover_models()
    with pytest.raises(CrusoeError, match="tools"):
        client.run_smoke()
    assert len(calls) == 1


def test_request_and_total_token_caps_are_enforced():
    bodies = []

    def handler(request):
        if request.method == "GET":
            return httpx.Response(200, json={"data": [{
                "id": "nvidia/test-nemotron", "supported_parameters": ["tools"]
            }]})
        bodies.append(json.loads(request.content))
        return httpx.Response(200, json=completion())

    client = client_for(handler, model_id="nvidia/test-nemotron", max_tokens=64,
                        max_total_tokens=30)
    client.discover_models()
    result = client.chat([], [{"type": "function", "function": {"name": "x"}}])
    assert bodies[0]["max_tokens"] == 64
    assert "reasoning_content" not in result.message
    with pytest.raises(CrusoeError, match="budget"):
        client.chat([], [])
    assert len(bodies) == 1


def test_503_retries_then_recovers_without_echoing_provider_error():
    calls = []

    def handler(request):
        calls.append(request)
        if len(calls) == 1:
            return httpx.Response(503, text="temporary")
        return httpx.Response(200, json={"data": [{"id": "Nemotron-Recovered"}]})

    client = client_for(handler)
    assert client.discover_models()[0].id == "Nemotron-Recovered"
    assert len(calls) == 2


def test_missing_usage_is_not_reported_as_free():
    def handler(request):
        if request.method == "GET":
            return httpx.Response(200, json={"data": [{
                "id": "nvidia/test-nemotron", "supported_parameters": ["tools"]
            }]})
        result = completion()
        del result["usage"]
        return httpx.Response(200, json=result)

    client = client_for(handler, model_id="nvidia/test-nemotron")
    client.discover_models()
    client.chat([], [])
    assert client.total_tokens > 0
    assert client.usage_records[0]["usage_estimated"] is True
    assert client.usage_records[0]["total_tokens"] is None


def test_malformed_usage_metadata_is_conservatively_accounted():
    def handler(request):
        if request.method == "GET":
            return httpx.Response(200, json={"data": [{
                "id": "nvidia/test-nemotron", "supported_parameters": ["tools"]
            }]})
        result = completion()
        result["usage"] = "invalid-provider-metadata"
        return httpx.Response(200, json=result)

    client = client_for(handler, model_id="nvidia/test-nemotron")
    client.discover_models()
    client.chat([], [])
    assert client.total_tokens > 0
    assert client.usage_records[0]["usage_estimated"] is True


def test_ambiguous_failed_completion_consumes_estimated_budget():
    attempts = []

    def handler(request):
        if request.method == "GET":
            return httpx.Response(200, json={"data": [{
                "id": "nvidia/test-nemotron", "supported_parameters": ["tools"]
            }]})
        attempts.append(request)
        return httpx.Response(200, text="malformed JSON after possibly paid execution")

    client = client_for(handler, model_id="nvidia/test-nemotron", max_total_tokens=1)
    client.discover_models()
    with pytest.raises(CrusoeError, match="invalid JSON"):
        client.chat([], [])
    assert client.usage_records[0]["usage_estimated"] is True
    assert client.usage_records[0]["request_status"] == "failed_or_unknown"
    with pytest.raises(CrusoeError, match="token budget"):
        client.chat([], [])
    assert len(attempts) == 1
