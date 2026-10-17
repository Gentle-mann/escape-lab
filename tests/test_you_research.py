import json

import httpx
import pytest

from escape_lab.you_research import (
    YouResearchClient,
    YouResearchError,
    read_api_key,
    run_contents_phase,
    run_search_phase,
    url_rejection_reason,
)

SECRET = "test-private-you-api-key"
URL = "https://pubmed.ncbi.nlm.nih.gov/123456/"


def plan_file(tmp_path, queries=None):
    path = tmp_path / "requested-plan.json"
    path.write_text(json.dumps({
        "queries": queries or ["FGFR kinase molecular brake structural mechanism"],
        "max_contents_pages": 2,
        "purpose": "evaluator-owned research quarantine",
    }))
    return path


def search_response(url=URL, **kwargs):
    return {"results": {"web": [{
        "url": url, "title": "Kinase structural mechanism", "description": "General evidence",
        "snippets": ["A structural mechanism."], **kwargs,
    }]}}


def create_client(handler):
    return YouResearchClient(SECRET, transport=httpx.MockTransport(handler))


def test_request_contract_plan_hash_and_quarantine(tmp_path, caplog):
    requests = []

    def handler(request):
        requests.append(request)
        assert request.method == "POST"
        assert str(request.url) == "https://ydc-index.io/v1/search"
        assert request.headers["X-API-Key"] == SECRET
        assert json.loads(request.content) == {
            "query": "FGFR kinase molecular brake structural mechanism", "count": 5,
            "include_domains": ["pubmed.ncbi.nlm.nih.gov", "pmc.ncbi.nlm.nih.gov", "rcsb.org"],
            "freshness": "1900-01-01to2023-12-31", "strict_operators": True, "language": "EN",
        }
        return httpx.Response(200, json=search_response())

    plan = plan_file(tmp_path)
    output = tmp_path / "research"
    summary = run_search_phase(plan, output, client=create_client(handler))
    assert summary["eligible_urls"] == [URL]
    assert summary["agent_exposure"] == "none; evaluator quarantine only"
    assert (output / "plan.json").read_bytes() == plan.read_bytes()
    assert len((output / "plan.sha256").read_text().strip()) == 64
    assert json.loads((output / "quarantine/search-01.json").read_text()) == search_response()
    usage = json.loads((output / "usage.json").read_text())
    assert usage["request_count"] == 1
    assert usage["estimated_cost_usd"] == 0.005
    assert usage["elapsed_seconds"] >= 0
    assert len(requests) == 1
    assert SECRET not in caplog.text
    assert SECRET not in "".join(path.read_text() for path in output.rglob("*.json"))


def test_credentials_use_process_then_dotenv_alias_without_printing(tmp_path, capsys):
    env = tmp_path / ".env"
    env.write_text("YDC_API_KEY=dotenv-secret\nYOU_API_KEY=secondary-secret\n")
    assert read_api_key(env, environ={"YDC_API_KEY": "process-secret"}) == "process-secret"
    assert read_api_key(env, environ={"YOU_API_KEY": "process-alias"}) == "process-alias"
    assert read_api_key(env, environ={}) == "dotenv-secret"
    env.write_text("YOU_API_KEY=alias-secret\n")
    assert read_api_key(env, environ={}) == "alias-secret"
    env.write_text("UNRELATED=not-a-key\n")
    with pytest.raises(YouResearchError, match="YDC_API_KEY"):
        read_api_key(env, environ={})
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("url", [
    "http://pubmed.ncbi.nlm.nih.gov/123/",
    "https://pubmed.ncbi.nlm.nih.gov.evil.org/123/",
    "https://evil.org/?url=https://pubmed.ncbi.nlm.nih.gov/123/",
    "https://user@pubmed.ncbi.nlm.nih.gov/123/",
    "https://pubmed.ncbi.nlm.nih.gov:8443/123/",
    "https://pubmed.ncbi.nlm.nih.gov./123/",
    "https://pmc.ncbi.nlm.nih.gov/articles/123/supplementary/",
    "https://pmc.ncbi.nlm.nih.gov/articles/123/bin/file.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/123/suppl_file/a.pdf",
    "https://pmc.ncbi.nlm.nih.gov/articles/123/%2573upplementary/a.pdf",
    "https://rcsb.org/structure/a?redirect=RLY-4008",
    "https://rcsb.org/structure/a#data",
    "https://rcsb.org/structure/s41587-026-03090-8",
    "https://rcsb.org/structure/s41588-025-02431-8",
    "https://rcsb.org/structure/saifudeen",
    "https://rcsb.org/structure/External%20S13",
])
def test_url_gate_rejects_scope_bypasses(url):
    assert url_rejection_reason(url)


@pytest.mark.parametrize("metadata", [
    {"title": "Saifudeen and colleagues"},
    {"doi": "10.1038/s41587-026-03090-8"},
    {"description": "RLY-4008 evidence"},
    {"metadata": {"doi": "s41588-025-02431-8"}},
    {"snippets": ["external S13"]},
])
def test_search_metadata_denied_before_contents(tmp_path, metadata):
    output = tmp_path / "research"
    client = create_client(lambda _: httpx.Response(200, json=search_response(**metadata)))
    summary = run_search_phase(plan_file(tmp_path), output, client=client)
    assert summary["eligible_urls"] == []
    assert summary["results"][0]["rejection_reason"]
    with pytest.raises(YouResearchError, match="eligible"):
        run_contents_phase(output, [URL], client=client)
    assert not (output / "contents").exists()


@pytest.mark.parametrize("queries", [
    ["query"] * 5, ["FGFR2 N549K activity"], ["FGFR2 p.Asn549Lys evidence"],
    ["measurements 0.42"], ["reports/revision-v2/results.json"], ["external S13"],
])
def test_plan_scope_is_rejected_before_network_or_directory_creation(tmp_path, queries):
    client = create_client(lambda _: pytest.fail("Network must not be used"))
    output = tmp_path / "research"
    with pytest.raises(YouResearchError):
        run_search_phase(plan_file(tmp_path, queries), output, client=client)
    assert not output.exists()


def test_http_failure_is_sanitized_recorded_and_not_retried(tmp_path):
    calls = []

    def handler(request):
        calls.append(request)
        return httpx.Response(503, text=f"provider diagnostic contains {SECRET}")

    output = tmp_path / "research"
    summary = run_search_phase(plan_file(tmp_path), output, client=create_client(handler))
    assert len(calls) == 1
    assert summary["errors"] == ["You.com returned HTTP 503."]
    assert summary["status"] == "failed"
    assert SECRET not in (output / "usage.json").read_text()
    assert SECRET not in (output / "search-summary.json").read_text()
    assert json.loads((output / "usage.json").read_text())["request_count"] == 1


def test_transport_failure_and_bad_json_never_echo_provider_data():
    def fail(request):
        raise httpx.ConnectError(SECRET, request=request)

    with pytest.raises(YouResearchError) as error:
        create_client(fail).search("FGFR mechanism")
    assert SECRET not in str(error.value)
    with pytest.raises(YouResearchError, match="invalid JSON"):
        create_client(lambda _: httpx.Response(200, text=SECRET)).search("FGFR mechanism")


def test_contents_only_fetches_observed_urls_and_preserves_nullable_page_failure(tmp_path):
    requests = []

    def handler(request):
        requests.append(request)
        if request.url.path.endswith("search"):
            return httpx.Response(200, json=search_response())
        assert json.loads(request.content) == {
            "urls": [URL], "formats": ["markdown", "metadata"], "crawl_timeout": 10,
        }
        return httpx.Response(200, json=[{"url": URL, "markdown": None}])

    client = create_client(handler)
    output = tmp_path / "research"
    run_search_phase(plan_file(tmp_path), output, client=client)
    with pytest.raises(YouResearchError, match="eligible"):
        run_contents_phase(output, ["https://pubmed.ncbi.nlm.nih.gov/999999/"], client=client)
    summary = run_contents_phase(output, [URL], client=client)
    assert summary["pages"][0]["status"] == "missing_markdown"
    assert summary["usable_page_count"] == 0
    assert summary["status"] == "failed"
    assert len(requests) == 2
    usage = json.loads((output / "contents/usage.json").read_text())
    assert usage["request_count"] == 1
    assert usage["estimated_cost_usd"] == 0.001


def test_contents_rechecks_returned_url_and_metadata(tmp_path):
    def handler(request):
        if request.url.path.endswith("search"):
            return httpx.Response(200, json=search_response())
        return httpx.Response(200, json=[{
            "url": "https://evil.org/redirected", "markdown": "Untrusted text",
        }])

    client = create_client(handler)
    output = tmp_path / "research"
    run_search_phase(plan_file(tmp_path), output, client=client)
    summary = run_contents_phase(output, [URL], client=client)
    assert summary["usable_page_count"] == 0
    assert summary["pages"][0]["status"] == "rejected"


def test_no_overwriting_attempt_artifacts_and_tampering_fails_closed(tmp_path):
    client = create_client(lambda _: httpx.Response(200, json=search_response()))
    plan = plan_file(tmp_path)
    output = tmp_path / "research"
    run_search_phase(plan, output, client=client)
    before = {str(p): p.read_bytes() for p in output.rglob("*") if p.is_file()}
    with pytest.raises(FileExistsError):
        run_search_phase(plan, output, client=client)
    assert before == {str(p): p.read_bytes() for p in output.rglob("*") if p.is_file()}
    (output / "plan.json").write_text('{"queries":["tampered"]}')
    with pytest.raises(YouResearchError, match="integrity"):
        run_contents_phase(output, [URL], client=client)


def test_successful_contents_is_still_evaluator_only_and_cannot_be_overwritten(tmp_path):
    def handler(request):
        if request.url.path.endswith("search"):
            return httpx.Response(200, json=search_response())
        return httpx.Response(200, json=[{"url": URL, "markdown": "General mechanism."}])

    client = create_client(handler)
    output = tmp_path / "research"
    run_search_phase(plan_file(tmp_path), output, client=client)
    summary = run_contents_phase(output, [URL], client=client)
    assert summary["usable_page_count"] == 1
    assert "markdown" not in summary["pages"][0]
    assert summary["agent_exposure"] == "none; evaluator quarantine only"
    with pytest.raises(FileExistsError):
        run_contents_phase(output, [URL], client=client)


def test_incorrect_search_shape_is_recorded_as_an_error(tmp_path):
    output = tmp_path / "research"
    summary = run_search_phase(plan_file(tmp_path), output, client=create_client(
        lambda _: httpx.Response(200, json={"results": {"web": "not a list"}})))
    assert summary["status"] == "failed"
    assert summary["eligible_urls"] == []
    assert (output / "quarantine/search-01.json").exists()


def test_provider_cannot_reflect_authentication_secret_into_saved_artifacts(tmp_path):
    output = tmp_path / "research"
    run_search_phase(plan_file(tmp_path), output, client=create_client(
        lambda _: httpx.Response(200, json=search_response(title=SECRET))))
    assert SECRET not in "".join(path.read_text() for path in output.rglob("*.json"))


def test_contents_scope_comes_from_raw_searches_not_an_edited_eligible_list(tmp_path):
    output = tmp_path / "research"
    client = create_client(lambda _: httpx.Response(200, json=search_response()))
    run_search_phase(plan_file(tmp_path), output, client=client)
    summary_path = output / "search-summary.json"
    summary = json.loads(summary_path.read_text())
    injected = "https://pubmed.ncbi.nlm.nih.gov/9999/"
    summary["eligible_urls"].append(injected)
    summary_path.write_text(json.dumps(summary))
    with pytest.raises(YouResearchError, match="eligible"):
        run_contents_phase(output, [injected], client=client)


def test_contents_denies_reserved_mentions_in_returned_page(tmp_path):
    def handler(request):
        if request.url.path.endswith("search"):
            return httpx.Response(200, json=search_response())
        return httpx.Response(200, json=[{
            "url": URL, "markdown": "Reserved RLY-4008 data must stay quarantined.",
        }])

    output = tmp_path / "research"
    client = create_client(handler)
    run_search_phase(plan_file(tmp_path), output, client=client)
    summary = run_contents_phase(output, [URL], client=client)
    assert summary["usable_page_count"] == 0
    assert summary["pages"][0]["status"] == "rejected"


def test_no_credential_and_missing_key_cannot_create_attempt(tmp_path, monkeypatch):
    monkeypatch.delenv("YDC_API_KEY", raising=False)
    monkeypatch.delenv("YOU_API_KEY", raising=False)
    output = tmp_path / "research"
    with pytest.raises(YouResearchError, match="requires"):
        run_search_phase(plan_file(tmp_path), output)
    assert not output.exists()
