"""Bounded You.com literature retrieval into an evaluator-owned quarantine.

Nothing in this module creates an agent evidence packet. Search snippets, page
contents, and normalized summaries remain untrusted evaluator-only material.
Domain/date controls and local deny checks are defense in depth, not a guarantee
that a page cannot contain outcomes. An evaluator must approve evidence separately.

The caller freezes one to four general queries in a JSON plan before retrieval.
The original plan bytes, their SHA256, raw responses, usage, and errors are saved
without overwriting earlier attempts. No retry, inline crawling, or import-time
network access occurs. Price estimates count attempted requests/pages; invoices
can differ, especially after transport failures or page-level failures.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import time
from collections.abc import Collection, Mapping
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Self
from urllib.parse import unquote, urlsplit

import httpx
from dotenv import dotenv_values

ALLOWED_HOSTS = ("pubmed.ncbi.nlm.nih.gov", "pmc.ncbi.nlm.nih.gov", "rcsb.org")
FRESHNESS = "1900-01-01to2023-12-31"
AGENT_EXPOSURE = "none; evaluator quarantine only"
PRICE_BASIS = {
    "search_usd_per_attempt": 0.005,
    "contents_usd_per_requested_page": 0.001,
    "source": "https://you.com/docs/administration/billing",
    "kind": "estimate; attempted requests/pages; not an invoice",
}
_DENIED_TEXT = re.compile(
    r"s41587[-\s]?026[-\s]?03090[-\s]?8|s41588[-\s]?025[-\s]?02431[-\s]?8"
    r"|saifudeen|rly[-\s]?4008|external[\s_-]*s13", re.IGNORECASE,
)
_SUPPLEMENT_PATH = re.compile(
    r"(?:supplement|suppl[/.\-_]|/bin/|\.(?:xlsx?|csv|tsv|zip)(?:[?#]|$))", re.IGNORECASE,
)


class YouResearchError(RuntimeError):
    """Sanitized operational or scope error without credentials/provider bodies."""


def read_api_key(
    env_path: Path | None = None, *, environ: Mapping[str, str] | None = None,
) -> str:
    """Read process YDC_API_KEY/YOU_API_KEY, then an explicitly selected .env."""
    environment = os.environ if environ is None else environ
    sources = [environment]
    if env_path is not None:
        # Do not expand other environment variables into a credential value.
        sources.append(dotenv_values(env_path, interpolate=False))
    for source in sources:
        for name in ("YDC_API_KEY", "YOU_API_KEY"):
            value = source.get(name)
            if isinstance(value, str) and value.strip():
                return value.strip()
    raise YouResearchError("You.com requires YDC_API_KEY or YOU_API_KEY.")


def _decoded(value: str) -> str:
    for _ in range(5):
        decoded = unquote(value)
        if decoded == value:
            return value
        value = decoded
    return value


def url_rejection_reason(url: object, metadata: Any = None) -> str | None:
    """Return a conservative local scope rejection, never a non-leakage claim."""
    if not isinstance(url, str) or not url or any(ord(c) < 33 for c in url):
        return "invalid_url"
    if "\\" in url:
        return "invalid_url"
    try:
        parsed = urlsplit(url)
        if parsed.scheme != "https" or parsed.hostname not in ALLOWED_HOSTS:
            return "outside_https_host_allowlist"
        if parsed.username is not None or parsed.password is not None or parsed.port is not None:
            return "userinfo_or_port_disallowed"
        if parsed.fragment:
            return "fragment_disallowed"
    except (TypeError, ValueError):
        return "invalid_url"
    decoded = _decoded(url)
    if "%" in decoded:
        return "ambiguous_url_encoding"
    if _DENIED_TEXT.search(decoded):
        return "blocked_study_or_reserved_panel"
    if _SUPPLEMENT_PATH.search(decoded):
        return "supplement_or_data_attachment"
    if metadata is not None and _DENIED_TEXT.search(json.dumps(metadata, ensure_ascii=False)):
        return "blocked_study_or_reserved_panel_metadata"
    return None


def _validate_query(query: object) -> str:
    if not isinstance(query, str) or not query.strip() or len(query) > 600:
        raise YouResearchError("Research queries must be nonempty general text of at most 600 chars.")
    if _DENIED_TEXT.search(query):
        raise YouResearchError("Research query refers to a blocked study or reserved panel.")
    without_gene_names = re.sub(r"\bFGFR[1-4]\b", "FGFR", query, flags=re.IGNORECASE)
    if re.search(r"\d|https?://|[/\\]|\b(?:holdout|outcomes?|results?|reports?|scores?|"
                 r"predictions?|test|measurements?)\b", without_gene_names, re.IGNORECASE):
        raise YouResearchError("Research queries cannot contain variants, values, reports or outcomes.")
    return query


def _validate_plan(plan: Any) -> list[str]:
    if not isinstance(plan, dict) or not isinstance(plan.get("queries"), list):
        raise YouResearchError("Research plan requires a queries list.")
    queries = plan["queries"]
    if not 1 <= len(queries) <= 4:
        raise YouResearchError("Research plan permits one to four predefined queries.")
    if len({str(query) for query in queries}) != len(queries):
        raise YouResearchError("Research queries must be distinct.")
    pages = plan.get("max_contents_pages", 4)
    if type(pages) is not int or not 1 <= pages <= 10:
        raise YouResearchError("max_contents_pages must be between one and ten.")
    return [_validate_query(query) for query in queries]


def _search_rows(payload: Any) -> list[dict]:
    if not isinstance(payload, dict) or not isinstance(payload.get("results"), dict):
        raise YouResearchError("You.com returned an invalid search response shape.")
    rows = payload["results"].get("web")
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise YouResearchError("You.com returned an invalid search result list.")
    return rows


class YouResearchClient:
    """Exact official API adapter with injectable transport and no retries."""

    def __init__(
        self, api_key: str, *, transport: httpx.BaseTransport | None = None, timeout: float = 30,
    ):
        if not isinstance(api_key, str) or not api_key.strip():
            raise YouResearchError("You.com requires a nonempty API key.")
        self._api_key = api_key
        self._eligible_urls: set[str] = set()
        self.requests: list[dict] = []
        self.last_payload: Any = None
        self._http = httpx.Client(
            base_url="https://ydc-index.io/v1/",
            headers={"X-API-Key": api_key, "Content-Type": "application/json"},
            transport=transport, timeout=timeout, follow_redirects=False,
        )

    def close(self) -> None:
        self._http.close()

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _request(self, endpoint: str, payload: dict) -> Any:
        self.last_payload = None
        start = time.monotonic()
        record = {
            "endpoint": endpoint, "status": "failed_or_unknown",
            "requested_pages": len(payload.get("urls", [])),
            "started_utc": datetime.now(UTC).isoformat(),
        }
        self.requests.append(record)
        try:
            try:
                response = self._http.post(endpoint, json=payload)
            except httpx.HTTPError as error:
                raise YouResearchError(
                    f"You.com transport failed ({type(error).__name__})."
                ) from None
            record["http_status"] = response.status_code
            if response.status_code >= 300:
                raise YouResearchError(f"You.com returned HTTP {response.status_code}.")
            try:
                result = response.json()
            except ValueError:
                raise YouResearchError("You.com returned invalid JSON.") from None
            # A provider must not reflect a credential into an otherwise raw archive.
            serialized = json.dumps(result, ensure_ascii=False)
            redacted = serialized.replace(self._api_key, "[REDACTED_API_KEY]")
            record["credential_reflection_redacted"] = redacted != serialized
            record["status"] = "completed"
            self.last_payload = json.loads(redacted)
            return self.last_payload
        except YouResearchError as error:
            record["error"] = str(error)
            raise
        finally:
            record["elapsed_seconds"] = time.monotonic() - start

    def search(self, query: str) -> dict:
        _validate_query(query)
        result = self._request("search", {
            "query": query, "count": 5, "include_domains": list(ALLOWED_HOSTS),
            "freshness": FRESHNESS, "strict_operators": True, "language": "EN",
        })
        for row in _search_rows(result):
            if url_rejection_reason(row.get("url"), row) is None:
                self._eligible_urls.add(row["url"])
        return result

    def contents(self, urls: list[str], *, eligible_urls: Collection[str] | None = None) -> list:
        """Fetch only returned eligible URLs, supplied from a verified search archive if new."""
        observed = self._eligible_urls if eligible_urls is None else set(eligible_urls)
        if not urls or len(urls) > 10 or len(set(urls)) != len(urls):
            raise YouResearchError("Contents requires one to ten distinct eligible URLs.")
        if any(url not in observed or url_rejection_reason(url) for url in urls):
            raise YouResearchError("Contents URL was not an eligible observed search result.")
        result = self._request("contents", {
            "urls": urls, "formats": ["markdown", "metadata"], "crawl_timeout": 10,
        })
        if not isinstance(result, list) or any(not isinstance(row, dict) for row in result):
            raise YouResearchError("You.com returned an invalid contents response shape.")
        return result


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_bytes(path: Path, data: bytes) -> None:
    with path.open("xb") as handle:
        handle.write(data)


def _write_json(path: Path, value: Any) -> str:
    data = (json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n").encode()
    _write_bytes(path, data)
    return _sha(data)


def _usage(client: YouResearchClient, start_index: int, start: float) -> dict:
    records = client.requests[start_index:]
    searches = sum(row["endpoint"] == "search" for row in records)
    pages = sum(row["requested_pages"] for row in records)
    return {
        "request_count": len(records), "search_requests": searches,
        "contents_requests": len(records) - searches, "requested_contents_pages": pages,
        "elapsed_seconds": time.monotonic() - start,
        "estimated_cost_usd": round(searches * 0.005 + pages * 0.001, 6),
        "price_basis": PRICE_BASIS, "automatic_retries": 0, "requests": records,
    }


def run_search_phase(
    plan_path: Path, output_dir: Path, *, client: YouResearchClient | None = None,
    env_path: Path | None = None,
) -> dict:
    """Freeze plan and quarantine bounded searches in a new exclusive directory."""
    plan_bytes = Path(plan_path).read_bytes()
    try:
        plan = json.loads(plan_bytes)
    except (ValueError, UnicodeDecodeError):
        raise YouResearchError("Research plan is not valid JSON.") from None
    queries = _validate_plan(plan)
    output_dir = Path(output_dir)
    owns_client = client is None
    client = client or YouResearchClient(read_api_key(env_path))
    start = time.monotonic()
    start_index = len(client.requests)
    try:
        output_dir.mkdir(parents=True, exist_ok=False)
        (output_dir / "quarantine").mkdir()
        _write_bytes(output_dir / "plan.json", plan_bytes)
        _write_bytes(output_dir / "plan.sha256", (_sha(plan_bytes) + "\n").encode())
        summary = {
            "agent_exposure": AGENT_EXPOSURE, "plan_sha256": _sha(plan_bytes),
            "queries": queries, "results": [], "eligible_urls": [], "errors": [],
            "raw_files": {}, "status": "completed",
            "limitation": "Local filters are defense in depth, not an outcome-leakage guarantee.",
        }
        successful_queries = 0
        for index, query in enumerate(queries, start=1):
            raw_name = f"quarantine/search-{index:02d}.json"
            try:
                raw = client.search(query)
                summary["raw_files"][raw_name] = _write_json(output_dir / raw_name, raw)
                successful_queries += 1
                for row in _search_rows(raw):
                    reason = url_rejection_reason(row.get("url"), row)
                    summary["results"].append({
                        "query_index": index, "url": row.get("url"),
                        "title": row.get("title"), "description": row.get("description"),
                        "snippets": row.get("snippets"), "metadata": row.get("metadata"),
                        "rejection_reason": reason,
                    })
                    if reason is None and row["url"] not in summary["eligible_urls"]:
                        summary["eligible_urls"].append(row["url"])
            except YouResearchError as error:
                summary["errors"].append(str(error))
                if client.last_payload is not None and raw_name not in summary["raw_files"]:
                    summary["raw_files"][raw_name] = _write_json(
                        output_dir / raw_name, client.last_payload,
                    )
        if summary["errors"]:
            summary["status"] = "partial" if successful_queries else "failed"
        _write_json(output_dir / "search-summary.json", summary)
        _write_json(output_dir / "usage.json", _usage(client, start_index, start))
        return summary
    finally:
        if owns_client:
            client.close()


def _verified_search_archive(research_dir: Path) -> tuple[dict, set[str]]:
    try:
        plan_bytes = (research_dir / "plan.json").read_bytes()
        expected_hash = (research_dir / "plan.sha256").read_text().strip()
        plan = json.loads(plan_bytes)
        summary = json.loads((research_dir / "search-summary.json").read_bytes())
        if _sha(plan_bytes) != expected_hash or summary["plan_sha256"] != expected_hash:
            raise ValueError
        _validate_plan(plan)
        eligible: set[str] = set()
        for filename, expected in summary["raw_files"].items():
            if re.fullmatch(r"quarantine/search-0[1-4]\.json", filename) is None:
                raise ValueError
            data = (research_dir / filename).read_bytes()
            if _sha(data) != expected:
                raise ValueError
            for row in _search_rows(json.loads(data)):
                if url_rejection_reason(row.get("url"), row) is None:
                    eligible.add(row["url"])
        return plan, eligible
    except (OSError, ValueError, KeyError, TypeError, YouResearchError):
        raise YouResearchError("Research archive integrity verification failed.") from None


def run_contents_phase(
    research_dir: Path, urls: list[str], *, client: YouResearchClient | None = None,
    env_path: Path | None = None,
) -> dict:
    """Quarantine chosen, verified search URLs once; page failure is not success."""
    research_dir = Path(research_dir)
    plan, eligible = _verified_search_archive(research_dir)
    if (not urls or len(urls) > plan.get("max_contents_pages", 4)
            or len(set(urls)) != len(urls)):
        raise YouResearchError("Contents selection exceeds the distinct eligible page budget.")
    if any(url not in eligible or url_rejection_reason(url) for url in urls):
        raise YouResearchError("Contents URL was not an eligible observed search result.")
    output_dir = research_dir / "contents"
    owns_client = client is None
    client = client or YouResearchClient(read_api_key(env_path))
    start = time.monotonic()
    start_index = len(client.requests)
    try:
        output_dir.mkdir(exist_ok=False)
        (output_dir / "quarantine").mkdir()
        _write_json(output_dir / "selected-urls.json", urls)
        summary = {
            "agent_exposure": AGENT_EXPOSURE, "pages": [], "errors": [],
            "usable_page_count": 0, "status": "completed", "raw_files": {},
        }
        try:
            raw = client.contents(urls, eligible_urls=eligible)
            raw_name = "quarantine/response.json"
            summary["raw_files"][raw_name] = _write_json(output_dir / raw_name, raw)
            returned: set[str] = set()
            for row in raw:
                url = row.get("url")
                reason = url_rejection_reason(url, row)
                if reason is None and (url not in urls or url in returned):
                    reason = "unexpected_or_duplicate_returned_url"
                if isinstance(url, str):
                    returned.add(url)
                markdown = row.get("markdown")
                if reason:
                    status = "rejected"
                elif not isinstance(markdown, str) or not markdown.strip():
                    status = "missing_markdown"
                else:
                    status = "quarantined_for_evaluator_review"
                    summary["usable_page_count"] += 1
                summary["pages"].append({
                    "url": url, "status": status, "rejection_reason": reason,
                    "markdown_characters": len(markdown) if isinstance(markdown, str) else 0,
                })
            for missing in urls:
                if missing not in returned:
                    summary["pages"].append({"url": missing, "status": "missing_response"})
        except YouResearchError as error:
            summary["errors"].append(str(error))
            if client.last_payload is not None and not summary["raw_files"]:
                summary["raw_files"]["quarantine/response.json"] = _write_json(
                    output_dir / "quarantine/response.json", client.last_payload,
                )
        if summary["usable_page_count"] != len(urls):
            summary["status"] = "partial" if summary["usable_page_count"] else "failed"
        _write_json(output_dir / "summary.json", summary)
        _write_json(output_dir / "usage.json", _usage(client, start_index, start))
        return summary
    finally:
        if owns_client:
            client.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    search = subparsers.add_parser("search")
    search.add_argument("--plan", required=True, type=Path)
    search.add_argument("--output", required=True, type=Path)
    search.add_argument("--env", type=Path, default=Path(".env"))
    contents = subparsers.add_parser("contents")
    contents.add_argument("--research-dir", required=True, type=Path)
    contents.add_argument("--urls", required=True, type=Path, help="JSON file containing URL list")
    contents.add_argument("--env", type=Path, default=Path(".env"))
    args = parser.parse_args()
    if args.command == "search":
        summary = run_search_phase(args.plan, args.output, env_path=args.env)
    else:
        summary = run_contents_phase(
            args.research_dir, json.loads(args.urls.read_text()), env_path=args.env,
        )
    # Only operational metadata goes to stdout; never queries, snippets, or bodies.
    print(json.dumps({"status": summary["status"], "agent_exposure": AGENT_EXPOSURE}))


if __name__ == "__main__":
    main()
