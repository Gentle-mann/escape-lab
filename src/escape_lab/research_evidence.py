"""Project reviewed, retrieved sources into a closed set of general biology cards.

Raw You.com snippets, paper titles, arbitrary claims and source instructions never
enter the proposing agent's context. A trusted curator verifies supporting text
and publication dates; this code binds that review to exact retrieved bytes.
This is a strict context gate, not a guarantee against model pretraining exposure.
"""

import argparse
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

from escape_lab.data import sha256
from escape_lab.improvement import write_json
from escape_lab.revision_cli import verify_hashes
from escape_lab.you_research import _verified_search_archive

CLAIMS = {
    "kinase_conformation": "Protein kinases can adopt multiple conformations; a single crystal "
                           "structure does not describe their full conformational ensemble.",
    "ligand_contacts": "Ligand recognition depends on local interaction geometry and chemical "
                       "environment, rather than sequence position alone.",
    "solvent_environment": "Solvent exposure and burial affect the physical environment of "
                           "amino acid side chains and ligand interactions.",
    "structure_kinetics": "Static structural proximity alone does not establish binding "
                         "kinetics or cellular response.",
}
CUTOFF_YEAR = 2023
HOSTS = {"pubmed.ncbi.nlm.nih.gov", "pmc.ncbi.nlm.nih.gov", "www.rcsb.org", "rcsb.org"}
BLOCKED = re.compile(r"s41587-026-03090-8|s41588-025-02431-8|saifudeen|rly[- ]?4008|"
                     r"supplement|/supp|table[ _-]?s13|external[-_]?records", re.IGNORECASE)
CARD_FIELDS = {"id", "claim", "source_url", "publication_year", "publication_year_support",
               "supporting_quote"}


def _read(path):
    return json.loads(Path(path).read_text())


def _url(value):
    parsed = urlsplit(value)
    if (parsed.scheme != "https" or parsed.hostname not in HOSTS or parsed.username
            or parsed.password or parsed.port not in (None, 443) or parsed.query or parsed.fragment
            or BLOCKED.search(unquote(value))):
        raise ValueError("Source URL is outside the reviewed evidence boundary")
    if parsed.hostname == "pubmed.ncbi.nlm.nih.gov":
        valid_path = re.fullmatch(r"/[1-9][0-9]*/?", parsed.path)
    elif parsed.hostname == "pmc.ncbi.nlm.nih.gov":
        valid_path = re.fullmatch(r"/articles/PMC[1-9][0-9]*/?", parsed.path)
    else:
        valid_path = re.fullmatch(r"/structure/[A-Za-z0-9]{4}", parsed.path)
    if valid_path is None:
        raise ValueError("Citation must use a canonical source identifier path")
    return value


def _normalize(value):
    return " ".join(value.split())


def _render_review(review):
    if review.get("schema_version") != "reviewed-research-claims-v1":
        raise ValueError("Unsupported evidence review schema")
    research = Path(review["research_dir"]).resolve()
    plan, content = research / "plan.json", research / "contents/quarantine/response.json"
    hashes = {str(plan): review["plan_sha256"], str(content): review["content_sha256"]}
    verify_hashes(hashes)
    _, eligible_urls = _verified_search_archive(research)
    search_summary = _read(research / "search-summary.json")
    contents_summary = _read(research / "contents/summary.json")
    selected_urls = _read(research / "contents/selected-urls.json")
    if not isinstance(selected_urls, list) or not set(selected_urls) <= eligible_urls:
        raise ValueError("Contents sources must originate in the planned Search results")
    if contents_summary["raw_files"].get("quarantine/response.json") != sha256(content):
        raise ValueError("Contents response differs from its retrieval receipt")
    reviewable_urls = {p["url"] for p in contents_summary["pages"]
                       if p["status"] == "quarantined_for_evaluator_review"}
    proof_paths = [research / "plan.sha256", research / "search-summary.json",
                   research / "usage.json", research / "contents/summary.json",
                   research / "contents/selected-urls.json", research / "contents/usage.json"]
    proof_paths.extend(research / filename for filename in search_summary["raw_files"])
    hashes.update({str(p): sha256(p) for p in proof_paths})
    sources = _read(content)
    if not isinstance(sources, list):
        raise TypeError("Contents response must be an array")
    indexed = {r["url"]: r for r in sources}
    if len(indexed) != len(sources):
        raise ValueError("Duplicate source responses are ambiguous")
    cards = review.get("cards")
    if not isinstance(cards, list) or not 1 <= len(cards) <= 4:
        raise ValueError("One to four reviewed evidence cards required")
    rendered, ids, urls = [], set(), set()
    for card in cards:
        if set(card) != CARD_FIELDS:
            raise ValueError("Unapproved evidence review fields")
        identity, claim = card["id"], card["claim"]
        if identity not in CLAIMS or identity in ids or claim != CLAIMS[identity]:
            raise ValueError("Claim is outside the fixed reviewed claim library")
        url = _url(card["source_url"])
        if (url in urls or url not in indexed or url not in selected_urls
                or url not in reviewable_urls):
            raise ValueError("Each card needs a distinct retrieved source")
        year = card["publication_year"]
        if type(year) is not int or not 1900 <= year <= CUTOFF_YEAR:
            raise ValueError("Evidence is outside the publication cutoff")
        source = indexed[url]
        text = source.get("markdown")
        if not isinstance(text, str) or not text.strip():
            raise ValueError("Retrieved full source text is unavailable")
        if BLOCKED.search(url + " " + str(source.get("title", ""))):
            raise ValueError("Excluded study source")
        quote, year_support = card["supporting_quote"], card["publication_year_support"]
        if (not isinstance(quote, str) or not quote.strip() or len(quote.split()) > 25
                or _normalize(quote) not in _normalize(text)):
            raise ValueError("Claim must have a short exact supporting source excerpt")
        if (not isinstance(year_support, str) or str(year) not in year_support
                or _normalize(year_support) not in _normalize(text)):
            raise ValueError("Publication date needs curator-verified source support")
        rendered.append({"id": identity, "claim": claim, "source_url": url,
                         "publication_year": year})
        ids.add(identity)
        urls.add(url)
    return rendered, hashes


def compile_dossier(review_path, output_path):
    review_path, output_path = Path(review_path).resolve(), Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError("Evidence dossier already exists")
    cards, hashes = _render_review(_read(review_path))
    dossier = {"schema_version": "curated-biological-evidence-v1", "cards": cards,
               "review_path": str(review_path), "review_sha256": sha256(review_path),
               "source_hashes": hashes, "gate_source_sha256": sha256(Path(__file__)),
               "boundary": "General mechanism context only. No raw snippets, target responses, "
                           "variant effects, prior scores, fitted models or selection advice.",
               "review_limit": "Publication date and claim support are curator-reviewed; "
                               "pretraining contamination cannot be excluded."}
    write_json(output_path, dossier)
    return dossier


def load_dossier(path):
    dossier = _read(path)
    if dossier.get("schema_version") != "curated-biological-evidence-v1":
        raise ValueError("Unsupported evidence dossier")
    if dossier.get("gate_source_sha256") != sha256(Path(__file__)):
        raise ValueError("Evidence gate source changed")
    review_path = Path(dossier["review_path"])
    verify_hashes({str(review_path): dossier["review_sha256"], **dossier["source_hashes"]})
    expected, hashes = _render_review(_read(review_path))
    if dossier["cards"] != expected or dossier["source_hashes"] != hashes:
        raise ValueError("Evidence dossier differs from the reviewed source claims")
    return dossier


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("review", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = compile_dossier(args.review, args.output)
    print(json.dumps({"dossier": str(args.output), "card_count": len(result["cards"])}))


if __name__ == "__main__":
    main()
