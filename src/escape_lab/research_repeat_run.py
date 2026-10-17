"""Frozen, repeated evidence and selector comparisons on grouped development folds.

Remote agents receive inner aggregates and fixed cards only. Every search finishes
before any outer predictions are scored. No external evaluator is invoked.
"""

import argparse
import json
import math
import random
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

from escape_lab.biological_features import BiologicalFeatureStore
from escape_lab.cli import credentials
from escape_lab.data import sha256
from escape_lab.experiment_planner import ExperimentPlanner
from escape_lab.improvement import load_labels, write_json
from escape_lab.model_swap import MODEL, REFERENCE_RUN
from escape_lab.research_evidence import load_dossier
from escape_lab.revision_cli import verify_hashes
from escape_lab.revision_evaluation import (
    evaluate_cv,
    grouped_folds,
    predict_selected,
    retained_candidates,
)
from escape_lab.revision_models import catalog
from escape_lab.you_research import _verified_search_archive

ARMS = ("control", "relevant", "placebo")
METHODS = (*ARMS, "fixed", "random", "numerical", "grid")
REPEATS = 10
WORKERS = 6
MAX_TOKENS = 16384
TOKEN_BUDGET = 120000
SECONDS = 600
BOOTSTRAPS = 10000
PLACEBO_CLAIMS = {
    "animal_communication": "Bee waggle dances communicate the direction and distance "
                            "of a food source to other bees in their colony.",
    "leaf_gas_exchange": "Leaf stomatal pores regulate carbon dioxide uptake and water "
                         "loss between plants and their surrounding atmosphere.",
}


def read(path):
    return json.loads(Path(path).read_text())


def make_jobs(folds):
    groups = [(repeat, fold) for repeat in range(REPEATS) for fold in sorted(folds)]
    random.Random(7300).shuffle(groups)
    jobs = []
    for repeat, fold in groups:
        index = int(fold.rsplit("-", 1)[1])
        rotation = (repeat + index) % 3
        for arm in ARMS[rotation:] + ARMS[:rotation]:
            jobs.append({"job_id": f"repeat-{repeat:02d}/{fold}/{arm}",
                         "repeat": repeat, "fold": fold, "arm": arm,
                         "selection_seed": 7100 + 100 * repeat + index})
    return jobs


def validate_blind_rows(rows):
    allowed = {"method", "candidate_id", "group_id", "gene", "drug", "probability"}
    for row in rows:
        if set(row) - allowed:
            raise ValueError("Blind prediction contains label or unapproved fields")
        score = row.get("probability")
        if type(score) not in (int, float) or not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError("Blind prediction probability is invalid")


def attach_labels(rows, labels, candidates):
    expected = {(c["candidate_id"], d) for c in candidates for d in ("pem", "fut")}
    indexed = {(r["candidate_id"], r["drug"]): r for r in rows}
    if len(indexed) != len(rows) or set(indexed) != expected:
        raise ValueError("Prediction coverage differs from the fixed test partition")
    if set(labels) != {c["candidate_id"] for c in candidates}:
        raise ValueError("Label coverage differs from the fixed test partition")
    return [{**row, "label": labels[row["candidate_id"]][row["drug"]]} for row in rows]


def load_placebo(path):
    from escape_lab.research_evidence import _normalize, _url

    dossier = read(path)
    verify_hashes(dossier["source_hashes"])
    review_path = Path(dossier["review_path"])
    verify_hashes({str(review_path): dossier["review_sha256"]})
    review = read(review_path)
    expected = []
    for card in review["cards"]:
        identity = card["id"]
        url = _url(card["source_url"])
        research = Path(review.get("source_research_dirs", {}).get(url, review["research_dir"]))
        _, eligible = _verified_search_archive(research)
        contents = research / "contents/quarantine/response.json"
        proof = [p for p in research.rglob("*") if p.is_file()]
        if any(dossier["source_hashes"].get(str(p)) != sha256(p) for p in proof):
            raise ValueError("Placebo retrieval provenance is missing")
        source_rows = read(contents)
        sources = {row["url"]: row for row in source_rows}
        if len(sources) != len(source_rows):
            raise ValueError("Duplicate placebo contents URLs are ambiguous")
        selected = read(research / "contents/selected-urls.json")
        summary = read(research / "contents/summary.json")
        if summary["raw_files"].get("quarantine/response.json") != sha256(contents):
            raise ValueError("Placebo contents receipt differs")
        usable = {p["url"] for p in summary["pages"]
                  if p["status"] == "quarantined_for_evaluator_review"}
        if identity not in PLACEBO_CLAIMS or card["claim"] != PLACEBO_CLAIMS[identity]:
            raise ValueError("Placebo claim is not the fixed canonical claim")
        if url not in eligible or url not in selected or url not in usable:
            raise ValueError("Placebo source was not retrieved through the planned search")
        text = _normalize(sources[url]["markdown"])
        quote, date = card["supporting_quote"], card["publication_year_support"]
        year = card["publication_year"]
        if (type(year) is not int or not 1900 <= year <= 2023 or str(year) not in date
                or not quote or len(quote.split()) > 25
                or _normalize(quote) not in text or _normalize(date) not in text):
            raise ValueError("Placebo source support or publication date failed")
        expected.append({k: card[k] for k in ("id", "claim", "source_url", "publication_year")})
    if (len(expected) != 2 or {c["id"] for c in expected} != set(PLACEBO_CLAIMS)
            or len({c["source_url"] for c in expected}) != 2
            or expected != dossier["cards"]):
        raise ValueError("Placebo cards differ from reviewed fixed claims")
    return dossier


def prepare(root, output, relevant_path, placebo_path):
    from escape_lab import research_repeat_analysis, research_repeat_selection

    original_path = root / "reports/revision-v2" / REFERENCE_RUN / "protocol.json"
    original = read(original_path)
    verify_hashes(original["hashes"])
    if {k: version(k) for k in original["dependencies"]} != original["dependencies"]:
        raise ValueError("Frozen numerical dependency versions changed")
    if original["catalog"] != {k: v.model_dump() for k, v in catalog().items()}:
        raise ValueError("Frozen catalog changed")
    relevant, placebo = load_dossier(relevant_path), load_placebo(placebo_path)
    lengths = {arm: len(json.dumps(packet["cards"], sort_keys=True))
               for arm, packet in (("relevant", relevant), ("placebo", placebo))}
    if abs(lengths["relevant"] - lengths["placebo"]) / lengths["relevant"] > .10:
        raise ValueError("Evidence packet lengths differ by more than ten percent")
    panel = read(root / "data/external_discovery/s13-panel.json")
    candidates = retained_candidates(read(root / "data/processed/candidates.json"), panel)
    if [c["candidate_id"] for c in candidates] != original["training_ids"]:
        raise ValueError("Original development cohort changed")
    tasks = grouped_folds(candidates, "outer-20261017", prefix="outer")
    jobs = make_jobs([t["task_id"] for t in tasks])
    protected = dict(original["hashes"])
    for folder in ("reports/revision-v2", "reports/model-swap-v2", "reports/research-pilot",
                   "research/you-evidence/you-pilot-20261017"):
        protected.update({str(p): sha256(p) for p in (root / folder).rglob("*") if p.is_file()})
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "preservation-before.json", protected)
    packets = {"control": [], "relevant": relevant["cards"], "placebo": placebo["cards"]}
    write_json(output / "evidence-packets.json", packets)
    write_json(output / "jobs.json", jobs)
    inputs = {**original["hashes"], str(original_path): sha256(original_path),
              str(relevant_path): sha256(relevant_path), str(placebo_path): sha256(placebo_path),
              **relevant["source_hashes"], relevant["review_path"]: relevant["review_sha256"],
              **placebo["source_hashes"], placebo["review_path"]: placebo["review_sha256"]}
    source_paths = [Path(__file__), Path(research_repeat_selection.__file__),
                    Path(research_repeat_analysis.__file__),
                    root / "docs/repeated-research-protocol.md",
                    *(root / "src/escape_lab" / name for name in
                      ("model_swap.py", "research_evidence.py", "you_research.py")),
                    *(root / "tests" / name for name in
                      ("test_research_repeat_run.py", "test_research_repeat_selection.py",
                       "test_research_repeat_analysis.py"))]
    inputs.update({str(p): sha256(p) for p in source_paths})
    protocol = {
        "schema_version": "repeated-research-v1", "created_at": datetime.now(UTC).isoformat(),
        "scientific_attempts": 1, "repeats": REPEATS, "outer_folds": 3,
        "searches": len(jobs), "arms": list(ARMS), "model": MODEL,
        "candidate_count": len(candidates), "group_count": len({c["group_id"] for c in candidates}),
        "training_ids": original["training_ids"], "catalog": original["catalog"],
        "dependencies": original["dependencies"], "source_hashes": inputs,
        "outer_groups": {t["task_id"]: sorted({c["group_id"] for c in t["pool"]}) for t in tasks},
        "numeric_seeds": [4100, 4101, 4102], "selection_seed_rule": "7100+100*repeat+fold_index",
        "schedule_seed": 7300, "schedule": jobs, "concurrent_searches": WORKERS,
        "max_tokens": MAX_TOKENS, "token_stop_threshold_per_search": TOKEN_BUDGET,
        "seconds_per_search": SECONDS, "http_timeout_seconds": 120, "http_retries": 0,
        "proposal_opportunities_per_search": 3, "format_repairs_per_opportunity": 1,
        "sampling_parameters": "Provider defaults; LLM sampling seed is not controlled.",
        "effective_system_prompt": research_repeat_selection.SYSTEM_PROMPT,
        "numerical_policy": research_repeat_selection.NUMERICAL_POLICY,
        "evidence_serialized_character_lengths": lengths,
        "evidence_length_note": "Approximate length match, not exact provider token matching.",
        "primary_metric": "mean of ten repeated pooled outer-fold macro average precisions",
        "primary_contrasts": ["relevant-control", "relevant-placebo", "control-fixed",
                              "control-random", "control-numerical"],
        "minimum_worthwhile_absolute_ap_gain": .002,
        "bootstrap": {"draws": BOOTSTRAPS, "seed": 8100, "units": "crossed aligned groups and repeat IDs",
                      "percentile_intervals": [95, 99], "primary_family_size": 5},
        "decision_rule": "Each required contrast needs delta>=.002, 99% lower>0, per-drug AP "
                         "decrease<=.01, each drug Brier<=comparator*1.05+.001. "
                         "Otherwise directional or inconclusive.",
        "failure_policy": "All searches retained with incumbent/reference fallback; no outcome-based rerun.",
        "comparison_costs": "Reference shared. Agent/random/numerical get 3 new evaluation opportunities; "
                            "fixed structure has zero search; full grid uses 11 and is not budget matched.",
        "data_boundary": "Only inner aggregates reach isolated remote planners. All outer predictions "
                         "committed after all 90 searches and before outer scoring. External panels unused.",
        "archive_boundary": "Existing load_labels parses a monolithic trusted archive but only returns "
                            "requested development partition labels. No final or external scoring occurs.",
        "limitations": ["Reused development data, not fresh blind external validation.",
                        "Ten repeats measure planner variability, not ten biological replications.",
                        "Related aligned blocks and conditional bootstrap limit inference.",
                        "Fixed source, prompt and catalog may make relevant/placebo context too weak.",
                        "Cohort remains original v2; future external transfer requires its own exclusions."]}
    write_json(output / "protocol.json", protocol)
    verify_hashes(protected)
    return candidates, tasks, jobs, packets, protocol, protected


def prepare_caches(root, output, tasks):
    store = BiologicalFeatureStore(root / "data/revision-v2/features")
    caches = {}
    for i, task in enumerate(tasks):
        training = task["training"]
        labels = load_labels(root / "data/processed", training)
        inner = grouped_folds(training, f"inner-seed-{4100+i}")
        metrics = {cid: evaluate_cv(cid, inner, labels, store, 4100 + i)["metrics"]
                   for cid in catalog()}
        caches[task["task_id"]] = metrics
        write_json(output / "inner-cache" / f"{task['task_id']}.json", {
            "training_ids": [c["candidate_id"] for c in training],
            "outer_pool_ids": [c["candidate_id"] for c in task["pool"]], "metrics": metrics})
    return caches


def run_job(job, metrics, cards, output, key, models):
    from escape_lab.research_repeat_selection import RepeatEvidenceClient, select_from_cache

    folder = output / "searches" / job["job_id"]
    folder.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    with RepeatEvidenceClient(key, model_id=MODEL, evidence_cards=cards,
                              max_tokens=MAX_TOKENS, max_retries=0, timeout=120,
                              max_total_tokens=TOKEN_BUDGET,
                              audit_path=folder / "requests.jsonl") as client:
        client.models = models.copy()
        client.configure_model(MODEL)
        result = select_from_cache(metrics, ExperimentPlanner(client),
                                   seed=job["selection_seed"], seconds=SECONDS)
        result.update(job=job, seconds=time.monotonic()-started,
                      budget_charge_tokens=client.total_tokens, usage_records=client.usage_records)
    write_json(folder / "search.json", result)
    return result


def fit_predictions(root, output, tasks, jobs, results, caches):
    from escape_lab.research_repeat_selection import numerical_selection

    choices = {str(r): {} for r in range(REPEATS)}
    for job in jobs:
        choices[str(job["repeat"])].setdefault(job["fold"], {})[job["arm"]] = (
            results[job["job_id"]]["selected"])
    controls = {}
    for repeat in range(REPEATS):
        for i, task in enumerate(tasks):
            fold = task["task_id"]
            control = numerical_selection(caches[fold], 7100 + 100 * repeat + i)
            choices[str(repeat)][fold].update(control["selections"])
            controls[f"{repeat}/{fold}"] = control
    write_json(output / "choices.json", choices)
    write_json(output / "controls.json", controls)
    store = BiologicalFeatureStore(root / "data/revision-v2/features")
    for i, task in enumerate(tasks):
        fold = task["task_id"]
        selected = sorted({config for r in choices.values() for config in r[fold].values()})
        labels = load_labels(root / "data/processed", task["training"])
        prediction = predict_selected({c: c for c in selected}, task["training"], labels,
                                      task["pool"], store, 4100+i, output / "models" / fold)
        validate_blind_rows(prediction["predictions"])
        write_json(output / "predictions" / f"{fold}.json", prediction)


def freeze_predictions(output, jobs, results, input_hashes):
    expected = {job["job_id"] for job in jobs}
    if set(results) != expected or len(expected) != len(jobs):
        raise ValueError("Completed search coverage differs from frozen jobs")
    hashes = dict(input_hashes)
    for path in output.rglob("*"):
        if path.is_file() and path.name not in ("run-status.json", "joint-freeze.json"):
            hashes[str(path.resolve())] = sha256(path)
    verify_hashes(hashes)
    write_json(output / "joint-freeze.json", {
        "schema_version": "repeated-research-joint-freeze-v1",
        "committed_at": datetime.now(UTC).isoformat(), "job_ids": sorted(expected), "hashes": hashes})
    return verify_freeze(output)


def verify_freeze(output):
    frozen = read(output / "joint-freeze.json")
    if (frozen.get("schema_version") != "repeated-research-joint-freeze-v1"
            or not frozen.get("hashes") or not frozen.get("job_ids")):
        raise ValueError("Invalid joint commitment")
    verify_hashes(frozen["hashes"])
    return frozen


def verified_scoring_inputs(output, candidates, tasks):
    """Validate complete committed coverage before requesting even one outer label."""
    frozen = verify_freeze(output)
    hashes = frozen["hashes"]

    def committed(path):
        if str(path.resolve()) not in hashes:
            raise ValueError("Required artifact is absent from joint commitment")
        return read(path)

    protocol = committed(output / "protocol.json")
    jobs = committed(output / "jobs.json")
    choices = committed(output / "choices.json")
    committed(output / "controls.json")
    committed(output / "evidence-packets.json")
    committed(output / "preservation-before.json")
    expected_jobs = make_jobs([t["task_id"] for t in tasks])
    if (len(jobs) != 90 or jobs != expected_jobs or protocol["schedule"] != jobs
            or sorted(j["job_id"] for j in jobs) != frozen["job_ids"]):
        raise ValueError("Frozen coverage must contain all ninety canonical search jobs")
    for filename, digest in protocol["source_hashes"].items():
        if hashes.get(filename) != digest:
            raise ValueError("Protocol source or evidence is not jointly committed")
    if set(choices) != {str(r) for r in range(REPEATS)}:
        raise ValueError("Frozen repeat coverage is incomplete")
    fold_ids = {t["task_id"] for t in tasks}
    if any(set(choice) != fold_ids for choice in choices.values()):
        raise ValueError("Frozen outer fold coverage is incomplete")
    for job in jobs:
        result = committed(output / "searches" / job["job_id"] / "search.json")
        selected = choices[str(job["repeat"])][job["fold"]][job["arm"]]
        if result["job"] != job or result["selected"] != selected:
            raise ValueError("Selection differs from committed search")
    banks = {}
    all_pool_ids = []
    for task in tasks:
        fold = task["task_id"]
        all_pool_ids.extend(c["candidate_id"] for c in task["pool"])
        if any(set(choices[str(r)][fold]) != set(METHODS) for r in range(REPEATS)):
            raise ValueError("Frozen method coverage is incomplete")
        configs = {config for r in choices.values() for config in r[fold].values()}
        if not configs <= set(catalog()):
            raise ValueError("Selected configuration is outside the frozen catalog")
        bank = committed(output / "predictions" / f"{fold}.json")
        rows = bank["predictions"]
        validate_blind_rows(rows)
        expected = {(config, c["candidate_id"], drug) for config in configs
                    for c in task["pool"] for drug in ("pem", "fut")}
        actual = {(r["method"], r["candidate_id"], r["drug"]) for r in rows}
        if actual != expected or len(actual) != len(rows):
            raise ValueError("Frozen prediction coverage is incomplete")
        metadata = {c["candidate_id"]: c for c in task["pool"]}
        if any(r["group_id"] != metadata[r["candidate_id"]]["group_id"]
               or r["gene"] != metadata[r["candidate_id"]]["gene"] for r in rows):
            raise ValueError("Prediction group or gene differs from frozen pool")
        models = {str((output / "models" / fold / config / f"{drug}.pkl").resolve())
                  for config in configs for drug in ("pem", "fut")}
        if set(bank["model_hashes"]) != models or any(
                hashes.get(path) != digest for path, digest in bank["model_hashes"].items()):
            raise ValueError("Selected fitted model is not jointly committed")
        banks[fold] = rows
    if (len(all_pool_ids) != len(set(all_pool_ids))
            or set(all_pool_ids) != {c["candidate_id"] for c in candidates}):
        raise ValueError("Outer pools must partition the complete development cohort")
    return choices, banks


def score_frozen(root, output, candidates, tasks):
    choices, banks = verified_scoring_inputs(output, candidates, tasks)
    # First outer scoring access, after every arm/fold/repeat prediction is committed.
    labels = load_labels(root / "data/processed", candidates)
    records = []
    for task in tasks:
        fold = task["task_id"]
        rows = banks[fold]
        validate_blind_rows(rows)
        cache = {}
        pool_labels = {c["candidate_id"]: labels[c["candidate_id"]] for c in task["pool"]}
        for config in sorted({r["method"] for r in rows}):
            cache[config] = attach_labels([r for r in rows if r["method"] == config],
                                          pool_labels, task["pool"])
        for repeat in range(REPEATS):
            if set(choices[str(repeat)][fold]) != set(METHODS):
                raise ValueError("Frozen method coverage is incomplete")
            for method, config in choices[str(repeat)][fold].items():
                records.extend({**row, "repeat": repeat, "method": method, "fold": fold}
                               for row in cache[config])
    verify_freeze(output)
    return records


def run(root, output, relevant_path, placebo_path):
    from escape_lab.research_repeat_analysis import analyze
    from escape_lab.research_repeat_selection import RepeatEvidenceClient

    root, output = root.resolve(), output.resolve()
    started = time.monotonic()
    candidates, tasks, jobs, packets, protocol, protected = prepare(
        root, output, relevant_path.resolve(), placebo_path.resolve())
    protocol_hash = sha256(output / "protocol.json")
    status = {"status": "started", "scientific_attempts": 1, "completed_searches": 0,
              "started_at": datetime.now(UTC).isoformat(), "protocol_sha256": protocol_hash}
    write_json(output / "run-status.json", status)
    try:
        caches = prepare_caches(root, output, tasks)
        print("Inner-only numerical caches prepared; no outer scores computed.", flush=True)
        key, _ = credentials(root)
        with RepeatEvidenceClient(key, evidence_cards=[], audit_path=output / "discovery-requests.jsonl") as client:
            client.discover_models(nemotron_only=False)
            model = client.configure_model(MODEL)
            if model.tools_supported is not True:
                raise ValueError("GLM must advertise tool support")
            models = client.models.copy()
            write_json(output / "model-discovery.json", model.model_dump())
        results = {}
        with ThreadPoolExecutor(max_workers=WORKERS) as executor:
            futures = {executor.submit(run_job, job, caches[job["fold"]], packets[job["arm"]],
                                       output, key, models): job for job in jobs}
            for future in as_completed(futures):
                job = futures[future]
                results[job["job_id"]] = future.result()
                status.update(completed_searches=len(results), last_completed=job["job_id"])
                write_json(output / "run-status.json", status)
                print(f"Completed search {len(results)}/90: {job['job_id']}; "
                      f"valid={results[job['job_id']]['successful_proposals']}/3", flush=True)
        verify_hashes(protected)
        verify_hashes(protocol["source_hashes"])
        verify_hashes({str(output / "protocol.json"): protocol_hash})
        fit_predictions(root, output, tasks, jobs, results, caches)
        freeze_predictions(output, jobs, results, protocol["source_hashes"])
        print("All ninety searches and all predictions committed. Starting evaluator scoring.", flush=True)
        records = score_frozen(root, output, candidates, tasks)
        analysis = analyze(records, replicates=BOOTSTRAPS, seed=8100)
        usage = {arm: {"searches": 0, "successful_proposals": 0, "failed_opportunities": 0,
                       "requests": 0, "tokens": 0, "prompt_tokens": 0, "completion_tokens": 0,
                       "unknown_usage_requests": 0, "seconds_sum": 0., "statuses": {}}
                 for arm in ARMS}
        for job in jobs:
            result, bucket = results[job["job_id"]], usage[job["arm"]]
            bucket["searches"] += 1
            for field in ("successful_proposals", "failed_opportunities"):
                bucket[field] += result[field]
            bucket["tokens"] += result["budget_charge_tokens"]
            bucket["seconds_sum"] += result["seconds"]
            for record in result["usage_records"]:
                bucket["requests"] += 1
                if record["prompt_tokens"] is None or record["completion_tokens"] is None:
                    bucket["unknown_usage_requests"] += 1
                bucket["prompt_tokens"] += record["prompt_tokens"] or 0
                bucket["completion_tokens"] += record["completion_tokens"] or 0
            for proposal in result["proposals"]:
                name = proposal["status"]
                bucket["statuses"][name] = bucket["statuses"].get(name, 0)+1
        for bucket in usage.values():
            bucket["estimated_uncached_usd"] = (bucket["prompt_tokens"]*1.4
                                                 + bucket["completion_tokens"]*4.4)/1e6
        analysis.update(run_status="completed", usage=usage, total_wall_seconds=time.monotonic()-started,
                        joint_freeze_sha256=sha256(output / "joint-freeze.json"),
                        protected_artifacts=len(protected), fresh_external_validation=False)
        write_json(output / "results.json", analysis)
        write_json(output / "private-scored-records.json", {"records": records})
        verify_hashes(protected)
        verify_freeze(output)
        status.update(status="completed", results_sha256=sha256(output / "results.json"))
        print(json.dumps({"status": "completed", "output": str(output)}), flush=True)
        return analysis
    except BaseException as error:
        status.update(status="failed", error_type=type(error).__name__)
        raise
    finally:
        status["seconds"] = time.monotonic()-started
        status["artifact_sha256"] = {str(p): sha256(p) for p in output.rglob("*")
                                     if p.is_file() and p != output / "run-status.json"}
        write_json(output / "run-status.json", status)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--relevant", type=Path, required=True)
    parser.add_argument("--placebo", type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.output, args.relevant, args.placebo)


if __name__ == "__main__":
    main()
