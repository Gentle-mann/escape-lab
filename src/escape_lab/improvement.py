"""Bounded, autonomous predictor revision with evaluator-owned data and scoring.

The proposer sees only development aggregates. It has no dataset, Python, shell,
or network-search tools. This is a capability boundary, not an OS sandbox claim.
"""

import csv
import json
import multiprocessing
import pickle
import platform
import random
import signal
import time
from contextlib import contextmanager
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from uuid import uuid4

import numpy as np
from sklearn.metrics import average_precision_score, brier_score_loss

from escape_lab.benchmark import build_tasks, validate_partition
from escape_lab.data import collapse_labels, sha256
from escape_lab.prediction import PredictionConfig, ResistancePredictor

DRUGS = ("pem", "fut")
PROTOCOL_VERSION = "autonomous-prediction-v1"


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def load_labels(data_dir, candidates, final_token=None):
    """Expose only the requested partition from the evaluator's private archive."""
    data_dir = Path(data_dir)
    if any(c["split"] == "final" for c in candidates):
        if final_token is None:
            raise ValueError("Final test remains sealed")
        receipt = json.loads((data_dir / "final-evaluation.json").read_text())
        if receipt.get("freeze_sha256") != final_token or receipt.get("status") != "started":
            raise ValueError("Invalid final evaluation authorization")
    # The archival JSON is monolithic. Logical access is enforced here; we do not
    # claim final bytes were never parsed by the trusted evaluator.
    outcomes = json.loads((data_dir / "outcomes.json").read_text())
    result = {}
    for candidate in candidates:
        cid = candidate["candidate_id"]
        result[cid] = {
            drug: collapse_labels([r[drug]["resistant"] for r in outcomes[cid]["records"]])
            for drug in DRUGS
        }
    return result


def aggregate_metrics(records):
    per_drug = {}
    for drug in DRUGS:
        all_rows = [r for r in records if r["drug"] == drug]
        rows = [r for r in all_rows if type(r["label"]) is bool]
        y = np.asarray([int(r["label"]) for r in rows])
        p = np.asarray([r["probability"] for r in rows])
        if p.size and (not np.isfinite(p).all() or np.any((p < 0) | (p > 1))):
            raise ValueError("Invalid prediction probability")
        positives = int(y.sum())
        ap = (float(average_precision_score(y, p)) if positives else 0.0) if len(rows) else None
        per_drug[drug] = {
            "average_precision": ap,
            "brier": float(brier_score_loss(y, p)) if len(rows) else None,
            "n": len(rows), "positive_count": positives,
            "inconclusive_count": len(all_rows) - len(rows),
        }
    scores = [per_drug[d]["average_precision"] for d in DRUGS]
    return {"macro_ap": float(np.mean(scores)) if all(s is not None for s in scores) else None,
            "per_drug": per_drug}


def accepts_revision(candidate, baseline, incumbent=None):
    if candidate.get("macro_ap") is None or baseline.get("macro_ap") is None:
        return False
    target = incumbent or baseline
    if candidate["macro_ap"] <= target["macro_ap"] + 1e-12:
        return False
    return all(candidate["per_drug"][d]["average_precision"] >=
               baseline["per_drug"][d]["average_precision"] - .02 - 1e-12 for d in DRUGS)


def fixed_configurations():
    return [PredictionConfig(model="logistic"), PredictionConfig(),
            PredictionConfig(model="hist_gradient_boosting", max_depth=3)]


def sample_configurations(count, seed, excluded):
    rng = random.Random(seed)
    selected, seen = [], set(excluded)
    for _ in range(10000):
        config = PredictionConfig(
            model=rng.choice(["logistic", "extra_trees", "hist_gradient_boosting"]),
            feature_set=rng.choice(["sequence", "position", "full"]),
            include_gene=rng.choice([True, False]), position_scale=rng.choice([10, 20, 40]),
            regularization=rng.choice([.1, 1., 10.]),
            min_samples_leaf=rng.choice([1, 3, 10]), max_features=rng.choice([.5, .8, 1.]),
            max_depth=rng.choice([3, 6, None]), learning_rate=rng.choice([.05, .1]),
            l2=rng.choice([0., 1., 10.]),
        )
        if config.config_id not in seen:
            selected.append(config)
            seen.add(config.config_id)
        if len(selected) == count:
            return selected
    raise ValueError("Could not sample distinct configurations")


def evaluate_configuration(config, tasks, labels, output_dir, seed):
    """Trusted numeric evaluation; no language model is invoked here."""
    output_dir = Path(output_dir)
    records, diagnostics, artifacts = [], [], {}
    start = time.monotonic()
    for task in tasks:
        training, pool = task["training"], task["pool"]
        validate_partition(training, pool)
        for drug in DRUGS:
            known = [c for c in training if type(labels[c["candidate_id"]][drug]) is bool]
            model = ResistancePredictor(config, seed=seed).fit(
                known, [labels[c["candidate_id"]][drug] for c in known])
            probabilities = model.predict(pool)
            if len(probabilities) != len(pool):
                raise ValueError("Prediction coverage does not match evaluation pool")
            model_path = output_dir / "models" / task["task_id"] / f"{drug}.pkl"
            model_path.parent.mkdir(parents=True, exist_ok=True)
            with model_path.open("wb") as handle:
                pickle.dump(model, handle)
            artifacts[str(model_path.resolve())] = sha256(model_path)
            diagnostics.append({"fold": task["task_id"], "drug": drug,
                                "training_n": len(known), "fallback_reason": model.fallback_reason})
            for candidate, probability in zip(pool, probabilities, strict=True):
                records.append({"candidate_id": candidate["candidate_id"], "drug": drug,
                                "gene": candidate["gene"], "group_id": candidate["group_id"],
                                "fold": task["task_id"], "probability": float(probability),
                                "label": labels[candidate["candidate_id"]][drug]})
    result = {"config": config.model_dump(), "config_id": config.config_id,
              "metrics": aggregate_metrics(records), "predictions": records,
              "artifacts": artifacts, "diagnostics": diagnostics,
              "seconds": time.monotonic() - start, "status": "completed"}
    result["error_analysis"] = {
        "by_gene": {gene: aggregate_metrics([r for r in records if r["gene"] == gene])["per_drug"]
                    for gene in sorted({r["gene"] for r in records})},
        "by_fold": {fold: aggregate_metrics([r for r in records if r["fold"] == fold])
                    for fold in sorted({r["fold"] for r in records})},
    }
    write_json(output_dir / "result.json", result)
    return result


def _worker(config_dict, tasks, labels, output_dir, seed):
    try:
        evaluate_configuration(PredictionConfig(**config_dict), tasks, labels, output_dir, seed)
    except Exception as error:  # noqa: BLE001 -- isolated worker reports all model failures.
        write_json(Path(output_dir) / "error.json", {"status": "failed", "error": type(error).__name__})


def bounded_evaluation(config, tasks, labels, output_dir, seed, seconds):
    if seconds <= 0:
        return {"status": "timed_out", "config": config.model_dump(), "config_id": config.config_id}
    context = multiprocessing.get_context("spawn")
    process = context.Process(target=_worker,
                              args=(config.model_dump(), tasks, labels, str(output_dir), seed))
    start = time.monotonic()
    process.start()
    try:
        process.join(timeout=seconds)
    finally:
        if process.is_alive():
            process.terminate()
            process.join(timeout=5)
            if process.is_alive():
                process.kill()
                process.join()
    if process.exitcode != 0 and time.monotonic() - start >= seconds:
        return {"status": "timed_out", "config": config.model_dump(),
                "config_id": config.config_id, "seconds": time.monotonic() - start}
    path = Path(output_dir) / "result.json"
    if not path.exists():
        return {"status": "failed", "config": config.model_dump(),
                "config_id": config.config_id, "seconds": time.monotonic() - start}
    result = json.loads(path.read_text())
    result["wall_seconds"] = time.monotonic() - start
    return result


def verify_freeze(freeze_path):
    freeze = json.loads(Path(freeze_path).read_text())
    for section in ("fingerprints", "artifacts"):
        for filename, expected in freeze[section].items():
            if sha256(Path(filename)) != expected:
                raise ValueError("Frozen source, data, or model artifact changed")
    return freeze


def reserve_final(data_dir, freeze_path, run_id):
    verify_freeze(freeze_path)
    receipt = {"status": "started", "run_id": run_id,
               "freeze_path": str(Path(freeze_path).resolve()),
               "freeze_sha256": sha256(Path(freeze_path)),
               "started_at": datetime.now(UTC).isoformat()}
    with (Path(data_dir) / "final-evaluation.json").open("x") as handle:
        json.dump(receipt, handle, indent=2)
    return receipt


def _bootstrap(final_records, seed, iterations=500):
    """Descriptive paired block-resampling, not independent biological replication."""
    groups = sorted({r["group_id"] for r in final_records["fixed"]})
    rng = random.Random(seed)
    by_group = {method: {g: [r for r in rows if r["group_id"] == g] for g in groups}
                for method, rows in final_records.items()}
    differences = {"agent_minus_fixed": [], "agent_minus_random": []}
    for _ in range(iterations):
        drawn = rng.choices(groups, k=len(groups))
        scores = {method: aggregate_metrics([r for g in drawn for r in grouped[g]])["macro_ap"]
                  for method, grouped in by_group.items()}
        if any(value is None for value in scores.values()):
            continue
        differences["agent_minus_fixed"].append(scores["agent"] - scores["fixed"])
        differences["agent_minus_random"].append(scores["agent"] - scores["random"])
    return {"method": "paired aligned-block bootstrap; descriptive sensitivity only",
            "groups": len(groups), "iterations": iterations,
            "percentiles_2_5_50_97_5": {
                key: np.percentile(values, [2.5, 50, 97.5]).tolist() if values else None
                for key, values in differences.items()}}


def evaluate_frozen_final(data_dir, freeze_path, output_dir, candidates, seed):
    freeze = verify_freeze(freeze_path)
    receipt = reserve_final(data_dir, freeze_path, freeze["run_id"])
    try:
        pool = [c for c in candidates if c["split"] == "final"]
        final_records = {}
        for method, entry in freeze["selected"].items():
            rows = []
            for drug in DRUGS:
                # These are our own checksum-verified local artifacts, never uploads.
                with Path(entry["model_paths"][drug]).open("rb") as handle:
                    model = pickle.load(handle)
                probabilities = model.predict(pool)
                for candidate, probability in zip(pool, probabilities, strict=True):
                    rows.append({"candidate_id": candidate["candidate_id"], "gene": candidate["gene"],
                                 "group_id": candidate["group_id"], "drug": drug,
                                 "probability": float(probability)})
            final_records[method] = rows
        # Persist all predictions before requesting the held-out labels.
        committed = Path(output_dir) / "final-predictions-committed.json"
        write_json(committed, final_records)
        receipt["predictions_sha256"] = sha256(committed)
        write_json(Path(data_dir) / "final-evaluation.json", receipt)
        if sha256(freeze_path) != receipt["freeze_sha256"]:
            raise ValueError("Frozen manifest changed after final reservation")
        verify_freeze(freeze_path)
        labels = load_labels(data_dir, pool, final_token=receipt["freeze_sha256"])
        for rows in final_records.values():
            for row in rows:
                row["label"] = labels[row["candidate_id"]][row["drug"]]
        metrics = {method: aggregate_metrics(rows) for method, rows in final_records.items()}
        agent = metrics["agent"]
        comparisons = {}
        for method in ("fixed", "random"):
            base = metrics[method]
            delta = agent["macro_ap"] - base["macro_ap"]
            drops = {d: agent["per_drug"][d]["average_precision"] -
                     base["per_drug"][d]["average_precision"] for d in DRUGS}
            comparisons[method] = {"macro_ap_delta": delta, "per_drug_delta": drops,
                                   "passes_margin": delta >= .05 and min(drops.values()) >= -.02}
        result = {"status": "completed", "metrics": metrics, "comparisons": comparisons,
                  "prototype_margin_passed": all(c["passes_margin"] for c in comparisons.values()),
                  "uncertainty": _bootstrap(final_records, seed),
                  "interpretation": "Retrospective held-out assay result; not biological superiority"}
        write_json(Path(output_dir) / "final-results.json", result)
        with (Path(output_dir) / "final-predictions.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["method", "candidate_id", "gene", "group_id",
                                                        "drug", "probability", "label"])
            writer.writeheader()
            for method, rows in final_records.items():
                writer.writerows({"method": method, **row} for row in rows)
        receipt.update(status="completed", report_path=str(Path(output_dir).resolve()),
                       completed_at=datetime.now(UTC).isoformat())
        write_json(Path(data_dir) / "final-evaluation.json", receipt)
        return result
    except Exception:
        receipt["status"] = "failed_after_reservation"
        write_json(Path(data_dir) / "final-evaluation.json", receipt)
        raise


def _compact(result):
    return {key: result[key] for key in ("config", "config_id", "metrics", "status", "seconds")
            if key in result}


def _macro_ap(result):
    return result["metrics"]["macro_ap"]


def _report_markdown(report):
    lines = ["# Autonomous FGFR prediction run", "", f"Run: `{report['run_id']}`",
             f"Outcome: **{report['outcome']}**", "", "## Development and validation", "",
             "| Method | Development macro AP | Validation macro AP | Promoted |",
             "|---|---:|---:|---|"]
    for method, row in report.get("selection", {}).items():
        lines.append(f"| {method} | {row['development_macro_ap']:.4f} | "
                     f"{row['validation_macro_ap']:.4f} | {row['promoted']} |")
    final = report.get("final", {})
    if final.get("status") == "completed":
        lines += ["", "## Frozen final evaluation", "",
                  "| Method | Pemigatinib AP | Futibatinib AP | Macro AP |",
                  "|---|---:|---:|---:|"]
        for method, metrics in final["metrics"].items():
            lines.append(f"| {method} | {metrics['per_drug']['pem']['average_precision']:.4f} | "
                         f"{metrics['per_drug']['fut']['average_precision']:.4f} | "
                         f"{metrics['macro_ap']:.4f} |")
    elif final.get("status") == "sealed":
        lines += ["", "Final test remains sealed. " + final.get("reason", "Evaluation did not qualify.")]
    else:
        lines += ["", ("Final evaluation was reserved but did not complete. Its labels must be "
                  "treated as potentially exposed; no fresh tuning or replacement run is permitted.")]
    lines += ["", "## Attempts", ""]
    for trial in report["trials"]:
        lines.append(f"- {trial['method']} {trial['attempt']}: {trial['status']}; "
                     f"accepted={trial.get('accepted', False)}. {trial.get('reason', '')}")
    lines += ["", (f"Elapsed: {report['elapsed_seconds']:.1f} seconds. "
              f"Within budget: {report['within_time_budget']}."), "",
              ("The ordinary search uses two random configurations, not an exhaustive AutoML search. "
              "All methods share training data, grouped folds, and evaluation rules. Matched trial "
              "budgets are not identical hardware time; actual timings are recorded."), "",
              ("Sequence blocks are not independent tumors or laboratories. Results are retrospective "
              "assay predictions, not treatment recommendations. Public-data pretraining exposure "
              "cannot be ruled out. No claim is made against unavailable historical CAGI submissions.")]
    return "\n".join(lines) + "\n"


@contextmanager
def _wallclock_limit(seconds):
    """POSIX main-process deadline; the supported local CLI runs on macOS/Linux."""
    previous_handler = signal.getsignal(signal.SIGALRM)
    previous_timer = signal.getitimer(signal.ITIMER_REAL)
    if previous_timer[0] != 0:
        raise RuntimeError("An existing process alarm prevents installing the run deadline")

    def expired(signum, frame):
        raise TimeoutError("Autonomous run reached its wall-clock limit")

    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)


def run_improvement(data_dir, output_dir, proposer, *, max_seconds=600, seed=0,
                    search_seed=20261017, evaluate_final=True):
    if not 60 <= max_seconds <= 600:
        raise ValueError("Run wall-clock budget must be between 60 and 600 seconds")
    with _wallclock_limit(max_seconds):
        return _run_improvement(data_dir, output_dir, proposer, max_seconds=max_seconds,
                                seed=seed, search_seed=search_seed, evaluate_final=evaluate_final)


def _run_improvement(data_dir, output_dir, proposer, *, max_seconds=600, seed=0,
                    search_seed=20261017, evaluate_final=True):
    data_dir, output_dir = Path(data_dir).resolve(), Path(output_dir).resolve()
    if (data_dir / "final-evaluation.json").exists():
        raise ValueError("Final evaluation already reserved; fresh tuning is disabled")
    if not 60 <= max_seconds <= 600:
        raise ValueError("Run wall-clock budget must be between 60 and 600 seconds")
    start = time.monotonic()
    deadline = start + max_seconds
    run_id = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    run_dir = output_dir / run_id
    run_dir.mkdir(parents=True)
    manifest = json.loads((data_dir / "manifest.json").read_text())
    for name in ("candidates", "outcomes"):
        if sha256(data_dir / f"{name}.json") != manifest[f"{name}_sha256"]:
            raise ValueError("Dataset checksum mismatch")
    sources = list(Path(__file__).parent.glob("*.py"))
    fingerprints = {str(p.resolve()): sha256(p) for p in sources}
    fingerprints.update({str(p): sha256(p) for p in [data_dir / "manifest.json",
                        data_dir / "candidates.json", data_dir / "outcomes.json"]})
    fixed = fixed_configurations()
    random_configs = sample_configurations(2, search_seed, {c.config_id for c in fixed})
    from escape_lab.improvement_agent import SOURCE_FINGERPRINTS
    protocol = {"version": PROTOCOL_VERSION, "seed": seed, "search_seed": search_seed,
                "max_seconds": max_seconds, "revisions_per_method": 2,
                "primary_metric": "mean of per-drug average precision on all grouped OOF predictions",
                "promotion": "strict macro AP increase; neither drug drops more than .02 from fixed",
                "final_success": "agent macro AP exceeds both controls by >=.05; each drug drop <=.02",
                "final_gate": "two completed agent proposals, complete controls, validation promotion",
                "fixed": [c.model_dump() for c in fixed],
                "random": [c.model_dump() for c in random_configs],
                "fingerprints": fingerprints,
                "runtime": {"python": platform.python_version(),
                            **{name: version(name) for name in ("numpy", "scikit-learn", "pydantic")}},
                "proposer_fingerprints": SOURCE_FINGERPRINTS,
                "generation": {"model_id": getattr(getattr(proposer, "client", None), "model_id", None),
                               "max_tokens": getattr(getattr(proposer, "client", None), "max_tokens", None),
                               "max_total_tokens": getattr(getattr(proposer, "client", None), "max_total_tokens", None),
                               "max_retries": getattr(getattr(proposer, "client", None), "max_retries", None),
                               "request_timeout_seconds": getattr(getattr(proposer, "client", None),
                                                                  "request_timeout_seconds", None),
                               "temperature": "provider default", "inference_seed": "provider default"}}
    write_json(run_dir / "protocol.json", protocol)
    report = {"run_id": run_id, "run_dir": str(run_dir), "protocol": protocol,
              "trials": [], "outcome": "incomplete", "final": {"status": "sealed"}}

    def checkpoint():
        report["elapsed_seconds"] = time.monotonic() - start
        report["within_time_budget"] = report["elapsed_seconds"] <= max_seconds
        write_json(run_dir / "report.json", report)
        (run_dir / "report.md").write_text(_report_markdown(report))

    candidates = json.loads((data_dir / "candidates.json").read_text())
    tasks = build_tasks(candidates, "development")
    dev = [c for c in candidates if c["split"] == "development"]
    labels = load_labels(data_dir, dev)
    evaluated = {}

    def evaluate(config, method, attempt, stage="development", stage_tasks=None, stage_labels=None):
        result = bounded_evaluation(config, stage_tasks or tasks, stage_labels or labels,
                                    run_dir / stage / config.config_id, seed,
                                    max(0, deadline - time.monotonic() - 30))
        trial = {"method": method, "attempt": attempt, **_compact(result), "stage": stage}
        report["trials"].append(trial)
        checkpoint()
        return result, trial

    checkpoint()
    try:
        for index, config in enumerate(fixed):
            result, _ = evaluate(config, "fixed", index)
            if result["status"] != "completed":
                raise RuntimeError("Fixed baseline evaluation failed or timed out")
            evaluated[config.config_id] = result
        baseline = max(evaluated.values(), key=_macro_ap)
        winners = {"fixed": baseline, "random": baseline, "agent": baseline}
        # Complete ordinary controls first so agent latency cannot weaken their budget.
        for index, config in enumerate(random_configs):
            result, trial = evaluate(config, "random", index + 1)
            if result["status"] == "completed":
                accepted = accepts_revision(result["metrics"], baseline["metrics"],
                                            winners["random"]["metrics"])
                trial["accepted"] = accepted
                if accepted:
                    winners["random"] = result
            evaluated[config.config_id] = result
            checkpoint()
        agent_history = []
        for index in range(2):
            # Reserve time for the initial tool call, one schema repair, and local evaluation.
            timeout = getattr(getattr(proposer, "client", None), "request_timeout_seconds", 60)
            if deadline - time.monotonic() < 2 * timeout + 35:
                report["trials"].append({"method": "agent", "attempt": index + 1,
                                         "status": "not_started_budget", "accepted": False})
                checkpoint()
                continue
            incumbent = winners["agent"]
            feedback = {"objective": "Improve separate-drug resistance prediction macro average precision",
                        "revision_index": index + 1, "incumbent": _compact(incumbent),
                        "history": [_compact(evaluated[c.config_id]) for c in fixed] + agent_history,
                        "error_analysis": incumbent["error_analysis"], "seed": seed,
                        "remaining_seconds": deadline - time.monotonic()}
            proposal = proposer.propose(feedback)
            proposal_record = {"rationale": proposal.rationale, "status": proposal.status,
                               "error": proposal.error, "model_id": proposal.model_id,
                               "prompt_version": proposal.prompt_version, "calls": proposal.calls}
            write_json(run_dir / f"proposal-{index+1}.json", proposal_record | {
                "config": proposal.config.model_dump() if proposal.config is not None else None})
            if proposal.config is None:
                report["trials"].append({"method": "agent", "attempt": index + 1,
                    "status": "proposal_failed", "accepted": False, "reason": proposal.error,
                    "proposal": proposal_record})
                checkpoint()
                continue
            config = proposal.config
            # Previously attempted random configurations are not agent-visible. Reuse their
            # deterministic evaluation if selected independently, without giving an extra trial.
            if config.config_id in {c.config_id for c in fixed} or any(
                    h.get("config_id") == config.config_id for h in agent_history):
                report["trials"].append({"method": "agent", "attempt": index + 1,
                    "status": "duplicate_proposal", "accepted": False, "proposal": proposal_record,
                    "config": config.model_dump(), "reason": "No new pipeline was proposed"})
                checkpoint()
                continue
            if config.config_id in evaluated:
                result = evaluated[config.config_id]
                trial = {"method": "agent", "attempt": index + 1, **_compact(result),
                         "stage": "development", "cached_shared_evaluation": True}
                report["trials"].append(trial)
            else:
                result, trial = evaluate(config, "agent", index + 1)
                evaluated[config.config_id] = result
            trial["proposal"] = proposal_record
            trial["reason"] = proposal.rationale
            accepted = result["status"] == "completed" and accepts_revision(
                result["metrics"], baseline["metrics"], winners["agent"]["metrics"])
            trial["accepted"] = accepted
            agent_history.append(_compact(result) | {"accepted": accepted, "revision_index": index+1})
            if accepted:
                winners["agent"] = result
            checkpoint()

        agent_trials = [t for t in report["trials"] if t["method"] == "agent" and
                        isinstance(t["attempt"], int)]
        control_trials = [t for t in report["trials"] if t["method"] in {"fixed", "random"}]
        complete = len(agent_trials) == 2 and all(t["status"] == "completed" for t in agent_trials)
        complete = complete and all(t["status"] == "completed" for t in control_trials)
        report["all_required_attempts_completed"] = complete
        if not complete:
            raise RuntimeError("Comparison incomplete: required proposals or evaluations failed; "
                               "validation and final evaluation were skipped")

        # No more model calls after this point; validation never enters agent feedback.
        validation_candidates = [c for c in candidates if c["split"] in {"development", "validation"}]
        val_labels = load_labels(data_dir, validation_candidates)
        val_tasks = build_tasks(candidates, "validation")
        validation = {}
        for method, winner in winners.items():
            cid = winner["config_id"]
            if cid not in validation:
                config = PredictionConfig(**winner["config"])
                result, _ = evaluate(config, method, "validation", "validation", val_tasks, val_labels)
                if result["status"] != "completed":
                    raise RuntimeError("Validation did not complete within the evaluation budget")
                validation[cid] = result
        fixed_val = validation[baseline["config_id"]]
        selected, selection = {}, {}
        for method, winner in winners.items():
            val = validation[winner["config_id"]]
            promoted = method != "fixed" and accepts_revision(val["metrics"], fixed_val["metrics"])
            selected[method] = val if promoted else fixed_val
            selection[method] = {"development_macro_ap": winner["metrics"]["macro_ap"],
                                 "validation_macro_ap": val["metrics"]["macro_ap"],
                                 "validation_metrics": val["metrics"], "promoted": promoted,
                                 "candidate_config": winner["config"],
                                 "selected_config": selected[method]["config"]}
        report["selection"] = selection
        report["outcome"] = ("validation_improvement" if selection["agent"]["promoted"]
                             else "no_agent_validation_improvement")
        artifacts = {}
        frozen_selection = {}
        for method, result in selected.items():
            artifacts.update(result["artifacts"])
            frozen_selection[method] = {"config": result["config"], "config_id": result["config_id"],
                "model_paths": {d: str(run_dir / "validation" / result["config_id"] /
                                      "models" / "validation-all-groups" / f"{d}.pkl") for d in DRUGS}}
        freeze_path = run_dir / "freeze.json"
        frozen_fingerprints = dict(fingerprints)
        frozen_fingerprints[str(run_dir / "protocol.json")] = sha256(run_dir / "protocol.json")
        freeze = {"run_id": run_id, "frozen_at": datetime.now(UTC).isoformat(),
                  "fingerprints": frozen_fingerprints, "artifacts": artifacts, "selected": frozen_selection,
                  "selection": selection, "criteria": protocol["final_success"]}
        write_json(freeze_path, freeze)
        verify_freeze(freeze_path)
        if evaluate_final and complete and selection["agent"]["promoted"] and deadline-time.monotonic()>30:
            report["final"] = evaluate_frozen_final(data_dir, freeze_path, run_dir, candidates, seed)
            report["outcome"] = ("prototype_margin_passed" if report["final"]["prototype_margin_passed"]
                                 else "prototype_margin_not_met")
        else:
            report["final"] = {"status": "sealed", "reason":
                "No complete, validation-improving agent revision qualified, or final evaluation was disabled/budget-limited."}
    except Exception as error:  # noqa: BLE001 -- persist failed attempts instead of losing the run.
        report["outcome"] = "incomplete"
        report["error"] = f"{type(error).__name__}: {str(error)[:240]}"
        receipt_path = data_dir / "final-evaluation.json"
        if receipt_path.exists():
            try:
                receipt = json.loads(receipt_path.read_text())
                final_status = receipt["status"]
            except (ValueError, KeyError):
                final_status = "reserved_unknown"
            report["final"] = {"status": final_status,
                               "reason": "One-shot final evaluation was reserved; see receipt."}
    client = getattr(proposer, "client", None)
    report["usage"] = getattr(client, "usage_records", [])
    report["proposal_attempts"] = getattr(proposer, "attempts", [])
    checkpoint()
    if not report["within_time_budget"]:
        report["outcome"] = "time_budget_exceeded"
        checkpoint()
    return report
