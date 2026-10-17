"""Evaluator-owned grouped selection; the planner receives inner aggregates only."""

import hashlib
import math
import pickle
import random
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

from escape_lab.benchmark import validate_partition
from escape_lab.data import sha256
from escape_lab.improvement import aggregate_metrics, write_json
from escape_lab.prediction import PredictionConfig, ResistancePredictor
from escape_lab.revision_models import BiologyPredictor, catalog

DRUGS = ("pem", "fut")
REFERENCE = "chemistry-logistic-c0.1"
BIOLOGY_FIXED = "structure-logistic-c0.1"
METHODS = ("v1_fixed", "biology_fixed", "random", "grid", "agent")


def check_deadline(deadline):
    if deadline is not None and time.monotonic() >= deadline:
        raise TimeoutError("Revision run deadline expired between bounded operations")


def retained_candidates(candidates, panel):
    excluded = set(panel["training_exclusion"]["external_aligned_groups"])
    kept = [c for c in candidates if c["split"] in ("development", "validation")
            and c["group_id"] not in excluded]
    if len({c["candidate_id"] for c in kept}) != len(kept):
        raise ValueError("Duplicate protein variants")
    return sorted(kept, key=lambda c: c["candidate_id"])


def grouped_folds(candidates, salt, prefix="inner"):
    groups = sorted({c["group_id"] for c in candidates},
                    key=lambda g: hashlib.sha256(f"fgfr-v2:{salt}:{g}".encode()).hexdigest())
    if len(groups) < 3:
        raise ValueError("At least three complete groups are required")
    assignment = {group: i % 3 for i, group in enumerate(groups)}
    rows = sorted(candidates, key=lambda c: c["candidate_id"])
    result = []
    for i in range(3):
        training = [c for c in rows if assignment[c["group_id"]] != i]
        pool = [c for c in rows if assignment[c["group_id"]] == i]
        validate_partition(training, pool)
        result.append({"task_id": f"{prefix}-{i}", "training": training, "pool": pool})
    return result


def summarize(records):
    result = aggregate_metrics(records)
    briers = []
    for drug in DRUGS:
        metric = result["per_drug"][drug]
        metric["prevalence"] = metric["positive_count"] / metric["n"] if metric["n"] else None
        briers.append(metric["brier"])
    result["macro_brier"] = float(np.mean(briers)) if all(x is not None for x in briers) else None
    folds = sorted({r["fold"] for r in records if "fold" in r})
    if folds:
        result["by_fold"] = {}
        for fold in folds:
            subset = [{k: v for k, v in row.items() if k != "fold"}
                      for row in records if row.get("fold") == fold]
            result["by_fold"][fold] = summarize(subset)
        aps = [m["macro_ap"] for m in result["by_fold"].values()]
        result["mean_fold_macro_ap"] = float(np.mean(aps)) if all(x is not None for x in aps) else None
        result["fold_ap_std"] = float(np.std(aps)) if all(x is not None for x in aps) else None
    return result


def acceptance_reasons(candidate, reference, incumbent=None):
    """Fixed engineering guards, not significance thresholds."""
    reasons = []
    incumbent = incumbent or reference
    for metric in (candidate, reference, incumbent):
        numbers = [metric.get("macro_ap")]
        numbers.extend(metric.get("per_drug", {}).get(d, {}).get(name)
                       for d in DRUGS for name in ("average_precision", "brier"))
        numbers.extend(f.get("macro_ap") for f in metric.get("by_fold", {}).values())
        if any(type(x) not in (float, int) or not math.isfinite(x) for x in numbers):
            return ["undefined_metric"]
    if candidate["macro_ap"] <= incumbent["macro_ap"] + 1e-12:
        reasons.append("no_pooled_ap_improvement")
    for drug in DRUGS:
        c, r = candidate["per_drug"][drug], reference["per_drug"][drug]
        if c["average_precision"] < r["average_precision"] - .01 - 1e-12:
            reasons.append(f"{drug}_ap_regression")
        if c["brier"] > r["brier"] * 1.05 + .001 + 1e-12:
            reasons.append(f"{drug}_brier_regression")
    folds = reference.get("by_fold", {})
    if len(folds) != 3 or set(folds) != set(candidate.get("by_fold", {})):
        return reasons + ["invalid_inner_folds"]
    deltas = [candidate["by_fold"][f]["macro_ap"] - folds[f]["macro_ap"] for f in folds]
    if sum(d >= -1e-12 for d in deltas) < 2:
        reasons.append("inconsistent_fold_gain")
    if min(deltas) < -.01 - 1e-12:
        reasons.append("worst_fold_regression")
    return reasons


def accepts(candidate, reference, incumbent=None):
    return not acceptance_reasons(candidate, reference, incumbent)


def _fit(config_id, training, labels, store, drug, seed):
    known = [c for c in training if type(labels[c["candidate_id"]][drug]) is bool]
    truth = [labels[c["candidate_id"]][drug] for c in known]
    if config_id == "v1_fixed":
        model = ResistancePredictor(PredictionConfig(model="logistic"), seed=seed)
    else:
        model = BiologyPredictor(catalog()[config_id], store, drug, seed=seed)
    return model.fit(known, truth)


def evaluate_cv(config_id, tasks, labels, store, seed, deadline=None):
    start = time.monotonic()
    predictions, diagnostics = [], []
    allowed = {c["candidate_id"] for t in tasks for c in t["training"] + t["pool"]}
    if set(labels) != allowed:
        raise ValueError("CV label boundary must exactly match the current search partition")
    for task in tasks:
        validate_partition(task["training"], task["pool"])
        for drug in DRUGS:
            check_deadline(deadline)
            model = _fit(config_id, task["training"], labels, store, drug, seed)
            probabilities = model.predict(task["pool"])
            diagnostics.append({"fold": task["task_id"], "drug": drug,
                                "fallback_reason": model.fallback_reason,
                                **getattr(model, "diagnostics", {})})
            for candidate, score in zip(task["pool"], probabilities, strict=True):
                predictions.append({"candidate_id": candidate["candidate_id"],
                                    "group_id": candidate["group_id"], "gene": candidate["gene"],
                                    "drug": drug, "fold": task["task_id"],
                                    "probability": float(score),
                                    "label": labels[candidate["candidate_id"]][drug]})
    metrics = summarize(predictions)
    norms = [d["coefficient_norm"] for d in diagnostics if "coefficient_norm" in d]
    if norms:
        metrics["coefficient_norm_mean"] = float(np.mean(norms))
    return {"candidate_id": config_id, "metrics": metrics, "diagnostics": diagnostics,
            "predictions": predictions, "seconds": time.monotonic() - start}


def search_candidates(training, labels, store, planner, seed, deadline, output_dir):
    """Only this partition's inner scores enter a fresh planner's feedback."""
    output_dir = Path(output_dir)
    tasks = grouped_folds(training, f"inner-seed-{seed}")
    if set(labels) != {c["candidate_id"] for c in training}:
        raise ValueError("Search received outcomes outside its training partition")
    check_deadline(deadline)
    configs = catalog()
    cache = {}

    def evaluate(candidate_id):
        if candidate_id not in cache:
            check_deadline(deadline)
            cache[candidate_id] = evaluate_cv(candidate_id, tasks, labels, store, seed, deadline)
            write_json(output_dir / "evaluations" / f"{candidate_id}.json", cache[candidate_id])
        return cache[candidate_id]["metrics"]

    reference = evaluate(REFERENCE)
    incumbent, history, proposals = REFERENCE, [], []
    # Only own history goes into prompts. Random/grid scores never enter feedback.
    history.append({"candidate_id": REFERENCE, "status": "baseline", "metrics": reference})
    for revision_index in range(3):
        check_deadline(deadline)
        feedback = {"incumbent": {"candidate_id": incumbent, "metrics": evaluate(incumbent)},
                    "history": history, "diagnostics": evaluate(incumbent), "seed": seed,
                    "revision_index": revision_index,
                    "remaining_seconds": max(0, deadline - time.monotonic())}
        proposal = planner.propose({key: config.model_dump() for key, config in configs.items()},
                                   feedback)
        record = asdict(proposal)
        candidate_id = proposal.candidate_id
        if proposal.status == "proposed" and candidate_id in configs:
            metrics = evaluate(candidate_id)
            reasons = acceptance_reasons(metrics, reference, evaluate(incumbent))
            accepted = not reasons
            record.update(metrics=metrics, accepted=accepted, rejection_reasons=reasons)
            history.append({"candidate_id": candidate_id, "parent_id": proposal.parent_id,
                            "status": "accepted" if accepted else "rejected",
                            "accepted": accepted, "metrics": metrics})
            if accepted:
                incumbent = candidate_id
        proposals.append(record)
        write_json(output_dir / "agent-progress.json", {"proposals": proposals,
                   "attempts": planner.attempts, "selected": incumbent})
    random_ids = random.Random(seed).sample([c for c in configs if c != REFERENCE], 3)

    def select(ids):
        selected = REFERENCE
        for candidate_id in ids:
            if accepts(evaluate(candidate_id), reference, evaluate(selected)):
                selected = candidate_id
        return selected

    selections = {"v1_fixed": "v1_fixed", "biology_fixed": BIOLOGY_FIXED,
                  "random": select(random_ids), "grid": select(list(configs)), "agent": incumbent}
    result = {"selections": selections, "random_candidates": random_ids,
              "agent_proposals": proposals, "agent_attempts": planner.attempts,
              "inner_metrics": {c: value["metrics"] for c, value in cache.items()},
              "inner_groups": {t["task_id"]: sorted({c["group_id"] for c in t["pool"]})
                               for t in tasks}}
    successes = sum(p["status"] == "proposed" and "metrics" in p for p in proposals)
    result.update(successful_proposals=successes, failed_opportunities=3-successes,
                  agent_status="completed" if successes == 3 else "degraded" if successes else
                  "unavailable")
    write_json(output_dir / "search.json", result)
    return result


def predict_selected(selections, training, labels, pool, store, seed, model_dir=None, deadline=None):
    validate_partition(training, pool)
    if set(labels) != {c["candidate_id"] for c in training}:
        raise ValueError("Prediction fit received outcomes outside training")
    predictions, hashes, diagnostics = [], {}, []
    fitted = {}
    for method, config_id in selections.items():
        for drug in DRUGS:
            check_deadline(deadline)
            key = (config_id, drug)
            if key not in fitted:
                fitted[key] = _fit(config_id, training, labels, store, drug, seed)
            model = fitted[key]
            scores = model.predict(pool)
            if model_dir is not None:
                path = Path(model_dir) / method / f"{drug}.pkl"
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("wb") as handle:
                    pickle.dump(model, handle)
                hashes[str(path.resolve())] = sha256(path)
            diagnostics.append({"method": method, "drug": drug,
                                "fallback_reason": model.fallback_reason,
                                **getattr(model, "diagnostics", {})})
            for c, score in zip(pool, scores, strict=True):
                predictions.append({"method": method, "candidate_id": c["candidate_id"],
                                    "group_id": c["group_id"], "gene": c["gene"],
                                    "drug": drug, "probability": float(score)})
    return {"predictions": predictions, "model_hashes": hashes, "diagnostics": diagnostics}


def nested_compare(candidates, labels, store, planner_factory, output_dir, deadline):
    output_dir = Path(output_dir)
    tasks = grouped_folds(candidates, "outer-20261017", prefix="outer")
    predictions, searches = [], {}
    for i, task in enumerate(tasks):
        check_deadline(deadline)
        fold = task["task_id"]
        # Deliberately give neither search nor fit any outer outcome dictionary.
        inner_labels = {c["candidate_id"]: labels[c["candidate_id"]] for c in task["training"]}
        search = search_candidates(task["training"], inner_labels, store, planner_factory(),
                                   4100 + i, deadline, output_dir / fold)
        searches[fold] = search
        result = predict_selected(search["selections"], task["training"], inner_labels,
                                  task["pool"], store, 4100 + i, deadline=deadline)
        # Only the independent scoring layer attaches outer labels after prediction.
        for row in result["predictions"]:
            row.update(fold=fold, label=labels[row["candidate_id"]][row["drug"]])
        predictions.extend(result["predictions"])
        write_json(output_dir / fold / "outer-result.json", result)
        print(f"Completed {fold}: {search['selections']}", flush=True)
    metrics = {method: summarize([r for r in predictions if r["method"] == method])
               for method in METHODS}
    result = {"metrics": metrics, "predictions": predictions,
              "selections": {fold: value["selections"] for fold, value in searches.items()},
              "agent_status": {fold: value["agent_status"] for fold, value in searches.items()},
              "outer_groups": {t["task_id"]: sorted({c["group_id"] for c in t["pool"]})
                               for t in tasks},
              "interpretation": "Exploratory nested CV on previously used development data."}
    write_json(output_dir / "nested-results.json", result)
    return result
