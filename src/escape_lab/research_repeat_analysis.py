"""Evaluator-only analysis of paired, committed research and selection predictions.

The experimental unit for predictive uncertainty is the complete aligned group.
Repeated agent runs measure procedural variation, not independent biology. No
model is refitted or selected here, and no result is sent back to a proposer.
"""

import hashlib
import math
import re

import numpy as np

from escape_lab.revision_evaluation import DRUGS, summarize

METHODS = ("control", "relevant", "placebo", "fixed", "random", "numerical", "grid")
REPEATS = tuple(range(10))
CONTRASTS = (
    ("relevant", "control"),
    ("relevant", "placebo"),
    ("control", "fixed"),
    ("control", "random"),
    ("control", "numerical"),
)
MINIMUM_AP_GAIN = .002


def _validate(records):
    indexed = {}
    for row in records:
        method, repeat = row["method"], row["repeat"]
        if method not in METHODS or type(repeat) is not int or repeat not in REPEATS:
            raise ValueError("Unexpected method or repeat; repeats must be 0 through 9")
        if row["drug"] not in DRUGS:
            raise ValueError("Unexpected drug")
        for field in ("candidate_id", "group_id", "gene"):
            if not isinstance(row[field], str) or not row[field]:
                raise ValueError(f"Invalid paired metadata: {field}")
        if not isinstance(row["fold"], str) or not re.fullmatch(r"outer-\d+", row["fold"]):
            raise ValueError("Invalid paired fold metadata")
        score = row["probability"]
        if type(score) not in (int, float) or not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError("Invalid probability")
        if row["label"] is not None and type(row["label"]) is not bool:
            raise ValueError("Invalid label; expected bool or None")
        key = (row["candidate_id"], row["drug"])
        cell = indexed.setdefault(method, {}).setdefault(repeat, {})
        if key in cell:
            raise ValueError("Duplicate prediction record")
        cell[key] = row
    required = set(METHODS) - {"grid"}
    if not required.issubset(indexed):
        raise ValueError("Incomplete method coverage")
    for cells in indexed.values():
        if set(cells) != set(REPEATS):
            raise ValueError("Incomplete repeat coverage")
    reference = indexed["control"][0]
    expected = set(reference)
    variants = {key[0] for key in expected}
    if not expected or expected != {(variant, drug) for variant in variants for drug in DRUGS}:
        raise ValueError("Both drugs required for every candidate in paired coverage")
    for cells in indexed.values():
        for cell in cells.values():
            if set(cell) != expected:
                raise ValueError("Unequal paired prediction coverage")
            for key, row in cell.items():
                if any(row[field] != reference[key][field]
                       for field in ("group_id", "gene", "fold", "label")):
                    raise ValueError("Disagreement on paired truth or metadata")
    for variant in variants:
        if any(reference[(variant, "pem")][field] != reference[(variant, "fut")][field]
               for field in ("group_id", "gene", "fold")):
            raise ValueError("Variant metadata differs across paired drugs")
    groups = sorted({row["group_id"] for row in reference.values()})
    if len(groups) < 2:
        raise ValueError("At least two aligned groups required")
    group_folds = {}
    for row in reference.values():
        previous = group_folds.setdefault(row["group_id"], row["fold"])
        if previous != row["fold"]:
            raise ValueError("An aligned group crosses outer folds")
    methods = tuple(method for method in METHODS if method in indexed)
    return indexed, methods, groups, sorted(expected)


def _mean_summaries(summaries):
    """Average repeat metrics; never pool duplicated records from ten runs."""
    first = summaries[0]
    if isinstance(first, dict):
        return {key: _mean_summaries([summary[key] for summary in summaries]) for key in first}
    if any(value is None for value in summaries):
        return None
    if all(type(value) is int for value in summaries) and len(set(summaries)) == 1:
        return first
    return float(np.mean(summaries))


def _bootstrap_counts(group_count, replicates, seed):
    rng = np.random.default_rng(seed)
    group_draws = np.empty((replicates, group_count), dtype=np.int64)
    repeat_draws = np.empty((replicates, len(REPEATS)), dtype=np.int64)
    for draw in range(replicates):
        group_draws[draw] = np.bincount(
            rng.integers(group_count, size=group_count), minlength=group_count)
        repeat_draws[draw] = np.bincount(
            rng.integers(len(REPEATS), size=len(REPEATS)), minlength=len(REPEATS))
    return group_draws, repeat_draws


def _weighted_ap_draws(truth, scores, groups, counts):
    """Weighted noninterpolated AP with sklearn-equivalent score-tie handling.

    Complete group multiplicities supply sample weights. Work in bounded batches
    and sort each score vector just once, rather than fitting or sorting once per
    bootstrap draw. A zero-positive sample has AP 0, matching summarize().
    """
    result = np.full(len(counts), np.nan)
    if not len(truth):
        return result
    order = np.argsort(-scores, kind="stable")
    sorted_scores, sorted_truth, sorted_groups = scores[order], truth[order], groups[order]
    tie_ends = np.r_[np.flatnonzero(np.diff(sorted_scores)), len(scores) - 1]
    for start in range(0, len(counts), 256):
        weights = counts[start:start + 256, sorted_groups]
        positive = np.cumsum(weights * sorted_truth, axis=1)[:, tie_ends]
        total = np.cumsum(weights, axis=1)[:, tie_ends]
        increments = np.diff(positive, axis=1, prepend=0)
        precision = np.divide(positive, total, out=np.zeros_like(positive), where=total > 0)
        denominator = positive[:, -1]
        values = np.divide(np.sum(increments * precision, axis=1), denominator,
                           out=np.zeros(len(weights)), where=denominator > 0)
        values[total[:, -1] == 0] = np.nan
        result[start:start + len(weights)] = values
    return result


def _bootstrap_metrics(indexed, methods, groups, keys, group_draws, repeat_draws):
    reference = indexed["control"][0]
    group_index = {group: index for index, group in enumerate(groups)}
    vector_ids, vectors = {}, []
    cell_ids = {(method, repeat): [] for method in methods for repeat in REPEATS}
    for drug in DRUGS:
        drug_keys = [key for key in keys if key[1] == drug and reference[key]["label"] is not None]
        truth = np.asarray([reference[key]["label"] for key in drug_keys], dtype=float)
        positions = np.asarray([group_index[reference[key]["group_id"]] for key in drug_keys],
                               dtype=int)
        for method in methods:
            for repeat in REPEATS:
                scores = np.asarray([indexed[method][repeat][key]["probability"]
                                     for key in drug_keys], dtype=float)
                identity = (drug, scores.tobytes())
                if identity not in vector_ids:
                    vector_ids[identity] = len(vectors)
                    vectors.append(_weighted_ap_draws(truth, scores, positions, group_draws))
                cell_ids[(method, repeat)].append(vector_ids[identity])
    pattern_cache, result = {}, {}
    for method in methods:
        pattern = tuple(tuple(cell_ids[(method, repeat)]) for repeat in REPEATS)
        if pattern not in pattern_cache:
            repeat_metrics = np.stack([
                (vectors[pem] + vectors[fut]) / 2 for pem, fut in pattern], axis=1)
            pattern_cache[pattern] = np.sum(repeat_metrics * repeat_draws, axis=1) / len(REPEATS)
        result[method] = pattern_cache[pattern]
    return result, len(vectors), len(pattern_cache)


def _defined_summary(metric):
    values = [metric["macro_ap"], metric["macro_brier"]]
    values.extend(metric["per_drug"][drug][field]
                  for drug in DRUGS for field in ("average_precision", "brier"))
    values.extend(fold["macro_ap"] for fold in metric.get("by_fold", {}).values())
    return all(value is not None and math.isfinite(value) for value in values)


def _comparison(method, comparator, summaries, draws, replicates):
    metric, baseline = summaries[method]["mean"], summaries[comparator]["mean"]
    defined = _defined_summary(metric) and _defined_summary(baseline)
    point_defined = metric["macro_ap"] is not None and baseline["macro_ap"] is not None
    delta = metric["macro_ap"] - baseline["macro_ap"] if point_defined else None
    differences = draws[method] - draws[comparator]
    valid = int(np.isfinite(differences).sum())
    complete = valid == replicates
    interval_95 = np.percentile(differences, [2.5, 97.5]).tolist() if complete else None
    interval_99 = np.percentile(differences, [.5, 99.5]).tolist() if complete else None
    reasons, guards = [], {}
    if not defined:
        reasons.append("undefined_metric")
    if not complete:
        reasons.append("undefined_bootstrap_replicate")
    if delta is not None and delta < MINIMUM_AP_GAIN:
        reasons.append("point_gain_below_0.002")
    if interval_99 is not None and interval_99[0] <= 0:
        reasons.append("99_percent_lower_bound_not_positive")
    for drug in DRUGS:
        candidate, reference = metric["per_drug"][drug], baseline["per_drug"][drug]
        numbers = [candidate[field] for field in ("average_precision", "brier")]
        numbers.extend(reference[field] for field in ("average_precision", "brier"))
        if any(value is None for value in numbers):
            guards[drug] = {"ap_noninferiority": None, "brier_noninferiority": None}
            continue
        ap_delta = candidate["average_precision"] - reference["average_precision"]
        brier_limit = reference["brier"] * 1.05 + .001
        ap_ok, brier_ok = ap_delta >= -.01 - 1e-12, candidate["brier"] <= brier_limit + 1e-12
        guards[drug] = {"ap_delta": ap_delta, "ap_noninferiority": ap_ok,
                        "brier": candidate["brier"], "brier_limit": brier_limit,
                        "brier_noninferiority": brier_ok}
        if not ap_ok:
            reasons.append(f"{drug}_ap_regression")
        if not brier_ok:
            reasons.append(f"{drug}_brier_regression")
    repeated_deltas = []
    for repeat in REPEATS:
        first = summaries[method]["per_repeat"][str(repeat)]["macro_ap"]
        second = summaries[comparator]["per_repeat"][str(repeat)]["macro_ap"]
        repeated_deltas.append(first - second if first is not None and second is not None else None)
    repeat_complete = all(value is not None for value in repeated_deltas)
    wins = sum(value > 1e-12 for value in repeated_deltas) if repeat_complete else None
    losses = sum(value < -1e-12 for value in repeated_deltas) if repeat_complete else None
    ties = len(REPEATS) - wins - losses if repeat_complete else None
    return {
        "method": method, "comparator_method": comparator,
        "estimate": metric["macro_ap"], "comparator_estimate": baseline["macro_ap"],
        "delta": delta, "direction": "undefined" if delta is None else (
            "higher" if delta > 1e-12 else "lower" if delta < -1e-12 else "equal"),
        "interval_95": interval_95, "interval_99_bonferroni": interval_99,
        "valid_replicates": valid, "invalid_replicates": replicates - valid,
        "repeat_only": {"delta_sd": float(np.std(repeated_deltas, ddof=1))
                        if repeat_complete else None,
                        "wins": wins, "ties": ties, "losses": losses,
                        "deltas": repeated_deltas,
                        "interpretation": "Run variability on the same biology, not independent "
                                          "biological replicates."},
        "guards": guards, "passes_improvement_rule": defined and complete and not reasons,
        "failed_rules": reasons,
    }


def _experiment(comparisons, names):
    selected = [comparisons[name] for name in names]
    if all(comparison["passes_improvement_rule"] for comparison in selected):
        status = "supported_development_improvement"
    elif any(comparison["delta"] is None or comparison["invalid_replicates"]
             or "undefined_metric" in comparison["failed_rules"] for comparison in selected):
        status = "inconclusive"
    elif all(comparison["delta"] <= 0 for comparison in selected):
        status = "no_observed_improvement"
    else:
        status = "inconclusive"
    return {"status": status, "required_contrasts": names,
            "contrast_directions": {name: comparisons[name]["direction"] for name in names},
            "interpretation": "Development evidence conditional on frozen predictions and this "
                              "cohort; not fresh biological generalization."}


def analyze(records, *, replicates=2000, seed=8100):
    """Analyze all ten paired repeats without treating runs as new biological data.

    Required methods are control, relevant, placebo, fixed, random, and numerical;
    a fully paired grid reference is optional. Both drugs must be represented for
    every candidate, including explicit None labels for inconclusive observations.
    Undefined metrics or bootstrap draws make the affected claim inconclusive;
    no valid-only interval silently drops those draws.
    """
    if type(replicates) is not int or replicates < 1:
        raise ValueError("Positive integer bootstrap replicate count required")
    if type(seed) is not int or seed < 0:
        raise ValueError("Nonnegative integer seed required")
    indexed, methods, groups, keys = _validate(records)
    summaries = {}
    for method in methods:
        per_repeat = {str(repeat): summarize([indexed[method][repeat][key] for key in keys])
                      for repeat in REPEATS}
        summaries[method] = {"mean": _mean_summaries(list(per_repeat.values())),
                             "per_repeat": per_repeat}
    group_draws, repeat_draws = _bootstrap_counts(len(groups), replicates, seed)
    draws, vector_count, pattern_count = _bootstrap_metrics(
        indexed, methods, groups, keys, group_draws, repeat_draws)
    comparisons = {f"{method}_minus_{comparator}": _comparison(
        method, comparator, summaries, draws, replicates) for method, comparator in CONTRASTS}
    return {
        "schema_version": "paired-research-repeat-analysis-v1",
        "metric": "mean_across_10_runs_of_pooled_outer_oof_macro_average_precision",
        "repeat_count": len(REPEATS), "aligned_group_count": len(groups),
        "aligned_groups": groups, "candidate_count": len(keys) // len(DRUGS),
        "methods": summaries, "comparisons": comparisons,
        "experiments": {
            "research": _experiment(comparisons, ["relevant_minus_control", "relevant_minus_placebo"]),
            "agent_selection": _experiment(comparisons, ["control_minus_fixed", "control_minus_random",
                                                         "control_minus_numerical"]),
        },
        "bootstrap": {
            "replicates": replicates, "seed": seed,
            "group_draws_shared_across_methods": True,
            "repeat_draws_shared_across_methods": True,
            "resampling": "Resample complete aligned groups and matched run indices independently; "
                          "use the same group and run multiplicities for every method in a draw.",
            "draw_sha256": hashlib.sha256(group_draws.tobytes() + repeat_draws.tobytes()).hexdigest(),
            "unique_drug_score_vectors": vector_count, "unique_method_repeat_patterns": pattern_count,
            "intervals": "Percentile 95% descriptive and 99% Bonferroni for five primary contrasts; "
                         "no interval is reported if any draw is undefined.",
        },
        "improvement_rule": {
            "minimum_absolute_macro_ap_gain": MINIMUM_AP_GAIN,
            "99_percent_lower_bound": "strictly greater than zero for every required contrast",
            "maximum_per_drug_ap_decrease": .01,
            "maximum_per_drug_brier": "comparator * 1.05 + 0.001",
        },
        "limitations": [
            ("Bootstrap intervals condition on committed predictions and do not refit models or "
             "repeat model selection, evidence discovery, or biological sampling."),
            ("Related aligned groups are not independent biological replicates; grouped "
             "resampling does not establish fresh-study generalization."),
            ("Ten runs measure agent and selection variability on the same cohort, not ten "
             "independent biological validations."),
            ("Mean fold scores and repeat-only variation are descriptive; primary estimates "
             "average pooled out-of-fold macro AP separately calculated for each run."),
        ],
    }
