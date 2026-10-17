"""Post-run descriptive uncertainty; never used for model or protocol selection.

Resamples complete aligned blocks with the same draw for every method. These
intervals condition on already committed predictions and do not include fitting,
selection, source-discovery, or independent biological sampling uncertainty.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score

from escape_lab.data import sha256
from escape_lab.improvement import write_json

METHODS = ("v1_fixed", "biology_fixed", "random", "grid", "agent")
CONTRASTS = (
    ("biology_fixed", "v1_fixed"),
    ("random", "v1_fixed"),
    ("grid", "v1_fixed"),
    ("agent", "v1_fixed"),
    ("agent", "biology_fixed"),
    ("agent", "random"),
    ("agent", "grid"),
)
DRUGS = ("pem", "fut")
REPLICATES = 2000
SEED = 6100


def _paired_arrays(records, endpoint):
    if endpoint not in ("internal", "external"):
        raise ValueError("Endpoint must be internal or external")
    identity = "candidate_id" if endpoint == "internal" else "variant_id"
    score_key = "probability" if endpoint == "internal" else "prediction_score"
    indexed = {method: {} for method in METHODS}
    for row in records:
        method, drug = row["method"], row["drug"]
        if method not in METHODS or drug not in DRUGS:
            raise ValueError("Unexpected method or drug")
        key = (row[identity], drug)
        if key in indexed[method]:
            raise ValueError("Duplicate prediction record")
        score = row[score_key]
        if type(score) not in (int, float) or not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError("Invalid prediction score")
        group = row["group_id"]
        if not isinstance(group, str) or not group:
            raise ValueError("Invalid aligned group")
        if endpoint == "internal":
            value = row["label"]
            if value is not None and type(value) is not bool:
                raise ValueError("Invalid internal label")
            status = "observed" if value is not None else "missing"
        else:
            status, value = row["measurement_status"], row["residual_activity"]
            if not isinstance(status, str) or not status:
                raise ValueError("Invalid measurement status")
            if status == "observed":
                if type(value) not in (int, float) or not math.isfinite(value):
                    raise ValueError("Invalid observed external measurement")
            elif value is not None:
                raise ValueError("Unobserved measurement must not contain a value")
        indexed[method][key] = {"score": score, "value": value, "status": status, "group": group}
    reference = indexed[METHODS[0]]
    if not reference:
        raise ValueError("No paired records")
    expected = set(reference)
    variants = {key[0] for key in expected}
    if expected != {(variant, drug) for variant in variants for drug in DRUGS}:
        raise ValueError("Both drugs are required for every variant")
    for method, rows in indexed.items():
        if set(rows) != expected:
            raise ValueError("Methods have unequal prediction coverage")
        for key, row in rows.items():
            baseline = reference[key]
            if any(row[field] != baseline[field] for field in ("value", "status", "group")):
                raise ValueError("Methods disagree on paired outcomes or groups")
    for variant in variants:
        if reference[(variant, "pem")]["group"] != reference[(variant, "fut")]["group"]:
            raise ValueError("Variant belongs to different groups across drugs")
    groups = sorted({row["group"] for row in reference.values()})
    if len(groups) < 2:
        raise ValueError("At least two aligned groups are required")
    group_index = {group: index for index, group in enumerate(groups)}
    arrays = {}
    for drug in DRUGS:
        keys = sorted(key for key in expected if key[1] == drug
                      and reference[key]["status"] == "observed")
        arrays[drug] = {
            "truth": np.asarray([reference[key]["value"] for key in keys], dtype=float),
            "groups": np.asarray([group_index[reference[key]["group"]] for key in keys], dtype=int),
            "scores": {method: np.asarray([indexed[method][key]["score"] for key in keys])
                       for method in METHODS},
        }
    return groups, arrays


def _metric(arrays, method, counts, endpoint):
    values = []
    for drug in DRUGS:
        rows = arrays[drug]
        weights = counts[rows["groups"]]
        keep = weights > 0
        truth, scores, weights = rows["truth"][keep], rows["scores"][method][keep], weights[keep]
        if not len(truth):
            return None
        if endpoint == "internal":
            # Identical to the frozen evaluator: a zero-positive sample has AP 0.
            value = (float(average_precision_score(truth, scores, sample_weight=weights))
                     if np.dot(truth, weights) > 0 else 0.0)
        else:
            if weights.sum() < 5 or len(set(truth)) < 2 or len(set(scores)) < 2:
                return None
            value = float(spearmanr(np.repeat(scores, weights),
                                    np.repeat(truth, weights)).statistic)
            if not math.isfinite(value):
                return None
        values.append(value)
    return float(np.mean(values))


def paired_group_bootstrap(records, endpoint, replicates=REPLICATES, seed=SEED):
    """Return seven predeclared paired descriptive percentile intervals."""
    if type(replicates) is not int or replicates < 1:
        raise ValueError("Positive replicate count required")
    groups, arrays = _paired_arrays(records, endpoint)
    ones = np.ones(len(groups), dtype=int)
    points = {method: _metric(arrays, method, ones, endpoint) for method in METHODS}
    draws = {contrast: [] for contrast in CONTRASTS}
    rng = np.random.default_rng(seed)
    for _ in range(replicates):
        counts = np.bincount(rng.integers(len(groups), size=len(groups)), minlength=len(groups))
        scores = {method: _metric(arrays, method, counts, endpoint) for method in METHODS}
        for (method, comparator), differences in draws.items():
            if scores[method] is not None and scores[comparator] is not None:
                differences.append(scores[method] - scores[comparator])
    comparisons = {}
    for (method, comparator), differences in draws.items():
        valid = points[method] is not None and points[comparator] is not None
        comparisons[f"{method}_minus_{comparator}"] = {
            "method": method, "comparator_method": comparator,
            "estimate": points[method], "comparator_estimate": points[comparator],
            "delta": points[method] - points[comparator] if valid else None,
            "percentile_2_5": float(np.percentile(differences, 2.5)) if differences else None,
            "percentile_97_5": float(np.percentile(differences, 97.5)) if differences else None,
            "valid_replicates": len(differences), "invalid_replicates": replicates - len(differences),
        }
    return {
        "endpoint": endpoint,
        "metric": "macro_average_precision" if endpoint == "internal" else "mean_per_drug_spearman",
        "aligned_group_count": len(groups), "aligned_groups": groups,
        "observed_pairs_per_drug": {drug: len(arrays[drug]["truth"]) for drug in DRUGS},
        "replicates": replicates, "seed": seed, "point_estimates": points,
        "comparisons": comparisons,
        "interpretation": "Descriptive paired block-bootstrap percentiles conditional on committed "
                          "predictions. Related sequence blocks are not independent biological "
                          "replicates. No refitting, model selection, p-values, or significance claim.",
    }


def _read_verified(path, digest):
    if not isinstance(digest, str) or sha256(path) != digest:
        raise ValueError(f"Released artifact checksum mismatch: {path.name}")
    return json.loads(path.read_text())


def run_uncertainty(run_dir, receipt_path):
    """Read only completed, released evaluation artifacts and write a separate report."""
    run_dir, receipt_path = Path(run_dir).resolve(), Path(receipt_path).resolve()
    output = run_dir / "uncertainty.json"
    if output.exists():
        raise FileExistsError("Uncertainty report already exists")
    status_path = run_dir / "run-status.json"
    status = json.loads(status_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    if status.get("status") != "completed" or receipt.get("status") != "completed":
        raise ValueError("Both run and external receipt must be completed before analysis")
    results_path = run_dir / "results.json"
    results = _read_verified(results_path, status.get("results_sha256"))
    report_path = run_dir / "external/external-report.json"
    if Path(receipt["report_path"]).resolve() != report_path:
        raise ValueError("External receipt belongs to another run")
    report = _read_verified(report_path, receipt.get("report_sha256"))
    records_path = run_dir / "external/external-records.json"
    if Path(report["records_path"]).resolve() != records_path:
        raise ValueError("External record path differs from the completed run")
    if report.get("records_sha256") != receipt.get("records_sha256"):
        raise ValueError("External record checksums disagree")
    records = _read_verified(records_path, receipt.get("records_sha256"))
    nested_path = run_dir / "nested/nested-results.json"
    nested = json.loads(nested_path.read_text())
    internal = paired_group_bootstrap(nested["predictions"], "internal")
    external = paired_group_bootstrap(records["records"], "external")
    for method in METHODS:
        if not math.isclose(internal["point_estimates"][method],
                            results["nested_metrics"][method]["macro_ap"], abs_tol=1e-12):
            raise ValueError("Internal bootstrap point estimate differs from frozen report")
        actual = external["point_estimates"][method]
        expected = report["methods"][method]["macro_spearman"]
        if actual != expected and (actual is None or expected is None
                                   or not math.isclose(actual, expected, abs_tol=1e-12)):
            raise ValueError("External bootstrap point estimate differs from frozen report")
    payload = {"analysis_only": True, "primary_results_unchanged": True,
               "source_sha256": sha256(Path(__file__)),
               "input_sha256": {str(path): sha256(path) for path in
                                (status_path, receipt_path, results_path, report_path,
                                 records_path, nested_path)},
               "internal": internal, "external": external}
    write_json(output, payload)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run_uncertainty(args.run_dir, args.receipt), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
