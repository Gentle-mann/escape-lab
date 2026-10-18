"""One fixed stability challenger and sequence ablation on approved development data.

This version deliberately supersedes source-code hashes, never historical input
hashes. Old protocols remain immutable. It uses the previous cohort allowlist and
fold mapping without opening external assay audits or outcomes.
"""

import argparse
import copy
import json
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from escape_lab.activation_data import load_development_activation
from escape_lab.activation_features import ActivationFeatureAugmenter
from escape_lab.activation_revision_run import SCORES
from escape_lab.biological_features import BiologicalFeatureStore
from escape_lab.data import sha256
from escape_lab.fgfr2_sequence_features import load_features
from escape_lab.fgfr2_sequence_models import (
    BIOLOGY_SCALE_FLOOR,
    MODEL_CONFIG,
    FGFR2SequenceModel,
    activation_label,
)
from escape_lab.fgfr2_sequence_run import REFERENCES, fold_partitions, summarize
from escape_lab.fgfr2_state_run import ESM_CACHE, verify_cached_predictions
from escape_lab.improvement import write_json
from escape_lab.mechanistic_pilot_run import WORKBOOK

VERSION = "fgfr2-sequence-revision-v1"
PRIOR = Path("reports/mechanistic-pilot/fgfr2-sequence-ap-v1")
ARM_CONFIGS = {
    "legacy_both": {"stability": "legacy", "sequence_branches": "both"},
    "floor_both": {"stability": "variance_floor", "sequence_branches": "both"},
    "legacy_wt": {"stability": "legacy", "sequence_branches": "wt"},
    "legacy_delta": {"stability": "legacy", "sequence_branches": "delta"},
}


def approved_cohort(candidates, previous_protocol, *, expected_n=668):
    allowed = previous_protocol["eligible_ids"]
    allowed_set = set(allowed)
    cohort = [candidate for candidate in candidates if candidate.get("candidate_id") in allowed_set]
    if (len(cohort) != expected_n or len(allowed) != expected_n
            or len(allowed_set) != expected_n
            or [candidate["candidate_id"] for candidate in cohort] != allowed):
        raise ValueError("Historical development cohort identities changed")
    reserved = set(previous_protocol["reserved_groups"])
    if any(candidate.get("group_id") in reserved for candidate in cohort):
        raise ValueError("Historical cohort includes a reserved group")
    fold_partitions(cohort, previous_protocol["fold_groups"])
    return cohort


def stability_gate(baseline, proposed, baseline_max, proposed_max):
    if baseline["mean_fold_high_ap"] is None or proposed["mean_fold_high_ap"] is None:
        return {"status": "inconclusive", "reason": "Undefined AP"}
    gain = proposed["mean_fold_high_ap"] - baseline["mean_fold_high_ap"]
    checks = {
        "mean_fold_ap_gain_at_least_001": gain >= .01 - 1e-12,
        "no_fold_ap_loss": all(proposed["by_fold"][fold]["high"]["ap"] >=
                               values["high"]["ap"] - 1e-12
                               for fold, values in baseline["by_fold"].items()),
        "no_high_hit_loss": proposed["high_hits"] >= baseline["high_hits"],
        "maximum_absolute_margin_reduced_tenfold": proposed_max <= baseline_max / 10.,
    }
    return {"status": "pass" if all(checks.values()) else "fail", "checks": checks,
            "mean_fold_ap_gain": gain,
            "high_hits_gain": proposed["high_hits"] - baseline["high_hits"],
            "margin_max_before": baseline_max, "margin_max_after": proposed_max}


def public_result(result):
    """Publish aggregate measurements using an explicit schema, never row records."""
    keys = ("version", "metrics", "gates", "fit_count", "protocol_sha256",
            "baseline_reproduction", "integrity", "fitting_seconds", "elapsed_seconds",
            "failures")
    public = {key: copy.deepcopy(result[key]) for key in keys}
    public["diagnostics"] = []
    for diagnostic in result["diagnostics"]:
        model = diagnostic["model"]
        public["diagnostics"].append({
            key: copy.deepcopy(value) for key, value in diagnostic.items() if key != "model"
        } | {"model": {key: copy.deepcopy(model[key]) for key in (
            "arm", "training_n", "high_count", "feature_count", "biology_feature_count",
            "sequence_components", "sequence_branches", "biology_scaling",
            "coefficient_norm", "iterations", "convergence_warning") if key in model}})
    return public


def report_text(result):
    lines = ["# FGFR2 model improvements", "",
             ("A fixed 12-fit experiment tests numerical stability and identifies which sequence "
             "representations contribute to activation ranking."), "",
             ("The same 668 FGFR2 variants, ten position blocks and three saved development folds "
             "are used throughout. Each fit uses only its training partition."), "",
             "| Arm | Mean-fold AP | High hits / 68 picks | Maximum absolute margin |",
             "|---|---:|---:|---:|"]
    for arm, metric in result["metrics"].items():
        maximum = max(item["margin_abs_max"] for item in result["diagnostics"] if item["arm"] == arm)
        lines.append(f"| {arm} | {metric['mean_fold_high_ap']:.6f} | {metric['high_hits']} | {maximum:.4g} |")
    lines += ["", f"Stability advancement gate: **{result['gates']['stability']['status']}**.",
              "", ("The historical probability and raw-margin predictions are reproduced before "
              "any challenger is fitted. Ranking uses raw margins in every arm."), "",
              ("The variance floor changes only biology feature scaling: fitted standard deviations "
              "below 0.01 use a denominator of 0.01. Sequence scaling and PCA stay unchanged. "
              "The WT-only and mutation-delta-only arms keep legacy scaling to isolate sequence content."),
              "", ("Representation arms are diagnostic ablations. They do not change the default "
              "model automatically. Reused development folds measure this experiment; an external "
              "biological validation remains the next test."), "",
              (f"Fits: {result['fit_count']}; fitting/scoring: {result['fitting_seconds']:.3f}s; "
              f"total: {result['elapsed_seconds']:.3f}s."), "",
              (f"Protected artifacts: {result['integrity']['protected_files']}; "
              f"changed: {len(result['integrity']['changed'])}. No final or external evaluation was run."),
              "", "## Predeclared decision", "", "```json",
              json.dumps(result["gates"], indent=2), "```", ""]
    if result["failures"]:
        lines += ["## Failures", "", json.dumps(result["failures"], indent=2), ""]
    return "\n".join(lines)


def _protected_inputs(root, prior):
    protected = {path: digest for path, digest in prior["protected_sha256"].items()
                 if Path(path).parts[0] in {"data", "reports", "research"}}
    changed = [path for path, digest in protected.items() if sha256(root / path) != digest]
    if changed:
        raise ValueError(f"Historical data or report artifacts changed: {changed}")
    extra = list((root / PRIOR).glob("*")) + list((root / ESM_CACHE).glob("*"))
    extra += [Path(__file__), root / "src/escape_lab/fgfr2_sequence_models.py",
              root / "tests/test_fgfr2_sequence_revision.py",
              root / "tests/test_fgfr2_sequence_models.py"]
    for path in extra:
        if path.is_file():
            protected[str(path.relative_to(root))] = sha256(path)
    return protected


def _scored_rows(test, predictions, activation, fold):
    rows = []
    for candidate, score in zip(test, predictions, strict=True):
        record = activation[candidate["candidate_id"]]
        rows.append({"candidate_id": candidate["candidate_id"], "group_id": candidate["group_id"],
                     "fold": fold, "truth": record["value"], "high": activation_label(record),
                     "strong": activation_label(record, 10.), "score": float(score)})
    return rows


def run(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError("Preserve attempts and choose a fresh run directory")
    # Both current raw outputs and all row-level data must stay in ignored local reports.
    if not output.is_relative_to(root / "reports"):
        raise ValueError("Run output must remain inside the local reports directory")
    started = time.monotonic()
    prior = json.loads((root / PRIOR / "protocol.json").read_text())
    protected = _protected_inputs(root, prior)
    candidates = json.loads((root / "data/processed/candidates.json").read_text())
    cohort = approved_cohort(candidates, prior)
    protocol = {
        "version": VERSION, "created_at": datetime.now(UTC).isoformat(),
        "supersedes_source_only": "The new model APIs and runner have new hashes. Historical "
        "protocols and their frozen data/report artifacts are preserved and verified.",
        "arms": ARM_CONFIGS, "fit_budget": 12, "model_config": MODEL_CONFIG,
        "biology_scale_floor": BIOLOGY_SCALE_FLOOR,
        "ranking": "decision_function margin for every arm",
        "eligible_ids": prior["eligible_ids"], "fold_groups": prior["fold_groups"],
        "reserved_groups": prior["reserved_groups"],
        "stability_gate": {"mean_fold_ap_gain": .01, "no_fold_ap_loss": True,
                           "no_high_hit_loss": True, "max_margin_ratio": .1},
        "representation_arms": "Diagnostic ablations only; both vs WT only vs mutation delta only. "
        "No post-score retuning or automatic promotion.",
        "target": "unchanged high activation: 3/4 replicates >=1.5 and median >5",
        "selection_budget": "ceil(0.1*n) in each saved fold; candidate ID tie breaking",
        "protected_sha256": protected,
    }
    output.mkdir(parents=True, mode=0o700)
    write_json(output / "protocol.json", protocol)
    protocol_hash = sha256(output / "protocol.json")
    write_json(output / "status.json", {"status": "started", "protocol_sha256": protocol_hash})
    rows, diagnostics, failures = {arm: [] for arm in ARM_CONFIGS}, [], []
    fit_count = 0
    try:
        # Outcomes are read only after the design and input hashes are committed to disk.
        activation, _ = load_development_activation(root / WORKBOOK, cohort)
        if set(activation) != {candidate["candidate_id"] for candidate in cohort}:
            raise ValueError("Incomplete approved activation cohort")
        cached = json.loads((root / PRIOR / "private-predictions.json").read_text())
        sequence = load_features(cohort, root / REFERENCES, root / ESM_CACHE)
        store = BiologicalFeatureStore(root / "data/revision-v2/features")
        augmenter = ActivationFeatureAugmenter(store, cohort, root / REFERENCES, root / SCORES)
        tasks = fold_partitions(cohort, prior["fold_groups"])
        reproduction = {}
        probability_rows = []
        fitting_started = time.monotonic()
        for arm, options in ARM_CONFIGS.items():
            if arm != "legacy_both" and not reproduction:
                raise ValueError("Historical baseline must reproduce before challenger fits")
            for fold, train, test in tasks:
                records = {candidate["candidate_id"]: activation[candidate["candidate_id"]]
                           for candidate in train}
                fit_count += 1
                model = FGFR2SequenceModel("esm_biology_logistic", store, augmenter, sequence,
                                          **options).fit(train, records)
                outputs = model.prediction_outputs(test)
                probability, margin = outputs["probability"], outputs["margin"]
                rows[arm] += _scored_rows(test, margin, activation, fold)
                if arm == "legacy_both":
                    probability_rows += _scored_rows(test, probability, activation, fold)
                diagnostics.append({"arm": arm, "fold": fold, "n": len(test),
                                    "probability_unique": len(np.unique(probability)),
                                    "probability_exact_one": int((probability == 1).sum()),
                                    "margin_abs_max": float(np.max(np.abs(margin))),
                                    "model": model.diagnostics})
            if arm == "legacy_both":
                reproduction = {
                    "probability": verify_cached_predictions(probability_rows, cached["probability"]),
                    "margin": verify_cached_predictions(rows[arm], cached["margin"]),
                }
            # Persist each completed arm so a later failure cannot erase its evidence.
            write_json(output / "private-predictions.json", rows)
            (output / "private-predictions.json").chmod(0o600)
        fitting_seconds = time.monotonic() - fitting_started
        metrics = {arm: summarize(values) for arm, values in rows.items()}
        maxima = {arm: max(item["margin_abs_max"] for item in diagnostics if item["arm"] == arm)
                  for arm in ARM_CONFIGS}
        gates = {"stability": stability_gate(metrics["legacy_both"], metrics["floor_both"],
                                             maxima["legacy_both"], maxima["floor_both"]),
                 "representation": {"status": "diagnostic_only", "default_changed": False,
                                    "arms": ["legacy_both", "legacy_wt", "legacy_delta"]}}
        changed = [path for path, digest in protected.items() if sha256(root / path) != digest]
        if changed or sha256(output / "protocol.json") != protocol_hash:
            raise ValueError(f"Protected artifacts or protocol changed: {changed}")
        result = {"version": VERSION, "metrics": metrics, "gates": gates,
                  "baseline_reproduction": reproduction, "diagnostics": diagnostics,
                  "fit_count": fit_count, "fitting_seconds": fitting_seconds,
                  "elapsed_seconds": time.monotonic() - started, "protocol_sha256": protocol_hash,
                  "integrity": {"protected_files": len(protected), "changed": changed},
                  "failures": failures}
        write_json(output / "results.json", public_result(result))
        (output / "report.md").write_text(report_text(result))
        write_json(output / "status.json", {"status": "completed", "protocol_sha256": protocol_hash,
                                           "results_sha256": sha256(output / "results.json")})
        return result
    except Exception as error:
        failures.append({"error_type": type(error).__name__, "error": str(error),
                         "attempted_fits": fit_count})
        write_json(output / "failures.json", failures)
        write_json(output / "private-predictions.json", rows)
        (output / "private-predictions.json").chmod(0o600)
        write_json(output / "status.json", {"status": "failed", "protocol_sha256": protocol_hash,
                                           "failures": failures})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.root, args.output)
    print(json.dumps({"metrics": {arm: {key: metric[key] for key in (
        "mean_fold_high_ap", "high_hits", "strong_hits", "budget")}
        for arm, metric in result["metrics"].items()}, "gates": result["gates"],
        "fit_count": result["fit_count"], "fitting_seconds": result["fitting_seconds"],
        "integrity": result["integrity"]}, indent=2))
