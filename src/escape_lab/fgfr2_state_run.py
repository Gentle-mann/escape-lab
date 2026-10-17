"""Frozen, matched continuous FGFR2 conformation-descriptor pilot."""

import argparse
import json
import math
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from escape_lab.activation_data import load_development_activation
from escape_lab.activation_features import ActivationFeatureAugmenter
from escape_lab.activation_revision_run import SCORES
from escape_lab.biological_features import BiologicalFeatureStore
from escape_lab.data import sha256
from escape_lab.fgfr2_sequence_features import load_features
from escape_lab.fgfr2_sequence_models import activation_label
from escape_lab.fgfr2_sequence_run import (
    REFERENCES,
    feature_candidates,
    fold_partitions,
    summarize,
)
from escape_lab.fgfr2_state_models import ARMS, MODEL_CONFIG, PRIMARY, FGFR2StateModel
from escape_lab.improvement import write_json
from escape_lab.mechanistic_pilot_run import WORKBOOK

PREVIOUS = Path("reports/mechanistic-pilot/fgfr2-sequence-v1")
RESEARCH = Path("research/fgfr2-state-pilot")
ESM_CACHE = Path("data/fgfr2-sequence-pilot/esm2-35m-v1")
HISTORICAL = "historical_biology_hgb"


def continuous_summary(rows):
    if len({r["candidate_id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate scored candidates")
    result = summarize(rows)
    for fold, metric in result["by_fold"].items():
        selected = [r for r in rows if r["fold"] == fold]
        y = np.asarray([r["truth"] for r in selected])
        p = np.asarray([r["score"] for r in selected])
        if not np.isfinite(y).all():
            raise ValueError("Nonfinite continuous targets")
        metric["mae"] = float(np.abs(y - p).mean())
        metric["rmse"] = float(np.sqrt(np.square(y - p).mean()))
    rhos = [r["spearman"] for r in result["by_fold"].values()]
    result["mean_fold_spearman"] = float(np.mean(rhos)) if all(r is not None for r in rhos) else None
    result["mean_fold_mae"] = float(np.mean([r["mae"] for r in result["by_fold"].values()]))
    result["pooled_mae"] = float(np.mean([abs(r["truth"] - r["score"]) for r in rows]))
    return result


def compare_gates(metrics):
    a, b, inactive, active = (metrics[name] for name in (PRIMARY, "baseline", "inactive", "active"))
    if any(m["mean_fold_spearman"] is None for m in (a, b, inactive, active)):
        return {name: {"status": "inconclusive", "reason": "Undefined fold correlation"}
                for name in ("incremental", "state_contrast")}
    gain = a["mean_fold_spearman"] - b["mean_fold_spearman"]
    wins = sum(a["by_fold"][f]["spearman"] > b["by_fold"][f]["spearman"]
               for f in a["by_fold"])
    checks = {"mean_fold_rho_gain_003": gain >= .03 - 1e-12,
              "rho_improves_two_folds": wins >= 2,
              "mae_within_5pct": a["mean_fold_mae"] <= b["mean_fold_mae"] * 1.05 + 1e-12,
              "no_high_hit_loss": a["high_hits"] >= b["high_hits"]}
    contrast = {"rho_at_least_inactive": a["mean_fold_spearman"] >= inactive["mean_fold_spearman"] - 1e-12,
                "rho_at_least_active": a["mean_fold_spearman"] >= active["mean_fold_spearman"] - 1e-12,
                "mae_within_5pct_best_single": a["mean_fold_mae"] <= 1.05 * min(
                    active["mean_fold_mae"], inactive["mean_fold_mae"]) + 1e-12}
    return {"incremental": {"status": "pass" if all(checks.values()) else "fail",
                            "checks": checks, "rho_gain": gain, "rho_fold_wins": wins,
                            "high_hits_gain": a["high_hits"] - b["high_hits"]},
            "state_contrast": {"status": "pass" if all(contrast.values()) else "fail",
                               "checks": contrast,
                               "interpretation": "Whether contrast beats the single-state controls; "
                                                 "not a causal mechanism test."}}


def verify_cached_predictions(rows, cached):
    first = {r["candidate_id"]: r for r in rows}
    second = {r["candidate_id"]: r for r in cached}
    if len(first) != len(rows) or len(second) != len(cached) or set(first) != set(second):
        raise ValueError("Cached prediction identities mismatch")
    if any(any(first[cid][key] != second[cid][key] for key in
               ("group_id", "fold", "truth", "high", "strong")) for cid in first):
        raise ValueError("Cached partition or target mismatch")
    differences = [abs(first[cid]["score"] - second[cid]["score"]) for cid in first]
    if not differences or not all(math.isfinite(d) for d in differences) or max(differences) > 1e-12:
        raise ValueError("Matched baseline did not reproduce cached predictions")
    return {"n": len(first), "max_abs_difference": max(differences), "tolerance": 1e-12}


def paired_uncertainty(rows, draws=1000):
    """Stratified block resampling conditional on the fixed out-of-fold scores."""
    rng = np.random.default_rng(17)
    folds = sorted({r["fold"] for r in rows["baseline"]})
    references = ("baseline", "inactive", "active")
    arrays = {arm: {fold: sorted([r for r in data if r["fold"] == fold],
                                key=lambda r: r["candidate_id"]) for fold in folds}
              for arm, data in rows.items()}
    group_rows = {}
    for fold in folds:
        baseline = arrays["baseline"][fold]
        for arm in (PRIMARY, *references):
            selected = arrays[arm][fold]
            if len(selected) != len(baseline) or any(any(a[key] != b[key] for key in
                    ("candidate_id", "group_id", "truth"))
                    for a, b in zip(selected, baseline, strict=True)):
                raise ValueError("Unpaired continuous bootstrap inputs")
        groups = sorted({r["group_id"] for r in baseline})
        group_rows[fold] = [np.asarray([i for i, r in enumerate(baseline) if r["group_id"] == g])
                            for g in groups]
    deltas = {ref: {"rho_gain": [], "mae_reduction": []} for ref in references}
    for _ in range(draws):
        fold_stats = {arm: {"rho": [], "mae": []} for arm in (PRIMARY, *references)}
        for fold in folds:
            groups = group_rows[fold]
            indices = np.concatenate([groups[i] for i in rng.integers(len(groups), size=len(groups))])
            y = np.asarray([r["truth"] for r in arrays["baseline"][fold]])[indices]
            for arm in fold_stats:
                p = np.asarray([r["score"] for r in arrays[arm][fold]])[indices]
                rho = float(spearmanr(y, p).statistic) if np.ptp(y) and np.ptp(p) else None
                fold_stats[arm]["rho"].append(rho)
                fold_stats[arm]["mae"].append(float(np.abs(y-p).mean()))
        for ref in references:
            a, b = fold_stats[PRIMARY], fold_stats[ref]
            if all(x is not None for x in a["rho"] + b["rho"]):
                deltas[ref]["rho_gain"].append(float(np.mean(a["rho"]) - np.mean(b["rho"])))
            deltas[ref]["mae_reduction"].append(float(np.mean(b["mae"]) - np.mean(a["mae"])))
    return {"draws": draws, "seed": 17, "unit": "aligned block within original outer fold",
            "contrasts": {ref: {name: {"valid_draws": len(v), "invalid_draws": draws-len(v),
                "percentile_2_5": float(np.percentile(v, 2.5)) if v else None,
                "percentile_97_5": float(np.percentile(v, 97.5)) if v else None}
                for name, v in contrast.items()} for ref, contrast in deltas.items()},
            "interpretation": "Descriptive fixed-prediction percentiles on reused development blocks; "
                              "no refitting or independent biological validation."}


def report_text(result):
    lines = ["# FGFR2 conformation-descriptor pilot", "",
             (f"Primary incremental gate: **{result['gates']['incremental']['status']}**. "
             f"State-contrast check: **{result['gates']['state_contrast']['status']}**."), "",
             ("Continuous target: log2(1 + median activation), unchanged four-replicate protein summaries. "
             "Same 668 variants, ten occupied position blocks, and three original folds."), "",
             "| Arm | Mean-fold Spearman | Mean-fold MAE | High hits /26 at 68 picks | Strong hits /7 |",
             "|---|---:|---:|---:|---:|"]
    for arm, m in result["metrics"].items():
        rho = "undefined" if m["mean_fold_spearman"] is None else f"{m['mean_fold_spearman']:.5f}"
        lines.append(f"| {arm} | {rho} | {m['mean_fold_mae']:.5f} | {m['high_hits']} | {m['strong_hits']} |")
    lines.extend(["", ("Higher Spearman and retrieval, lower MAE are better. These are not percentage "
                  "accuracy. Historical HGB is an unchanged cached reference, not a newly trained arm."),
                  "", "## Fixed decisions", "", "```json", json.dumps(result["gates"], indent=2),
                  "```", "", "## Fold results", "",
                  "| Fold | Arm | Spearman | MAE | High hits |", "|---|---|---:|---:|---:|"])
    for arm, m in result["metrics"].items():
        for fold, value in m["by_fold"].items():
            rho = "undefined" if value["spearman"] is None else f"{value['spearman']:.5f}"
            lines.append(f"| {fold} | {arm} | {rho} | {value['mae']:.5f} | {value['high']['hits']} |")
    lines.extend(["", "## Scope and verification", "",
                  ("- All four fitted arms use the same continuous ridge head and preprocessing. "
                  "Structural arms have equal feature count and a common geometry mask."),
                  ("- Inactive coordinates are a published FGFR1-templated FGFR2 homology model; "
                  "active coordinates derive from prepared 2PVF. Both retain C491A and lack phosphotyrosine."),
                  ("- Engineered position 491 and declared uncertain loops are excluded from queries and "
                  "contact environments. Only common residue/atom identities contribute."),
                  ("- Static contact/packing differences are mechanistic proxies, not energy, state "
                  "population, relaxed mutant geometry or simulated catalytic activity."),
                  ("- Reused development data only; external assay outcomes and original final evaluation "
                  "remain excluded. No post-score parameter tuning."),
                  f"- Exact previous ridge reproduction: {result['baseline_reproduction']}.",
                  f"- Protected artifacts unchanged: {result['integrity']['protected_files']}.",
                  (f"- Fitting/scoring: {result['fitting_seconds']:.3f}s; evaluator total: "
                  f"{result['elapsed_seconds']:.3f}s."), ""])
    return "\n".join(lines)


def run(root, output, feature_dir):
    from escape_lab.fgfr2_state_features import FGFR2StateFeatureStore

    root, output, feature_dir = (Path(p).resolve() for p in (root, output, feature_dir))
    if output.exists():
        raise FileExistsError("Use a fresh run directory; never overwrite a previous attempt")
    started = time.monotonic()
    cohort, old, _ = feature_candidates(root)
    previous_protocol = json.loads((root / PREVIOUS / "protocol.json").read_text())
    if previous_protocol["eligible_ids"] != [c["candidate_id"] for c in cohort]:
        raise ValueError("Cohort differs from previous FGFR2 pilot")
    if previous_protocol["fold_groups"] != old["fold_groups"]:
        raise ValueError("Original fold mapping changed")
    protected = dict(previous_protocol["protected_sha256"])
    changed = [p for p, h in protected.items() if sha256(root / p) != h]
    if changed:
        raise ValueError(f"Previous protected artifacts changed: {changed}")
    added = list((root / PREVIOUS).glob("*")) + list((root / RESEARCH).rglob("*"))
    added += list(feature_dir.glob("*"))
    added += list((root / "src/escape_lab").glob("fgfr2_state_*.py"))
    added += list((root / "tests").glob("test_fgfr2_state_*.py"))
    for p in added:
        if p.is_file() and "__pycache__" not in p.parts:
            protected[str(p.relative_to(root))] = sha256(p)
    protocol = {"created_at": datetime.now(UTC).isoformat(), "version": "fgfr2-state-v1",
                "arms": ARMS, "primary": PRIMARY, "model_config": MODEL_CONFIG,
                "historical_reference_only": HISTORICAL, "fold_groups": old["fold_groups"],
                "eligible_ids": [c["candidate_id"] for c in cohort],
                "reserved_groups": old["reserved_groups"],
                "target": "unchanged log2(1 + median of four protein replicate summaries)",
                "incremental_gate": {"mean_fold_spearman_gain": .03, "fold_wins": 2,
                    "max_mean_fold_mae_ratio": 1.05, "high_hits_no_loss": True},
                "state_contrast_check": "rho at least both single states; MAE <=1.05 best single state",
                "tail_secondary": {"high": "3/4 replicates>=1.5 AND median>5",
                    "strong": "3/4 replicates>=1.5 AND median>10",
                    "budget": "ceil(.1*n) per original fold", "ties": "candidate_id ascending"},
                "bootstrap": {"draws": 1000, "seed": 17, "unit": "block within fold"},
                "protected_sha256": protected,
                "scope": "Reused development only; structural proxies, not equilibrium simulation"}
    output.mkdir(parents=True)
    write_json(output / "protocol.json", protocol)
    protocol_hash = sha256(output / "protocol.json")
    write_json(output / "status.json", {"status": "started", "protocol_sha256": protocol_hash})
    try:
        activation, coverage = load_development_activation(root / WORKBOOK, cohort)
        if set(activation) != {c["candidate_id"] for c in cohort}:
            raise ValueError("Incomplete activation cohort")
        previous = json.loads((root / PREVIOUS / "private-predictions.json").read_text())
        cached = previous["esm_biology_ridge"]
        group_fold = {g: f for f, groups in old["fold_groups"].items() for g in groups}
        expected = []
        for c in cohort:
            r = activation[c["candidate_id"]]
            expected.append({"candidate_id": c["candidate_id"], "group_id": c["group_id"],
                             "fold": group_fold[c["group_id"]], "truth": r["value"],
                             "high": activation_label(r), "strong": activation_label(r, 10.),
                             "score": 0.})
        reference_map = {r["candidate_id"]: r for r in cached}
        for r in expected:
            r["score"] = reference_map[r["candidate_id"]]["score"]
        verify_cached_predictions(expected, cached)
        # Check the historical comparator against the same evaluator-owned targets.
        historical = previous["cached_v2_biology"]
        historical_by_id = {r["candidate_id"]: r for r in historical}
        historical_expected = [{**r, "score": historical_by_id[r["candidate_id"]]["score"]}
                               for r in expected]
        verify_cached_predictions(historical_expected, historical)
        sequence = load_features(cohort, root / REFERENCES, root / ESM_CACHE)
        state = FGFR2StateFeatureStore(feature_dir)
        store = BiologicalFeatureStore(root / "data/revision-v2/features")
        augmenter = ActivationFeatureAugmenter(store, cohort, root / REFERENCES, root / SCORES)
        rows = {arm: [] for arm in ARMS}
        diagnostics = []
        tasks = fold_partitions(cohort, old["fold_groups"])
        fitting_started = time.monotonic()
        reproduction = None
        for arm in ARMS:
            if arm != "baseline" and reproduction is None:
                raise ValueError("Baseline must reproduce before new structural model fits")
            for fold, train, test in tasks:
                records = {c["candidate_id"]: activation[c["candidate_id"]] for c in train}
                model = FGFR2StateModel(arm, store, augmenter, sequence, state).fit(train, records)
                scores = model.predict(test)
                diagnostics.append({"fold": fold, **model.diagnostics})
                for c, score in zip(test, scores, strict=True):
                    r = activation[c["candidate_id"]]
                    rows[arm].append({"candidate_id": c["candidate_id"], "group_id": c["group_id"],
                                      "fold": fold, "score": float(score), "truth": r["value"],
                                      "high": activation_label(r), "strong": activation_label(r, 10.)})
            if arm == "baseline":
                reproduction = verify_cached_predictions(rows[arm], cached)
        fitting_seconds = time.monotonic() - fitting_started
        metrics = {arm: continuous_summary(data) for arm, data in rows.items()}
        metrics[HISTORICAL] = continuous_summary(historical)
        uncertainty = paired_uncertainty(rows)
        changed = [p for p, h in protected.items() if sha256(root / p) != h]
        if changed or sha256(output / "protocol.json") != protocol_hash:
            raise ValueError(f"Protected artifacts changed: {changed}")
        result = {"metrics": metrics, "gates": compare_gates(metrics), "coverage": coverage,
                  "baseline_reproduction": reproduction, "uncertainty": uncertainty,
                  "model_diagnostics": diagnostics, "feature_manifest": state.manifest,
                  "fitting_seconds": fitting_seconds, "elapsed_seconds": time.monotonic()-started,
                  "integrity": {"protected_files": len(protected), "changed": changed},
                  "protocol_sha256": protocol_hash}
        write_json(output / "private-predictions.json", rows)
        (output / "private-predictions.json").chmod(0o600)
        write_json(output / "results.json", result)
        (output / "report.md").write_text(report_text(result))
        write_json(output / "status.json", {"status": "completed", "protocol_sha256": protocol_hash,
                                           "results_sha256": sha256(output / "results.json")})
        return result
    except Exception as error:
        write_json(output / "status.json", {"status": "failed", "protocol_sha256": protocol_hash,
                                           "error_type": type(error).__name__, "error": str(error)})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.root, args.output, args.features)
    print(json.dumps({"gates": result["gates"], "baseline_reproduction": result["baseline_reproduction"],
                      "fitting_seconds": result["fitting_seconds"]}, indent=2))
