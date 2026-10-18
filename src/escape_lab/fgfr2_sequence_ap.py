"""One fixed AP sensitivity: preserve classifier margins before sigmoid rounding."""

import argparse
import json
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.special import expit
from sklearn.metrics import average_precision_score

from escape_lab.activation_data import load_development_activation
from escape_lab.activation_features import ActivationFeatureAugmenter
from escape_lab.activation_revision_run import SCORES
from escape_lab.biological_features import BiologicalFeatureStore
from escape_lab.data import sha256
from escape_lab.fgfr2_sequence_features import load_features
from escape_lab.fgfr2_sequence_models import FGFR2SequenceModel, activation_label
from escape_lab.fgfr2_sequence_run import REFERENCES, feature_candidates, fold_partitions, summarize
from escape_lab.fgfr2_state_run import ESM_CACHE, PREVIOUS, verify_cached_predictions
from escape_lab.improvement import write_json
from escape_lab.mechanistic_pilot_run import WORKBOOK

PRIOR = Path("reports/mechanistic-pilot/fgfr2-state-tree-v1")
RESEARCH = Path("research/fgfr2-sequence-ap-pilot")


def ranking_outputs(model, candidates):
    """Keep parent prediction guards, then expose the same head's raw margin."""
    if model.arm != "esm_biology_logistic":
        raise ValueError("Only the frozen sequence classifier is permitted")
    outputs = model.prediction_outputs(candidates)
    probability, margin = outputs["probability"], outputs["margin"]
    if margin.shape != (len(candidates),) or not np.isfinite(margin).all():
        raise ValueError("Finite one-dimensional classifier margins required")
    if not np.allclose(expit(margin), probability, rtol=1e-12, atol=1e-14):
        raise ValueError("Margin and cached probability head disagree")
    return outputs


def screening_gate(baseline, proposed):
    if baseline["mean_fold_high_ap"] is None or proposed["mean_fold_high_ap"] is None:
        return {"status": "inconclusive", "reason": "Undefined AP"}
    gain = proposed["mean_fold_high_ap"] - baseline["mean_fold_high_ap"]
    checks = {"mean_fold_ap_gain_at_least_001": gain >= .01 - 1e-12,
              "no_fold_ap_loss": all(proposed["by_fold"][f]["high"]["ap"] >=
                  b["high"]["ap"] - 1e-12 for f, b in baseline["by_fold"].items()),
              "no_high_hit_loss": proposed["high_hits"] >= baseline["high_hits"]}
    return {"status": "pass" if all(checks.values()) else "fail", "checks": checks,
            "mean_fold_ap_gain": gain,
            "high_hits_gain": proposed["high_hits"] - baseline["high_hits"],
            "interpretation": "Technical ranking-recovery screen, not independent model validation"}


def paired_ap_interval(rows, draws=1000):
    baseline = sorted(rows["probability"], key=lambda r: r["candidate_id"])
    proposed = sorted(rows["margin"], key=lambda r: r["candidate_id"])
    if len(baseline) != len(proposed) or any(any(a[k] != b[k] for k in
        ("candidate_id", "group_id", "fold", "truth", "high", "strong"))
        for a, b in zip(baseline, proposed, strict=True)):
        raise ValueError("Unpaired AP bootstrap")
    folds = []
    for fold in sorted({r["fold"] for r in baseline}):
        indices = [i for i, r in enumerate(baseline) if r["fold"] == fold]
        records = [baseline[i] for i in indices]
        groups = sorted({r["group_id"] for r in records})
        blocks = [np.array([i for i, r in enumerate(records) if r["group_id"] == g]) for g in groups]
        folds.append((blocks, np.array([r["high"] for r in records], dtype=int),
                      np.array([r["score"] for r in records]),
                      np.array([proposed[i]["score"] for i in indices])))
    rng, differences = np.random.default_rng(17), []
    for _ in range(draws):
        gains = []
        for blocks, y, b, p in folds:
            indices = np.concatenate([blocks[i] for i in rng.integers(len(blocks), size=len(blocks))])
            if not y[indices].sum():
                break
            gains.append(float(average_precision_score(y[indices], p[indices]) -
                               average_precision_score(y[indices], b[indices])))
        if len(gains) == len(folds):
            differences.append(float(np.mean(gains)))
    return {"draws": draws, "valid_draws": len(differences),
            "invalid_draws": draws-len(differences), "seed": 17,
            "percentile_2_5": float(np.percentile(differences, 2.5)) if differences else None,
            "percentile_97_5": float(np.percentile(differences, 97.5)) if differences else None,
            "interpretation": "Descriptive paired-block mean-fold AP gains on reused development data"}


def make_report(result):
    lines = ["# Tiny sequence-classifier AP test", "",
             f"Technical ranking screen: **{result['gate']['status']}**.", "",
             ("Only change: rank by the same fitted classifier's decision margin instead of its "
             "sigmoid probability. Three original classifier fits, no parameter search."), "",
             "| Ranking | Mean-fold high-activation AP ↑ | High hits /26 at 68 picks | Strong hits /7 |",
             "|---|---:|---:|---:|"]
    for arm, m in result["metrics"].items():
        lines.append(f"| {arm} | {m['mean_fold_high_ap']:.5f} | {m['high_hits']} | {m['strong_hits']} |")
    lines += ["", "| Fold | Probability AP | Margin AP | Probability hits | Margin hits |",
              "|---|---:|---:|---:|---:|"]
    for fold in result["metrics"]["probability"]["by_fold"]:
        b, a = [result["metrics"][arm]["by_fold"][fold]["high"] for arm in ("probability", "margin")]
        lines.append(f"| {fold} | {b['ap']:.5f} | {a['ap']:.5f} | {b['hits']} | {a['hits']} |")
    lines += ["", "## Scope", "",
              "- Same 668 FGFR2 variants, original folds, target, embeddings, PCA, scaler, head and class weights.",
              ("- Raw margins resolve sigmoid rounding ties; they do not fix unsafe feature extrapolation, "
              "improve calibration or add biological information."),
              "- Sequence AP superiority was concentrated in pilot-1; it was not consistent across folds.",
              "- Development outcomes are reused. External assay outcomes and final evaluation remain reserved.",
              f"- Exact original probability reproduction: {result['baseline_reproduction']}.",
              f"- Three fits/scoring: {result['fitting_seconds']:.3f}s; total: {result['elapsed_seconds']:.3f}s.",
              f"- Protected files unchanged: {result['integrity']['protected_files']}.", "",
              "```json", json.dumps({"gate": result["gate"], "uncertainty": result["uncertainty"],
                                       "diagnostics": result["diagnostics"]}, indent=2), "```", ""]
    return "\n".join(lines)


def run(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError("Preserve previous attempts; use a new directory")
    start = time.monotonic()
    cohort, old, _ = feature_candidates(root)
    prior = json.loads((root / PRIOR / "protocol.json").read_text())
    protected = dict(prior["protected_sha256"])
    changed = [p for p, h in protected.items() if sha256(root / p) != h]
    if changed:
        raise ValueError(f"Protected files changed: {changed}")
    added = list((root / PRIOR).glob("*")) + list((root / RESEARCH).glob("*"))
    added += [Path(__file__), root / "tests/test_fgfr2_sequence_ap.py"]
    for path in added:
        if path.is_file():
            protected[str(path.relative_to(root))] = sha256(path)
    protocol = {"created_at": datetime.now(UTC).isoformat(), "version": "fgfr2-sequence-ap-v1",
                "change": "Raw decision margin instead of rounded sigmoid probability; identical head",
                "eligible_ids": [c["candidate_id"] for c in cohort], "fold_groups": old["fold_groups"],
                "reserved_groups": old["reserved_groups"], "fit_budget": 3,
                "screening_gate": {"mean_fold_ap_gain": .01, "no_fold_ap_loss": True,
                                   "no_high_hit_loss": True},
                "bootstrap": {"draws": 1000, "seed": 17}, "protected_sha256": protected,
                "scope": "One technical ranking sensitivity on reused development data"}
    output.mkdir(parents=True)
    write_json(output / "protocol.json", protocol)
    protocol_hash = sha256(output / "protocol.json")
    write_json(output / "status.json", {"status": "started", "protocol_sha256": protocol_hash})
    try:
        activation, _ = load_development_activation(root / WORKBOOK, cohort)
        if set(activation) != {c["candidate_id"] for c in cohort}:
            raise ValueError("Incomplete approved activation cohort")
        previous = json.loads((root / PREVIOUS / "private-predictions.json").read_text())
        cached = previous["esm_biology_logistic"]
        states = load_features(cohort, root / REFERENCES, root / ESM_CACHE)
        store = BiologicalFeatureStore(root / "data/revision-v2/features")
        augmenter = ActivationFeatureAugmenter(store, cohort, root / REFERENCES, root / SCORES)
        rows, diagnostics = {"probability": [], "margin": []}, []
        fitting_start = time.monotonic()
        for fold, train, test in fold_partitions(cohort, old["fold_groups"]):
            records = {c["candidate_id"]: activation[c["candidate_id"]] for c in train}
            model = FGFR2SequenceModel("esm_biology_logistic", store, augmenter, states).fit(train, records)
            predictions = ranking_outputs(model, test)
            p, m = predictions["probability"], predictions["margin"]
            diagnostics.append({"fold": fold, "n": len(test), "probability_unique": len(np.unique(p)),
                "margin_unique": len(np.unique(m)), "probability_exact_one": int((p == 1).sum()),
                "margin_min": float(m.min()), "margin_max": float(m.max()),
                "model": model.diagnostics})
            for arm, scores in predictions.items():
                for c, score in zip(test, scores, strict=True):
                    r = activation[c["candidate_id"]]
                    rows[arm].append({"candidate_id": c["candidate_id"], "group_id": c["group_id"],
                        "fold": fold, "truth": r["value"], "high": activation_label(r),
                        "strong": activation_label(r, 10.), "score": float(score)})
        fitting_seconds = time.monotonic() - fitting_start
        reproduction = verify_cached_predictions(rows["probability"], cached)
        metrics = {arm: summarize(values) for arm, values in rows.items()}
        metrics["historical_tree_reference"] = summarize(previous["cached_v2_biology"])
        uncertainty = paired_ap_interval(rows)
        changed = [p for p, h in protected.items() if sha256(root / p) != h]
        if changed or sha256(output / "protocol.json") != protocol_hash:
            raise ValueError(f"Protected artifacts changed: {changed}")
        result = {"metrics": metrics, "gate": screening_gate(metrics["probability"], metrics["margin"]),
                  "uncertainty": uncertainty, "baseline_reproduction": reproduction,
                  "diagnostics": diagnostics, "fit_count": 3, "fitting_seconds": fitting_seconds,
                  "elapsed_seconds": time.monotonic() - start, "protocol_sha256": protocol_hash,
                  "integrity": {"protected_files": len(protected), "changed": changed}}
        write_json(output / "private-predictions.json", rows)
        (output / "private-predictions.json").chmod(0o600)
        write_json(output / "results.json", result)
        (output / "report.md").write_text(make_report(result))
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
    args = parser.parse_args()
    result = run(args.root, args.output)
    print(json.dumps({k: result[k] for k in ("gate", "baseline_reproduction", "fit_count",
                                          "fitting_seconds", "elapsed_seconds")}, indent=2))
