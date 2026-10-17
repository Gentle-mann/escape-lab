"""Six-fit development smoke test of frozen FGFR2 state features in pooled HGB."""

import argparse
import json
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from escape_lab.activation_features import ActivationFeatureAugmenter
from escape_lab.activation_revision import CAPACITY
from escape_lab.activation_revision_run import SCORES
from escape_lab.benchmark import validate_partition
from escape_lab.biological_features import BiologicalFeatureStore
from escape_lab.data import sha256
from escape_lab.fgfr2_sequence_run import REFERENCES, feature_candidates
from escape_lab.fgfr2_state_features import FGFR2StateFeatureStore
from escape_lab.fgfr2_state_run import PREVIOUS, continuous_summary, verify_cached_predictions
from escape_lab.improvement import write_json

PRIOR = Path("reports/mechanistic-pilot/fgfr2-state-bounded-v1")
FEATURES = Path("data/fgfr2-state-pilot/paired-v1")
RESEARCH = Path("research/fgfr2-state-tree-pilot")
ARMS = ("baseline", "contrast")


def fixed_partitions(cohort, mapping):
    if (not cohort or len({c["candidate_id"] for c in cohort}) != len(cohort)
            or any(c["split"] not in ("development", "validation") for c in cohort)):
        raise ValueError("Unique development cohort required")
    groups = [g for values in mapping.values() for g in values]
    if (len(mapping) != 3 or len(groups) != len(set(groups))
            or {c["group_id"] for c in cohort} != set(groups)):
        raise ValueError("Original complete three-fold mapping required")
    tasks = []
    for fold, held in sorted(mapping.items()):
        train = [c for c in cohort if c["group_id"] not in held]
        test = [c for c in cohort if c["group_id"] in held and c["gene"] == "FGFR2"]
        validate_partition(train, test)
        if not train or not test:
            raise ValueError("Empty fold")
        tasks.append((fold, train, test))
    return tasks


def screening_gate(metrics):
    b, a = metrics["baseline"], metrics["contrast"]
    if any(m["mean_fold_spearman"] is None for m in (a, b)):
        return {"status": "inconclusive", "reason": "Undefined correlation"}
    gain = a["mean_fold_spearman"] - b["mean_fold_spearman"]
    wins = sum(a["by_fold"][f]["spearman"] > b["by_fold"][f]["spearman"] for f in b["by_fold"])
    checks = {"rho_gain_at_least_003": gain >= .03 - 1e-12,
              "rho_improves_two_folds": wins >= 2,
              "mae_within_5pct": a["mean_fold_mae"] <= 1.05 * b["mean_fold_mae"] + 1e-12,
              "no_high_hit_loss": a["high_hits"] >= b["high_hits"]}
    return {"status": "pass" if all(checks.values()) else "fail", "checks": checks,
            "rho_gain": gain, "rho_fold_wins": wins,
            "mae_reduction": b["mean_fold_mae"] - a["mean_fold_mae"],
            "high_hits_gain": a["high_hits"] - b["high_hits"]}


def paired_uncertainty(rows, draws=1000):
    baseline = sorted(rows["baseline"], key=lambda r: r["candidate_id"])
    contrast = sorted(rows["contrast"], key=lambda r: r["candidate_id"])
    if (not baseline or len(baseline) != len(contrast)
            or len({r["candidate_id"] for r in baseline}) != len(baseline)
            or any(any(a[k] != b[k] for k in ("candidate_id", "group_id", "fold", "truth",
                                             "high", "strong"))
                   for a, b in zip(baseline, contrast, strict=True))):
        raise ValueError("Unpaired bootstrap rows")
    arrays = []
    for fold in sorted({r["fold"] for r in baseline}):
        idx = [i for i, r in enumerate(baseline) if r["fold"] == fold]
        records = [baseline[i] for i in idx]
        groups = sorted({r["group_id"] for r in records})
        blocks = [np.array([i for i, r in enumerate(records) if r["group_id"] == g]) for g in groups]
        arrays.append((blocks, np.array([r["truth"] for r in records]),
                       np.array([r["score"] for r in records]),
                       np.array([contrast[i]["score"] for i in idx])))
    rng = np.random.default_rng(17)
    values = {"rho_gain": [], "mae_reduction": []}
    for _ in range(draws):
        rhos, errors = [], []
        for blocks, y, b, a in arrays:
            indices = np.concatenate([blocks[i] for i in rng.integers(len(blocks), size=len(blocks))])
            truth, ref, new = y[indices], b[indices], a[indices]
            rhos.append(float(spearmanr(truth, new).statistic - spearmanr(truth, ref).statistic)
                        if np.ptp(truth) and np.ptp(ref) and np.ptp(new) else None)
            errors.append(float(np.abs(truth-ref).mean() - np.abs(truth-new).mean()))
        if all(v is not None for v in rhos):
            values["rho_gain"].append(float(np.mean(rhos)))
        values["mae_reduction"].append(float(np.mean(errors)))
    return {"draws": draws, "seed": 17, "unit": "aligned block within original fold",
            "contrasts": {k: {"valid_draws": len(v), "invalid_draws": draws-len(v),
                "percentile_2_5": float(np.percentile(v, 2.5)) if v else None,
                "percentile_97_5": float(np.percentile(v, 97.5)) if v else None}
                for k, v in values.items()},
            "interpretation": "Descriptive fixed-prediction resampling of reused development data; "
                              "no refitting or external validation."}


def report_text(result):
    lines = ["# Tiny FGFR2 structure-contrast tree pilot", "",
             f"Fixed development screening gate: **{result['gate']['status']}**.", "",
             ("Six fits only: pooled FGFR1–4 biology HGB baseline and the same HGB with cached FGFR2 "
             "active-minus-inactive descriptors, each on the three original folds. No tuning, agent "
             "loop, nested crossfits, resistance evaluation or external assay reads."), "",
             "| Model | Mean-fold Spearman ↑ | Mean-fold MAE ↓ | High hits /26 | Strong hits /7 |",
             "|---|---:|---:|---:|---:|"]
    for arm, m in result["metrics"].items():
        lines.append(f"| {arm} | {m['mean_fold_spearman']:.5f} | {m['mean_fold_mae']:.5f} | "
                     f"{m['high_hits']} | {m['strong_hits']} |")
    lines.extend(["", ("668 FGFR2 predictions; 68 total nominations. MAE uses the unchanged "
                  "log2(1 + activation) target. Spearman is not percentage accuracy."), "",
                  "| Fold | Baseline rho | Contrast rho | Baseline MAE | Contrast MAE |",
                  "|---|---:|---:|---:|---:|"])
    for fold in result["metrics"]["baseline"]["by_fold"]:
        b, a = [result["metrics"][arm]["by_fold"][fold] for arm in ARMS]
        lines.append(f"| {fold} | {b['spearman']:.5f} | {a['spearman']:.5f} | "
                     f"{b['mae']:.5f} | {a['mae']:.5f} |")
    lines.extend(["", "## Gate and descriptive uncertainty", "", "```json",
                  json.dumps({"gate": result["gate"], "uncertainty": result["uncertainty"]}, indent=2),
                  "```", "", "## Verification and limits", "",
                  f"- Baseline reproduction: {result['baseline_reproduction']}.",
                  f"- {result['integrity']['protected_files']} protected files unchanged.",
                  f"- Six fits/scoring: {result['fitting_seconds']:.3f}s; total: {result['elapsed_seconds']:.3f}s.",
                  "- Cached development labels only; no workbook ingestion or new outcome datasets.",
                  "- Same training rows/order, HGB hyperparameters and train-only imputation in both arms.",
                  ("- Additional features are available only for FGFR2; other genes have zero geometry "
                  "and availability=0. No new cross-gene structural mapping was introduced."),
                  "- No ESM, standardization, PCA, range clamp, output clipping or new biological features.",
                  "- Reused development folds cannot establish generalization; all external assays remain reserved.",
                  "- Structural features remain static prepared-coordinate proxies, not mutant energy or activity simulations.", ""])
    return "\n".join(lines)


def run(root, output):
    from escape_lab.fgfr2_state_tree import FGFR2StateTreeModel

    root, output = Path(root).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError("Use a fresh pilot directory")
    start = time.monotonic()
    fgfr2, old, saved = feature_candidates(root)
    saved_by_id = {r["candidate_id"]: r for r in saved}
    if len(saved_by_id) != len(saved) or len(saved) != 2719:
        raise ValueError("Historical complete pooled cohort mismatch")
    all_candidates = {c["candidate_id"]: c for c in json.loads(
        (root / "data/processed/candidates.json").read_text())}
    cohort = [all_candidates[cid] for cid in old["eligible_ids"] if cid in saved_by_id]
    if {c["candidate_id"] for c in cohort} != set(saved_by_id):
        raise ValueError("Historical training eligibility mismatch")
    tasks = fixed_partitions(cohort, old["fold_groups"])
    if {c["candidate_id"] for _, _, test in tasks for c in test} != {c["candidate_id"] for c in fgfr2}:
        raise ValueError("FGFR2 test cohort changed")
    group_fold = {g: f for f, groups in old["fold_groups"].items() for g in groups}
    for c in cohort:
        r = saved_by_id[c["candidate_id"]]
        if any(r[k] != c[k] for k in ("gene", "group_id")) or r["fold"] != group_fold[c["group_id"]]:
            raise ValueError("Cached training identity mismatch")
    prior_protocol = json.loads((root / PRIOR / "protocol.json").read_text())
    protected = dict(prior_protocol["protected_sha256"])
    changed = [p for p, h in protected.items() if sha256(root / p) != h]
    if changed:
        raise ValueError(f"Protected files changed: {changed}")
    added = list((root / PRIOR).glob("*")) + list((root / RESEARCH).glob("*"))
    added += list((root / "src/escape_lab").glob("fgfr2_state_tree*.py"))
    added += list((root / "tests").glob("test_fgfr2_state_tree*.py"))
    for path in added:
        if path.is_file():
            protected[str(path.relative_to(root))] = sha256(path)
    protocol = {"version": "fgfr2-state-tree-v1", "created_at": datetime.now(UTC).isoformat(),
                "arms": ARMS, "primary": "contrast", "hgb_capacity": CAPACITY,
                "pooled_training_cohort": [c["candidate_id"] for c in cohort],
                "scoring_cohort": [c["candidate_id"] for c in fgfr2],
                "fold_groups": old["fold_groups"], "reserved_groups": old["reserved_groups"],
                "screening_gate": {"rho_gain": .03, "fold_wins": 2, "max_mae_ratio": 1.05,
                                   "no_high_hit_loss": True},
                "fit_budget": 6, "bootstrap": {"draws": 1000, "seed": 17},
                "source": "Previously evaluated cached development activation labels and structural features",
                "protected_sha256": protected}
    output.mkdir(parents=True)
    write_json(output / "protocol.json", protocol)
    protocol_hash = sha256(output / "protocol.json")
    write_json(output / "status.json", {"status": "started", "protocol_sha256": protocol_hash})
    try:
        cached = json.loads((root / PREVIOUS / "private-predictions.json").read_text())["cached_v2_biology"]
        truth = {r["candidate_id"]: r for r in cached}
        if len(truth) != len(fgfr2) or any(truth[c["candidate_id"]]["truth"] !=
                                          saved_by_id[c["candidate_id"]]["truth"] for c in fgfr2):
            raise ValueError("Cached FGFR2 targets mismatch")
        store = BiologicalFeatureStore(root / "data/revision-v2/features")
        augmenter = ActivationFeatureAugmenter(store, cohort, root / REFERENCES, root / SCORES)
        states = FGFR2StateFeatureStore(root / FEATURES)
        rows, diagnostics, reproduction = {arm: [] for arm in ARMS}, [], None
        fitting_start = time.monotonic()
        for arm in ARMS:
            if arm != "baseline" and reproduction is None:
                raise ValueError("Baseline must reproduce before challenger fits")
            for fold, train, test in tasks:
                labels = {c["candidate_id"]: {"value": saved_by_id[c["candidate_id"]]["truth"],
                    "active": saved_by_id[c["candidate_id"]]["active"]} for c in train}
                model = FGFR2StateTreeModel(store, augmenter, states, add_contrast=arm == "contrast")
                model.fit(train, labels)
                scores = model.predict(test)["value"]
                diagnostics.append({"arm": arm, "fold": fold, "training_n": len(train),
                                    "test_n": len(test),
                                    "features": len(model.vectorizer.get_feature_names_out()),
                                    "hgb_params": model.regressor.get_params()})
                for c, score in zip(test, scores, strict=True):
                    r = truth[c["candidate_id"]]
                    if r["fold"] != fold or r["group_id"] != c["group_id"]:
                        raise ValueError("Scored cached fold mismatch")
                    rows[arm].append({**r, "score": float(score)})
            if arm == "baseline":
                reproduction = verify_cached_predictions(rows[arm], cached)
        fitting_seconds = time.monotonic() - fitting_start
        metrics = {arm: continuous_summary(values) for arm, values in rows.items()}
        uncertainty = paired_uncertainty(rows)
        changed = [p for p, h in protected.items() if sha256(root / p) != h]
        if changed or sha256(output / "protocol.json") != protocol_hash:
            raise ValueError(f"Protected artifacts changed: {changed}")
        result = {"metrics": metrics, "gate": screening_gate(metrics), "uncertainty": uncertainty,
                  "baseline_reproduction": reproduction, "model_diagnostics": diagnostics,
                  "fit_count": len(diagnostics), "fitting_seconds": fitting_seconds,
                  "elapsed_seconds": time.monotonic() - start, "protocol_sha256": protocol_hash,
                  "integrity": {"protected_files": len(protected), "changed": changed}}
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
    args = parser.parse_args()
    result = run(args.root, args.output)
    print(json.dumps({k: result[k] for k in ("gate", "baseline_reproduction", "fit_count",
                                          "fitting_seconds", "elapsed_seconds")}, indent=2))
