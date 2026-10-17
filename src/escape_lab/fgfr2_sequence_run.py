"""Fixed FGFR2 activation retrieval pilot; no external or resistance evaluation."""

import argparse
import json
import math
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score

from escape_lab.activation_data import load_development_activation
from escape_lab.activation_features import ActivationFeatureAugmenter
from escape_lab.activation_revision_run import SCORES
from escape_lab.benchmark import validate_partition
from escape_lab.biological_features import BiologicalFeatureStore
from escape_lab.data import sha256
from escape_lab.improvement import write_json
from escape_lab.mechanistic_pilot import eligible_candidates
from escape_lab.mechanistic_pilot_run import MANO_AUDIT, RLY_AUDIT, WORKBOOK

V2 = Path("reports/mechanistic-pilot/activation-revision-v2")
RESEARCH = Path("research/fgfr2-sequence-pilot")
REFERENCES = Path("data/revision-v2/features/references.json")
ARMS = ("cached_v2_biology", "biology_logistic", "esm_biology_logistic", "esm_biology_ridge")
PRIMARY = "esm_biology_logistic"


def fold_partitions(candidates, mapping):
    """Subset the saved global folds, never recalculate them on the gene subset."""
    if not candidates or any(c["gene"] != "FGFR2" or c["split"] not in
                             ("development", "validation") for c in candidates):
        raise ValueError("Only FGFR2 development candidates are allowed")
    if len({c["candidate_id"] for c in candidates}) != len(candidates):
        raise ValueError("Duplicate candidates")
    groups = [g for values in mapping.values() for g in values]
    if len(set(groups)) != len(groups) or len(mapping) != 3:
        raise ValueError("Saved mapping must contain three disjoint folds")
    if {c["group_id"] for c in candidates} - set(groups):
        raise ValueError("Candidate group missing from saved folds")
    tasks = []
    for fold, held in sorted(mapping.items()):
        train = [c for c in candidates if c["group_id"] not in held]
        test = [c for c in candidates if c["group_id"] in held]
        if not train or not test:
            raise ValueError("Empty training or scoring fold")
        validate_partition(train, test)
        tasks.append((fold, train, test))
    return tasks


def ranking_metrics(rows):
    rows = sorted(rows, key=lambda r: r["candidate_id"])
    if not rows:
        raise ValueError("Cannot score an empty fold")
    scores = np.asarray([r["score"] for r in rows])
    if not np.isfinite(scores).all():
        raise ValueError("Nonfinite predictions")
    k = math.ceil(.1 * len(rows))
    selected = np.argsort(-scores, kind="stable")[:k]
    result = {"n": len(rows), "budget": k}
    for endpoint in ("high", "strong"):
        y = np.asarray([r[endpoint] for r in rows], dtype=int)
        hits = int(y[selected].sum())
        result[endpoint] = {
            "positives": int(y.sum()), "prevalence": float(y.mean()),
            "ap": float(average_precision_score(y, scores)) if y.sum() else None,
            "hits": hits, "recall": hits / int(y.sum()) if y.sum() else None,
            "precision": hits / k,
        }
    truth = np.asarray([r["truth"] for r in rows])
    result["spearman"] = float(spearmanr(truth, scores).statistic) if (
        np.ptp(truth) and np.ptp(scores)) else None
    return result


def summarize(rows):
    by_fold = {fold: ranking_metrics([r for r in rows if r["fold"] == fold])
               for fold in sorted({r["fold"] for r in rows})}
    aps = [m["high"]["ap"] for m in by_fold.values()]
    hits = sum(m["high"]["hits"] for m in by_fold.values())
    strong_hits = sum(m["strong"]["hits"] for m in by_fold.values())
    budget = sum(m["budget"] for m in by_fold.values())
    positives = sum(m["high"]["positives"] for m in by_fold.values())
    return {"by_fold": by_fold, "mean_fold_high_ap": float(np.mean(aps))
            if all(a is not None for a in aps) else None,
            "high_hits": hits, "strong_hits": strong_hits, "budget": budget,
            "high_recall": hits / positives, "high_precision": hits / budget,
            "random_expected_hits": sum(m["budget"] * m["high"]["prevalence"]
                                        for m in by_fold.values()),
            "pooled_diagnostic": ranking_metrics(rows)}


def advancement_gate(metrics):
    a, b, c = (metrics[name] for name in (PRIMARY, "cached_v2_biology", "biology_logistic"))
    aps = [m["mean_fold_high_ap"] for m in (a, b, c)]
    if any(x is None for x in aps):
        return {"status": "inconclusive", "reason": "Undefined primary metric"}
    fold_wins = sum(a["by_fold"][f]["high"]["ap"] > b["by_fold"][f]["high"]["ap"]
                    for f in a["by_fold"])
    checks = {"three_extra_hits_vs_v2": a["high_hits"] >= b["high_hits"] + 3,
              "ap_gain_003_vs_v2": aps[0] - aps[1] >= .03 - 1e-12,
              "ap_gain_003_vs_biology": aps[0] - aps[2] >= .03 - 1e-12,
              "no_hit_loss_vs_biology": a["high_hits"] >= c["high_hits"],
              "ap_improves_two_folds_vs_v2": fold_wins >= 2}
    return {"status": "pass" if all(checks.values()) else "fail", "checks": checks,
            "high_hits_gain_vs_v2": a["high_hits"] - b["high_hits"],
            "mean_fold_ap_gain_vs_v2": aps[0] - aps[1],
            "mean_fold_ap_gain_vs_biology": aps[0] - aps[2], "fold_wins": fold_wins}


def validate_reference(candidates, activation, saved, group_to_fold):
    ids = {c["candidate_id"] for c in candidates}
    rows = {r["candidate_id"]: r for r in saved if r["candidate_id"] in ids}
    if set(rows) != ids or len(rows) != len([r for r in saved if r["candidate_id"] in ids]):
        raise ValueError("Cached reference cohort mismatch")
    for c in candidates:
        cid = c["candidate_id"]
        r, target = rows[cid], activation[cid]
        if (r["gene"] != "FGFR2" or r["group_id"] != c["group_id"] or
                r["fold"] != group_to_fold[c["group_id"]] or
                r["truth"] != target["value"] or r["active"] != target["active"] or
                not math.isfinite(r["prediction"])):
            raise ValueError("Cached reference identity, fold, or target mismatch")
    return {cid: r["prediction"] for cid, r in rows.items()}


def feature_candidates(root):
    """IDs only from prior complete cohort; no outcome-derived feature selection."""
    root = Path(root)
    old = json.loads((root / V2 / "protocol.json").read_text())
    candidates = json.loads((root / "data/processed/candidates.json").read_text())
    mano, rly = (json.loads((root / p).read_text()) for p in (MANO_AUDIT, RLY_AUDIT))
    reserved = set(mano["mapped_groups"]) | {r["group_id"] for r in rly["panel_variants"]["rly2023"]}
    eligible = eligible_candidates(candidates, reserved)
    if [c["candidate_id"] for c in eligible] != old["eligible_ids"]:
        raise ValueError("Development eligibility changed")
    previous = json.loads((root / V2 / "private-predictions.json").read_text())
    saved = previous["activation_records"]["biology"]
    complete = {r["candidate_id"] for r in saved if r["gene"] == "FGFR2"}
    cohort = [c for c in eligible if c["gene"] == "FGFR2" and c["candidate_id"] in complete]
    if len(cohort) != 668 or len({c["group_id"] for c in cohort}) != 10:
        raise ValueError("Expected 668 complete FGFR2 variants across ten occupied blocks")
    fold_partitions(cohort, old["fold_groups"])
    return cohort, old, saved


def paired_uncertainty(rows, draws=1000):
    """Resample blocks within each fixed outer fold; no refitting or model selection."""
    rng = np.random.default_rng(17)
    reference = {r["candidate_id"]: r for r in rows["cached_v2_biology"]}
    contrasts = {name: [] for name in ("cached_v2_biology", "biology_logistic")}
    folds = sorted({r["fold"] for r in reference.values()})
    grouped = {fold: sorted({r["group_id"] for r in reference.values() if r["fold"] == fold})
               for fold in folds}
    arrays = {}
    for arm, records in rows.items():
        if (len(records) != len(reference) or
                {r["candidate_id"] for r in records} != set(reference)):
            raise ValueError("Unpaired bootstrap")
        if any(any(r[key] != reference[r["candidate_id"]][key]
                   for key in ("fold", "group_id", "high", "strong", "truth")) for r in records):
            raise ValueError("Bootstrap target or partition mismatch")
        arrays[arm] = {fold: sorted([r for r in records if r["fold"] == fold],
                                    key=lambda r: r["candidate_id"]) for fold in folds}
    for _ in range(draws):
        weights = {fold: dict(zip(grouped[fold], np.bincount(
            rng.integers(len(grouped[fold]), size=len(grouped[fold])),
            minlength=len(grouped[fold])), strict=True)) for fold in folds}
        scores = {}
        for arm in (PRIMARY, *contrasts):
            values = []
            for fold in folds:
                data = arrays[arm][fold]
                w = np.asarray([weights[fold][r["group_id"]] for r in data])
                y = np.asarray([r["high"] for r in data])
                if np.dot(w, y) == 0:
                    break
                values.append(float(average_precision_score(y, [r["score"] for r in data],
                                                             sample_weight=w)))
            scores[arm] = float(np.mean(values)) if len(values) == len(folds) else None
        for ref, values in contrasts.items():
            if scores[PRIMARY] is not None and scores[ref] is not None:
                values.append(scores[PRIMARY] - scores[ref])
    return {"draws": draws, "seed": 17, "unit": "aligned block within original fold",
            "mean_fold_ap_gain": {ref: {"valid_draws": len(v), "invalid_draws": draws-len(v),
                "percentile_2_5": float(np.percentile(v, 2.5)) if v else None,
                "percentile_97_5": float(np.percentile(v, 97.5)) if v else None}
                for ref, v in contrasts.items()},
            "interpretation": "Descriptive fixed-prediction uncertainty on ten reused blocks; "
                              "not independent validation or a model-selection correction."}


def report_text(result):
    lines = ["# FGFR2 sequence representation pilot", "",
             f"Primary development gate: **{result['gate']['status']}**.", "",
             "668 FGFR2 variants; 26 high-activation positives; ten occupied position blocks.",
             ("High means >=3/4 replicate medians >=1.5 AND median >5. This is a protein-level "
             "adaptation of the study categories. Strong (>10) is a secondary endpoint."), "",
             "| Model | Mean fold AP | High hits / 68 picks | Recall of 26 | Strong hits / 7 |",
             "|---|---:|---:|---:|---:|"]
    for arm in ARMS:
        m = result["metrics"][arm]
        lines.append(f"| {arm} | {m['mean_fold_high_ap']:.4f} | {m['high_hits']} | "
                     f"{m['high_recall']:.1%} | {m['strong_hits']} |")
    lines.extend(["", ("AP is ranking quality, not accuracy. Selections are top 10% within each "
                  "fold (28+18+22=68), with candidate-ID tie breaking. Pooled scores are secondary."),
                  "", "```json", json.dumps(result["gate"], indent=2), "```", "",
                  "## Per-fold results", "",
                  "| Fold | Model | AP | High hits | Positive / N |",
                  "|---|---|---:|---:|---:|"])
    for arm in ARMS:
        for fold, m in result["metrics"][arm]["by_fold"].items():
            lines.append(f"| {fold} | {arm} | {m['high']['ap']:.4f} | {m['high']['hits']} | "
                         f"{m['high']['positives']} / {m['n']} |")
    lines.extend(["", (f"Head fitting/scoring: {result['fitting_seconds']:.2f}s; "
                  f"feature extraction: {result['feature_extraction_seconds']:.2f}s."), "",
                  "- ESM2-35M remains frozen. No LLM calls or supervised PLM fine-tuning.",
                  "- Scalers, imputers and separate PCA16 branches fit only each training fold.",
                  "- Class-balanced logistic outputs are ranking scores, not calibrated probabilities.",
                  ("- Cached V2 reference was trained across genes. The matched biology-only classifier "
                  "isolates the representation change in FGFR2-only training."),
                  ("- No external assay outcomes, original final evaluation, resistance labels or "
                  "post-score hyperparameter tuning used. These development folds have been reused."),
                  "- Strong activation has no positives in pilot-1; its AP there is undefined.",
                  f"- All {result['integrity']['protected_files']} protected artifacts hash-identical.",
                  "", "See protocol.json, results.json and feature manifest for full provenance.", ""])
    return "\n".join(lines)


def run(root, output, feature_dir):
    from escape_lab.fgfr2_sequence_features import load_features
    from escape_lab.fgfr2_sequence_models import (
        MODEL_CONFIG,
        TRAINABLE_ARMS,
        FGFR2SequenceModel,
        activation_label,
    )
    root, output, feature_dir = (Path(p).resolve() for p in (root, output, feature_dir))
    if output.exists():
        raise FileExistsError("Use a new run directory; preserve failed attempts")
    started = time.monotonic()
    cohort, old, saved = feature_candidates(root)
    protected = dict(old["protected_sha256"])
    changed = [p for p, h in protected.items() if sha256(root / p) != h]
    if changed:
        raise ValueError(f"Previous protected artifacts changed: {changed}")
    added = list((root / V2).glob("*")) + list((root / RESEARCH).glob("*"))
    added += list(feature_dir.glob("*"))
    added += list((root / "src/escape_lab").glob("fgfr2_sequence_*.py"))
    added += list((root / "tests").glob("test_fgfr2_sequence_*.py"))
    added += [root / "pyproject.toml", root / "uv.lock"]
    for p in added:
        if p.is_file():
            protected[str(p.relative_to(root))] = sha256(p)
    protocol = {"created_at": datetime.now(UTC).isoformat(), "version": "fgfr2-sequence-v1",
                "arms": ARMS, "primary": PRIMARY, "model_config": MODEL_CONFIG,
                "fold_groups": old["fold_groups"], "eligible_ids": [c["candidate_id"] for c in cohort],
                "reserved_groups": old["reserved_groups"], "high": "3/4 reps>=1.5 AND median>5",
                "strong_secondary": "3/4 reps>=1.5 AND median>10", "budget": "ceil(0.1*n) per fold",
                "ties": "candidate_id ascending", "gate": {"hits_gain_vs_v2": 3,
                "mean_fold_ap_gain_vs_both": .03, "hits_no_loss_vs_biology": True,
                "ap_fold_wins_vs_v2": 2}, "protected_sha256": protected,
                "scope": "Reused development only; no external or resistance testing"}
    output.mkdir(parents=True)
    write_json(output / "protocol.json", protocol)
    protocol_hash = sha256(output / "protocol.json")
    write_json(output / "status.json", {"status": "started", "protocol_sha256": protocol_hash})
    try:
        activation, coverage = load_development_activation(root / WORKBOOK, cohort)
        if set(activation) != {c["candidate_id"] for c in cohort}:
            raise ValueError("Missing activation targets")
        if sum(activation_label(r) for r in activation.values()) != 26:
            raise ValueError("Unexpected high-activation target count")
        group_map = {g: fold for fold, groups in old["fold_groups"].items() for g in groups}
        reference = validate_reference(cohort, activation, saved, group_map)
        sequence = load_features(cohort, root / REFERENCES, feature_dir)
        store = BiologicalFeatureStore(root / "data/revision-v2/features")
        augmenter = ActivationFeatureAugmenter(store, cohort, root / REFERENCES, root / SCORES)
        rows = {arm: [] for arm in ARMS}
        diagnostics = []
        fitting_started = time.monotonic()
        for fold, train, test in fold_partitions(cohort, old["fold_groups"]):
            records = {c["candidate_id"]: activation[c["candidate_id"]] for c in train}
            predictions = {"cached_v2_biology": np.asarray([reference[c["candidate_id"]] for c in test])}
            for arm in TRAINABLE_ARMS:
                model = FGFR2SequenceModel(arm, store, augmenter, sequence).fit(train, records)
                predictions[arm] = model.predict(test)
                diagnostics.append({"arm": arm, "fold": fold, **model.diagnostics})
            for arm, scores in predictions.items():
                for c, score in zip(test, scores, strict=True):
                    r = activation[c["candidate_id"]]
                    rows[arm].append({"candidate_id": c["candidate_id"], "group_id": c["group_id"],
                                      "fold": fold, "score": float(score), "truth": r["value"],
                                      "high": activation_label(r),
                                      "strong": activation_label(r, threshold=10.)})
        fitting_seconds = time.monotonic() - fitting_started
        metrics = {arm: summarize(r) for arm, r in rows.items()}
        uncertainty = paired_uncertainty(rows)
        changed = [p for p, h in protected.items() if sha256(root / p) != h]
        if changed or sha256(output / "protocol.json") != protocol_hash:
            raise ValueError(f"Frozen artifacts changed: {changed}")
        manifest = sequence.manifest
        result = {"metrics": metrics, "gate": advancement_gate(metrics), "coverage": coverage,
                  "uncertainty": uncertainty, "model_diagnostics": diagnostics,
                  "feature_extraction_seconds": manifest.get("elapsed_seconds", 0.),
                  "fitting_seconds": fitting_seconds, "elapsed_seconds": time.monotonic()-started,
                  "integrity": {"protected_files": len(protected), "changed": changed},
                  "protocol_sha256": protocol_hash, "feature_manifest": str(feature_dir / "manifest.json"),
                  "reference_validation": "Exact cached IDs, folds, continuous and active targets match"}
        write_json(output / "private-predictions.json", rows)
        (output / "private-predictions.json").chmod(0o600)
        write_json(output / "results.json", result)
        (output / "report.md").write_text(report_text(result))
        write_json(output / "status.json", {"status": "completed", "protocol_sha256": protocol_hash,
                                           "results_sha256": sha256(output / "results.json")})
        return result
    except Exception as error:
        write_json(output / "status.json", {"status": "failed", "error_type": type(error).__name__,
                                           "error": str(error), "protocol_sha256": protocol_hash})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.root, args.output, args.features)
    print(json.dumps({"gate": result["gate"], "fitting_seconds": result["fitting_seconds"]}, indent=2))
