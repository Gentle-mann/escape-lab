"""One explicitly post-result, matched feature-extrapolation sensitivity check."""

import argparse
import json
import time
from datetime import UTC, datetime
from pathlib import Path

from escape_lab.activation_data import load_development_activation
from escape_lab.activation_features import ActivationFeatureAugmenter
from escape_lab.activation_revision_run import SCORES
from escape_lab.biological_features import BiologicalFeatureStore
from escape_lab.data import sha256
from escape_lab.fgfr2_sequence_features import load_features
from escape_lab.fgfr2_sequence_models import activation_label
from escape_lab.fgfr2_sequence_run import REFERENCES, feature_candidates, fold_partitions
from escape_lab.fgfr2_state_bounded_model import BoundedFGFR2StateModel
from escape_lab.fgfr2_state_features import FGFR2StateFeatureStore
from escape_lab.fgfr2_state_models import ARMS, MODEL_CONFIG, PRIMARY
from escape_lab.fgfr2_state_run import (
    ESM_CACHE,
    HISTORICAL,
    PREVIOUS,
    compare_gates,
    continuous_summary,
    paired_uncertainty,
    verify_cached_predictions,
)
from escape_lab.improvement import write_json
from escape_lab.mechanistic_pilot_run import WORKBOOK

FIRST = Path("reports/mechanistic-pilot/fgfr2-state-v1")
FEATURES = Path("data/fgfr2-state-pilot/paired-v1")
RESEARCH = Path("research/fgfr2-state-bounded")


def report_text(result):
    lines = ["# FGFR2 bounded-feature sensitivity", "",
             "**Post-result development sensitivity. The original failed run remains preserved.**", "",
             ("Biology features are limited to their training ranges at inference. Training matrices, "
             "coefficients, sequence PCA and all four feature arms are unchanged."), "",
             (f"Incremental gate: **{result['gates']['incremental']['status']}**; "
             f"state-contrast check: **{result['gates']['state_contrast']['status']}**."), "",
             "| Arm | Mean-fold Spearman | Mean-fold MAE | High hits /26 at 68 picks | Strong hits /7 |",
             "|---|---:|---:|---:|---:|"]
    for arm, m in result["metrics"].items():
        rho = "undefined" if m["mean_fold_spearman"] is None else f"{m['mean_fold_spearman']:.5f}"
        lines.append(f"| {arm} | {rho} | {m['mean_fold_mae']:.5f} | {m['high_hits']} | {m['strong_hits']} |")
    lines.extend(["", ("Spearman measures ranking, not percentage accuracy. MAE is on the unchanged "
                  "log2(1 + activation) target. Historical HGB is an unchanged cached reference."), "",
                  "## Fixed gates", "", "```json", json.dumps(result["gates"], indent=2), "```", "",
                  "## Per-fold results", "",
                  "| Fold | Arm | Spearman | MAE | High hits |", "|---|---|---:|---:|---:|"])
    for arm, m in result["metrics"].items():
        for fold, v in m["by_fold"].items():
            rho = "undefined" if v["spearman"] is None else f"{v['spearman']:.5f}"
            lines.append(f"| {fold} | {arm} | {rho} | {v['mae']:.5f} | {v['high']['hits']} |")
    lines.extend(["", "## Verification and limits", "",
                  f"- Unbounded score reproduction for every arm: {result['unbounded_reproduction']}.",
                  f"- Protected files unchanged: {result['integrity']['protected_files']}.",
                  f"- Fitting/scoring: {result['fitting_seconds']:.3f}s; total: {result['elapsed_seconds']:.3f}s.",
                  "- Exactly the same 668 variants and three original folds; no new assay outcomes.",
                  "- Training-range bounds are fitted without labels and never updated using test rows.",
                  "- No parameter search, output clipping, external evaluation or selection of a new primary arm.",
                  "- Structural inputs remain prepared C491A-background proxies; inactive is a homology model.",
                  ("- Static contacts are not free energies or mutant relaxation. Masking C491A cannot remove "
                  "possible long-range effects of that shared engineered background."),
                  "- Block-resampling intervals are descriptive on reused development predictions.", ""])
    return "\n".join(lines)


def run(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError("Use a fresh sensitivity directory")
    started = time.monotonic()
    cohort, old, _ = feature_candidates(root)
    first_protocol = json.loads((root / FIRST / "protocol.json").read_text())
    if (first_protocol["eligible_ids"] != [c["candidate_id"] for c in cohort]
            or first_protocol["fold_groups"] != old["fold_groups"]):
        raise ValueError("First-run cohort or fold mismatch")
    protected = dict(first_protocol["protected_sha256"])
    changed = [p for p, h in protected.items() if sha256(root / p) != h]
    if changed:
        raise ValueError(f"Prior protected artifacts changed: {changed}")
    added = list((root / FIRST).glob("*")) + list((root / RESEARCH).glob("*"))
    added += list((root / "src/escape_lab").glob("fgfr2_state_bounded_*.py"))
    added += list((root / "tests").glob("test_fgfr2_state_bounded_*.py"))
    for path in added:
        if path.is_file():
            protected[str(path.relative_to(root))] = sha256(path)
    protocol = {"created_at": datetime.now(UTC).isoformat(), "version": "fgfr2-state-bounded-v1",
                "post_result_sensitivity": True, "first_run": str(FIRST), "arms": ARMS,
                "primary": PRIMARY, "model_config": MODEL_CONFIG,
                "correction": "Clip standardized biology columns to training-only min/max at prediction; "
                              "training matrices, ridge coefficients and sequence PCA unchanged",
                "fold_groups": old["fold_groups"], "eligible_ids": [c["candidate_id"] for c in cohort],
                "reserved_groups": old["reserved_groups"], "protected_sha256": protected,
                **{key: first_protocol[key] for key in ("target", "incremental_gate",
                    "state_contrast_check", "tail_secondary", "bootstrap")}}
    output.mkdir(parents=True)
    write_json(output / "protocol.json", protocol)
    protocol_hash = sha256(output / "protocol.json")
    write_json(output / "status.json", {"status": "started", "protocol_sha256": protocol_hash})
    try:
        activation, coverage = load_development_activation(root / WORKBOOK, cohort)
        if set(activation) != {c["candidate_id"] for c in cohort}:
            raise ValueError("Incomplete activation cohort")
        first = json.loads((root / FIRST / "private-predictions.json").read_text())
        historical = json.loads((root / PREVIOUS / "private-predictions.json").read_text())["cached_v2_biology"]
        store = BiologicalFeatureStore(root / "data/revision-v2/features")
        augmenter = ActivationFeatureAugmenter(store, cohort, root / REFERENCES, root / SCORES)
        sequence = load_features(cohort, root / REFERENCES, root / ESM_CACHE)
        state = FGFR2StateFeatureStore(root / FEATURES)
        tasks = fold_partitions(cohort, old["fold_groups"])
        rows, reproduction, diagnostics = {arm: [] for arm in ARMS}, {}, []
        fitting_started = time.monotonic()
        for arm in ARMS:
            unbounded = []
            for fold, train, test in tasks:
                records = {c["candidate_id"]: activation[c["candidate_id"]] for c in train}
                model = BoundedFGFR2StateModel(arm, store, augmenter, sequence, state).fit(train, records)
                scores = model.predict(test)
                diagnostics.append({"fold": fold, **model.diagnostics})
                for c, score, raw in zip(test, scores, model.unbounded_predictions, strict=True):
                    r = activation[c["candidate_id"]]
                    common = {"candidate_id": c["candidate_id"], "group_id": c["group_id"],
                              "fold": fold, "truth": r["value"], "high": activation_label(r),
                              "strong": activation_label(r, 10.)}
                    rows[arm].append({**common, "score": float(score)})
                    unbounded.append({**common, "score": float(raw)})
            reproduction[arm] = verify_cached_predictions(unbounded, first[arm])
        fitting_seconds = time.monotonic() - fitting_started
        metrics = {arm: continuous_summary(data) for arm, data in rows.items()}
        metrics[HISTORICAL] = continuous_summary(historical)
        uncertainty = paired_uncertainty(rows)
        changed = [p for p, h in protected.items() if sha256(root / p) != h]
        if changed or sha256(output / "protocol.json") != protocol_hash:
            raise ValueError(f"Protected artifacts changed: {changed}")
        result = {"post_result_sensitivity": True, "metrics": metrics, "gates": compare_gates(metrics),
                  "uncertainty": uncertainty, "coverage": coverage, "model_diagnostics": diagnostics,
                  "unbounded_reproduction": reproduction, "fitting_seconds": fitting_seconds,
                  "elapsed_seconds": time.monotonic() - started, "protocol_sha256": protocol_hash,
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
    print(json.dumps({"gates": result["gates"], "unbounded_reproduction": result["unbounded_reproduction"],
                      "fitting_seconds": result["fitting_seconds"]}, indent=2))
