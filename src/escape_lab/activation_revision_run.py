"""Freeze and run the second developmental activation pilot without external evaluation."""

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
from escape_lab.activation_revision import (
    ARMS,
    CAPACITY,
    GENES,
    PRIMARY,
    component_gate,
    evaluate_revision,
    resistance_gate,
)
from escape_lab.biological_features import BiologicalFeatureStore
from escape_lab.data import sha256
from escape_lab.improvement import load_labels, write_json
from escape_lab.mechanistic_pilot import SALT, SEED, eligible_candidates
from escape_lab.mechanistic_pilot_run import MANO_AUDIT, RLY_AUDIT, WORKBOOK
from escape_lab.revision_evaluation import grouped_folds

V1 = Path("reports/mechanistic-pilot/activation-transfer-v1")
SCORES = Path("research/external-predictions/2026-10-17-aggregate-sources/"
              "existing_model_scores_by_protein.csv")


def replicate_diagnostics(candidates, activation):
    result = {}
    for gene in ("pooled", *GENES):
        rows = [activation[c["candidate_id"]] for c in candidates if gene == "pooled" or c["gene"] == gene]
        matrix = np.log2(1 + np.asarray([r["replicates"] for r in rows]))
        first, second = np.median(matrix[:, :2], axis=1), np.median(matrix[:, 2:], axis=1)
        active = np.asarray([r["active"] for r in rows])
        sd = matrix.std(axis=1)
        pairwise = [float(spearmanr(matrix[:, i], matrix[:, j]).statistic)
                    for i in range(4) for j in range(i + 1, 4)]
        result[gene] = {"n": len(rows), "split_half_spearman": float(spearmanr(first, second).statistic),
                        "median_pairwise_replicate_spearman": float(np.median(pairwise)),
                        "median_sd_active": float(np.median(sd[active])) if active.any() else None,
                        "median_sd_inactive": float(np.median(sd[~active])) if (~active).any() else None}
    return {"metrics": result, "interpretation": "Reproducibility of summary replicates, not a "
            "rigorous model-performance ceiling; source count uncertainty is unavailable."}


def paired_uncertainty(raw, draws=1000):
    """Fixed-prediction, paired aligned-block resampling; no model selection/refitting."""
    rows = raw["activation_records"]
    ref = sorted(rows["original"], key=lambda r: r["candidate_id"])
    primary = sorted(rows[PRIMARY], key=lambda r: r["candidate_id"])
    if [r["candidate_id"] for r in ref] != [r["candidate_id"] for r in primary]:
        raise ValueError("Unpaired activation predictions")
    groups = sorted({r["group_id"] for r in ref})
    group_index = {g: i for i, g in enumerate(groups)}
    by_gene = {}
    for gene in GENES:
        selected = [i for i, r in enumerate(ref) if r["gene"] == gene]
        by_gene[gene] = {
            "truth": np.asarray([ref[i]["truth"] for i in selected]),
            "active": np.asarray([ref[i]["active"] for i in selected], dtype=int),
            "original": np.asarray([ref[i]["prediction"] for i in selected]),
            "primary": np.asarray([primary[i]["prediction"] for i in selected]),
            "original_score": np.asarray([ref[i]["active_score"] for i in selected]),
            "primary_score": np.asarray([primary[i]["active_score"] for i in selected]),
            "groups": np.asarray([group_index[ref[i]["group_id"]] for i in selected]),
        }
    resistance = {}
    for drug in ("pem", "fut"):
        ordered = {}
        for arm in ("structure", PRIMARY):
            ordered[arm] = sorted([r for r in raw["resistance_records"][arm]
                                   if r["drug"] == drug and type(r["label"]) is bool],
                                  key=lambda r: r["candidate_id"])
        resistance[drug] = {
            "truth": np.asarray([r["label"] for r in ordered["structure"]], dtype=int),
            "groups": np.asarray([group_index[r["group_id"]] for r in ordered["structure"]]),
            "structure": np.asarray([r["probability"] for r in ordered["structure"]]),
            "primary": np.asarray([r["probability"] for r in ordered[PRIMARY]]),
        }
    differences = {"macro_gene_spearman": [], "macro_gene_active_ap": [], "resistance_macro_ap": []}
    rng = np.random.default_rng(SEED)
    for _ in range(draws):
        counts = np.bincount(rng.integers(len(groups), size=len(groups)), minlength=len(groups))
        rho, act_ap = [], []
        for array in by_gene.values():
            indices = np.repeat(np.arange(len(array["truth"])), counts[array["groups"]])
            y, p, b = (array[key][indices] for key in ("truth", "primary", "original"))
            if len(indices) < 2 or not np.ptp(y) or not np.ptp(p) or not np.ptp(b):
                rho.append(None)
            else:
                rho.append(float(spearmanr(y, p).statistic - spearmanr(y, b).statistic))
            active = array["active"][indices]
            if len(np.unique(active)) < 2:
                act_ap.append(None)
            else:
                act_ap.append(float(average_precision_score(active, array["primary_score"][indices])
                                    - average_precision_score(active, array["original_score"][indices])))
        drug_ap = []
        for array in resistance.values():
            weights = counts[array["groups"]]
            if np.dot(array["truth"], weights) == 0:
                drug_ap.append(0.)
            else:
                drug_ap.append(float(average_precision_score(array["truth"], array["primary"],
                                                              sample_weight=weights)
                                     - average_precision_score(array["truth"], array["structure"],
                                                               sample_weight=weights)))
        for name, values in (("macro_gene_spearman", rho), ("macro_gene_active_ap", act_ap),
                             ("resistance_macro_ap", drug_ap)):
            if all(v is not None and math.isfinite(v) for v in values):
                differences[name].append(float(np.mean(values)))
    return {"draws": draws, "seed": SEED, "groups": len(groups),
            "contrasts": {name: {"percentile_2_5": float(np.percentile(v, 2.5)) if v else None,
                                 "percentile_97_5": float(np.percentile(v, 97.5)) if v else None,
                                 "valid_draws": len(v), "invalid_draws": draws - len(v)}
                          for name, v in differences.items()},
            "interpretation": "Conditional descriptive percentiles on 11 related development "
                              "blocks, not independent biological confirmation. No refitting. "
                              "Activation contrasts versus original; resistance versus structure."}


def verify_v1_reproduction(raw, previous):
    old_activation = {r["candidate_id"]: r for r in previous["activation_records"]
                      if r["arm"] == "predicted_activation"}
    new_activation = {r["candidate_id"]: r for r in raw["activation_records"]["original"]}
    if set(old_activation) != set(new_activation):
        raise ValueError("V1 reproduction candidate mismatch")
    delta = max(abs(r["prediction"] - old_activation[cid]["prediction"])
                for cid, r in new_activation.items())
    if delta > 1e-12:
        raise ValueError("Original activation predictions do not reproduce")
    for arm, old_arm in (("structure", "structure"), ("original", "predicted_activation")):
        old = {(r["candidate_id"], r["drug"]): r for r in previous["records"][old_arm]}
        new = {(r["candidate_id"], r["drug"]): r for r in raw["resistance_records"][arm]}
        if set(new) != set(old) or any(abs(r["probability"] - old[k]["probability"]) > 1e-12
                                     for k, r in new.items()):
            raise ValueError("Original resistance predictions do not reproduce")
    return {"matched_variants": len(new_activation), "max_activation_prediction_difference": delta,
            "resistance_predictions_reproduced": True}


def make_report(result):
    def number(x):
        return "undefined" if x is None else f"{x:.5f}"

    lines = ["# Second activation pilot", "", (f"Primary challenger: `{PRIMARY}`. "
             f"Component: **{result['component_gate']['status']}**. "
             f"Resistance: **{result['resistance_gate']['status']}**."), "",
             (f"Same {result['n']} variants and 11 aligned blocks as pilot v1. "
             f"Fitting: {result['fitting_seconds']:.2f}s; total: {result['elapsed_seconds']:.2f}s."), "",
             ("| Arm | Macro within-gene rho | MAE | Activation-class macro AP | "
             "Activation-class macro recall at 10% | Resistance macro AP | Macro Brier |"),
             "|---|---:|---:|---:|---:|---:|---:|"]
    for arm in ARMS:
        a, r = result["activation_metrics"][arm], result["resistance_metrics"][arm]
        values = [a["macro_gene_spearman"], a["pooled"]["mae"], a["macro_gene_active_ap"],
                  a["macro_gene_active_recall_at_10pct"], r["macro_ap"], r["macro_brier"]]
        lines.append(f"| {arm} | " + " | ".join(number(v) for v in values) + " |")
    base = result["resistance_metrics"]["structure"]
    lines.extend(["", (f"Structural-only resistance baseline: AP {base['macro_ap']:.5f}; "
                  f"Brier {base['macro_brier']:.5f}."), "",
                  ("Higher correlation, AP and recall are better; lower MAE and Brier are better. "
                  "AP is not percentage accuracy. The two enriched arms share the same continuous "
                  "regressor. Their difference is the additional activation classifier."), "",
                  "## Per-gene activation results", "",
                  "| Gene | Arm | Rho | Activation-class AP | Recall at 10% | Active / N |",
                  "|---|---|---:|---:|---:|---:|"])
    for gene in GENES:
        for arm in ARMS:
            m = result["activation_metrics"][arm]["per_gene"][gene]
            lines.append(f"| {gene} | {arm} | {number(m['spearman'])} | "
                         f"{number(m['active_ap'])} | {number(m['active_recall_at_10pct'])} | "
                         f"{m['active_count']} / {m['n']} |")
    lines.extend(["", "## Decisions and uncertainty", "", "```json", json.dumps({
        "component_gate": result["component_gate"], "resistance_gate": result["resistance_gate"],
        "uncertainty": result["uncertainty"]}, indent=2), "```", "",
        "## Interpretation limits", "",
        ("- Original development outcomes and folds have been reused. This is a hypothesis-directed "
        "development comparison, not a new final test."),
        ("- The binary activation target is our protein-level adaptation of the study's replicate "
        "rule; it is not asserted to reproduce its nucleotide-level annotations."),
        ("- Added motif features are sequence distances, not mutant binding/free-energy predictions. "
        "AlphaMissense/EVE are auxiliary predictions, not activation labels. No curated outcome "
        "annotations or pretrained FGFR-specific supervised models were imported."),
        ("- Missing EVE scores are explicitly flagged and imputed only in training folds. "
        "Both activation heads are generated out of fold for downstream training."),
        "- No LLM calls, hyperparameter search, external-test reveal or post-score model tuning occurred.",
        (f"- {result['integrity']['protected_files']} protected files stayed hash-identical; "
        "original activation and resistance predictions reproduce to 1e-12."),
        "", "See results.json for all per-drug, fold, prior-reference, replicate and fitting diagnostics.", ""])
    return "\n".join(lines)


def run(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError("Use a new run directory; never replace a prior attempt")
    start = time.monotonic()
    old_protocol = json.loads((root / V1 / "protocol.json").read_text())
    # Another chat is independently adding You.com research modules. Freeze the
    # established original inventory and this experiment's complete dependency
    # set, rather than unrelated files that its owner is still changing.
    protected = dict(old_protocol["protected_sha256"])
    changed = [p for p, h in old_protocol["protected_sha256"].items() if sha256(root / p) != h]
    if changed:
        raise ValueError(f"Original pilot artifacts changed: {changed}")
    added = [SCORES, Path("research/activation-v2/research-and-protocol.md"),
             Path("research/activation-v2/diagnostics.json"), Path("research/activation-v2/diagnostics.md")]
    added.extend(p.relative_to(root) for p in (root / V1).iterdir() if p.is_file())
    added.extend(Path("src/escape_lab") / name for name in (
        "activation_data.py", "activation_features.py", "activation_revision.py", "activation_revision_run.py"))
    added.extend(Path("tests") / name for name in (
        "test_activation_data.py", "test_activation_features.py", "test_activation_revision.py"))
    for p in added:
        protected[str(p)] = sha256(root / p)
    candidates = json.loads((root / "data/processed/candidates.json").read_text())
    mano, rly = (json.loads((root / path).read_text()) for path in (MANO_AUDIT, RLY_AUDIT))
    reserved = set(mano["mapped_groups"]) | {r["group_id"] for r in rly["panel_variants"]["rly2023"]}
    eligible = eligible_candidates(candidates, reserved)
    if [c["candidate_id"] for c in eligible] != old_protocol["eligible_candidate_ids"]:
        raise ValueError("Eligibility differs from pilot v1")
    protocol = {"version": "activation-revision-v2", "created_at": datetime.now(UTC).isoformat(),
                "arms": ARMS, "primary": PRIMARY, "capacity": CAPACITY, "seed": SEED, "salt": SALT,
                "eligible_ids": [c["candidate_id"] for c in eligible], "reserved_groups": sorted(reserved),
                "fold_groups": {t["task_id"]: sorted({c["group_id"] for c in t["pool"]})
                                for t in grouped_folds(eligible, SALT, prefix="pilot")},
                "required_complete_n": 2719, "required_groups": 11,
                "active_definition": "3/4 per-replicate codon medians >=1.5; protein-level adaptation",
                "continuous_target": "exact pilot-v1 log2(1 + median of four replicate summaries)",
                "component_gate": {"macro_gene_rho_gain": .05, "positive_fold_count": 2,
                                   "max_gene_rho_loss": .05, "max_mae_ratio": 1.05,
                                   "macro_gene_active_ap_gain": .02, "max_gene_active_ap_loss": .02,
                                   "undefined_class_or_fallback": "inconclusive"},
                "resistance_gate": {"macro_ap_gain": .005, "positive_fold_count": 2,
                                    "max_drug_ap_loss": .01, "brier_max": "baseline*1.05+.001"},
                "bootstrap": {"draws": 1000, "seed": SEED, "unit": "aligned block"},
                "protected_sha256": protected,
                "scope": "Development only, no new independent confirmation or clinical claims"}
    output.mkdir(parents=True)
    write_json(output / "protocol.json", protocol)
    protocol_hash = sha256(output / "protocol.json")
    write_json(output / "status.json", {"status": "started", "protocol_sha256": protocol_hash})
    try:
        activation, coverage = load_development_activation(root / WORKBOOK, eligible)
        cohort = [c for c in eligible if c["candidate_id"] in activation]
        previous = json.loads((root / V1 / "private-predictions.json").read_text())
        old_truth = {r["candidate_id"]: r["truth"] for r in previous["activation_records"]
                     if r["arm"] == "predicted_activation"}
        if len(cohort) != 2719 or len({c["group_id"] for c in cohort}) != 11:
            raise ValueError("Complete cohort differs from pilot v1")
        if set(activation) != set(old_truth) or any(v["value"] != old_truth[cid]
                                                   for cid, v in activation.items()):
            raise ValueError("Activation target changed from pilot v1")
        resistance = load_labels(root / "data/processed", cohort)
        store = BiologicalFeatureStore(root / "data/revision-v2/features")
        augmenter = ActivationFeatureAugmenter(store, cohort,
                                               root / "data/revision-v2/features/references.json",
                                               root / SCORES)
        fitting_start = time.monotonic()
        raw = evaluate_revision(cohort, activation, resistance, store, augmenter)
        fitting_seconds = time.monotonic() - fitting_start
        reproduced = verify_v1_reproduction(raw, previous)
        uncertainty = paired_uncertainty(raw)
        changed = [p for p, h in protected.items() if sha256(root / p) != h]
        if changed or sha256(output / "protocol.json") != protocol_hash:
            raise ValueError(f"Protected artifacts changed: {changed}")
        result = {"n": len(cohort), "coverage": coverage, "feature_score_coverage": augmenter.score_coverage,
                  "activation_metrics": raw["activation_metrics"],
                  "resistance_metrics": raw["resistance_metrics"],
                  "gene_prior_metrics": raw["gene_prior_metrics"],
                  "replicate_diagnostics": replicate_diagnostics(cohort, activation),
                  "component_gate": component_gate(raw), "resistance_gate": resistance_gate(raw),
                  "uncertainty": uncertainty, "v1_reproduction": reproduced,
                  "model_diagnostics": raw["resistance_diagnostics"],
                  "tail_fit_fallbacks": [a for a in raw["split_audit"] if a["outer_tail_fallback"]
                                        or any(i["tail_fallback"] for i in a["inner"])],
                  "fitting_seconds": fitting_seconds, "elapsed_seconds": time.monotonic() - start,
                  "protocol_sha256": protocol_hash,
                  "integrity": {"protected_files": len(protected), "changed": changed}}
        write_json(output / "private-predictions.json", raw)
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
    print(json.dumps({k: result[k] for k in ("component_gate", "resistance_gate", "n",
                                           "fitting_seconds", "elapsed_seconds", "v1_reproduction")},
                     indent=2))
