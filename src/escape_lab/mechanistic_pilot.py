"""Fixed developmental activation-transfer pilot. Never authorizes final evaluation."""

import hashlib
import math
import re
import statistics
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.feature_extraction import DictVectorizer
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error

from escape_lab.benchmark import validate_partition
from escape_lab.revision_evaluation import grouped_folds, summarize
from escape_lab.revision_models import MISSING_MEASUREMENTS, BioConfig, BiologyPredictor

ARMS = ("structure", "predicted_activation", "measured_activation", "shuffled_activation")
SEED = 17
SALT = "activation-transfer-pilot-v1"
ACTIVATION_PARAMS = {"max_iter": 100, "max_depth": 3, "max_leaf_nodes": 7,
                     "min_samples_leaf": 30, "l2_regularization": 10.0,
                     "early_stopping": False, "random_state": SEED}


def eligible_candidates(candidates, reserved_groups):
    kept = [c for c in candidates if c["split"] in ("development", "validation")
            and c["group_id"] not in reserved_groups]
    if len({c["candidate_id"] for c in kept}) != len(kept):
        raise ValueError("Duplicate candidate identity")
    return sorted(kept, key=lambda c: c["candidate_id"])


def activation_from_rows(rows, candidates):
    """Read only eligible activation cells, never resistance or discovery columns.

    Aggregate codon representations within each replicate, then take the median
    across four complete replicate summaries. Missing measurements are not zeros.
    """
    if any(c["split"] not in ("development", "validation") for c in candidates):
        raise ValueError("Activation pilot only permits development candidates")
    lookup = {(c["gene"], c["wt"], c["position"], c["mut"]): c["candidate_id"]
              for c in candidates}
    replicates = defaultdict(lambda: [[], [], [], []])
    seen = set()
    count = 0
    for row in rows:
        gene, aa = row[0], row[4]
        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", str(aa))
        if match is None:
            continue
        wt, position, mut = match.groups()
        cid = lookup.get((gene, wt, int(position), mut))
        if cid is None:
            continue
        identity = (gene, row[3])
        if identity in seen:
            raise ValueError("Duplicate nucleotide identity in activation source")
        seen.add(identity)
        count += 1
        for i in range(4):
            value = row[5 + i]
            if value is None or str(value).strip().lower() in ("", "na", "nan", "n/a"):
                continue
            value = float(value)
            if not math.isfinite(value) or value < 0:
                raise ValueError("Activation requires finite nonnegative enrichment")
            replicates[cid][i].append(value)
    values, spreads = {}, []
    for c in candidates:
        cid = c["candidate_id"]
        reps = replicates.get(cid, [[], [], [], []])
        if all(reps):
            medians = [statistics.median(rep) for rep in reps]
            values[cid] = math.log2(1 + statistics.median(medians))
            spreads.append(float(np.std(np.log2(1 + np.asarray(medians)))))
    return values, {
        "eligible_nucleotide_rows": count,
        "eligible_protein_variants": len(candidates),
        "complete_protein_variants": len(values),
        "incomplete_protein_variants": len(candidates) - len(values),
        "median_within_variant_replicate_sd_log2_1p": statistics.median(spreads)
        if spreads else None,
        "aggregation": "median across codons within replicate, median of four replicates",
        "transform": "log2(1 + enrichment)",
    }


class ActivationRegressor:
    """Nonlinear auxiliary predictor with fixed capacity and train-only preprocessing."""

    def __init__(self, store):
        self.store = store
        self.vectorizer = DictVectorizer(sparse=False)
        self.imputer = SimpleImputer(strategy="median", keep_empty_features=True)
        self.model = HistGradientBoostingRegressor(**ACTIVATION_PARAMS)

    def features(self, candidates):
        result = []
        for c in candidates:
            raw = self.store.features(c, "pem", "structure")
            row = {k: np.nan if v == 999 or k in MISSING_MEASUREMENTS and v == -1 else v
                   for k, v in raw.items() if "ligand" not in k and "covalent_cys" not in k}
            row.update({f"gene_{gene}": float(c["gene"] == gene)
                        for gene in ("FGFR1", "FGFR2", "FGFR3", "FGFR4")})
            result.append(row)
        return result

    def fit(self, candidates, values):
        allowed = {c["candidate_id"] for c in candidates}
        if set(values) != allowed:
            raise ValueError("Activation labels must exactly match the training partition")
        target = np.asarray([values[c["candidate_id"]] for c in candidates])
        if not len(target) or not np.isfinite(target).all():
            raise ValueError("Finite activation training values required")
        matrix = self.vectorizer.fit_transform(self.features(candidates))
        matrix = self.imputer.fit_transform(matrix)
        self.model.fit(matrix, target)
        return self

    def predict(self, candidates):
        matrix = self.vectorizer.transform(self.features(candidates))
        result = self.model.predict(self.imputer.transform(matrix))
        if not np.isfinite(result).all():
            raise ValueError("Nonfinite activation prediction")
        return result


def training_values(candidates, values, shuffled, salt):
    result = {c["candidate_id"]: values[c["candidate_id"]] for c in candidates}
    if not shuffled:
        return result
    # Shuffle separately inside each fit, so inner-held-out outcomes cannot leak
    # into training through an earlier outer-partition permutation.
    seed = int(hashlib.sha256(salt.encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    for gene in sorted({c["gene"] for c in candidates}):
        ids = sorted(c["candidate_id"] for c in candidates if c["gene"] == gene)
        shuffled_values = rng.permutation([result[cid] for cid in ids])
        result.update(zip(ids, map(float, shuffled_values), strict=True))
    return result


def crossfit_activation(candidates, values, factory, salt, shuffled):
    predictions, audit = {}, []
    for task in grouped_folds(candidates, salt, prefix="activation-inner"):
        training, pool = task["training"], task["pool"]
        validate_partition(training, pool)
        subset = training_values(training, values, shuffled, salt + task["task_id"])
        model = factory().fit(training, subset)
        scores = model.predict(pool)
        predictions.update((c["candidate_id"], float(p))
                           for c, p in zip(pool, scores, strict=True))
        audit.append({"training_groups": sorted({c["group_id"] for c in training}),
                      "prediction_groups": sorted({c["group_id"] for c in pool})})
    if set(predictions) != {c["candidate_id"] for c in candidates}:
        raise ValueError("Crossfit did not cover training candidates exactly")
    return predictions, audit


class ActivationStore:
    def __init__(self, store, activation):
        self.store, self.activation = store, activation

    def features(self, candidate, drug, family):
        row = dict(self.store.features(candidate, drug, family))
        activation = self.activation[candidate["candidate_id"]]
        if not math.isfinite(activation):
            raise ValueError("Finite activation feature required")
        distance = row["structure_ligand_ca_distance"]
        proximity = math.exp(-distance / 5) if row["structure_ca_available"] else 0.0
        row["activation_log2_1p"] = activation
        row["activation_x_ligand_proximity"] = activation * proximity
        return row


def evaluate_pilot(candidates, labels, activation, store):
    ids = {c["candidate_id"] for c in candidates}
    if set(labels) != ids or set(activation) != ids:
        raise ValueError("Pilot labels must exactly match the common cohort")
    if any(c["split"] not in ("development", "validation") for c in candidates):
        raise ValueError("Final candidates forbidden in pilot")
    records = {arm: [] for arm in ARMS}
    diagnostics, activation_records, audit = [], [], []
    tasks = grouped_folds(candidates, SALT, prefix="pilot")
    for task in tasks:
        training, pool, fold = task["training"], task["pool"], task["task_id"]
        validate_partition(training, pool)
        allowed_activation = {c["candidate_id"]: activation[c["candidate_id"]] for c in training}
        arm_stores = {"structure": store, "measured_activation": ActivationStore(store, activation)}
        for arm, shuffled in (("predicted_activation", False), ("shuffled_activation", True)):
            def factory():
                return ActivationRegressor(store)

            inner, inner_audit = crossfit_activation(
                training, allowed_activation, factory, SALT + fold, shuffled)
            fit_values = training_values(training, allowed_activation, shuffled, SALT + fold)
            model = factory().fit(training, fit_values)
            predicted = model.predict(pool)
            feature_values = {**inner, **{c["candidate_id"]: float(p)
                                         for c, p in zip(pool, predicted, strict=True)}}
            arm_stores[arm] = ActivationStore(store, feature_values)
            audit.append({"fold": fold, "arm": arm, "inner_splits": inner_audit,
                          "outer_training_groups": sorted({c["group_id"] for c in training}),
                          "outer_prediction_groups": sorted({c["group_id"] for c in pool})})
            for c, prediction in zip(pool, predicted, strict=True):
                activation_records.append({"candidate_id": c["candidate_id"], "gene": c["gene"],
                                           "group_id": c["group_id"], "fold": fold, "arm": arm,
                                           "truth": activation[c["candidate_id"]],
                                           "prediction": float(prediction)})
        for arm in ARMS:
            for drug in ("pem", "fut"):
                known = [c for c in training if type(labels[c["candidate_id"]][drug]) is bool]
                truth = [labels[c["candidate_id"]][drug] for c in known]
                model = BiologyPredictor(BioConfig(representation="structure", C=.1),
                                         arm_stores[arm], drug, seed=SEED).fit(known, truth)
                probabilities = model.predict(pool)
                diagnostics.append({"fold": fold, "arm": arm, "drug": drug,
                                    "fallback": model.fallback_reason, **model.diagnostics})
                for c, p in zip(pool, probabilities, strict=True):
                    records[arm].append({"candidate_id": c["candidate_id"], "gene": c["gene"],
                                         "group_id": c["group_id"], "drug": drug, "fold": fold,
                                         "label": labels[c["candidate_id"]][drug],
                                         "probability": float(p)})
    return {"records": records, "metrics": {arm: summarize(r) for arm, r in records.items()},
            "activation_records": activation_records, "diagnostics": diagnostics,
            "split_audit": audit}


def activation_metrics(records):
    def metrics(rows):
        y = np.asarray([r["truth"] for r in rows])
        p = np.asarray([r["prediction"] for r in rows])
        rho = float(spearmanr(y, p).statistic) if np.std(y) and np.std(p) else None
        return {"n": len(rows), "spearman": rho, "mae_log2_1p": float(mean_absolute_error(y, p))}

    result = {}
    for arm in ("predicted_activation", "shuffled_activation"):
        rows = [r for r in records if r["arm"] == arm]
        result[arm] = {"pooled": metrics(rows), "per_gene": {
            gene: metrics([r for r in rows if r["gene"] == gene])
            for gene in sorted({r["gene"] for r in rows})}}
    return result


def advancement_gate(result, activation_result):
    base = result["metrics"]["structure"]
    candidate = result["metrics"]["predicted_activation"]
    shuffled = result["metrics"]["shuffled_activation"]
    required = []
    for metric in (base, candidate, shuffled):
        required.append(metric.get("macro_ap"))
        required.extend(metric["per_drug"][drug].get(name)
                        for drug in ("pem", "fut") for name in ("average_precision", "brier"))
        required.extend(f.get("macro_ap") for f in metric.get("by_fold", {}).values())
    if any(type(x) not in (int, float) or not math.isfinite(x) for x in required):
        return {"status": "inconclusive", "failed_guards": [],
                "inconclusive_reasons": ["undefined_metric"], "fold_macro_ap_deltas": {},
                "interpretation": "Insufficient conclusive labels for the predeclared gate."}
    reasons, inconclusive = [], []
    if candidate["macro_ap"] - base["macro_ap"] < .005:
        reasons.append("macro_ap_gain_below_0.005")
    if candidate["macro_ap"] <= shuffled["macro_ap"]:
        reasons.append("does_not_beat_shuffled_activation")
    fold_deltas = {fold: candidate["by_fold"][fold]["macro_ap"] - m["macro_ap"]
                   for fold, m in base["by_fold"].items()}
    if sum(d > 0 for d in fold_deltas.values()) < 2:
        reasons.append("gain_in_fewer_than_two_folds")
    for drug in ("pem", "fut"):
        a, b = candidate["per_drug"][drug], base["per_drug"][drug]
        if a["average_precision"] < b["average_precision"] - .01:
            reasons.append(f"{drug}_ap_regression")
        if a["brier"] > b["brier"] * 1.05 + .001:
            reasons.append(f"{drug}_brier_regression")
        for fold, m in base["by_fold"].items():
            if m["per_drug"][drug]["positive_count"] == 0:
                inconclusive.append(f"{fold}_{drug}_no_positive_examples")
    if any(d["fallback"] or d["convergence_warning"] for d in result["diagnostics"]):
        inconclusive.append("fit_fallback_or_convergence_warning")
    rho = activation_result["predicted_activation"]["pooled"]["spearman"]
    if rho is None or rho <= 0:
        reasons.append("activation_rank_correlation_not_positive")
    return {"status": "inconclusive" if inconclusive else "advance" if not reasons else "do_not_advance",
            "failed_guards": reasons, "inconclusive_reasons": inconclusive,
            "fold_macro_ap_deltas": fold_deltas,
            "interpretation": "Development triage only; thresholds are not significance tests."}
