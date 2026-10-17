"""Small fixed activation-feature/capacity/tail ablation, with nested group isolation."""

import math

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.feature_extraction import DictVectorizer
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, mean_absolute_error, mean_squared_error

from escape_lab.benchmark import validate_partition
from escape_lab.mechanistic_pilot import SALT, SEED, ActivationRegressor, ActivationStore
from escape_lab.revision_evaluation import grouped_folds, summarize
from escape_lab.revision_models import BioConfig, BiologyPredictor

ARMS = ("original", "capacity", "biology", "biology_tail")
PRIMARY = "biology_tail"
CAPACITY = {"max_iter": 200, "max_depth": 4, "max_leaf_nodes": 15,
            "min_samples_leaf": 10, "l2_regularization": 1.0,
            "early_stopping": False, "random_state": SEED}
GENES = ("FGFR1", "FGFR2", "FGFR3", "FGFR4")


class ActivationRevisionModel:
    def __init__(self, store, augmenter, arm):
        if arm not in ARMS:
            raise ValueError("Unknown activation revision arm")
        self.store, self.augmenter, self.arm = store, augmenter, arm
        self.original = ActivationRegressor(store)
        self.vectorizer = DictVectorizer(sparse=False)
        self.imputer = SimpleImputer(strategy="median", keep_empty_features=True)
        self.regressor = HistGradientBoostingRegressor(**CAPACITY)
        self.classifier = None
        self.tail_fallback = None
        self.tail_prior = None

    def features(self, candidates):
        rows = self.original.features(candidates)
        if self.arm in ("biology", "biology_tail"):
            for c, row in zip(candidates, rows, strict=True):
                extra = self.augmenter.features(c)
                if set(extra) & set(row):
                    raise ValueError("Augmentation must not overwrite original features")
                row.update(extra)
        return rows

    def fit(self, candidates, records):
        if set(records) != {c["candidate_id"] for c in candidates}:
            raise ValueError("Activation labels must exactly match training candidates")
        values = {cid: r["value"] for cid, r in records.items()}
        if self.arm == "original":
            self.original.fit(candidates, values)
            return self
        target = np.asarray([values[c["candidate_id"]] for c in candidates])
        if not len(target) or not np.isfinite(target).all():
            raise ValueError("Finite nonempty targets required")
        matrix = self.vectorizer.fit_transform(self.features(candidates))
        matrix = self.imputer.fit_transform(matrix)
        self.regressor.fit(matrix, target)
        if self.arm == "biology_tail":
            active = [records[c["candidate_id"]]["active"] for c in candidates]
            if any(type(a) is not bool for a in active):
                raise ValueError("Boolean activation class required")
            self.tail_prior = float(np.mean(active))
            if len(set(active)) < 2:
                self.tail_fallback = "single_class_training"
            else:
                self.classifier = HistGradientBoostingClassifier(**CAPACITY)
                self.classifier.fit(matrix, np.asarray(active, dtype=int))
        return self

    def predict(self, candidates):
        if self.arm == "original":
            values = self.original.predict(candidates)
            return {"value": values, "active_score": values.copy()}
        matrix = self.vectorizer.transform(self.features(candidates))
        matrix = self.imputer.transform(matrix)
        values = self.regressor.predict(matrix)
        scores = values.copy()
        if self.arm == "biology_tail":
            if self.tail_fallback:
                scores = np.full(len(candidates), self.tail_prior)
            else:
                index = int(np.flatnonzero(self.classifier.classes_ == 1)[0])
                scores = self.classifier.predict_proba(matrix)[:, index]
        if not np.isfinite(values).all() or not np.isfinite(scores).all():
            raise ValueError("Nonfinite prediction")
        return {"value": values, "active_score": scores}


def crossfit_outputs(candidates, records, factory, salt):
    predictions, audit = {}, []
    for task in grouped_folds(candidates, salt, prefix="activation-inner"):
        train, pool = task["training"], task["pool"]
        validate_partition(train, pool)
        labels = {c["candidate_id"]: records[c["candidate_id"]] for c in train}
        model = factory().fit(train, labels)
        outputs = model.predict(pool)
        for i, c in enumerate(pool):
            predictions[c["candidate_id"]] = {k: float(v[i]) for k, v in outputs.items()}
        audit.append({"training_groups": sorted({c["group_id"] for c in train}),
                      "prediction_groups": sorted({c["group_id"] for c in pool}),
                      "tail_fallback": model.tail_fallback})
    if set(predictions) != {c["candidate_id"] for c in candidates}:
        raise ValueError("Incomplete crossfit coverage")
    return predictions, audit


class DualActivationStore:
    def __init__(self, store, predictions, use_tail):
        self.base = ActivationStore(store, {cid: r["value"] for cid, r in predictions.items()})
        self.predictions, self.use_tail = predictions, use_tail

    def features(self, candidate, drug, family):
        row = self.base.features(candidate, drug, family)
        if self.use_tail:
            score = self.predictions[candidate["candidate_id"]]["active_score"]
            if not math.isfinite(score) or not 0 <= score <= 1:
                raise ValueError("Invalid tail-head score")
            proximity = (math.exp(-row["structure_ligand_ca_distance"] / 5)
                         if row["structure_ca_available"] else 0.)
            row["activation_active_score"] = score
            row["activation_active_x_ligand_proximity"] = score * proximity
        return row


def _metrics(rows):
    truth = np.asarray([r["truth"] for r in rows])
    prediction = np.asarray([r["prediction"] for r in rows])
    active = np.asarray([r["active"] for r in rows], dtype=int)
    scores = np.asarray([r["active_score"] for r in rows])
    rho = (float(spearmanr(truth, prediction).statistic)
           if len(rows) > 1 and np.ptp(truth) > 0 and np.ptp(prediction) > 0 else None)
    classes = len(np.unique(active)) == 2
    ap = float(average_precision_score(active, scores)) if classes else None
    k = max(1, math.ceil(len(rows) * .1))
    selected = np.argsort(-scores, kind="stable")[:k]
    recall = float(active[selected].sum() / active.sum()) if classes else None
    top = np.argsort(-truth, kind="stable")[:k]
    predicted_top = np.argsort(-prediction, kind="stable")[:k]
    return {"n": len(rows), "spearman": rho,
            "mae": float(mean_absolute_error(truth, prediction)),
            "rmse": float(np.sqrt(mean_squared_error(truth, prediction))),
            "bias": float(np.mean(prediction - truth)), "truth_sd": float(np.std(truth)),
            "prediction_sd": float(np.std(prediction)), "active_count": int(active.sum()),
            "active_prevalence": float(np.mean(active)), "active_ap": ap,
            "active_recall_at_10pct": recall,
            "top_decile_overlap": len(set(top) & set(predicted_top)) / k,
            "active_subset_bias": float(np.mean(prediction[active == 1] - truth[active == 1]))
            if active.sum() else None}


def _macro(values):
    return float(np.mean(values)) if all(v is not None for v in values) else None


def summarize_activation(rows, include_folds=True):
    ordered = sorted(rows, key=lambda r: r["candidate_id"])
    genes = {gene: _metrics([r for r in ordered if r["gene"] == gene])
             for gene in GENES if any(r["gene"] == gene for r in ordered)}
    result = {"pooled": _metrics(ordered), "per_gene": genes,
              "macro_gene_spearman": _macro([genes[g]["spearman"] if g in genes else None
                                             for g in GENES]),
              "macro_gene_active_ap": _macro([genes[g]["active_ap"] if g in genes else None
                                              for g in GENES]),
              "macro_gene_active_recall_at_10pct": _macro([
                  genes[g]["active_recall_at_10pct"] if g in genes else None for g in GENES])}
    if include_folds:
        result["by_fold"] = {fold: summarize_activation([r for r in rows if r["fold"] == fold], False)
                             for fold in sorted({r["fold"] for r in rows})}
    return result


def evaluate_revision(candidates, activation, resistance, store, augmenter):
    ids = {c["candidate_id"] for c in candidates}
    if set(activation) != ids or set(resistance) != ids:
        raise ValueError("Every label map must exactly match the cohort")
    if any(c["split"] not in ("development", "validation") for c in candidates):
        raise ValueError("Final candidates forbidden")
    activation_rows = {arm: [] for arm in ARMS}
    resistance_rows = {arm: [] for arm in ("structure", *ARMS)}
    audit, diagnostics, gene_prior = [], [], []
    for task in grouped_folds(candidates, SALT, prefix="pilot"):
        train, pool, fold = task["training"], task["pool"], task["task_id"]
        validate_partition(train, pool)
        records = {c["candidate_id"]: activation[c["candidate_id"]] for c in train}
        stores = {"structure": store}
        for gene in GENES:
            train_gene = [records[c["candidate_id"]] for c in train if c["gene"] == gene]
            if not train_gene:
                raise ValueError("Training fold lacks a gene")
            mean = float(np.mean([r["value"] for r in train_gene]))
            prior = float(np.mean([r["active"] for r in train_gene]))
            for c in pool:
                if c["gene"] == gene:
                    target = activation[c["candidate_id"]]
                    gene_prior.append({"candidate_id": c["candidate_id"], "gene": gene,
                                       "group_id": c["group_id"], "fold": fold,
                                       "truth": target["value"], "active": target["active"],
                                       "prediction": mean, "active_score": prior})
        for arm in ARMS:
            def factory(selected_arm=arm):
                return ActivationRevisionModel(store, augmenter, selected_arm)

            inner, inner_audit = crossfit_outputs(train, records, factory, SALT + fold)
            model = factory().fit(train, records)
            outputs = model.predict(pool)
            outer = {c["candidate_id"]: {k: float(v[i]) for k, v in outputs.items()}
                     for i, c in enumerate(pool)}
            stores[arm] = DualActivationStore(store, {**inner, **outer}, arm == "biology_tail")
            audit.append({"fold": fold, "arm": arm, "inner": inner_audit,
                          "outer_training_groups": sorted({c["group_id"] for c in train}),
                          "outer_prediction_groups": sorted({c["group_id"] for c in pool}),
                          "outer_tail_fallback": model.tail_fallback})
            for c in pool:
                cid = c["candidate_id"]
                activation_rows[arm].append({"candidate_id": cid, "gene": c["gene"],
                                            "group_id": c["group_id"], "fold": fold,
                                            "truth": activation[cid]["value"],
                                            "active": activation[cid]["active"],
                                            "prediction": outer[cid]["value"],
                                            "active_score": outer[cid]["active_score"]})
        for arm, features in stores.items():
            for drug in ("pem", "fut"):
                known = [c for c in train if type(resistance[c["candidate_id"]][drug]) is bool]
                labels = [resistance[c["candidate_id"]][drug] for c in known]
                model = BiologyPredictor(BioConfig(representation="structure", C=.1),
                                         features, drug, seed=SEED).fit(known, labels)
                for c, p in zip(pool, model.predict(pool), strict=True):
                    resistance_rows[arm].append({"candidate_id": c["candidate_id"],
                                                "gene": c["gene"], "group_id": c["group_id"],
                                                "fold": fold, "drug": drug,
                                                "label": resistance[c["candidate_id"]][drug],
                                                "probability": float(p)})
                diagnostics.append({"fold": fold, "arm": arm, "drug": drug,
                                    "fallback": model.fallback_reason, **model.diagnostics})
    return {"activation_records": activation_rows, "resistance_records": resistance_rows,
            "activation_metrics": {arm: summarize_activation(rows)
                                   for arm, rows in activation_rows.items()},
            "resistance_metrics": {arm: summarize(rows) for arm, rows in resistance_rows.items()},
            "gene_prior_metrics": summarize_activation(gene_prior), "gene_prior_records": gene_prior,
            "split_audit": audit, "resistance_diagnostics": diagnostics}


def component_gate(result):
    a, b = result["activation_metrics"][PRIMARY], result["activation_metrics"]["original"]
    required = [m[name] for m in (a, b) for name in ("macro_gene_spearman", "macro_gene_active_ap")]
    required.extend(m["pooled"]["mae"] for m in (a, b))
    for m in (a, b):
        required.extend(m["per_gene"][g][name] for g in GENES for name in ("spearman", "active_ap"))
        required.extend(f["macro_gene_spearman"] for f in m["by_fold"].values())
        # No silently dropping a scarce gene/fold from the activation classification evaluation.
        required.extend(f["macro_gene_active_ap"] for f in m["by_fold"].values())
    inconclusive = []
    if any(v is None or not math.isfinite(v) for v in required):
        inconclusive.append("undefined_gene_or_fold_metric")
    if any(t["outer_tail_fallback"] or any(i["tail_fallback"] for i in t["inner"])
           for t in result["split_audit"] if t["arm"] == PRIMARY):
        inconclusive.append("tail_training_fallback")
    if inconclusive:
        return {"status": "inconclusive", "reasons": inconclusive}
    reasons = []
    if a["macro_gene_spearman"] - b["macro_gene_spearman"] < .05:
        reasons.append("macro_gene_spearman_gain_below_0.05")
    if a["macro_gene_active_ap"] - b["macro_gene_active_ap"] < .02:
        reasons.append("active_macro_ap_gain_below_0.02")
    gains = [a["by_fold"][f]["macro_gene_spearman"] - b["by_fold"][f]["macro_gene_spearman"]
             for f in b["by_fold"]]
    if sum(g > 0 for g in gains) < 2:
        reasons.append("spearman_gain_in_fewer_than_two_folds")
    for gene in GENES:
        if a["per_gene"][gene]["spearman"] < b["per_gene"][gene]["spearman"] - .05:
            reasons.append(f"{gene}_spearman_drop")
        if a["per_gene"][gene]["active_ap"] < b["per_gene"][gene]["active_ap"] - .02:
            reasons.append(f"{gene}_active_ap_drop")
    if a["pooled"]["mae"] > b["pooled"]["mae"] * 1.05:
        reasons.append("mae_worse_than_five_percent_guard")
    return {"status": "advance_component" if not reasons else "do_not_advance",
            "reasons": reasons, "spearman_fold_deltas": gains,
            "interpretation": "Predeclared engineering triage; not a significance test."}


def resistance_gate(result):
    a, b = result["resistance_metrics"][PRIMARY], result["resistance_metrics"]["structure"]
    required = []
    inconclusive = []
    for metric in (a, b):
        required.append(metric.get("macro_ap"))
        required.extend(metric["per_drug"][drug].get(name)
                        for drug in ("pem", "fut") for name in ("average_precision", "brier"))
        required.extend(f.get("macro_ap") for f in metric.get("by_fold", {}).values())
        if any(f["per_drug"][drug]["positive_count"] == 0
               for f in metric.get("by_fold", {}).values() for drug in ("pem", "fut")):
            inconclusive.append("zero_positive_drug_fold")
    if any(type(v) not in (float, int) or not math.isfinite(v) for v in required):
        inconclusive.append("undefined_metric")
    if inconclusive:
        return {"status": "inconclusive", "reasons": sorted(set(inconclusive))}
    reasons = []
    if a["macro_ap"] - b["macro_ap"] < .005:
        reasons.append("macro_ap_gain_below_0.005")
    deltas = [a["by_fold"][f]["macro_ap"] - b["by_fold"][f]["macro_ap"] for f in b["by_fold"]]
    if sum(d > 0 for d in deltas) < 2:
        reasons.append("gain_in_fewer_than_two_folds")
    for drug in ("pem", "fut"):
        if a["per_drug"][drug]["average_precision"] < b["per_drug"][drug]["average_precision"] - .01:
            reasons.append(f"{drug}_ap_drop")
        if a["per_drug"][drug]["brier"] > b["per_drug"][drug]["brier"] * 1.05 + .001:
            reasons.append(f"{drug}_brier_drop")
    if any(d["fallback"] or d["convergence_warning"] for d in result["resistance_diagnostics"]):
        return {"status": "inconclusive", "reasons": ["fit_fallback_or_convergence_warning"]}
    return {"status": "advance_resistance" if not reasons else "do_not_advance",
            "reasons": reasons, "fold_macro_ap_deltas": deltas}
