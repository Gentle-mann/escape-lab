"""Fixed low-capacity heads for the development-only FGFR2 sequence pilot.

Sequence embeddings are label-free. Every fitted transform is trained solely
on the supplied training partition; scoring a training group is forbidden.
Cached historical predictions and fold membership remain evaluator-owned.
"""

import math
import warnings

import numpy as np
from scipy.special import expit
from sklearn.decomposition import PCA
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_extraction import DictVectorizer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler

from escape_lab.mechanistic_pilot import ActivationRegressor

ARMS = ("cached_v2_biology", "biology_logistic", "esm_biology_logistic", "esm_biology_ridge")
TRAINABLE_ARMS = ARMS[1:]
PRIMARY = "esm_biology_logistic"
SEED = 17
EMBEDDING_DIMENSION = 480
PCA_COMPONENTS = 16
BIOLOGY_SCALE_FLOOR = 0.01
STABILITY_MODES = ("legacy", "variance_floor")
SEQUENCE_BRANCHES = ("both", "wt", "delta")
MODEL_CONFIG = {
    "primary": PRIMARY,
    "logistic": {"C": 0.1, "l1_ratio": 0.0, "class_weight": "balanced",
                 "solver": "lbfgs", "max_iter": 2000, "random_state": SEED},
    "ridge": {"alpha": 100.0, "solver": "svd"},
    "sequence": {"branch_dimensions": EMBEDDING_DIMENSION, "pca_components_per_branch": 16,
                 "pca_solver": "full", "random_state": SEED,
                 "branches": ["wt", "mutant_minus_wt"], "train_only_scaling_and_pca": True},
    "target": {"consistent_replicates": 3, "consistency_threshold_inclusive": 1.5,
               "high_median_threshold_strict": 5.0, "strong_median_threshold_strict": 10.0},
}


def _validated_record(record):
    try:
        reps = np.asarray(record["replicates"], dtype=float)
        value = float(record["value"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Activation record requires four replicates and a target") from exc
    if reps.shape != (4,) or not np.isfinite(reps).all() or (reps < 0).any():
        raise ValueError("Activation record requires four finite nonnegative replicates")
    median = float(np.median(reps))
    if not math.isfinite(value) or not math.isclose(
        value, math.log2(1 + median), rel_tol=1e-10, abs_tol=1e-10
    ):
        raise ValueError("Activation target disagrees with replicate median")
    consistent = bool(np.count_nonzero(reps >= 1.5) >= 3)
    if "active" in record and (
        type(record["active"]) is not bool or record["active"] != consistent
    ):
        raise ValueError("Activation consistency label disagrees with replicates")
    return median, value, consistent


def activation_label(record, threshold=5.0):
    """Use the frozen strict median threshold and inclusive replicate criterion."""
    if not math.isfinite(threshold) or threshold < 1.5:
        raise ValueError("Invalid activation threshold")
    median, _, consistent = _validated_record(record)
    return bool(consistent and median > threshold)


def _candidate_scope(candidates):
    candidates = list(candidates)
    if not candidates:
        raise ValueError("Nonempty candidate partition required")
    required = ("candidate_id", "group_id", "gene", "split")
    if any(any(key not in candidate for key in required) for candidate in candidates):
        raise ValueError("Incomplete candidate identity")
    if any(candidate["gene"] != "FGFR2"
           or candidate["split"] not in ("development", "validation")
           for candidate in candidates):
        raise ValueError("Only FGFR2 development candidates are permitted")
    ids = [candidate["candidate_id"] for candidate in candidates]
    groups = [candidate["group_id"] for candidate in candidates]
    if any(not isinstance(item, str) or not item for item in [*ids, *groups]):
        raise ValueError("Candidate and group identities must be nonempty strings")
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate candidate identity")
    return candidates, ids, set(groups)


class FGFR2SequenceModel:
    """One fixed supervised head; the caller supplies only approved candidates."""

    def __init__(self, arm, store, augmenter, sequence_store=None, *, stability="legacy",
                 sequence_branches="both"):
        if arm == "cached_v2_biology":
            raise ValueError("The cached reference arm is evaluator-owned")
        if arm not in TRAINABLE_ARMS:
            raise ValueError("Unknown FGFR2 sequence arm")
        if arm.startswith("esm_") and sequence_store is None:
            raise ValueError("Sequence arm requires frozen sequence features")
        if stability not in STABILITY_MODES:
            raise ValueError("Unknown biology scaling stability mode")
        if sequence_branches not in SEQUENCE_BRANCHES:
            raise ValueError("Unknown sequence branch selection")
        self.arm = arm
        self.stability = stability
        self.sequence_branches = sequence_branches
        self.store = store
        self.augmenter = augmenter
        self.sequence_store = sequence_store
        self.original = ActivationRegressor(store)
        self.vectorizer = DictVectorizer(sparse=False)
        self.imputer = SimpleImputer(strategy="median", keep_empty_features=True)
        self.scaler = StandardScaler()
        self.wt_scaler, self.delta_scaler = StandardScaler(), StandardScaler()
        self.wt_pca = PCA(n_components=PCA_COMPONENTS, svd_solver="full", random_state=SEED)
        self.delta_pca = PCA(n_components=PCA_COMPONENTS, svd_solver="full", random_state=SEED)
        self.model = (Ridge(**MODEL_CONFIG["ridge"]) if arm.endswith("ridge")
                      else LogisticRegression(**MODEL_CONFIG["logistic"]))
        self.fitted = False
        self.training_ids = set()
        self.training_groups = set()
        self.diagnostics = {}
        self.biology_scaling = {}
        self._prediction_diagnostics = {}

    @property
    def prediction_diagnostics(self):
        """Identity-free aggregates for the last guarded prediction, returned by value."""
        return dict(self._prediction_diagnostics)

    def _biology_features(self, candidates):
        rows = self.original.features(candidates)
        for candidate, row in zip(candidates, rows, strict=True):
            extra = self.augmenter.features(candidate)
            if set(extra) & set(row):
                raise ValueError("Biology augmentation must not overwrite original features")
            row.update(extra)
        return rows

    def _sequence_matrix(self, ids):
        values = np.asarray(self.sequence_store.matrix(ids, "wt_delta"), dtype=float)
        if values.shape != (len(ids), EMBEDDING_DIMENSION * 2):
            raise ValueError("Sequence matrix must have exactly 960 columns")
        if not np.isfinite(values).all():
            raise ValueError("Sequence features must be finite")
        return values[:, :EMBEDDING_DIMENSION], values[:, EMBEDDING_DIMENSION:]

    def _fit_matrix(self, candidates, ids):
        matrix = self.vectorizer.fit_transform(self._biology_features(candidates))
        matrix = self.imputer.fit_transform(matrix)
        self.scaler.fit(matrix)
        original_scale = self.scaler.scale_.copy()
        floored = original_scale < BIOLOGY_SCALE_FLOOR
        if self.stability == "variance_floor":
            self.scaler.scale_ = np.maximum(original_scale, BIOLOGY_SCALE_FLOOR)
        self.biology_scaling = {
            "mode": self.stability, "floor": BIOLOGY_SCALE_FLOOR,
            "minimum_training_scale": float(original_scale.min()),
            "minimum_effective_scale": float(self.scaler.scale_.min()),
            "below_floor_count": int(floored.sum()),
            "floored_feature_names": self.vectorizer.get_feature_names_out()[floored].tolist()
            if self.stability == "variance_floor" else [],
        }
        matrix = self.scaler.transform(matrix)
        if self.arm.startswith("esm_"):
            if len(candidates) < PCA_COMPONENTS + 1:
                raise ValueError("Sequence fitting requires at least 17 training variants")
            wt, delta = self._sequence_matrix(ids)
            branches = [matrix]
            if self.sequence_branches in ("both", "wt"):
                branches.append(self.wt_pca.fit_transform(self.wt_scaler.fit_transform(wt)))
            if self.sequence_branches in ("both", "delta"):
                branches.append(self.delta_pca.fit_transform(self.delta_scaler.fit_transform(delta)))
            matrix = np.column_stack(branches)
        if not np.isfinite(matrix).all():
            raise ValueError("Nonfinite fitted model features")
        return matrix

    def _predict_matrix(self, candidates, ids):
        matrix = self.vectorizer.transform(self._biology_features(candidates))
        matrix = self.scaler.transform(self.imputer.transform(matrix))
        if self.arm.startswith("esm_"):
            wt, delta = self._sequence_matrix(ids)
            branches = [matrix]
            if self.sequence_branches in ("both", "wt"):
                branches.append(self.wt_pca.transform(self.wt_scaler.transform(wt)))
            if self.sequence_branches in ("both", "delta"):
                branches.append(self.delta_pca.transform(self.delta_scaler.transform(delta)))
            matrix = np.column_stack(branches)
        if not np.isfinite(matrix).all():
            raise ValueError("Nonfinite prediction features")
        return matrix

    def fit(self, candidates, records):
        self.fitted = False
        self.diagnostics = {}
        self._prediction_diagnostics = {}
        candidates, ids, groups = _candidate_scope(candidates)
        if set(records) != set(ids):
            raise ValueError("Activation records must exactly match training candidates")
        labels = np.asarray([activation_label(records[cid]) for cid in ids], dtype=int)
        values = np.asarray([_validated_record(records[cid])[1] for cid in ids])
        if not self.arm.endswith("ridge") and len(np.unique(labels)) != 2:
            raise ValueError("Classifier training requires two classes")
        matrix = self._fit_matrix(candidates, ids)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", ConvergenceWarning)
                self.model.fit(matrix, values if self.arm.endswith("ridge") else labels)
        except ConvergenceWarning as exc:
            raise ValueError("FGFR2 classifier failed to converge") from exc
        if not np.isfinite(self.model.coef_).all():
            raise ValueError("Nonfinite fitted model coefficients")
        self.training_ids, self.training_groups = set(ids), groups
        self.diagnostics = {
            "arm": self.arm, "training_n": len(ids), "high_count": int(labels.sum()),
            "training_groups": sorted(groups), "feature_count": int(matrix.shape[1]),
            "biology_feature_count": len(self.vectorizer.get_feature_names_out()),
            "sequence_components": {branch: PCA_COMPONENTS for branch in ("wt", "delta")
                                    if self.sequence_branches in ("both", branch)}
            if self.arm.startswith("esm_") else {},
            "sequence_branches": self.sequence_branches,
            "biology_scaling": self.biology_scaling,
            "coefficient_norm": float(np.linalg.norm(self.model.coef_)),
            "iterations": int(np.max(self.model.n_iter_))
            if hasattr(self.model, "n_iter_") and self.model.n_iter_ is not None else None,
            "convergence_warning": False,
        }
        self.fitted = True
        return self

    def _guarded_prediction_matrix(self, candidates):
        self._prediction_diagnostics = {}
        if not self.fitted:
            raise ValueError("Model must be successfully fitted before prediction")
        candidates, ids, groups = _candidate_scope(candidates)
        if set(ids) & self.training_ids or groups & self.training_groups:
            raise ValueError("Prediction candidates overlap the training partition")
        matrix = self._predict_matrix(candidates, ids)
        # Biology columns are standardized. Appended PCA scores have a different
        # scale and must not be mistaken for standardized feature-shift evidence.
        biology = matrix[:, :len(self.scaler.scale_)]
        self._prediction_diagnostics = {
            "scaled_feature_abs_max": float(min(np.max(np.abs(biology)), 1e12)),
            "feature_shift_fraction": float(np.mean(np.abs(biology) > 8.)),
            "near_constant_feature_count": int(np.count_nonzero(self.scaler.var_ < 1e-12)),
        }
        return matrix

    @staticmethod
    def _validate_scores(scores, n):
        scores = np.asarray(scores, dtype=float)
        if scores.shape != (n,) or not np.isfinite(scores).all():
            raise ValueError("Nonfinite or malformed prediction scores")
        return scores

    def ranking_scores(self, candidates):
        """Return continuous scores for AP and selection without sigmoid rounding."""
        matrix = self._guarded_prediction_matrix(candidates)
        scores = (self.model.predict(matrix) if self.arm.endswith("ridge")
                  else self.model.decision_function(matrix))
        scores = self._validate_scores(scores, len(matrix))
        self._prediction_diagnostics.update({
            "ranking_score_abs_max": float(min(np.max(np.abs(scores)), 1e12)),
            "score_count": len(scores), "nonfinite_score_fraction": 0.,
        })
        return scores

    def prediction_outputs(self, candidates):
        """Produce a classifier's probabilities and margins from one guarded transform.

        Class balancing means these probabilities are not calibration estimates.
        The logistic transform exactly matches the historical binary sklearn head.
        """
        if self.arm.endswith("ridge"):
            raise ValueError("Probability outputs require a classifier")
        matrix = self._guarded_prediction_matrix(candidates)
        if not np.array_equal(self.model.classes_, [0, 1]):
            raise ValueError("Classifier must use the frozen binary activation labels")
        margin = self._validate_scores(self.model.decision_function(matrix), len(matrix))
        probability = self._validate_scores(expit(margin), len(matrix))
        self._prediction_diagnostics.update({
            "ranking_score_abs_max": float(min(np.max(np.abs(margin)), 1e12)),
            "probability_saturated_fraction": float(np.mean((probability == 0) | (probability == 1))),
            "probability_unique_fraction": float(len(np.unique(probability)) / len(probability)),
            "score_count": len(probability), "nonfinite_score_fraction": 0.,
        })
        return {"probability": probability, "margin": margin}

    def predict_probability(self, candidates):
        """Return the historical sigmoid output; use ranking_scores for retrieval."""
        return self.prediction_outputs(candidates)["probability"]

    def predict(self, candidates):
        """Backward-compatible probability output, or continuous ridge output."""
        if self.arm.endswith("ridge"):
            return self.ranking_scores(candidates)
        return self.predict_probability(candidates)
