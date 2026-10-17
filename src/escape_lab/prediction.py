"""Bounded resistance predictors using public sequence features and disclosed labels.

No split assignment, outcome lookup, evaluation, or network access occurs here.
The trusted evaluator supplies training labels and keeps test outcomes separate.
"""

import hashlib
import json
import math
from typing import Literal

import numpy as np
from Bio.Align import substitution_matrices
from Bio.SeqUtils.ProtParamData import kd
from pydantic import BaseModel, ConfigDict, model_validator
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
BLOSUM = substitution_matrices.load("BLOSUM62")
PUBLIC_FEATURE_FIELDS = frozenset({"wt", "mut", "gene", "aligned_position", "context"})


class PredictionConfig(BaseModel):
    """Strict proposal surface. ``regularization`` is sklearn's inverse-penalty C."""

    model_config = ConfigDict(frozen=True, strict=True, extra="forbid")

    model: Literal["logistic", "extra_trees", "hist_gradient_boosting"] = "extra_trees"
    feature_set: Literal["sequence", "position", "full"] = "full"
    include_gene: bool = True
    position_scale: Literal[10, 20, 40] = 20
    regularization: Literal[0.1, 1.0, 10.0] = 1.0
    min_samples_leaf: Literal[1, 3, 10] = 3
    max_features: Literal[0.5, 0.8, 1.0] = 0.8
    max_depth: Literal[3, 6, None] = None  # noqa: PYI061 -- keep the frozen single-enum schema.
    learning_rate: Literal[0.05, 0.1] = 0.1
    l2: Literal[0.0, 1.0, 10.0] = 1.0

    @model_validator(mode="before")
    @classmethod
    def reject_numeric_coercion(cls, value):
        # Literal comparisons otherwise accept True == 1 and False == 0 even
        # under strict mode. Normalize real JSON numbers for semantic hashing.
        if not isinstance(value, dict):
            return value
        result = value.copy()
        for field in ("position_scale", "min_samples_leaf", "max_depth"):
            if field not in result or field == "max_depth" and result[field] is None:
                continue
            if type(result[field]) is not int:
                raise ValueError(f"{field} requires an integer literal")
        for field in ("regularization", "max_features", "learning_rate", "l2"):
            if field not in result:
                continue
            if type(result[field]) not in (int, float) or not math.isfinite(result[field]):
                raise ValueError(f"{field} requires a finite numerical literal")
            result[field] = float(result[field])
        return result

    def effective_config(self):
        """Canonical active settings for deduplication, not a full input schema.

        Use model_dump() to reconstruct a complete PredictionConfig. Constants
        are included here so a model's fixed computational bounds are recorded.
        """
        result = {
            "model": self.model,
            "feature_set": self.feature_set,
            "include_gene": self.include_gene,
        }
        if self.feature_set in {"position", "full"}:
            result["position_scale"] = self.position_scale
        if self.model == "logistic":
            result.update(regularization=self.regularization, max_iter=1000, solver="lbfgs")
        else:
            result.update(
                min_samples_leaf=self.min_samples_leaf,
                max_features=self.max_features,
                max_depth=self.max_depth,
            )
            if self.model == "extra_trees":
                result.update(n_estimators=100)
            else:
                result.update(
                    learning_rate=self.learning_rate,
                    l2=self.l2,
                    max_iter=100,
                    early_stopping=False,
                    max_leaf_nodes=31,
                )
        return result

    @property
    def config_id(self):
        encoded = json.dumps(
            self.effective_config(), sort_keys=True, separators=(",", ":"), allow_nan=False
        )
        return hashlib.sha256(encoded.encode()).hexdigest()


def candidate_features(candidate, config):
    """Read only named sequence inputs; never candidate IDs, labels or split IDs."""
    public = {key: candidate[key] for key in PUBLIC_FEATURE_FIELDS if key in candidate}
    wt, mut = public["wt"], public["mut"]
    if wt not in AMINO_ACIDS or mut not in AMINO_ACIDS:
        raise ValueError("Wild-type and mutant residues must be standard amino acids")
    features = {
        "wt": wt,
        "mut": mut,
        "substitution": wt + mut,
        "blosum": float(BLOSUM[wt, mut]) / 10,
        "hydropathy_delta": (kd[mut] - kd[wt]) / 10,
        "charge_delta": int(mut in "KR") - int(mut in "DE") - int(wt in "KR") + int(wt in "DE"),
        "aromatic_delta": int(mut in "FWY") - int(wt in "FWY"),
        "proline_delta": int(mut == "P") - int(wt == "P"),
        "glycine_delta": int(mut == "G") - int(wt == "G"),
    }
    if config.include_gene:
        gene = public["gene"]
        if not isinstance(gene, str) or not gene:
            raise ValueError("Gene identity must be a nonempty string")
        features["gene"] = gene
    if config.feature_set in {"position", "full"}:
        position = public["aligned_position"]
        if type(position) is not int or position < 1:
            raise ValueError("Aligned position must be a positive integer")
        features["position"] = position / 1000
        for center in range(450, 811, 20):
            features[f"position_rbf_{center}"] = math.exp(
                -(((position - center) / config.position_scale) ** 2)
            )
    if config.feature_set in {"sequence", "full"}:
        context = public.get("context", "")
        if (
            not isinstance(context, str)
            or len(context) > 7
            or any(residue not in AMINO_ACIDS | {"X"} for residue in context)
        ):
            raise ValueError("Sequence context must contain at most seven amino-acid residues")
        for offset, residue in enumerate(context):
            features[f"context_{offset}"] = residue
    return features


class ResistancePredictor:
    """A pickleable numerical model with explicit no-data/single-class behavior."""

    def __init__(self, config, seed=0):
        if not isinstance(config, PredictionConfig):
            raise TypeError("config must be a validated PredictionConfig")
        if type(seed) is not int or seed < 0:
            raise ValueError("seed must be a nonnegative integer")
        self.config = config
        self.seed = seed
        self.vectorizer = DictVectorizer(sparse=False)
        self.model = None
        self.prior = 0.5
        self.training_size = 0
        self.fallback_reason = None
        self._fitted = False

    @property
    def using_prior_fallback(self):
        return self._fitted and self.fallback_reason is not None

    def fit(self, candidates, labels):
        if len(candidates) != len(labels) or any(type(label) is not bool for label in labels):
            raise ValueError("Training requires one defined boolean label per candidate")
        self._fitted = False
        self.model = None
        self.vectorizer = DictVectorizer(sparse=False)
        self.training_size = len(candidates)
        self.prior = (sum(labels) + 1) / (len(labels) + 2)
        self.fallback_reason = None
        if not candidates:
            self.fallback_reason = "empty_training"
            self._fitted = True
            return self
        matrix = self.vectorizer.fit_transform(
            [candidate_features(candidate, self.config) for candidate in candidates]
        )
        if len(set(labels)) == 1:
            self.fallback_reason = "single_class_training"
            self._fitted = True
            return self
        config = self.config
        if config.model == "logistic":
            model = LogisticRegression(
                C=config.regularization, max_iter=1000, solver="lbfgs", random_state=self.seed
            )
        elif config.model == "extra_trees":
            model = ExtraTreesClassifier(
                n_estimators=100,
                min_samples_leaf=config.min_samples_leaf,
                max_features=config.max_features,
                max_depth=config.max_depth,
                random_state=self.seed,
                n_jobs=1,
            )
        else:
            model = HistGradientBoostingClassifier(
                max_iter=100,
                early_stopping=False,
                min_samples_leaf=config.min_samples_leaf,
                max_features=config.max_features,
                max_depth=config.max_depth,
                max_leaf_nodes=31,
                learning_rate=config.learning_rate,
                l2_regularization=config.l2,
                random_state=self.seed,
            )
        model.fit(matrix, np.asarray(labels, dtype=int))
        self.model = model
        self._fitted = True
        return self

    def predict(self, candidates):
        if not self._fitted:
            raise RuntimeError("Call fit before requesting predictions")
        if not candidates:
            return np.asarray([], dtype=float)
        if self.using_prior_fallback:
            return np.full(len(candidates), self.prior)
        matrix = self.vectorizer.transform(
            [candidate_features(candidate, self.config) for candidate in candidates]
        )
        positive_class = int(np.flatnonzero(self.model.classes_ == 1)[0])
        probability = self.model.predict_proba(matrix)[:, positive_class]
        if not np.isfinite(probability).all() or np.any((probability < 0) | (probability > 1)):
            raise RuntimeError("Predictor returned invalid probabilities")
        return probability


def fit_predict(config, training, labels, pool, seed=0):
    return ResistancePredictor(config, seed=seed).fit(training, labels).predict(pool)
