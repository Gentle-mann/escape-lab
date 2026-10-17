"""Compact, drug-specific biological predictors with train-only preprocessing."""

import math
import warnings
from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_extraction import DictVectorizer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

MISSING_MEASUREMENTS = frozenset({"structure_sasa", "structure_ligand_atom_pairs_4a",
                                 "structure_ligand_atom_pairs_5a", "structure_ca_neighbors_8a"})


class BioConfig(BaseModel):
    model_config = ConfigDict(frozen=True, strict=True, extra="forbid")
    representation: Literal["chemistry", "structure", "interactions"] = "chemistry"
    model: Literal["logistic", "hist_gradient_boosting"] = "logistic"
    C: Literal[0.01, 0.1, 1.0] = 0.1
    max_depth: Literal[2, 3] = 2

    @model_validator(mode="before")
    @classmethod
    def reject_coercion(cls, value):
        if not isinstance(value, dict):
            return value
        value = value.copy()
        if "C" in value:
            if type(value["C"]) not in (int, float) or not math.isfinite(value["C"]):
                raise ValueError("C requires a finite numerical literal")
            value["C"] = float(value["C"])
        if "max_depth" in value and type(value["max_depth"]) is not int:
            raise ValueError("max_depth requires an integer literal")
        return value


def catalog():
    result = {}
    for representation in ("chemistry", "structure", "interactions"):
        for penalty in (.01, .1, 1.):
            result[f"{representation}-logistic-c{penalty:g}"] = BioConfig(
                representation=representation, C=penalty)
    for representation, depth in (("structure", 2), ("structure", 3), ("interactions", 2)):
        result[f"{representation}-hgb-depth{depth}"] = BioConfig(
            representation=representation, model="hist_gradient_boosting", max_depth=depth)
    return result


class BiologyPredictor:
    def __init__(self, config, feature_store, drug, seed=0):
        if not isinstance(config, BioConfig) or drug not in ("pem", "fut"):
            raise ValueError("Validated biological configuration and drug required")
        self.config, self.feature_store, self.drug, self.seed = config, feature_store, drug, seed
        self.vectorizer = DictVectorizer(sparse=False)
        self.imputer = SimpleImputer(strategy="median", keep_empty_features=True)
        self.scaler = StandardScaler()
        self.model = None
        self.fallback_reason = None
        self.prior = .5
        self.fitted = False
        self.diagnostics = {}

    def _features(self, candidates):
        rows = [self.feature_store.features(c, self.drug, self.config.representation)
                for c in candidates]
        # The feature contract uses 999 for unavailable geometry. Preserve explicit
        # availability columns and impute only from the fitted training partition.
        result = []
        for candidate, row in zip(candidates, rows, strict=True):
            gene = candidate.get("gene")
            if gene not in ("FGFR1", "FGFR2", "FGFR3", "FGFR4"):
                raise ValueError("Explicit FGFR gene identity is required")
            features = {key: np.nan if value == 999 or key in MISSING_MEASUREMENTS and value == -1
                        else value for key, value in row.items()}
            features.update({f"gene_{g}": float(gene == g)
                             for g in ("FGFR1", "FGFR2", "FGFR3", "FGFR4")})
            result.append(features)
        return result

    def fit(self, candidates, labels):
        if len(candidates) != len(labels) or any(type(y) is not bool for y in labels):
            raise ValueError("Training requires one defined boolean per example")
        self.fitted = False
        self.model = None
        self.fallback_reason = None
        self.prior = (sum(labels) + 1) / (len(labels) + 2)
        self.diagnostics = {"training_n": len(labels), "positive_count": sum(labels),
                            "convergence_warning": False}
        features = self._features(candidates)
        if not labels or len(set(labels)) < 2:
            self.fallback_reason = "empty_training" if not labels else "single_class_training"
            self.fitted = True
            return self
        matrix = self.vectorizer.fit_transform(features)
        matrix = self.imputer.fit_transform(matrix)
        matrix = self.scaler.fit_transform(matrix)
        if self.config.model == "logistic":
            self.model = LogisticRegression(C=self.config.C, max_iter=2000, solver="lbfgs",
                                            random_state=self.seed)
        else:
            self.model = HistGradientBoostingClassifier(
                max_iter=80, max_depth=self.config.max_depth, max_leaf_nodes=7,
                min_samples_leaf=30, l2_regularization=10., early_stopping=False,
                random_state=self.seed)
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always", ConvergenceWarning)
            self.model.fit(matrix, np.asarray(labels, dtype=int))
        self.diagnostics["convergence_warning"] = any(
            issubclass(w.category, ConvergenceWarning) for w in captured)
        if hasattr(self.model, "coef_"):
            self.diagnostics["coefficient_norm"] = float(np.linalg.norm(self.model.coef_))
        self.fitted = True
        return self

    def predict(self, candidates):
        if not self.fitted:
            raise RuntimeError("Fit before prediction")
        if not candidates:
            return np.asarray([], dtype=float)
        features = self._features(candidates)
        if self.fallback_reason:
            return np.full(len(candidates), self.prior)
        matrix = self.vectorizer.transform(features)
        matrix = self.scaler.transform(self.imputer.transform(matrix))
        positive = int(np.flatnonzero(self.model.classes_ == 1)[0])
        scores = self.model.predict_proba(matrix)[:, positive]
        if not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 1)):
            raise RuntimeError("Invalid probability")
        return scores
