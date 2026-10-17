"""Development sensitivity: bound only standardized biology extrapolation.

Training matrices, fitted coefficients, sequence PCA, and ridge settings match
the frozen state pilot exactly. Prediction biology columns are clamped to their
own training minima/maxima; this is not an output-range or label-derived bound.
"""

import numpy as np

from escape_lab.fgfr2_state_models import FGFR2StateModel


class BoundedFGFR2StateModel(FGFR2StateModel):
    """Apply one parameter-free prediction-only sensitivity to all matched arms."""

    def _fit_matrix(self, candidates, ids):
        matrix = super()._fit_matrix(candidates, ids)
        self.biology_feature_count_ = len(self.vectorizer.get_feature_names_out())
        biology = matrix[:, :self.biology_feature_count_]
        self.biology_minimum_ = biology.min(axis=0).copy()
        self.biology_maximum_ = biology.max(axis=0).copy()
        self.biology_minimum_.setflags(write=False)
        self.biology_maximum_.setflags(write=False)
        return matrix

    def fit(self, candidates, records):
        self.unbounded_predictions = None
        super().fit(candidates, records)
        self.diagnostics["bounded_extrapolation"] = {
            "rule": "clamp standardized biology to per-column training minima and maxima",
            "biology_feature_count": self.biology_feature_count_,
            "feature_names": list(self.vectorizer.get_feature_names_out()),
            "training_minimum": self.biology_minimum_.tolist(),
            "training_maximum": self.biology_maximum_.tolist(),
            "training_matrix_changed": False,
            "sequence_branch_changed": False,
            "bounds_are_fixed_at_fit": True,
        }
        return self

    def _predict_matrix(self, candidates, ids):
        matrix = super()._predict_matrix(candidates, ids)
        self.unbounded_predictions = np.asarray(self.model.predict(matrix), dtype=float).copy()
        biology = matrix[:, :self.biology_feature_count_]
        bounded_biology = np.clip(biology, self.biology_minimum_, self.biology_maximum_)
        changed = bounded_biology != biology
        diagnostics = {
            "prediction_rows": len(candidates),
            "clipped_value_count": int(changed.sum()),
            "clipped_row_count": int(np.any(changed, axis=1).sum()),
            "clipped_row_fraction": float(np.any(changed, axis=1).mean()),
            "maximum_excess_beyond_training_bound": float(np.abs(biology - bounded_biology).max()),
            "raw_max_abs_biology": float(np.abs(biology).max()),
            "bounded_max_abs_biology": float(np.abs(bounded_biology).max()),
        }
        self.diagnostics["bounded_extrapolation"].update(diagnostics)
        result = matrix.copy()
        result[:, :self.biology_feature_count_] = bounded_biology
        return result
