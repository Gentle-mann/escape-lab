"""Matched continuous FGFR2 regressors with frozen conformational features.

The inherited regression head, sequence representation, and all preprocessing
are unchanged. Only the additional label-free structural columns differ.
"""

import copy
import math
from collections.abc import Mapping
from numbers import Real

from escape_lab.fgfr2_sequence_models import MODEL_CONFIG as SEQUENCE_MODEL_CONFIG
from escape_lab.fgfr2_sequence_models import FGFR2SequenceModel

ARMS = ("baseline", "inactive", "active", "contrast")
PRIMARY = "contrast"
BASELINE_ARM = "esm_biology_ridge"
MAX_STATE_FEATURES = 32
MODEL_CONFIG = {
    "arms": ARMS,
    "primary": PRIMARY,
    "baseline_arm": BASELINE_ARM,
    "ridge": copy.deepcopy(SEQUENCE_MODEL_CONFIG["ridge"]),
    "sequence": copy.deepcopy(SEQUENCE_MODEL_CONFIG["sequence"]),
    "target": "log2(1 + median of four codon-collapsed activation replicate summaries)",
    "maximum_added_features": MAX_STATE_FEATURES,
    "state_feature_scaling": "same train-only biology scaler used for all added columns",
    "negative_predictions": "retained without clipping, matching the previous ridge model",
}


class FGFR2StateModel(FGFR2SequenceModel):
    """Append one fixed state representation to the exact previous ridge model.

``arm`` remains the parent's internal sequence-regression selector. The public
pilot arm is ``state_arm`` and is recorded as ``diagnostics['arm']``.
"""

    def __init__(self, arm, store, augmenter, sequence_store, state_store):
        if arm not in ARMS:
            raise ValueError("Unknown FGFR2 state arm")
        if arm != "baseline" and state_store is None:
            raise ValueError("State features are required for structural arms")
        super().__init__(BASELINE_ARM, store, augmenter, sequence_store)
        self.state_arm = arm
        self.state_store = state_store
        self.state_feature_names = None

    def _biology_features(self, candidates):
        rows = super()._biology_features(candidates)
        if self.state_arm == "baseline":
            return rows
        for candidate, row in zip(candidates, rows, strict=True):
            extra = self.state_store.features(candidate, self.state_arm)
            if not isinstance(extra, Mapping) or not 1 <= len(extra) <= MAX_STATE_FEATURES:
                raise ValueError("State features must contain between 1 and 32 named values")
            if any(not isinstance(key, str) or not key for key in extra):
                raise ValueError("State feature names must be nonempty strings")
            if set(extra) & set(row):
                raise ValueError("State features must not overwrite existing biology features")
            if any(not isinstance(value, Real) or math.isinf(value) for value in extra.values()):
                raise ValueError("State features must be numeric and cannot contain infinity")
            names = tuple(sorted(extra))
            if self.state_feature_names is None:
                self.state_feature_names = names
            elif names != self.state_feature_names:
                raise ValueError("State feature schema must remain fixed across candidates")
            row.update(extra)
        return rows

    def fit(self, candidates, records):
        self.state_feature_names = None
        super().fit(candidates, records)
        self.diagnostics.update({
            "arm": self.state_arm,
            "baseline_arm": BASELINE_ARM,
            "state_feature_count": len(self.state_feature_names or ()),
            "state_feature_names": list(self.state_feature_names or ()),
            "matched_ridge_settings": copy.deepcopy(MODEL_CONFIG["ridge"]),
        })
        return self
