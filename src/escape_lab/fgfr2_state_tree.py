"""Tiny continuous tree sensitivity using only frozen FGFR2 state contrasts.

The inherited biology HGB model and preprocessing are unchanged. Pooled FGFR
training is permitted, but only FGFR2 receives cached contrast features; all
other genes receive zeros, without mapping or opening additional structures.
"""

import math
from numbers import Real

from escape_lab.activation_revision import GENES, ActivationRevisionModel
from escape_lab.fgfr2_state_features import FEATURE_NAMES


def _scope(candidates):
    candidates = list(candidates)
    if not candidates:
        raise ValueError("nonempty development candidates required")
    ids, groups = [], []
    for candidate in candidates:
        if (candidate.get("gene") not in GENES
                or candidate.get("split") not in {"development", "validation"}):
            raise ValueError("only known FGFR development/validation candidates permitted")
        candidate_id, group = candidate.get("candidate_id"), candidate.get("group_id")
        if not isinstance(candidate_id, str) or not candidate_id:
            raise ValueError("candidate identity must be a nonempty string")
        if not isinstance(group, str) or not group:
            raise ValueError("group identity must be a nonempty string")
        ids.append(candidate_id)
        groups.append(group)
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate candidate identity")
    return candidates, set(ids), set(groups)


class FGFR2StateTreeModel(ActivationRevisionModel):
    def __init__(self, store, augmenter, state_store, add_contrast=False):
        if type(add_contrast) is not bool:
            raise ValueError("add_contrast must be boolean")
        if add_contrast and state_store is None:
            raise ValueError("contrast arm requires a state feature store")
        super().__init__(store, augmenter, "biology")
        self.state_store = state_store
        self.add_contrast = add_contrast
        self.fitted = False
        self.training_ids, self.training_groups = set(), set()
        self.diagnostics = {}

    def features(self, candidates):
        candidates, _, _ = _scope(candidates)
        rows = super().features(candidates)
        if not self.add_contrast:
            return rows
        for candidate, row in zip(candidates, rows, strict=True):
            if set(row) & set(FEATURE_NAMES):
                raise ValueError("state features must not overwrite existing biology")
            extra = (self.state_store.features(candidate, "contrast")
                     if candidate["gene"] == "FGFR2" else dict.fromkeys(FEATURE_NAMES, 0.0))
            if not isinstance(extra, dict) or set(extra) != set(FEATURE_NAMES):
                raise ValueError("state contrast schema must equal the frozen 31 features")
            if any(not isinstance(value, Real) or not math.isfinite(value)
                   for value in extra.values()):
                raise ValueError("state contrasts must be finite numeric values")
            row.update(extra)
        return rows

    def fit(self, candidates, records):
        self.fitted = False
        candidates, ids, groups = _scope(candidates)
        super().fit(candidates, records)
        self.training_ids, self.training_groups = ids, groups
        self.diagnostics = {"training_n": len(ids), "training_groups": sorted(groups),
                            "feature_count": len(self.vectorizer.get_feature_names_out()),
                            "state_feature_count": len(FEATURE_NAMES) if self.add_contrast else 0,
                            "add_contrast": self.add_contrast,
                            "state_gene": "FGFR2", "parent_arm": "biology"}
        self.fitted = True
        return self

    def predict(self, candidates):
        if not self.fitted:
            raise ValueError("model must be successfully fitted before prediction")
        candidates, ids, groups = _scope(candidates)
        if ids & self.training_ids or groups & self.training_groups:
            raise ValueError("prediction candidates overlap training identities or groups")
        return super().predict(candidates)
