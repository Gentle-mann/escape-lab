"""The tiny tree sensitivity adds only FGFR2's frozen contrast features."""

import copy

import numpy as np
import pytest

from escape_lab.activation_revision import CAPACITY, ActivationRevisionModel
from escape_lab.fgfr2_state_features import FEATURE_NAMES
from escape_lab.fgfr2_state_tree import FGFR2StateTreeModel


class Biology:
    def features(self, candidate, drug, family):
        return {"position": float(candidate["position"]), "missing": np.nan}


class Augmenter:
    def features(self, candidate):
        return {"context": candidate["position"] % 7}


class StateStore:
    def __init__(self):
        self.calls = []

    def features(self, candidate, arm):
        self.calls.append((candidate["candidate_id"], arm))
        return dict.fromkeys(FEATURE_NAMES, float(candidate["position"] % 3))


def candidates(prefix="train", n=64):
    return [{"candidate_id": f"{prefix}-{i}", "group_id": f"{prefix}-group-{i // 16}",
             "gene": f"FGFR{1 + i % 4}", "split": "development", "position": 500 + i}
            for i in range(n)]


def records(rows):
    return {row["candidate_id"]: {"value": float(np.log2(1 + row["position"] % 13))}
            for row in rows}


def model(add_contrast=False):
    return FGFR2StateTreeModel(Biology(), Augmenter(), StateStore(), add_contrast=add_contrast)


def test_baseline_exactly_matches_parent_and_never_calls_state_store():
    train, test = candidates(), candidates("test", 8)
    old = ActivationRevisionModel(Biology(), Augmenter(), "biology").fit(train, records(train))
    new = model().fit(train, records(train))
    old_prediction, new_prediction = old.predict(test), new.predict(test)
    for key in ("value", "active_score"):
        np.testing.assert_array_equal(old_prediction[key], new_prediction[key])
    assert new.state_store.calls == []
    assert new.features(test) == old.features(test)
    for name, value in CAPACITY.items():
        assert new.regressor.get_params()[name] == value


def test_contrast_uses_only_fgfr2_and_zeros_other_genes():
    rows = candidates(n=4)
    current = model(add_contrast=True)
    original_rows = ActivationRevisionModel(Biology(), Augmenter(), "biology").features(rows)
    combined = current.features(rows)
    assert current.state_store.calls == [(rows[1]["candidate_id"], "contrast")]
    for candidate, before, after in zip(rows, original_rows, combined, strict=True):
        assert set(after) == set(before) | set(FEATURE_NAMES)
        for name in FEATURE_NAMES:
            assert after[name] == (float(candidate["position"] % 3)
                                   if candidate["gene"] == "FGFR2" else 0.0)


@pytest.mark.parametrize("change", [{"split": "final"}, {"gene": "BRAF"},
                                     {"group_id": ""}, {"candidate_id": ""}])
def test_invalid_scope_rejected_before_state_access(change):
    train = candidates()
    train[0].update(change)
    current = model(add_contrast=True)
    with pytest.raises(ValueError):
        current.fit(train, records(train))
    assert current.state_store.calls == []


def test_duplicate_identity_rejected():
    train = candidates()
    train.append(copy.deepcopy(train[0]))
    with pytest.raises(ValueError, match="duplicate"):
        model().fit(train, records(train))


@pytest.mark.parametrize("change", [{"group_id": "train-group-0"}, {"split": "final"}])
def test_prediction_keeps_group_and_final_guards(change):
    train, test = candidates(), candidates("test", 4)
    current = model(add_contrast=True).fit(train, records(train))
    test[0].update(change)
    with pytest.raises(ValueError):
        current.predict(test)


def test_state_schema_and_collisions_rejected():
    current = model(add_contrast=True)
    rows = candidates(n=4)
    current.state_store.features = lambda candidate, arm: {"unexpected": 1.0}
    with pytest.raises(ValueError, match="schema"):
        current.features(rows)
    current = model(add_contrast=True)
    current.augmenter.features = lambda candidate: {FEATURE_NAMES[0]: 1.0}
    with pytest.raises(ValueError, match="overwrite"):
        current.features(rows)


def test_prefit_predict_and_missing_contrast_store_fail():
    with pytest.raises(ValueError, match="fitted"):
        model().predict(candidates("test", 4))
    with pytest.raises(ValueError, match="state"):
        FGFR2StateTreeModel(Biology(), Augmenter(), None, add_contrast=True)
