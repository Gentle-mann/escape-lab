"""A bounded sensitivity changes only out-of-range prediction biology features."""

import math

import numpy as np
import pytest

from escape_lab.fgfr2_state_bounded_model import BoundedFGFR2StateModel
from escape_lab.fgfr2_state_models import FGFR2StateModel


class Biology:
    def features(self, candidate, drug, family):
        return {"near_constant": candidate["near_constant"]}


class Augmenter:
    def features(self, candidate):
        return {"other": candidate["other"]}


class Sequence:
    def __init__(self, rows):
        rng = np.random.default_rng(103)
        self.values = {row["candidate_id"]: rng.normal(size=960) for row in rows}

    def matrix(self, ids, representation):
        assert representation == "wt_delta"
        return np.asarray([self.values[cid] for cid in ids])


class States:
    def features(self, candidate, arm):
        return {"state_fixed": candidate["other"] / 2}


@pytest.fixture
def experiment():
    train = [{"candidate_id": f"train-{i}", "gene": "FGFR2", "split": "development",
              "group_id": f"train-group-{i // 8}", "near_constant": 1e-8 + i * 1e-10,
              "other": float(i % 3)} for i in range(24)]
    test = [{**train[i], "candidate_id": f"test-{i}", "group_id": "test-group"}
            for i in range(3)]
    labels = {}
    for i, row in enumerate(train):
        value = 1 + i / 2
        labels[row["candidate_id"]] = {"value": math.log2(1 + value),
                                        "replicates": [value] * 4, "active": value >= 1.5}
    sequence = Sequence(train + test)
    return train, test, labels, sequence


def models(experiment, arm="baseline"):
    train, _, labels, sequence = experiment
    previous = FGFR2StateModel(arm, Biology(), Augmenter(), sequence, States())
    bounded = BoundedFGFR2StateModel(arm, Biology(), Augmenter(), sequence, States())
    previous.fit(train, labels)
    bounded.fit(train, labels)
    return previous, bounded


@pytest.mark.parametrize("arm", ["baseline", "inactive", "active", "contrast"])
def test_training_head_and_transforms_remain_identical(experiment, arm):
    previous, bounded = models(experiment, arm)
    np.testing.assert_array_equal(previous.model.coef_, bounded.model.coef_)
    assert previous.model.intercept_ == bounded.model.intercept_
    np.testing.assert_array_equal(previous.scaler.mean_, bounded.scaler.mean_)
    np.testing.assert_array_equal(previous.scaler.scale_, bounded.scaler.scale_)
    np.testing.assert_array_equal(previous.wt_pca.components_, bounded.wt_pca.components_)
    np.testing.assert_array_equal(previous.delta_pca.components_, bounded.delta_pca.components_)
    train = experiment[0]
    ids = [row["candidate_id"] for row in train]
    old_matrix = previous._predict_matrix(train, ids)
    new_matrix = bounded._predict_matrix(train, ids)
    np.testing.assert_allclose(old_matrix, new_matrix, atol=1e-14)
    assert bounded.model.alpha == 100


def test_near_constant_training_feature_cannot_amplify_out_of_range_test_value(experiment):
    previous, bounded = models(experiment)
    test = experiment[1]
    test[0]["near_constant"] = 1.0
    unbounded = previous.predict(test)
    corrected = bounded.predict(test)
    assert abs(unbounded[0]) > 1e5
    assert abs(corrected[0]) < 10
    np.testing.assert_array_equal(unbounded, bounded.unbounded_predictions)
    diagnostics = bounded.diagnostics["bounded_extrapolation"]
    assert diagnostics["clipped_value_count"] == 1
    assert diagnostics["clipped_row_count"] == 1
    assert diagnostics["clipped_row_fraction"] == pytest.approx(1 / 3)
    assert diagnostics["maximum_excess_beyond_training_bound"] > 1e8
    assert diagnostics["raw_max_abs_biology"] > 1e8
    assert diagnostics["bounded_max_abs_biology"] < 3


def test_in_range_predictions_are_unchanged(experiment):
    previous, bounded = models(experiment)
    test = experiment[1]
    np.testing.assert_array_equal(previous.predict(test), bounded.predict(test))
    assert bounded.diagnostics["bounded_extrapolation"]["clipped_value_count"] == 0


def test_bounds_never_use_or_update_from_test_features(experiment):
    _, bounded = models(experiment)
    minima, maxima = bounded.biology_minimum_.copy(), bounded.biology_maximum_.copy()
    test = experiment[1]
    test[0]["near_constant"] = 1e5
    bounded.predict(test)
    test[0]["near_constant"] = -1e5
    bounded.predict(test)
    np.testing.assert_array_equal(minima, bounded.biology_minimum_)
    np.testing.assert_array_equal(maxima, bounded.biology_maximum_)


def test_sequence_pca_columns_are_never_clipped(experiment):
    previous, bounded = models(experiment)
    test, sequence = experiment[1], experiment[3]
    test[0]["near_constant"] = 1e3
    sequence.values[test[0]["candidate_id"]] *= 1000
    ids = [row["candidate_id"] for row in test]
    old = previous._predict_matrix(test, ids)
    new = bounded._predict_matrix(test, ids)
    start = bounded.biology_feature_count_
    np.testing.assert_array_equal(old[:, start:], new[:, start:])
    assert np.abs(new[:, start:]).max() > 100


@pytest.mark.parametrize("change", [{"split": "final"}, {"group_id": "train-group-0"}])
def test_parent_group_and_development_guards_remain_in_force(experiment, change):
    _, bounded = models(experiment)
    test = experiment[1]
    test[0].update(change)
    with pytest.raises(ValueError):
        bounded.predict(test)


def test_predict_before_fit_still_fails(experiment):
    _, test, _, sequence = experiment
    bounded = BoundedFGFR2StateModel("baseline", Biology(), Augmenter(), sequence, None)
    with pytest.raises(ValueError, match="fitted"):
        bounded.predict(test)
