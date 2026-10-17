import copy
import math

import numpy as np
import pytest

from escape_lab.fgfr2_sequence_models import FGFR2SequenceModel
from escape_lab.fgfr2_state_models import ARMS, PRIMARY, FGFR2StateModel


class BiologyStore:
    def features(self, candidate, drug, family):
        return {"position": candidate["position"], "geometry_missing": np.nan}


class Augmenter:
    def features(self, candidate):
        return {"context": candidate["context"], "charge": candidate["position"] % 3}


class SequenceStore:
    def __init__(self, rows):
        rng = np.random.default_rng(101)
        self.values = {c["candidate_id"]: rng.normal(size=960) for c in rows}
        self.requests = []

    def matrix(self, ids, representation="wt_delta"):
        assert representation == "wt_delta"
        self.requests.append(tuple(ids))
        return np.asarray([self.values[cid] for cid in ids])


class StateStore:
    def __init__(self, rows):
        self.requests = []
        self.values = {c["candidate_id"]: {
            "state_contacts": float(c["position"] % 5),
            "state_distance": float(c["position"] % 7),
        } for c in rows}

    def features(self, candidate, arm):
        self.requests.append((candidate["candidate_id"], arm))
        scale = {"inactive": 1, "active": 2, "contrast": -1}[arm]
        return {key: value * scale for key, value in self.values[candidate["candidate_id"]].items()}


def candidates(n=24, prefix="train"):
    return [{"candidate_id": f"{prefix}-{i}", "gene": "FGFR2",
             "group_id": f"{prefix}-{i // 8}", "split": "development",
             "position": i + 500, "context": "A" if i % 2 else "G"} for i in range(n)]


def records(rows):
    values = {}
    for i, row in enumerate(rows):
        value = 1.0 + (i % 8)
        values[row["candidate_id"]] = {
            "replicates": [value] * 4,
            "value": math.log2(1 + value),
            "active": value >= 1.5,
        }
    return values


def make_model(arm, rows):
    return FGFR2StateModel(arm, BiologyStore(), Augmenter(), SequenceStore(rows), StateStore(rows))


def test_arms_and_primary_are_fixed():
    assert ARMS == ("baseline", "inactive", "active", "contrast")
    assert PRIMARY == "contrast"


def test_baseline_exactly_reproduces_previous_sequence_regressor_without_state_access():
    train, heldout = candidates(), candidates(5, "heldout")
    sequence = SequenceStore(train + heldout)
    previous = FGFR2SequenceModel("esm_biology_ridge", BiologyStore(), Augmenter(), sequence)
    current = FGFR2StateModel("baseline", BiologyStore(), Augmenter(), sequence, None)
    previous.fit(train, records(train))
    current.fit(train, records(train))
    np.testing.assert_array_equal(previous.predict(heldout), current.predict(heldout))
    np.testing.assert_array_equal(previous.model.coef_, current.model.coef_)
    np.testing.assert_array_equal(previous.scaler.mean_, current.scaler.mean_)
    assert current.diagnostics["arm"] == "baseline"
    assert current.diagnostics["baseline_arm"] == "esm_biology_ridge"
    assert current.diagnostics["state_feature_count"] == 0


def test_all_structural_arms_have_identical_columns_and_fixed_head():
    train, heldout = candidates(), candidates(5, "heldout")
    models = [make_model(arm, train + heldout).fit(train, records(train)) for arm in ARMS]
    reference = models[0]
    assert reference.state_store.requests == []
    for model in models[1:]:
        assert model.model.alpha == 100.0
        assert model.model.solver == "svd"
        assert model.diagnostics["state_feature_count"] == 2
        assert model.diagnostics["feature_count"] == reference.diagnostics["feature_count"] + 2
        assert model.diagnostics["arm"] == model.state_arm
        assert model.diagnostics["sequence_components"] == {"wt": 16, "delta": 16}
        assert np.isfinite(model.predict(heldout)).all()
        np.testing.assert_array_equal(model.vectorizer.get_feature_names_out(),
                                      models[1].vectorizer.get_feature_names_out())


def test_heldout_features_cannot_change_any_fitted_transform_or_coefficients():
    train, heldout = candidates(), candidates(5, "heldout")
    first, second = [make_model("contrast", train + heldout) for _ in range(2)]
    for row in heldout:
        cid = row["candidate_id"]
        second.sequence_store.values[cid] *= 1e9
        second.state_store.values[cid] = {"state_contacts": 1e9, "state_distance": -1e9}
    first.fit(train, records(train))
    second.fit(train, records(train))
    for attribute in ("scaler", "wt_scaler", "delta_scaler"):
        np.testing.assert_array_equal(getattr(first, attribute).mean_,
                                      getattr(second, attribute).mean_)
    np.testing.assert_array_equal(first.imputer.statistics_, second.imputer.statistics_)
    np.testing.assert_array_equal(first.wt_pca.components_, second.wt_pca.components_)
    np.testing.assert_array_equal(first.delta_pca.components_, second.delta_pca.components_)
    np.testing.assert_array_equal(first.model.coef_, second.model.coef_)
    assert first.sequence_store.requests == [tuple(c["candidate_id"] for c in train)]
    assert first.state_store.requests == [(c["candidate_id"], "contrast") for c in train]


def test_missing_geometry_is_imputed_from_training_only():
    train, heldout = candidates(), candidates(5, "heldout")
    model = make_model("active", train + heldout)
    model.state_store.values[train[0]["candidate_id"]]["state_distance"] = np.nan
    model.state_store.values[heldout[0]["candidate_id"]]["state_distance"] = np.nan
    model.fit(train, records(train))
    assert np.isfinite(model.predict(heldout)).all()


@pytest.mark.parametrize("arm", ["unknown", "esm_biology_ridge", "cached_v2_biology"])
def test_unknown_arms_rejected(arm):
    with pytest.raises(ValueError, match="Unknown"):
        make_model(arm, candidates())


def test_state_store_required_for_structural_arms():
    with pytest.raises(ValueError, match="State features"):
        FGFR2StateModel("active", BiologyStore(), Augmenter(), SequenceStore(candidates()), None)


@pytest.mark.parametrize("bad", [None, {}, {"position": 1.0}, {"state_distance": "categorical"},
                                  {"state_distance": np.inf},
                                  {f"state_{i}": 1.0 for i in range(33)}])
def test_invalid_or_conflicting_state_features_rejected(bad):
    train = candidates()
    model = make_model("contrast", train)

    def bad_features(candidate, arm):
        return bad

    model.state_store.features = bad_features
    with pytest.raises(ValueError, match="State feature"):
        model.fit(train, records(train))


def test_state_schema_cannot_change_between_candidates_or_at_prediction():
    train, heldout = candidates(), candidates(5, "heldout")
    model = make_model("contrast", train + heldout)
    model.state_store.values[train[1]["candidate_id"]]["unexpected"] = 1.0
    with pytest.raises(ValueError, match="schema"):
        model.fit(train, records(train))
    model = make_model("contrast", train + heldout).fit(train, records(train))
    model.state_store.values[heldout[0]["candidate_id"]]["unexpected"] = 1.0
    with pytest.raises(ValueError, match="schema"):
        model.predict(heldout)


def test_parent_partition_and_exact_label_guards_remain_in_force():
    train, heldout = candidates(), candidates(5, "heldout")
    with pytest.raises(ValueError, match="exactly"):
        make_model("contrast", train + heldout).fit(train, records(train + heldout))
    model = make_model("contrast", train + heldout).fit(train, records(train))
    heldout[0]["group_id"] = train[0]["group_id"]
    with pytest.raises(ValueError, match="training"):
        model.predict(heldout)


@pytest.mark.parametrize("change", [{"gene": "FGFR1"}, {"split": "final"}])
def test_parent_development_scope_guard_remains_in_force(change):
    train = candidates()
    train[0].update(change)
    with pytest.raises(ValueError, match="FGFR2 development"):
        make_model("contrast", train).fit(train, records(train))


def test_features_and_labels_remain_unchanged():
    train = candidates()
    labels = records(train)
    model = make_model("contrast", train)
    before = copy.deepcopy((train, labels, model.state_store.values))
    model.fit(train, labels)
    assert (train, labels, model.state_store.values) == before
