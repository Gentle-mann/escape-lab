import copy
import math
import warnings

import numpy as np
import pytest
from sklearn.exceptions import ConvergenceWarning

from escape_lab.fgfr2_sequence_models import (
    TRAINABLE_ARMS,
    FGFR2SequenceModel,
    activation_label,
)


class Store:
    def features(self, candidate, drug, family):
        return {
            "position": candidate["position"],
            "structure_missing": 999,
            "structure_ligand_distance": 1234,
        }


class Augmenter:
    def features(self, candidate):
        return {"context": candidate["context"], "property": candidate["position"] % 3}


class SequenceStore:
    def __init__(self, candidates):
        rng = np.random.default_rng(101)
        self.values = {c["candidate_id"]: rng.normal(size=960) for c in candidates}
        self.requests = []

    def matrix(self, ids, representation="wt_delta"):
        assert representation == "wt_delta"
        self.requests.append(tuple(ids))
        return np.asarray([self.values[cid] for cid in ids])


def candidates(n=24, group_prefix="training"):
    return [{"candidate_id": f"{group_prefix}-{i}", "gene": "FGFR2",
             "group_id": f"{group_prefix}-{i // 8}", "split": "development",
             "position": i + 500, "context": "A" if i % 2 else "G"} for i in range(n)]


def records(rows):
    result = {}
    for i, candidate in enumerate(rows):
        enrichment = 8.0 if i % 4 == 0 else 1.0
        result[candidate["candidate_id"]] = {
            "replicates": [enrichment] * 4,
            "value": math.log2(1 + enrichment),
            "active": enrichment >= 1.5,
        }
    return result


def make_model(arm, rows):
    return FGFR2SequenceModel(arm, Store(), Augmenter(), SequenceStore(rows))


def test_activation_thresholds_are_strict_and_require_replicate_consistency():
    def record(reps):
        return {"replicates": reps, "value": math.log2(1 + np.median(reps)),
                "active": sum(x >= 1.5 for x in reps) >= 3}

    assert not activation_label(record([5, 5, 5, 5]))
    assert activation_label(record([6, 6, 6, 6]))
    assert not activation_label(record([0, 0, 20, 20]))
    assert not activation_label(record([10, 10, 10, 10]), threshold=10)
    assert activation_label(record([11, 11, 11, 11]), threshold=10)
    with pytest.raises(ValueError, match="target"):
        activation_label({"replicates": [8] * 4, "value": 0, "active": True})


@pytest.mark.parametrize("arm", TRAINABLE_ARMS)
def test_fit_and_predict_have_exact_partition_and_finite_output(arm):
    train, heldout = candidates(), candidates(5, "heldout")
    model = make_model(arm, train + heldout).fit(train, records(train))
    output = model.predict(heldout)
    assert output.shape == (5,)
    assert np.isfinite(output).all()
    assert model.diagnostics["training_n"] == 24
    assert model.diagnostics["high_count"] == 6
    assert "structure_ligand_distance" not in model.vectorizer.get_feature_names_out()
    if "logistic" in arm:
        assert ((output >= 0) & (output <= 1)).all()
    if arm.startswith("esm"):
        assert model.diagnostics["sequence_components"] == {"wt": 16, "delta": 16}


def test_heldout_poison_does_not_change_training_preprocessors_or_fit():
    train, heldout = candidates(), candidates(5, "heldout")
    first = make_model("esm_biology_logistic", train + heldout)
    second = make_model("esm_biology_logistic", train + heldout)
    for candidate in heldout:
        second.sequence_store.values[candidate["candidate_id"]] *= 1e8
    first.fit(train, records(train))
    second.fit(train, records(train))
    np.testing.assert_array_equal(first.scaler.mean_, second.scaler.mean_)
    np.testing.assert_array_equal(first.wt_scaler.mean_, second.wt_scaler.mean_)
    np.testing.assert_array_equal(first.delta_scaler.mean_, second.delta_scaler.mean_)
    np.testing.assert_array_equal(first.wt_pca.components_, second.wt_pca.components_)
    np.testing.assert_array_equal(first.model.coef_, second.model.coef_)
    assert first.sequence_store.requests == [tuple(c["candidate_id"] for c in train)]


def test_labels_must_exactly_match_training_scope():
    train, heldout = candidates(), candidates(5, "heldout")
    with pytest.raises(ValueError, match="exactly"):
        make_model("biology_logistic", train).fit(train, records(train + heldout))
    with pytest.raises(ValueError, match="Duplicate"):
        make_model("biology_logistic", train).fit(train + [train[0]], records(train))


@pytest.mark.parametrize("mutation", [{"gene": "FGFR1"}, {"split": "final"}])
def test_unapproved_candidates_rejected(mutation):
    train = candidates()
    train[0].update(mutation)
    with pytest.raises(ValueError, match="FGFR2 development"):
        make_model("biology_logistic", train).fit(train, records(train))


def test_predictions_reject_shared_groups_and_training_identity():
    train = candidates()
    model = make_model("biology_logistic", train).fit(train, records(train))
    with pytest.raises(ValueError, match="training"):
        model.predict([train[0]])
    heldout = candidates(1, "heldout")
    heldout[0]["group_id"] = train[0]["group_id"]
    with pytest.raises(ValueError, match="training"):
        model.predict(heldout)


def test_no_single_class_classifier_or_inconsistent_activation():
    train = candidates()
    data = records(train)
    for record in data.values():
        record.update(replicates=[8] * 4, value=math.log2(9), active=True)
    with pytest.raises(ValueError, match="two classes"):
        make_model("biology_logistic", train).fit(train, data)
    data = records(train)
    data[train[0]["candidate_id"]]["active"] = False
    with pytest.raises(ValueError, match="consistency"):
        make_model("biology_logistic", train).fit(train, data)


def test_bad_sequence_values_and_dimensions_rejected():
    train = candidates()
    model = make_model("esm_biology_logistic", train)
    model.sequence_store.values[train[0]["candidate_id"]][0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        model.fit(train, records(train))
    model = make_model("esm_biology_logistic", train)
    model.sequence_store.values = {k: v[:12] for k, v in model.sequence_store.values.items()}
    with pytest.raises(ValueError, match="960"):
        model.fit(train, records(train))


def test_nonconvergence_is_fatal(monkeypatch):
    train = candidates()
    model = make_model("biology_logistic", train)

    def warning_fit(*args, **kwargs):
        warnings.warn("did not converge", ConvergenceWarning)

    monkeypatch.setattr(model.model, "fit", warning_fit)
    with pytest.raises(ValueError, match="converge"):
        model.fit(train, records(train))
    with pytest.raises(ValueError, match="fitted"):
        model.predict(candidates(2, "heldout"))


def test_features_and_labels_are_not_mutated():
    train = candidates()
    data, original = records(train), copy.deepcopy(train)
    previous = copy.deepcopy(data)
    make_model("biology_logistic", train).fit(train, data)
    assert train == original
    assert data == previous


def test_cached_arm_is_evaluator_owned():
    with pytest.raises(ValueError, match="cached"):
        make_model("cached_v2_biology", candidates())
