import numpy as np
import pytest
from pydantic import ValidationError

from escape_lab.revision_models import BioConfig, BiologyPredictor, catalog


class FakeStore:
    def features(self, candidate, drug, family="chemistry"):
        return {"x": candidate["x"], "missing": float(candidate["x"] == 999)}


def test_catalog_is_twelve_distinct_bounded_configs():
    configs = catalog()
    assert len(configs) == 12
    assert len({tuple(c.model_dump().items()) for c in configs.values()}) == 12
    assert "chemistry-logistic-c0.1" in configs


@pytest.mark.parametrize("kwargs", [{"C": True}, {"C": "0.1"}, {"C": 10},
                                   {"max_depth": True}, {"unknown": 1}])
def test_invalid_config_is_rejected(kwargs):
    with pytest.raises(ValidationError):
        BioConfig(**kwargs)


def test_imputer_and_scaler_fit_only_training_and_missing_is_not_zero():
    train = [{"x": x, "gene": "FGFR1"} for x in (1., 3., 999., 5., 7., 9.)]
    model = BiologyPredictor(BioConfig(), FakeStore(), "pem").fit(
        train, [False, False, False, True, True, True])
    means = model.scaler.mean_.copy()
    statistics = model.imputer.statistics_.copy()
    probability = model.predict([{"x": 100000., "gene": "FGFR1"},
                                 {"x": 999., "gene": "FGFR1"}])
    np.testing.assert_equal(model.scaler.mean_, means)
    np.testing.assert_equal(model.imputer.statistics_, statistics)
    assert statistics[model.vectorizer.feature_names_.index("x")] == 5
    assert np.isfinite(probability).all()
    assert ((probability >= 0) & (probability <= 1)).all()


def test_prior_fallback_is_explicit_and_predict_requires_fit():
    model = BiologyPredictor(BioConfig(), FakeStore(), "fut")
    with pytest.raises(RuntimeError):
        model.predict([{"x": 1}])
    model.fit([{"x": 1, "gene": "FGFR1"}, {"x": 2, "gene": "FGFR1"}], [False, False])
    assert model.fallback_reason == "single_class_training"
    np.testing.assert_equal(model.predict([{"x": 3, "gene": "FGFR1"}]), [.25])
    with pytest.raises(ValueError, match="gene"):
        model.predict([{"x": 3}])


def test_tree_runs_and_predictions_ignore_candidate_outcome_fields():
    train = [{"x": float(i), "gene": "FGFR1"} for i in range(80)]
    model = BiologyPredictor(BioConfig(model="hist_gradient_boosting"),
                             FakeStore(), "pem").fit(train, [i > 40 for i in range(80)])
    first = model.predict([{"x": 30., "gene": "FGFR1"}])
    second = model.predict([{"x": 30., "gene": "FGFR1", "label": True, "split": "final"}])
    np.testing.assert_equal(first, second)


def test_minus_one_only_means_missing_in_declared_structure_columns():
    class MissingStore:
        def features(self, candidate, drug, family):
            return {"structure_sasa": -1, "chemistry_charge_delta": -1}

    model = BiologyPredictor(BioConfig(), MissingStore(), "pem")
    features = model._features([{"gene": "FGFR2"}])[0]
    assert np.isnan(features["structure_sasa"])
    assert features["chemistry_charge_delta"] == -1
    assert features["gene_FGFR2"] == 1
    assert features["gene_FGFR3"] == 0
