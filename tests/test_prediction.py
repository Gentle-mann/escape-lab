import pickle

import numpy as np
import pytest
from pydantic import ValidationError

from escape_lab.prediction import (
    PredictionConfig,
    ResistancePredictor,
    candidate_features,
    fit_predict,
)


def candidate(index, gene="FGFR1", mutant="K"):
    position = 450 + index
    return {
        "candidate_id": f"{gene}:A{position}{mutant}",
        "gene": gene,
        "wt": "A",
        "mut": mutant,
        "aligned_position": position,
        "context": "AAAAGAA",
        "group_id": f"group-{index // 12}",
        "split": "development",
        "position": position,
    }


def training_data():
    candidates = [candidate(index, mutant="K" if index % 2 else "E") for index in range(40)]
    return candidates, [bool(index % 2) for index in range(40)]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"model": "unsupported"},
        {"feature_set": "labels"},
        {"include_gene": "true"},
        {"position_scale": 21},
        {"regularization": "1.0"},
        {"regularization": True},
        {"min_samples_leaf": True},
        {"max_features": False},
        {"max_depth": 9},
        {"learning_rate": 0.2},
        {"l2": -1},
        {"unknown_parameter": True},
    ],
)
def test_configuration_rejects_unsupported_coerced_or_unknown_settings(kwargs):
    with pytest.raises(ValidationError):
        PredictionConfig(**kwargs)


def test_configuration_is_frozen_json_serializable_and_has_a_schema():
    config = PredictionConfig()
    assert PredictionConfig.model_validate_json(config.model_dump_json()) == config
    assert "model" in PredictionConfig.model_json_schema()["properties"]
    with pytest.raises(ValidationError):
        config.include_gene = False
    assert len(config.config_id) == 64
    assert int(config.config_id, 16) >= 0


def test_semantic_ids_ignore_inactive_model_and_feature_parameters():
    extra = PredictionConfig(model="extra_trees", feature_set="sequence")
    equivalent = PredictionConfig(
        model="extra_trees",
        feature_set="sequence",
        regularization=10.0,
        learning_rate=0.05,
        l2=10.0,
        position_scale=40,
    )
    assert extra.effective_config() == equivalent.effective_config()
    assert extra.config_id == equivalent.config_id
    assert "position_scale" not in extra.effective_config()
    assert "regularization" not in extra.effective_config()
    assert extra.effective_config()["n_estimators"] == 100
    assert (
        PredictionConfig(model="logistic", regularization=1.0).config_id
        != PredictionConfig(model="logistic", regularization=10.0).config_id
    )
    assert PredictionConfig(model="hist_gradient_boosting").effective_config()["max_iter"] == 100
    assert (
        PredictionConfig(include_gene=True).config_id
        != PredictionConfig(include_gene=False).config_id
    )


@pytest.mark.parametrize("feature_set", ["sequence", "position", "full"])
def test_features_ignore_outcome_fields_identifiers_and_split(feature_set):
    config = PredictionConfig(feature_set=feature_set)
    clean = candidate(3)
    poisoned = {
        **clean,
        "candidate_id": "answer-is-positive",
        "score": 999,
        "dual_resistant": True,
        "group_id": "positive-only",
        "split": "final",
        "evidence_ids": ["secret"],
        "pem_median": 999,
        "position": -9000,
    }
    assert candidate_features(clean, config) == candidate_features(poisoned, config)
    assert "gene" not in candidate_features(clean, PredictionConfig(include_gene=False))


def test_feature_sets_respect_sequence_position_and_gene_switches():
    original = candidate(1)
    different_position = {**original, "aligned_position": 780}
    different_context = {**original, "context": "WWWWWWW"}
    sequence = PredictionConfig(feature_set="sequence")
    position = PredictionConfig(feature_set="position")
    full = PredictionConfig(feature_set="full")
    assert candidate_features(original, sequence) == candidate_features(
        different_position, sequence
    )
    assert candidate_features(original, sequence) != candidate_features(different_context, sequence)
    assert candidate_features(original, position) == candidate_features(different_context, position)
    assert candidate_features(original, position) != candidate_features(
        different_position, position
    )
    assert candidate_features(original, full) != candidate_features(different_context, full)
    assert candidate_features(original, full) != candidate_features(different_position, full)
    without_gene = PredictionConfig(include_gene=False)
    assert candidate_features(original, without_gene) == candidate_features(
        {**original, "gene": "FGFR4"}, without_gene
    )


@pytest.mark.parametrize("model", ["logistic", "extra_trees", "hist_gradient_boosting"])
def test_models_produce_deterministic_finite_probabilities_and_pickle_roundtrip(model):
    training, labels = training_data()
    pool = [candidate(50 + index, mutant="K" if index % 2 else "E") for index in range(8)]
    config = PredictionConfig(model=model)
    first = ResistancePredictor(config, seed=7).fit(training, labels)
    probabilities = first.predict(pool)
    second = fit_predict(config, training, labels, pool, seed=7)
    assert probabilities.shape == (8,)
    assert np.isfinite(probabilities).all()
    assert ((0 <= probabilities) & (probabilities <= 1)).all()
    np.testing.assert_allclose(probabilities, second)
    np.testing.assert_allclose(probabilities, pickle.loads(pickle.dumps(first)).predict(pool))
    assert first.fallback_reason is None
    assert first.using_prior_fallback is False
    if model == "extra_trees":
        assert len(first.model.estimators_) == 100
    elif model == "hist_gradient_boosting":
        assert first.model.n_iter_ == 100


def test_vectorizer_fits_training_only_and_predict_never_adds_new_features():
    training, labels = training_data()
    predictor = ResistancePredictor(PredictionConfig(), seed=0).fit(training, labels)
    vocabulary = dict(predictor.vectorizer.vocabulary_)
    pool = [{**candidate(85, gene="FGFR4", mutant="W"), "context": "WWWWWWW"}]
    predictor.predict(pool)
    assert predictor.vectorizer.vocabulary_ == vocabulary
    assert "gene=FGFR4" not in vocabulary
    assert "context_0=W" not in vocabulary
    assert "mut=W" not in vocabulary


def test_empty_single_class_and_refit_reset_are_explicit_laplace_fallbacks():
    predictor = ResistancePredictor(PredictionConfig(), seed=0)
    with pytest.raises(RuntimeError, match="fit"):
        predictor.predict([candidate(1)])
    predictor.fit([], [])
    assert predictor.fallback_reason == "empty_training"
    assert predictor.using_prior_fallback is True
    np.testing.assert_allclose(predictor.predict([candidate(1)]), [0.5])
    predictor.fit([candidate(1), candidate(2)], [False, False])
    assert predictor.fallback_reason == "single_class_training"
    np.testing.assert_allclose(predictor.predict([candidate(5)]), [0.25])
    predictor.fit([candidate(1), candidate(2)], [True, True])
    np.testing.assert_allclose(predictor.predict([candidate(5)]), [0.75])
    training, labels = training_data()
    predictor.fit(training, labels)
    assert predictor.using_prior_fallback is False
    predictor.fit([], [])
    assert predictor.model is None
    assert not hasattr(predictor.vectorizer, "vocabulary_")
    assert predictor.predict([]).shape == (0,)


@pytest.mark.parametrize("labels", [[1, 0], [True, None], [True]])
def test_training_requires_matching_defined_boolean_labels(labels):
    with pytest.raises(ValueError, match="boolean"):
        ResistancePredictor(PredictionConfig()).fit([candidate(1), candidate(2)], labels)
