import json
import time

import pytest

from escape_lab.data import sha256
from escape_lab.improvement import (
    _wallclock_limit,
    accepts_revision,
    aggregate_metrics,
    load_labels,
    reserve_final,
    sample_configurations,
    verify_freeze,
)
from escape_lab.prediction import PredictionConfig


def metrics(pem, fut):
    return {"macro_ap": (pem + fut) / 2, "per_drug": {
        "pem": {"average_precision": pem}, "fut": {"average_precision": fut}}}


def test_promotion_requires_gain_and_each_drug_guardrail():
    baseline = metrics(.4, .5)
    assert accepts_revision(metrics(.5, .49), baseline)
    assert not accepts_revision(metrics(.6, .47), baseline)
    assert not accepts_revision(metrics(.39, .49), baseline)
    assert not accepts_revision(baseline, baseline)


def test_metrics_score_each_drug_and_keep_missing_labels_out_of_scores():
    rows = [
        {"drug": "pem", "label": True, "probability": .9},
        {"drug": "pem", "label": False, "probability": .1},
        {"drug": "pem", "label": None, "probability": .99},
        {"drug": "fut", "label": True, "probability": .2},
        {"drug": "fut", "label": False, "probability": .8},
    ]
    score = aggregate_metrics(rows)
    assert score["macro_ap"] == pytest.approx(.75)
    assert score["per_drug"]["pem"]["n"] == 2
    assert score["per_drug"]["pem"]["inconclusive_count"] == 1
    assert score["per_drug"]["pem"]["brier"] == pytest.approx(.01)


def test_random_search_is_reproducible_and_excludes_semantic_duplicates():
    initial = PredictionConfig()
    first = sample_configurations(2, seed=7, excluded={initial.config_id})
    second = sample_configurations(2, seed=7, excluded={initial.config_id})
    assert first == second
    assert len({c.config_id for c in first}) == 2
    assert initial.config_id not in {c.config_id for c in first}


def test_final_labels_refused_before_file_access(tmp_path):
    with pytest.raises(ValueError, match="sealed"):
        load_labels(tmp_path, [{"candidate_id": "x", "split": "final"}])


def test_labels_collapse_each_drug_independently(tmp_path):
    outcome = {"x": {"records": [
        {"pem": {"resistant": True}, "fut": {"resistant": False}},
        {"pem": {"resistant": True}, "fut": {"resistant": None}},
    ]}}
    (tmp_path / "outcomes.json").write_text(json.dumps(outcome))
    labels = load_labels(tmp_path, [{"candidate_id": "x", "split": "development"}])
    assert labels == {"x": {"pem": True, "fut": None}}


def test_final_reservation_is_atomic_and_freeze_detects_tampering(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    source = tmp_path / "source.py"
    source.write_text("original")
    artifact = tmp_path / "model.pkl"
    artifact.write_bytes(b"model")
    freeze = tmp_path / "freeze.json"
    freeze.write_text(json.dumps({"fingerprints": {str(source): sha256(source)},
                                 "artifacts": {str(artifact): sha256(artifact)}}))
    verify_freeze(freeze)
    receipt = reserve_final(data, freeze, "test-run")
    assert receipt["freeze_sha256"] == sha256(freeze)
    with pytest.raises(FileExistsError):
        reserve_final(data, freeze, "another-run")
    source.write_text("modified")
    with pytest.raises(ValueError, match="changed"):
        verify_freeze(freeze)


def test_final_read_requires_matching_reserved_freeze(tmp_path):
    (tmp_path / "final-evaluation.json").write_text(json.dumps({
        "freeze_sha256": "a" * 64, "status": "started"}))
    with pytest.raises(ValueError, match="authorization"):
        load_labels(tmp_path, [{"candidate_id": "x", "split": "final"}],
                    final_token="b" * 64)


def test_wallclock_limit_interrupts_work_and_restores_alarm():
    import signal

    previous = signal.getsignal(signal.SIGALRM)
    with pytest.raises(TimeoutError, match="wall-clock"), _wallclock_limit(.02):
        time.sleep(.2)
    assert signal.getsignal(signal.SIGALRM) == previous
    assert signal.getitimer(signal.ITIMER_REAL)[0] == 0
