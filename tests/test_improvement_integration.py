"""Synthetic end-to-end boundary tests; no network, real data, or real final set."""

import json
import pickle
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import pytest

from escape_lab import improvement
from escape_lab.data import sha256
from escape_lab.improvement_agent import Proposal
from escape_lab.prediction import PredictionConfig, ResistancePredictor


def make_candidate(index, split, mutant="K"):
    position = 450 + index
    return {
        "candidate_id": f"FGFR1:A{position}{mutant}",
        "gene": "FGFR1",
        "wt": "A",
        "mut": mutant,
        "aligned_position": position,
        "position": position,
        "protein_variant": f"A{position}{mutant}",
        "context": "AAAAGAA",
        "group_id": f"group-{index // 2}",
        "split": split,
    }


@pytest.fixture
def synthetic_data(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    candidates = []
    for index in range(16):
        split = "development" if index < 8 else "validation" if index < 12 else "final"
        candidates.append(make_candidate(index, split, "K" if index % 2 else "E"))
    outcomes = {}
    for index, candidate in enumerate(candidates):
        outcomes[candidate["candidate_id"]] = {
            "records": [
                {
                    "pem": {"resistant": candidate["mut"] == "K"},
                    "fut": {"resistant": None if index == 12 else candidate["mut"] != "K"},
                }
            ]
        }
    (data / "candidates.json").write_text(json.dumps(candidates))
    (data / "outcomes.json").write_text(json.dumps(outcomes))
    (data / "manifest.json").write_text(
        json.dumps(
            {
                "candidates_sha256": sha256(data / "candidates.json"),
                "outcomes_sha256": sha256(data / "outcomes.json"),
            }
        )
    )
    return data, candidates


@pytest.fixture
def frozen_run(synthetic_data, tmp_path):
    data, candidates = synthetic_data
    output = tmp_path / "frozen-output"
    output.mkdir()
    training = [candidate for candidate in candidates if candidate["split"] == "development"]
    config = PredictionConfig(model="logistic", feature_set="sequence", include_gene=False)
    model = ResistancePredictor(config).fit(
        training, [candidate["mut"] == "K" for candidate in training]
    )
    artifacts, paths = {}, {}
    for drug in improvement.DRUGS:
        path = output / f"own-{drug}.pkl"
        with path.open("wb") as handle:
            pickle.dump(model, handle)
        paths[drug] = str(path.resolve())
        artifacts[str(path.resolve())] = sha256(path)
    source = tmp_path / "frozen-source.py"
    source.write_text("# Synthetic source snapshot\n")
    freeze_path = output / "freeze.json"
    freeze = {
        "run_id": "synthetic-run",
        "fingerprints": {
            str(source.resolve()): sha256(source),
            str((data / "candidates.json").resolve()): sha256(data / "candidates.json"),
            str((data / "outcomes.json").resolve()): sha256(data / "outcomes.json"),
        },
        "artifacts": artifacts,
        "selected": {
            method: {
                "config": config.model_dump(),
                "config_id": config.config_id,
                "model_paths": paths.copy(),
            }
            for method in ("fixed", "random", "agent")
        },
        "criteria": "Synthetic test criteria",
    }
    improvement.write_json(freeze_path, freeze)
    return data, candidates, output, freeze_path, source


def quick_bootstrap(records, seed):
    assert set(records) == {"fixed", "random", "agent"}
    return {"method": "synthetic bootstrap stub", "seed": seed}


def test_final_predictions_are_committed_before_labels_and_receipt_is_one_shot(
    frozen_run, monkeypatch
):
    data, candidates, output, freeze_path, _ = frozen_run
    original = improvement.load_labels
    reads = []

    def guarded_read(data_dir, requested, final_token=None):
        assert all(candidate["split"] == "final" for candidate in requested)
        receipt = json.loads((data / "final-evaluation.json").read_text())
        assert receipt["status"] == "started"
        assert receipt["freeze_sha256"] == final_token == sha256(freeze_path)
        committed = output / "final-predictions-committed.json"
        assert receipt["predictions_sha256"] == sha256(committed)
        predictions = json.loads(committed.read_text())
        expected_ids = {candidate["candidate_id"] for candidate in requested}
        assert set(predictions) == {"fixed", "random", "agent"}
        for rows in predictions.values():
            assert len(rows) == len(requested) * 2
            assert "label" not in json.dumps(rows)
            for drug in improvement.DRUGS:
                assert {row["candidate_id"] for row in rows if row["drug"] == drug} == expected_ids
        reads.append(final_token)
        return original(data_dir, requested, final_token=final_token)

    monkeypatch.setattr(improvement, "load_labels", guarded_read)
    monkeypatch.setattr(improvement, "_bootstrap", quick_bootstrap)
    result = improvement.evaluate_frozen_final(data, freeze_path, output, candidates, seed=3)
    assert len(reads) == 1
    assert result["status"] == "completed"
    for metrics in result["metrics"].values():
        assert metrics["per_drug"]["pem"]["n"] == 4
        assert metrics["per_drug"]["pem"]["positive_count"] == 2
        assert metrics["per_drug"]["pem"]["average_precision"] == pytest.approx(1.0)
        assert metrics["per_drug"]["fut"]["n"] == 3
        assert metrics["per_drug"]["fut"]["inconclusive_count"] == 1
        assert metrics["per_drug"]["fut"]["positive_count"] == 1
        assert metrics["per_drug"]["fut"]["average_precision"] < 1.0
    assert result["prototype_margin_passed"] is False
    assert (output / "final-predictions.csv").exists()
    assert (output / "final-results.json").exists()
    assert json.loads((data / "final-evaluation.json").read_text())["status"] == "completed"
    with pytest.raises(FileExistsError):
        improvement.evaluate_frozen_final(data, freeze_path, output, candidates, seed=3)
    assert len(reads) == 1


def test_two_concurrent_reservations_allow_exactly_one_owner(frozen_run):
    data, _, _, freeze_path, _ = frozen_run

    def reserve(index):
        try:
            return improvement.reserve_final(data, freeze_path, f"contender-{index}")
        except FileExistsError:
            return None

    with ThreadPoolExecutor(max_workers=2) as workers:
        results = list(workers.map(reserve, (0, 1)))
    successful = [result for result in results if result is not None]
    assert len(successful) == 1
    assert json.loads((data / "final-evaluation.json").read_text()) == successful[0]


def test_final_failure_after_reservation_is_recorded_and_rerun_is_blocked(frozen_run, monkeypatch):
    data, candidates, output, freeze_path, _ = frozen_run
    reads = []

    def interrupted_read(*args, **kwargs):
        assert (output / "final-predictions-committed.json").exists()
        reads.append("attempted")
        raise RuntimeError("Synthetic interruption after final reservation")

    monkeypatch.setattr(improvement, "load_labels", interrupted_read)
    with pytest.raises(RuntimeError, match="Synthetic interruption"):
        improvement.evaluate_frozen_final(data, freeze_path, output, candidates, seed=0)
    receipt = json.loads((data / "final-evaluation.json").read_text())
    assert receipt["status"] == "failed_after_reservation"
    assert receipt["predictions_sha256"] == sha256(output / "final-predictions-committed.json")
    with pytest.raises(FileExistsError):
        improvement.evaluate_frozen_final(data, freeze_path, output, candidates, seed=0)
    assert len(reads) == 1


def test_source_tampering_is_rejected_before_reserving_or_reading_final(frozen_run, monkeypatch):
    data, candidates, output, freeze_path, source = frozen_run
    source.write_text("# Changed source\n")

    def forbidden_read(*args, **kwargs):
        pytest.fail("Final labels must not be read after source tampering")

    monkeypatch.setattr(improvement, "load_labels", forbidden_read)
    with pytest.raises(ValueError, match="changed"):
        improvement.evaluate_frozen_final(data, freeze_path, output, candidates, seed=0)
    assert not (data / "final-evaluation.json").exists()


def test_freeze_manifest_tampering_after_reservation_prevents_label_read(frozen_run, monkeypatch):
    data, candidates, output, freeze_path, _ = frozen_run
    original_write = improvement.write_json

    def tampering_write(path, value):
        original_write(path, value)
        if Path(path).name == "final-predictions-committed.json":
            altered = json.loads(freeze_path.read_text())
            altered["criteria"] = "Illegitimate changed criteria"
            original_write(freeze_path, altered)

    def forbidden_read(*args, **kwargs):
        pytest.fail("The frozen manifest changed after reservation; do not read labels")

    monkeypatch.setattr(improvement, "write_json", tampering_write)
    monkeypatch.setattr(improvement, "load_labels", forbidden_read)
    with pytest.raises(ValueError, match="[Ff]roz|[Ff]reez|changed"):
        improvement.evaluate_frozen_final(data, freeze_path, output, candidates, seed=0)
    assert (
        json.loads((data / "final-evaluation.json").read_text())["status"]
        == "failed_after_reservation"
    )


def synthetic_metrics(score):
    return {
        "macro_ap": score,
        "per_drug": {
            drug: {
                "average_precision": score,
                "brier": 0.2,
                "n": 8,
                "positive_count": 4,
                "inconclusive_count": 0,
            }
            for drug in improvement.DRUGS
        },
    }


@pytest.mark.parametrize("validation_mode", ["no_gain", "evaluation_failed"])
def test_two_autonomous_revisions_never_receive_validation_feedback_and_failure_keeps_final_sealed(
    synthetic_data, tmp_path, monkeypatch, validation_mode
):
    data, _ = synthetic_data
    agent_configs = [
        PredictionConfig(model="logistic", feature_set="sequence", regularization=0.1),
        PredictionConfig(feature_set="sequence", min_samples_leaf=10, include_gene=False),
    ]
    random_configs = [
        PredictionConfig(model="logistic", feature_set="position", regularization=10.0),
        PredictionConfig(model="hist_gradient_boosting", feature_set="sequence", l2=10.0),
    ]
    evaluations, feedbacks, reads = [], [], []
    original_labels = improvement.load_labels

    def selected_random_configs(count, seed, excluded):
        assert count == 2
        assert not {config.config_id for config in random_configs} & set(excluded)
        return random_configs

    def guarded_labels(data_dir, requested, final_token=None):
        splits = {candidate["split"] for candidate in requested}
        assert "final" not in splits
        reads.append(splits)
        return original_labels(data_dir, requested, final_token=final_token)

    def mock_evaluation(config, tasks, labels, output_dir, seed, seconds):
        stage = "validation" if tasks[0]["task_id"].startswith("validation") else "development"
        evaluations.append((stage, config.config_id))
        if stage == "validation":
            assert len(feedbacks) == 2
            if validation_mode == "evaluation_failed":
                return {
                    "status": "failed",
                    "config": config.model_dump(),
                    "config_id": config.config_id,
                    "seconds": 0.001,
                }
            score = (
                0.5 if config.config_id not in {item.config_id for item in agent_configs} else 0.4
            )
        elif config.config_id == agent_configs[0].config_id:
            score = 0.7
        elif config.config_id == agent_configs[1].config_id:
            score = 0.8
        else:
            score = 0.5
        artifacts = {}
        for task in tasks:
            for drug in improvement.DRUGS:
                path = Path(output_dir) / "models" / task["task_id"] / f"{drug}.pkl"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(pickle.dumps({"synthetic_fixture": True}))
                artifacts[str(path.resolve())] = sha256(path)
        return {
            "status": "completed",
            "config": config.model_dump(),
            "config_id": config.config_id,
            "seconds": 0.001,
            "metrics": synthetic_metrics(score),
            "artifacts": artifacts,
            "predictions": [{"label": True, "never_send_this_record": stage}],
            "error_analysis": {"by_fold": {tasks[0]["task_id"]: synthetic_metrics(score)}},
        }

    class FakeProposer:
        def __init__(self):
            self.client = SimpleNamespace(
                model_id="synthetic-no-network", max_tokens=100, usage_records=[]
            )
            self.attempts = []

        def propose(self, feedback):
            assert all(stage == "development" for stage, _ in evaluations)
            assert "validation" not in json.dumps(feedback)
            assert "never_send_this_record" not in json.dumps(feedback)
            feedbacks.append(json.loads(json.dumps(feedback)))
            config = agent_configs[len(feedbacks) - 1]
            return Proposal(
                config,
                "Synthetic distinct revision",
                "proposed",
                None,
                self.client.model_id,
                "synthetic-prompt-v1",
                1,
            )

    monkeypatch.setattr(improvement, "sample_configurations", selected_random_configs)
    monkeypatch.setattr(improvement, "bounded_evaluation", mock_evaluation)
    monkeypatch.setattr(improvement, "load_labels", guarded_labels)
    report = improvement.run_improvement(data, tmp_path / "runs", FakeProposer(), max_seconds=600)
    attempts = [
        trial
        for trial in report["trials"]
        if trial["method"] == "agent" and isinstance(trial["attempt"], int)
    ]
    assert len(feedbacks) == len(attempts) == 2
    assert [trial["config_id"] for trial in attempts] == [
        config.config_id for config in agent_configs
    ]
    assert all(trial["status"] == "completed" and trial["accepted"] for trial in attempts)
    assert report["final"]["status"] == "sealed"
    assert not (data / "final-evaluation.json").exists()
    assert reads == [{"development"}, {"development", "validation"}]
    if validation_mode == "no_gain":
        assert report["outcome"] == "no_agent_validation_improvement"
        assert report["selection"]["agent"]["promoted"] is False
    else:
        assert report["outcome"] == "incomplete"
    assert (Path(report["run_dir"]) / "report.json").exists()


def test_failed_proposals_consume_both_attempts_without_numerical_replacement(
    synthetic_data, tmp_path, monkeypatch
):
    data, _ = synthetic_data
    proposal_calls, evaluated = [], []

    class FailedProposer:
        def __init__(self):
            self.client = SimpleNamespace(
                model_id="synthetic-no-network", max_tokens=100, usage_records=[]
            )
            self.attempts = []

        def propose(self, feedback):
            proposal_calls.append(feedback)
            return Proposal(
                None,
                "",
                "failed",
                "Synthetic provider failure",
                self.client.model_id,
                "synthetic-prompt-v1",
                1,
            )

    def evaluation(config, tasks, labels, output_dir, seed, seconds):
        evaluated.append(config.config_id)
        if tasks[0]["task_id"].startswith("validation"):
            return {
                "status": "failed",
                "config": config.model_dump(),
                "config_id": config.config_id,
            }
        return {
            "status": "completed",
            "config": config.model_dump(),
            "config_id": config.config_id,
            "metrics": synthetic_metrics(0.5),
            "error_analysis": {},
            "seconds": 0.001,
            "artifacts": {},
        }

    monkeypatch.setattr(improvement, "bounded_evaluation", evaluation)
    report = improvement.run_improvement(
        data, tmp_path / "failed-proposals", FailedProposer(), max_seconds=600
    )
    attempts = [
        trial
        for trial in report["trials"]
        if trial["method"] == "agent" and isinstance(trial["attempt"], int)
    ]
    assert len(proposal_calls) == len(attempts) == 2
    assert all(trial["status"] == "proposal_failed" for trial in attempts)
    assert all(trial["accepted"] is False for trial in attempts)
    assert len(evaluated) == 5  # Failed proposals must not trigger validation.
    assert report["outcome"] == "incomplete"
    assert report["all_required_attempts_completed"] is False
    assert "required proposals" in report["error"]
    assert report["final"]["status"] == "sealed"
    assert not (data / "final-evaluation.json").exists()
