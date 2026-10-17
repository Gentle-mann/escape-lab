import json
from types import SimpleNamespace

import pytest

from escape_lab.benchmark import build_tasks, run_benchmark, validate_partition
from escape_lab.data import sha256


@pytest.fixture
def prepared(tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    candidates, outcomes = [], {}
    for group in range(15):
        split = "development" if group < 10 else "validation" if group < 14 else "final"
        for gene in ("FGFR1", "FGFR2"):
            for mut in ("K", "E"):
                position = 450 + group * 12
                candidate_id = f"{gene}:p.Ala{position}{'Lys' if mut == 'K' else 'Glu'}"
                candidates.append(
                    {
                        "candidate_id": candidate_id,
                        "gene": gene,
                        "protein_variant": f"A{position}{mut}",
                        "position": position,
                        "wt": "A",
                        "mut": mut,
                        "aligned_position": position,
                        "group_id": f"block-{group}",
                        "split": split,
                        "context": "AAAAGAA",
                        "reference_sequence_sha256": "a" * 64,
                    }
                )
                outcomes[candidate_id] = {
                    "dual_resistant": mut == "K",
                    "records": [
                        {
                            "nucleotide_variant": f"c.{position}A>T",
                            "pem": {"replicates": [2, 3, 4], "median": 3},
                            "fut": {"replicates": [3, 4, 5], "median": 4},
                        }
                    ],
                    "evidence_ids": [f"evidence:{candidate_id}"],
                }
    # Undefined labels are excluded from fitting, but remain queryable in their pool.
    outcomes[candidates[0]["candidate_id"]]["dual_resistant"] = None
    for name, content in [("candidates", candidates), ("outcomes", outcomes)]:
        (data_dir / f"{name}.json").write_text(json.dumps(content))
    manifest = {
        "candidates_sha256": sha256(data_dir / "candidates.json"),
        "outcomes_sha256": sha256(data_dir / "outcomes.json"),
        "label_rule": "Synthetic boolean fixtures only",
        "grouping": "Fixed aligned blocks",
        "final_test_status": "sealed",
    }
    (data_dir / "manifest.json").write_text(json.dumps(manifest))
    return data_dir, candidates, outcomes


def test_final_stage_is_refused_before_data_access(tmp_path):
    with pytest.raises(ValueError, match="sealed"):
        run_benchmark(tmp_path / "does-not-exist", tmp_path / "out", stage="final")


def test_partitions_hold_out_whole_groups_and_validation_uses_only_development(prepared):
    _, candidates, _ = prepared
    tasks = build_tasks(candidates, "development")
    assert len(tasks) == 3
    all_test_ids = []
    for task in tasks:
        train, pool = task["training"], task["pool"]
        validate_partition(train, pool)
        assert {c["split"] for c in train + pool} == {"development"}
        all_test_ids.extend(c["candidate_id"] for c in pool)
    expected = [c["candidate_id"] for c in candidates if c["split"] == "development"]
    assert sorted(all_test_ids) == sorted(expected)
    validation = build_tasks(candidates, "validation")
    assert len(validation) == 1
    assert {c["split"] for c in validation[0]["training"]} == {"development"}
    assert {c["split"] for c in validation[0]["pool"]} == {"validation"}


def test_cross_manifest_group_leakage_is_rejected(prepared):
    _, candidates, _ = prepared
    copied = [dict(c) for c in candidates]
    copied[0]["split"] = "validation"
    with pytest.raises(ValueError, match="group"):
        build_tasks(copied, "validation")
    with pytest.raises(ValueError, match="group"):
        validate_partition([candidates[0]], [candidates[1]])


def test_benchmark_uses_vault_commits_and_reveals_and_reproducible_random_seeds(prepared, tmp_path):
    data_dir, _, outcomes = prepared
    first = run_benchmark(
        data_dir,
        tmp_path / "report-a",
        stage="validation",
        policies=("random", "gene_prior"),
        seeds=(3, 4),
        budget=3,
    )
    second = run_benchmark(
        data_dir,
        tmp_path / "report-b",
        stage="validation",
        policies=("random", "gene_prior"),
        seeds=(3, 4),
        budget=3,
    )
    assert len(first["trials"]) == 4
    for a, b in zip(first["trials"], second["trials"], strict=True):
        assert a["selected_candidate_ids"] == b["selected_candidate_ids"]
        assert a["queries"] == a["committed_count"] == a["revealed_count"] == 3
        assert a["hits"] == sum(
            outcomes[cid]["dual_resistant"] for cid in a["selected_candidate_ids"]
        )
        assert a["precision"] == a["hits"] / 3
        assert a["valid"] is True
        assert set(a["training_groups"]).isdisjoint(a["pool_groups"])
        assert 0 <= a["initial_average_precision"] <= 1
        assert 0 <= a["initial_brier"] <= 1
    random_trials = [t for t in first["trials"] if t["policy"] == "random"]
    assert random_trials[0]["selected_candidate_ids"] != random_trials[1]["selected_candidate_ids"]
    assert (tmp_path / "report-a" / "validation.json").exists()
    assert (tmp_path / "report-a" / "validation.md").exists()
    assert list((tmp_path / "report-a" / "private" / "vaults").glob("*.sqlite"))
    public_report = json.dumps(first)
    assert "run_token" not in public_report
    assert "block-14" not in public_report
    assert "evidence:FGFR1:p.Ala618Lys" not in public_report


def test_numerical_refits_on_disclosed_observations_and_records_paired_comparisons(
    prepared, tmp_path
):
    data_dir, _, _ = prepared
    report = run_benchmark(
        data_dir,
        tmp_path / "numerical",
        stage="validation",
        policies=("logistic", "extra_trees", "extra_explore", "extra_diverse"),
        seeds=(0,),
        budget=6,
    )
    assert all(trial["fit_count"] == 2 for trial in report["trials"])
    assert all(trial["queries"] == 6 for trial in report["trials"])
    assert all(trial["initial_training_count"] == 39 for trial in report["trials"])
    assert report["paired_comparisons"]
    assert all(
        trial["verdict"]["status"] in {"refuted", "unresolved"} for trial in report["trials"]
    )


def test_agent_sees_only_authorized_observations_and_fallback_is_reported(prepared, tmp_path):
    data_dir, _, _ = prepared
    states = []

    class Investigator:
        decisions_made = 0
        client = SimpleNamespace(model_id="synthetic-nemotron", total_tokens=0, usage_records=[])

        def choose_action(self, state, shortlists):
            states.append(state)
            assert "initial_average_precision" not in state
            assert "outcomes" not in state
            assert len(state["observations"]) == state["step"]
            for observation in state["observations"]:
                assert observation["evidence_id"].startswith("evidence:")
                assert type(observation["dual_resistant"]) is bool
            return SimpleNamespace(
                action="select_candidate",
                candidate_id=shortlists["exploit"][0]["candidate_id"],
                reason="Synthetic test decision",
                fallback_used=True,
                fallback_reason="Synthetic failure",
            )

    report = run_benchmark(
        data_dir,
        tmp_path / "agent",
        stage="validation",
        policies=("agent",),
        seeds=(0,),
        budget=2,
        investigator=Investigator(),
    )
    trial = report["trials"][0]
    assert len(states) == 2
    assert trial["fallback_count"] == 2
    assert trial["agent_unassisted"] is False
    assert trial["queries"] == 2
    assert report["summary"]["agent"]["mean_hits"] == trial["hits"]
    assert report["summary"]["agent"]["unassisted_mean_hits"] is None
    assert report["summary"]["agent"]["fallback_assisted_mean_hits"] == trial["hits"]
    assert trial["trace"][0]["fallback_reason"] == "Synthetic failure"
    assert trial["decision_attempts"][0]["fallback_reason"] == "Synthetic failure"
    assert "agent.py" in report["source_code_fingerprints"]
    assert "crusoe.py" in report["source_code_fingerprints"]


def test_premature_agent_finish_is_an_explicit_failed_trial(prepared, tmp_path):
    data_dir, _, _ = prepared

    class Investigator:
        client = SimpleNamespace(model_id="synthetic", total_tokens=0, usage_records=[])

        def choose_action(self, state, shortlists):
            return SimpleNamespace(action="finish", reason="Pretend done", fallback_used=False)

    report = run_benchmark(
        data_dir,
        tmp_path / "bad-agent",
        stage="validation",
        policies=("agent",),
        seeds=(0,),
        budget=2,
        investigator=Investigator(),
    )
    trial = report["trials"][0]
    assert trial["valid"] is False
    assert trial["queries"] == 0
    assert trial["failure"] == "Agent stopped before the coverage audit budget was exhausted."
    assert report["summary"]["agent"]["valid_trials"] == 0
    assert report["summary"]["agent"]["completion_rate"] == 0
    assert report["summary"]["agent"]["mean_hits"] == 0
    assert report["summary"]["agent"]["completed_only_mean_hits"] is None
    assert trial["decision_attempts"][0]["action"] == "finish"


def test_inconclusive_results_consume_budget_without_leaking_eligibility(prepared, tmp_path):
    data_dir, candidates, outcomes = prepared
    for candidate in candidates:
        if candidate["split"] == "validation":
            outcomes[candidate["candidate_id"]]["dual_resistant"] = None
    path = data_dir / "outcomes.json"
    path.write_text(json.dumps(outcomes))
    manifest_path = data_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["outcomes_sha256"] = sha256(path)
    manifest_path.write_text(json.dumps(manifest))
    report = run_benchmark(
        data_dir,
        tmp_path / "inconclusive",
        stage="validation",
        policies=("logistic",),
        seeds=(0,),
        budget=6,
    )
    trial = report["trials"][0]
    assert trial["queries"] == trial["inconclusive_reveals"] == 6
    assert trial["committed_count"] == trial["revealed_count"] == 6
    assert trial["hits"] == 0
    assert trial["valid"] is True
    assert trial["fit_count"] == 2
    assert trial["pool_count"] == 16
    assert trial["initial_audit_count"] == 0
    assert trial["initial_average_precision"] is None
    assert trial["initial_brier"] is None
    assert report["summary"]["logistic"]["mean_initial_brier"] is None


@pytest.mark.parametrize(
    "kwargs", [{"budget": 0}, {"seeds": ()}, {"policies": ("invented",)}, {"policies": ("agent",)}]
)
def test_invalid_configuration_is_rejected(prepared, tmp_path, kwargs):
    data_dir, _, _ = prepared
    with pytest.raises(ValueError):
        run_benchmark(data_dir, tmp_path / "invalid", **kwargs)


def test_source_changes_fail_checksum_verification(prepared, tmp_path):
    data_dir, _, _ = prepared
    path = data_dir / "candidates.json"
    path.write_text(path.read_text() + " ")
    with pytest.raises(ValueError, match="checksum"):
        run_benchmark(data_dir, tmp_path / "mismatch", policies=("random",), seeds=(0,), budget=2)
