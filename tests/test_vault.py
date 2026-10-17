import json
import sqlite3
from concurrent.futures import ThreadPoolExecutor

import pytest

from escape_lab.vault import Vault, VaultError


def dataset():
    candidates = [
        {
            "candidate_id": f"FGFR2:{variant}",
            "gene": "FGFR2",
            "protein_variant": variant,
            "position": position,
            "wt": variant[0],
            "mut": variant[-1],
            "aligned_position": position,
            "group_id": f"position:{position}",
            "split": "development",
        }
        for variant, position in [("N549K", 549), ("V564G", 564), ("A600V", 600)]
    ]
    outcomes = {
        item["candidate_id"]: {
            "dual_resistant": index == 0,
            "records": [
                {
                    "source_row": index + 11,
                    "pem_replicates": [2, 3, 4],
                    "fut_replicates": [3, 4, 5],
                    "hidden_marker": "DO_NOT_LEAK",
                }
            ],
            "evidence_ids": [f"source:row:{index + 11}"],
        }
        for index, item in enumerate(candidates)
    }
    return candidates, outcomes


@pytest.fixture
def lab(tmp_path):
    candidates, outcomes = dataset()
    vault = Vault(tmp_path / "vault.sqlite", candidates, outcomes)
    run = vault.create_run([c["candidate_id"] for c in candidates], 2, "v1", "development", 7)
    return vault, run, candidates


def commit(vault, run, candidate, key="commit-1", probability=0.8):
    return vault.commit(
        run["run_id"],
        run["run_token"],
        candidate["candidate_id"],
        probability,
        "Test a falsifiable resistance hypothesis",
        "v1",
        key,
    )


def test_unrevealed_data_hidden_and_measurements_require_prior_commit(lab):
    vault, run, candidates = lab
    args = (run["run_id"], run["run_token"])
    receipt = commit(vault, run, candidates[0])
    public = json.dumps([vault.candidates(*args), vault.state(*args), vault.events(*args), receipt])
    assert "DO_NOT_LEAK" not in public
    assert "dual_resistant" not in public
    with pytest.raises(VaultError):
        vault.reveal(*args, "unknown-commit", "unknown-reveal")
    result = vault.reveal(*args, receipt["commit_id"], "reveal-1")
    assert result["outcome"]["dual_resistant"] is True
    assert result["outcome"]["evidence_ids"] == ["source:row:11"]
    assert result["outcome"]["records"][0]["hidden_marker"] == "DO_NOT_LEAK"
    assert vault.state(*args)["revealed_count"] == 1


def test_budget_reserved_at_commit_and_retry_is_idempotent(lab):
    vault, run, candidates = lab
    args = (run["run_id"], run["run_token"])
    first = commit(vault, run, candidates[0])
    assert commit(vault, run, candidates[0]) == first
    with pytest.raises(VaultError, match="idempotency"):
        commit(vault, run, candidates[1])
    with pytest.raises(VaultError, match="already committed"):
        commit(vault, run, candidates[0], "different-key")
    second = commit(vault, run, candidates[1], "commit-2")
    with pytest.raises(VaultError, match="budget"):
        commit(vault, run, candidates[2], "commit-3")
    result = vault.reveal(*args, first["commit_id"], "reveal-1")
    assert vault.reveal(*args, first["commit_id"], "reveal-1") == result
    assert vault.reveal(*args, first["commit_id"], "another-retry") == result
    with pytest.raises(VaultError, match="idempotency"):
        vault.reveal(*args, second["commit_id"], "reveal-1")
    state = vault.state(*args)
    assert state["committed_count"] == 2
    assert state["revealed_count"] == 1
    assert state["remaining_budget"] == 0


def test_authentication_and_cross_run_commit_isolation(lab):
    vault, run, candidates = lab
    other = vault.create_run([c["candidate_id"] for c in candidates], 1, "v1", "development")
    receipt = commit(vault, run, candidates[0])
    for method in (vault.candidates, vault.state, vault.events):
        with pytest.raises(VaultError, match="Unauthorized"):
            method(run["run_id"], other["run_token"])
    with pytest.raises(VaultError, match="not found"):
        vault.reveal(other["run_id"], other["run_token"], receipt["commit_id"], "cross-run")
    assert vault.state(other["run_id"], other["run_token"])["revealed_count"] == 0


@pytest.mark.parametrize("probability", [-0.1, 1.1, float("nan"), float("inf"), True, "0.5"])
def test_invalid_probabilities_rejected_without_spending_budget(lab, probability):
    vault, run, candidates = lab
    with pytest.raises(VaultError, match="probability"):
        commit(vault, run, candidates[0], probability=probability)
    assert vault.state(run["run_id"], run["run_token"])["committed_count"] == 0


def test_out_of_pool_unknown_and_wrong_policy_rejected(lab):
    vault, _, candidates = lab
    run = vault.create_run([candidates[0]["candidate_id"]], 1, "v1", "development")
    with pytest.raises(VaultError, match="pool"):
        commit(vault, run, candidates[1])
    with pytest.raises(VaultError, match="pool"):
        commit(vault, run, {"candidate_id": "invented"})
    with pytest.raises(VaultError, match="policy"):
        vault.commit(
            run["run_id"],
            run["run_token"],
            candidates[0]["candidate_id"],
            0.5,
            "reason",
            "invented",
            "c1",
        )


def test_concurrent_commits_cannot_oversubscribe_and_reveals_spend_once(lab):
    vault, _, candidates = lab
    run = vault.create_run([c["candidate_id"] for c in candidates], 1, "v1", "development")

    def attempt(index):
        try:
            return commit(vault, run, candidates[index], f"c-{index}")
        except VaultError:
            return None

    with ThreadPoolExecutor(max_workers=3) as workers:
        receipts = list(workers.map(attempt, range(3)))
    receipt = next(item for item in receipts if item)
    assert sum(item is not None for item in receipts) == 1
    with ThreadPoolExecutor(max_workers=8) as workers:
        reveals = list(
            workers.map(
                lambda _: vault.reveal(
                    run["run_id"], run["run_token"], receipt["commit_id"], "same"
                ),
                range(8),
            )
        )
    assert all(item == reveals[0] for item in reveals)
    state = vault.state(run["run_id"], run["run_token"])
    assert state["revealed_count"] == state["committed_count"] == 1
    assert (
        len(
            [
                e
                for e in vault.events(run["run_id"], run["run_token"])
                if e["event_type"] == "revealed"
            ]
        )
        == 1
    )


def test_verdict_requires_revealed_positive_evidence(lab):
    vault, run, candidates = lab
    args = (run["run_id"], run["run_token"])
    with pytest.raises(VaultError):
        vault.verdict(*args, "refuted", ["source:row:11"])
    with pytest.raises(VaultError):
        vault.verdict(*args, "safe", [])
    assert vault.verdict(*args, "unresolved", [])["status"] == "unresolved"
    negative = commit(vault, run, candidates[1], "negative")
    vault.reveal(*args, negative["commit_id"], "negative-reveal")
    with pytest.raises(VaultError):
        vault.verdict(*args, "refuted", ["source:row:12"])
    positive = commit(vault, run, candidates[0], "positive")
    vault.reveal(*args, positive["commit_id"], "positive-reveal")
    result = vault.verdict(*args, "refuted", ["source:row:11"])
    assert result["verified_candidate_ids"] == [candidates[0]["candidate_id"]]


def test_reopen_fingerprint_validation_and_append_only_ledger(lab, tmp_path):
    vault, run, candidates = lab
    commit(vault, run, candidates[0])
    reopened = Vault(tmp_path / "vault.sqlite")
    assert reopened.state(run["run_id"], run["run_token"])["committed_count"] == 1
    original_candidates, outcomes = dataset()
    Vault(tmp_path / "vault.sqlite", original_candidates, outcomes)
    outcomes[candidates[0]["candidate_id"]]["dual_resistant"] = False
    with pytest.raises(VaultError, match="fingerprint"):
        Vault(tmp_path / "vault.sqlite", original_candidates, outcomes)
    with sqlite3.connect(tmp_path / "vault.sqlite") as connection:
        for statement in ("DELETE FROM events", "UPDATE commits SET reason = 'tampered'"):
            with pytest.raises(sqlite3.DatabaseError, match="append-only"):
                connection.execute(statement)


def test_import_rejects_hidden_public_fields_and_split_crossing(tmp_path):
    candidates, outcomes = dataset()
    candidates[0]["dual_resistant"] = True
    with pytest.raises(VaultError, match="public"):
        Vault(tmp_path / "bad.sqlite", candidates, outcomes)
    del candidates[0]["dual_resistant"]
    candidates[1]["split"] = "validation"
    candidates[2]["split"] = "final"
    vault = Vault(tmp_path / "good.sqlite", candidates, outcomes)
    with pytest.raises(VaultError, match="split"):
        vault.create_run([c["candidate_id"] for c in candidates], 1, "v1", "development")
    with pytest.raises(VaultError, match="sealed"):
        vault.create_run([candidates[2]["candidate_id"]], 1, "v1", "final")


def test_optional_sequence_features_are_public_but_unknown_fields_are_rejected(tmp_path):
    candidates, outcomes = dataset()
    candidates[0]["context"] = "AAVNQDF"
    candidates[0]["reference_sequence_sha256"] = "a" * 64
    vault = Vault(tmp_path / "sequence.sqlite", candidates, outcomes)
    run = vault.create_run([candidates[0]["candidate_id"]], 1, "v1", "development")
    public = vault.candidates(run["run_id"], run["run_token"])
    assert public[0]["context"] == "AAVNQDF"
    assert public[0]["reference_sequence_sha256"] == "a" * 64
    candidates[1]["resistance_prediction_from_source"] = 0.99
    with pytest.raises(VaultError, match="public"):
        Vault(tmp_path / "unknown.sqlite", candidates, outcomes)
