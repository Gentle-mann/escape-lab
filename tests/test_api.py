import pytest
from fastapi.testclient import TestClient

from escape_lab.api import create_app
from escape_lab.vault import Vault


@pytest.fixture
def client(tmp_path):
    candidate = {
        "candidate_id": "FGFR2:N549K",
        "gene": "FGFR2",
        "protein_variant": "N549K",
        "position": 549,
        "wt": "N",
        "mut": "K",
        "aligned_position": 549,
        "group_id": "position:549",
        "split": "development",
    }
    outcome = {
        "dual_resistant": True,
        "records": [{"hidden": "secret measurement"}],
        "evidence_ids": ["workbook:11"],
    }
    vault = Vault(tmp_path / "vault.sqlite", [candidate], {candidate["candidate_id"]: outcome})
    return TestClient(create_app(vault, "test-admin-token"))


def create_run(client):
    return client.post(
        "/runs",
        headers={"X-Admin-Token": "test-admin-token"},
        json={
            "candidate_ids": ["FGFR2:N549K"],
            "budget": 1,
            "policy_version": "v1",
            "split": "development",
            "seed": 1,
        },
    ).json()


def test_admin_and_run_auth_are_distinct_and_fail_closed(client):
    response = client.post(
        "/runs",
        json={
            "candidate_ids": ["FGFR2:N549K"],
            "budget": 1,
            "policy_version": "v1",
            "split": "development",
        },
    )
    assert response.status_code == 401
    run = create_run(client)
    url = f"/runs/{run['run_id']}/state"
    assert client.get(url).status_code == 401
    assert client.get(url, headers={"X-Run-Token": "test-admin-token"}).status_code == 401
    assert client.get(url, headers={"X-Run-Token": run["run_token"]}).status_code == 200


def test_api_end_to_end_hidden_outcomes_and_validated_verdict(client):
    run = create_run(client)
    base = f"/runs/{run['run_id']}"
    headers = {"X-Run-Token": run["run_token"]}
    candidates = client.get(f"{base}/candidates", headers=headers)
    assert "secret measurement" not in candidates.text
    receipt = client.post(
        f"{base}/commits",
        headers=headers,
        json={
            "candidate_id": "FGFR2:N549K",
            "predicted_probability": 0.8,
            "reason": "Hypothesis",
            "policy_version": "v1",
            "idempotency_key": "c1",
        },
    )
    assert receipt.status_code == 200
    assert "dual_resistant" not in receipt.text
    response = client.post(
        f"{base}/reveals",
        headers=headers,
        json={"commit_id": receipt.json()["commit_id"], "idempotency_key": "r1"},
    )
    assert response.status_code == 200
    assert response.json()["outcome"]["dual_resistant"] is True
    verdict = client.post(
        f"{base}/verdict",
        headers=headers,
        json={"status": "refuted", "evidence_ids": ["workbook:11"]},
    )
    assert verdict.status_code == 200
    assert verdict.json()["status"] == "refuted"
    assert (
        client.post(
            f"{base}/verdict", headers=headers, json={"status": "safe", "evidence_ids": []}
        ).status_code
        == 422
    )


def test_api_rejects_unrecognized_fields_and_invalid_probabilities(client):
    run = create_run(client)
    body = {
        "candidate_id": "FGFR2:N549K",
        "predicted_probability": 0.8,
        "reason": "Hypothesis",
        "policy_version": "v1",
        "idempotency_key": "c1",
        "outcome_override": True,
    }
    url = f"/runs/{run['run_id']}/commits"
    headers = {"X-Run-Token": run["run_token"]}
    assert client.post(url, headers=headers, json=body).status_code == 422
    del body["outcome_override"]
    body["predicted_probability"] = True
    assert client.post(url, headers=headers, json=body).status_code == 422


def test_api_must_not_start_with_an_empty_admin_token(tmp_path):
    with pytest.raises(ValueError, match="admin"):
        create_app(Vault(tmp_path / "empty.sqlite", [], {}), "")
