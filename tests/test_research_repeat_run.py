import json

import pytest

from escape_lab import research_repeat_run as run


def test_jobs_have_ten_matched_triplets_for_every_fold():
    jobs = run.make_jobs(["outer-0", "outer-1", "outer-2"])
    assert len(jobs) == 90
    assert len({job["job_id"] for job in jobs}) == 90
    for repeat in range(10):
        for fold in ("outer-0", "outer-1", "outer-2"):
            matched = [j for j in jobs if j["repeat"] == repeat and j["fold"] == fold]
            assert {j["arm"] for j in matched} == {"control", "relevant", "placebo"}
            assert len({j["selection_seed"] for j in matched}) == 1
    assert jobs == run.make_jobs(["outer-0", "outer-1", "outer-2"])


def test_joint_freeze_refuses_incomplete_job_coverage(tmp_path):
    (tmp_path / "protocol.json").write_text("{}")
    with pytest.raises(ValueError, match="coverage"):
        run.freeze_predictions(tmp_path, [{"job_id": "missing"}], {}, {})


def test_scoring_does_not_read_labels_before_valid_commitment(tmp_path, monkeypatch):
    accessed = []
    monkeypatch.setattr(run, "load_labels", lambda *a, **kw: accessed.append(True))
    with pytest.raises((FileNotFoundError, ValueError)):
        run.score_frozen(tmp_path, tmp_path, [], [])
    assert not accessed


def test_ninety_jobs_without_committed_predictions_cannot_read_labels(tmp_path, monkeypatch):
    from escape_lab.data import sha256

    jobs = run.make_jobs(["outer-0", "outer-1", "outer-2"])
    paths = {"protocol.json": {"schedule": jobs, "source_hashes": {}},
             "jobs.json": jobs, "choices.json": {}, "controls.json": {},
             "evidence-packets.json": {}, "preservation-before.json": {}}
    hashes = {}
    for name, value in paths.items():
        path = tmp_path / name
        path.write_text(json.dumps(value))
        hashes[str(path)] = sha256(path)
    manifest = {"schema_version": "repeated-research-joint-freeze-v1",
                "job_ids": sorted(j["job_id"] for j in jobs), "hashes": hashes}
    (tmp_path / "joint-freeze.json").write_text(json.dumps(manifest))
    accessed = []
    monkeypatch.setattr(run, "load_labels", lambda *a, **kw: accessed.append(True))
    tasks = [{"task_id": f"outer-{i}", "pool": []} for i in range(3)]
    with pytest.raises(ValueError, match="coverage"):
        run.score_frozen(tmp_path, tmp_path, [], tasks)
    assert not accessed


def test_label_attachment_requires_exact_prediction_identity():
    rows = [{"candidate_id": "a", "group_id": "g", "drug": "pem", "probability": 0.5}]
    labels = {"a": {"pem": True, "fut": False}}
    candidates = [{"candidate_id": "a", "group_id": "g", "gene": "FGFR1"}]
    with pytest.raises(ValueError, match="coverage"):
        run.attach_labels(rows, labels, candidates)


def test_prediction_rows_never_contain_truth_before_commitment():
    rows = [{"candidate_id": "a", "group_id": "g", "gene": "FGFR1", "drug": "pem",
             "probability": 0.2, "label": True}]
    with pytest.raises(ValueError, match="label"):
        run.validate_blind_rows(rows)


def test_frozen_manifest_detects_modified_prediction(tmp_path):
    pred = tmp_path / "predictions.json"
    pred.write_text("[]")
    from escape_lab.data import sha256
    manifest = {"schema_version": "repeated-research-joint-freeze-v1",
                "job_ids": ["one"], "hashes": {str(pred): sha256(pred)}}
    (tmp_path / "joint-freeze.json").write_text(json.dumps(manifest))
    pred.write_text("[1]")
    with pytest.raises(ValueError):
        run.verify_freeze(tmp_path)
