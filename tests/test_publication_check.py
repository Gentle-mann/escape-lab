import hashlib
import json

import pytest

from escape_lab.publication_check import PublicationError, check_publication


def fixture_repo(tmp_path):
    archive = tmp_path / "docs/archive"
    archive.mkdir(parents=True)
    result = archive / "result.json"
    result.write_text('{"aggregate": true}\n')
    entry = {"public_path": "docs/archive/result.json",
             "public_sha256": hashlib.sha256(result.read_bytes()).hexdigest()}
    (archive / "manifest.json").write_text(json.dumps({"artifacts": [entry]}))
    (tmp_path / "README.md").write_text("# Example\n[Evidence](docs/archive/result.json)\n")
    return tmp_path, result


def test_checks_only_public_files_and_validates_hashes(tmp_path):
    root, _ = fixture_repo(tmp_path)
    (root / ".env").write_text("not public")
    result = check_publication(root, check_headline=False)
    assert result["artifacts"] == 1
    assert result["relative_links"] == 1


def test_changed_archive_fails(tmp_path):
    root, result = fixture_repo(tmp_path)
    result.write_text("changed")
    with pytest.raises(PublicationError, match="hash"):
        check_publication(root, check_headline=False)


def test_missing_link_fails(tmp_path):
    root, _ = fixture_repo(tmp_path)
    (root / "README.md").write_text("[Missing](docs/missing.md)\n")
    with pytest.raises(PublicationError, match="link"):
        check_publication(root, check_headline=False)


def test_manifest_cannot_read_outside_public_archive(tmp_path):
    root, _ = fixture_repo(tmp_path)
    path = root / "docs/archive/manifest.json"
    path.write_text(json.dumps({"artifacts": [{"public_path": ".env", "public_sha256": "0"}]}))
    with pytest.raises(PublicationError, match="archive"):
        check_publication(root, check_headline=False)


def test_unknown_archive_file_is_detected(tmp_path):
    root, _ = fixture_repo(tmp_path)
    (root / "docs/archive/unregistered.json").write_text("{}")
    with pytest.raises(PublicationError, match="unregistered"):
        check_publication(root, check_headline=False)


def test_inline_and_fenced_code_are_not_link_targets(tmp_path):
    root, _ = fixture_repo(tmp_path)
    (root / "README.md").write_text("`[example](missing.md)`\n```md\n[x](missing.md)\n```\n")
    assert check_publication(root, check_headline=False)["relative_links"] == 0


def test_headline_drift_is_rejected(tmp_path):
    root, _ = fixture_repo(tmp_path)
    from escape_lab.publication_check import verify_headline
    source = root / "docs/archive/reports/mechanistic-pilot/fgfr2-sequence-ap-v1/results.json"
    source.parent.mkdir(parents=True)
    metrics = {
        "probability": {"mean_fold_high_ap": .09732091648, "high_hits": 1, "budget": 68},
        "margin": {"mean_fold_high_ap": .11134574838, "high_hits": 4, "budget": 68},
    }
    for arm in metrics.values():
        arm.update(strong_hits=0, by_fold={"fold": {"high": {"positives": 26},
                                                   "strong": {"positives": 7}}})
    source.write_text(json.dumps({"metrics": metrics}))
    for path in [root / "README.md", root / "docs/RESULTS.md"]:
        path.write_text("| Original classifier probabilities | 0.09732 | 1 / 26 | 0 / 7 |\n"
                        "| Same classifier, raw decision margin | 0.21135 | 4 / 26 | 0 / 7 |\n"
                        "14.4% relative gain\n")
    with pytest.raises(PublicationError, match="headline"):
        verify_headline(root)


@pytest.mark.parametrize("old,new", [("4 / 26", "4 / 20"), ("0 / 7", "1 / 7")])
def test_headline_retrieval_denominators_and_strong_hits_are_checked(tmp_path, old, new):
    root, _ = fixture_repo(tmp_path)
    from escape_lab.publication_check import verify_headline
    source = root / "docs/archive/reports/mechanistic-pilot/fgfr2-sequence-ap-v1/results.json"
    source.parent.mkdir(parents=True)
    metrics = {
        "probability": {"mean_fold_high_ap": .09732091648, "high_hits": 1, "budget": 68},
        "margin": {"mean_fold_high_ap": .11134574838, "high_hits": 4, "budget": 68},
    }
    for arm in metrics.values():
        arm.update(strong_hits=0, by_fold={"fold": {"high": {"positives": 26},
                                                   "strong": {"positives": 7}}})
    source.write_text(json.dumps({"metrics": metrics}))
    valid = ("| Original classifier probabilities | 0.09732 | 1 / 26 | 0 / 7 |\n"
             "| Same classifier, raw decision margin | 0.11135 | 4 / 26 | 0 / 7 |\n"
             "14.4% relative gain within 68 picks\n")
    for path in [root / "README.md", root / "docs/RESULTS.md"]:
        path.write_text(valid)
    verify_headline(root)
    (root / "README.md").write_text(valid.replace(old, new))
    with pytest.raises(PublicationError, match="headline"):
        verify_headline(root)
