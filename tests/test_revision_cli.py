import json

import pytest

from escape_lab.data import sha256
from escape_lab.revision_cli import verify_hashes


def test_frozen_artifact_change_is_detected(tmp_path):
    path = tmp_path / "artifact.json"
    path.write_text(json.dumps({"a": 1}))
    hashes = {str(path): sha256(path)}
    verify_hashes(hashes)
    path.write_text(json.dumps({"a": 2}))
    with pytest.raises(ValueError, match="changed"):
        verify_hashes(hashes)
