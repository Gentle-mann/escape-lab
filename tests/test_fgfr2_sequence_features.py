"""Outcome-free invariants for the pinned, frozen sequence representation."""

import hashlib
import json
import sys
from contextlib import contextmanager
from types import SimpleNamespace

import numpy as np
import pytest

from escape_lab import fgfr2_sequence_features as module

SEQUENCE = "ACDEFG"


def candidates(sequence=SEQUENCE):
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    return [
        {"candidate_id": "FGFR2:p.Cys2Trp", "protein_variant": "p.Cys2Trp",
         "gene": "FGFR2", "position": 2, "wt": "C", "mut": "W",
         "split": "development", "group_id": "block-a", "reference_sequence_sha256": digest},
        {"candidate_id": "FGFR2:p.Gly6Ala", "protein_variant": "p.Gly6Ala",
         "gene": "FGFR2", "position": 6, "wt": "G", "mut": "A",
         "split": "validation", "group_id": "block-b", "reference_sequence_sha256": digest},
    ]


@pytest.fixture
def setup(tmp_path, monkeypatch):
    references = tmp_path / "references.json"
    references.write_text(json.dumps({
        "sequences": {"FGFR2": SEQUENCE},
        "sequence_sha256": {"FGFR2": hashlib.sha256(SEQUENCE.encode()).hexdigest()},
    }))
    files = {}
    for name in module.CHECKPOINT_FILENAMES:
        path = tmp_path / name
        path.write_bytes(b"synthetic-checkpoint-for-invariant-tests")
        files[name] = {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    weight = tmp_path / "model.safetensors"

    class FakeEncoder:
        def __init__(self):
            self.calls = []

        def encode(self, sequences):
            self.calls.append(list(sequences))
            output = np.zeros((len(sequences), len(SEQUENCE) + 2, 480), dtype=np.float32)
            output[:, 0] = -1000  # BOS must never be mistaken for residue one.
            output[:, -1] = 1000
            for row, sequence in enumerate(sequences):
                for index, residue in enumerate(sequence, start=1):
                    output[row, index] = ord(residue) + np.arange(480) / 1000
            return output

    encoder = FakeEncoder()
    provenance = {
        "files": files,
        "packages": {"torch": "synthetic", "transformers": "synthetic"},
        "device": "cpu",
    }
    monkeypatch.setattr(module, "_load_encoder", lambda **kwargs: (encoder, provenance))
    return references, encoder, tmp_path / "cache", weight


def test_bos_position_signed_delta_and_single_mutation(setup):
    references, encoder, output, _ = setup
    cache = module.extract_features(candidates(), references, output, batch_size=2)
    assert encoder.calls == [[SEQUENCE], ["AWDEFG", "ACDEFA"]]
    np.testing.assert_allclose(cache.wt[:, 0], [ord("C"), ord("G")])
    np.testing.assert_allclose(cache.delta[:, 0], [ord("W") - ord("C"), ord("A") - ord("G")])
    assert cache.delta[1, 0] < 0  # Direction must not be erased by an absolute value.
    assert cache.features("FGFR2:p.Cys2Trp")["wt"].shape == (480,)
    assert cache.matrix(list(reversed(cache.candidate_ids)), "wt_delta").shape == (2, 960)
    assert cache.manifest["outcomes_used"] is False


@pytest.mark.parametrize("field,value", [
    ("gene", "FGFR1"), ("split", "final"), ("wt", "A"), ("mut", "C"),
    ("position", 0), ("position", 7), ("reference_sequence_sha256", "bad"),
    ("candidate_id", "FGFR2:p.Cys2Ala"), ("protein_variant", "p.Cys3Trp"),
])
def test_invalid_input_rejected_before_any_model_execution(setup, field, value):
    references, encoder, output, _ = setup
    rows = candidates()
    rows[0][field] = value
    with pytest.raises(ValueError):
        module.extract_features(rows, references, output)
    assert encoder.calls == []
    assert not output.exists()


def test_duplicate_candidate_rejected(setup):
    references, encoder, output, _ = setup
    with pytest.raises(ValueError, match="duplicate"):
        module.extract_features([candidates()[0]] * 2, references, output)
    assert encoder.calls == []


def test_reference_hash_mismatch_rejected(setup):
    references, _, output, _ = setup
    contents = json.loads(references.read_text())
    contents["sequence_sha256"]["FGFR2"] = "bad"
    references.write_text(json.dumps(contents))
    with pytest.raises(ValueError, match="sequence hash"):
        module.extract_features(candidates(), references, output)


@pytest.mark.parametrize("bad_output", [np.zeros((1, 8, 479)), np.full((1, 8, 480), np.nan)])
def test_invalid_model_output_is_rejected_without_cache(setup, bad_output):
    references, encoder, output, _ = setup
    encoder.encode = lambda sequences: bad_output
    with pytest.raises(ValueError):
        module.extract_features(candidates(), references, output)
    assert not (output / "manifest.json").exists()


def test_cache_roundtrip_reorders_identical_candidate_set(setup):
    references, encoder, output, _ = setup
    original = module.extract_features(candidates(), references, output)
    loaded = module.load_features(list(reversed(candidates())), references, output)
    np.testing.assert_array_equal(loaded.features(candidates()[0]["candidate_id"])["delta"],
                                  original.delta[0])
    module.extract_features(candidates(), references, output)
    assert len(encoder.calls) == 2  # Valid cache reused, no second model run.


def test_changed_input_identity_or_subset_rejected(setup):
    references, _, output, _ = setup
    module.extract_features(candidates(), references, output)
    altered = candidates()
    altered[0]["group_id"] = "different-block"
    with pytest.raises(ValueError, match="candidate"):
        module.load_features(altered, references, output)
    with pytest.raises(ValueError, match="candidate"):
        module.load_features(candidates()[:1], references, output)


@pytest.mark.parametrize("field,value", [
    ("model_revision", "unfrozen-main"), ("embedding_dim", 640),
    ("representation", "absolute-difference"), ("outcomes_used", True),
])
def test_cache_provenance_tampering_rejected(setup, field, value):
    references, _, output, _ = setup
    module.extract_features(candidates(), references, output)
    path = output / "manifest.json"
    manifest = json.loads(path.read_text())
    manifest[field] = value
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="provenance"):
        module.load_features(candidates(), references, output)


def test_feature_file_corruption_rejected(setup):
    references, _, output, _ = setup
    module.extract_features(candidates(), references, output)
    with (output / "features.npz").open("ab") as file:
        file.write(b"changed")
    with pytest.raises(ValueError, match="feature file hash"):
        module.load_features(candidates(), references, output)


def test_checkpoint_file_corruption_rejected(setup):
    references, _, output, weight = setup
    module.extract_features(candidates(), references, output)
    weight.write_bytes(b"changed")
    with pytest.raises(ValueError, match="checkpoint"):
        module.load_features(candidates(), references, output)


def test_partial_cache_is_not_silently_overwritten(setup):
    references, _, output, _ = setup
    output.mkdir()
    (output / "features.npz").write_bytes(b"partial")
    with pytest.raises(ValueError, match="incomplete"):
        module.extract_features(candidates(), references, output)


def test_cache_feature_accessor_does_not_expose_mutable_backing_array(setup):
    references, _, output, _ = setup
    cache = module.extract_features(candidates(), references, output)
    values = cache.features(cache.candidate_ids[0])
    values["wt"][0] = -1
    assert cache.features(cache.candidate_ids[0])["wt"][0] == ord("C")


def test_loader_pins_revision_safe_weights_and_frozen_evaluation(tmp_path, monkeypatch):
    calls = {}
    for name in module.CHECKPOINT_FILENAMES:
        (tmp_path / name).write_text("synthetic")

    def snapshot_download(**kwargs):
        calls["download"] = kwargs
        return str(tmp_path)

    class FakeModel:
        config = SimpleNamespace(hidden_size=480, num_hidden_layers=12)

        def requires_grad_(self, value):
            calls["requires_grad"] = value

        def to(self, device):
            calls["device"] = device

        def eval(self):
            calls["eval"] = True

    class FakeAutoModel:
        @staticmethod
        def from_pretrained(path, **kwargs):
            calls["model"] = kwargs
            return FakeModel()

    class FakeAutoTokenizer:
        @staticmethod
        def from_pretrained(path, **kwargs):
            calls["tokenizer"] = kwargs
            return object()

    monkeypatch.setitem(sys.modules, "torch", SimpleNamespace())
    monkeypatch.setitem(sys.modules, "huggingface_hub", SimpleNamespace(
        snapshot_download=snapshot_download,
    ))
    monkeypatch.setitem(sys.modules, "transformers", SimpleNamespace(
        AutoModel=FakeAutoModel, AutoTokenizer=FakeAutoTokenizer,
    ))
    monkeypatch.setattr(module.importlib.metadata, "version", lambda name: "synthetic")
    _, provenance = module._load_encoder(device="cpu", cache_dir=None)
    assert calls["download"]["repo_id"] == module.MODEL_ID
    assert calls["download"]["revision"] == module.MODEL_REVISION
    assert calls["model"]["use_safetensors"] is True
    assert calls["model"]["trust_remote_code"] is False
    assert calls["model"]["local_files_only"] is True
    assert calls["tokenizer"]["trust_remote_code"] is False
    assert calls["requires_grad"] is False
    assert calls["eval"] is True
    assert set(provenance["files"]) == set(module.CHECKPOINT_FILENAMES)


@pytest.mark.parametrize("wrong_bos", [False, True])
def test_forward_pass_has_inference_mode_and_exact_token_alignment(wrong_bos):
    state = {"inference": False, "calls": 0}

    class Tensor:
        def __init__(self, values):
            self.values = np.asarray(values)

        def cpu(self):
            return self

        def numpy(self):
            return self.values

        def to(self, *args, **kwargs):
            return self

        def detach(self):
            return self

        def __getitem__(self, index):
            return Tensor(self.values[index])

        def sum(self):
            return Tensor(self.values.sum())

        def item(self):
            return self.values.item()

    class Tokenizer:
        cls_token_id, eos_token_id = 0, 1

        def convert_tokens_to_ids(self, residues):
            return [ord(residue) for residue in residues]

        def __call__(self, sequences, **kwargs):
            assert kwargs["truncation"] is False
            assert kwargs["add_special_tokens"] is True
            ids = [[99 if wrong_bos else 0, *self.convert_tokens_to_ids(sequence), 1]
                   for sequence in sequences]
            return {"input_ids": Tensor(ids), "attention_mask": Tensor(np.ones_like(ids))}

    @contextmanager
    def inference_mode():
        state["inference"] = True
        try:
            yield
        finally:
            state["inference"] = False

    def model(**encoded):
        assert state["inference"] is True
        state["calls"] += 1
        return SimpleNamespace(last_hidden_state=Tensor(np.zeros((1, 8, 480))))

    torch = SimpleNamespace(inference_mode=inference_mode, float32=np.float32)
    encoder = module._FrozenEncoder(Tokenizer(), model, torch, "cpu")
    if wrong_bos:
        with pytest.raises(ValueError, match="alignment"):
            encoder.encode([SEQUENCE])
        assert state["calls"] == 0
    else:
        assert encoder.encode([SEQUENCE]).shape == (1, 8, 480)
        assert state["calls"] == 1
        assert state["inference"] is False
