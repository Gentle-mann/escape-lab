"""Frozen, outcome-free ESM2 representations for the FGFR2 development pilot.

This module never opens an assay outcome file. The evaluator supplies an eligible
candidate list. Protein positions are one-based: residue ``position - 1`` in the
sequence occupies token ``position`` after the leading BOS token.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import sys
import time
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from Bio.Data.IUPACData import protein_letters_1to3

MODEL_ID = "facebook/esm2_t12_35M_UR50D"
MODEL_REVISION = "6fbf070e65b0b7291e7bbcd451118c216cff79d8"
EMBEDDING_DIM = 480
SCHEMA_VERSION = "fgfr2-frozen-sequence-v1"
REPRESENTATION = "final-layer-wt-residue-and-signed-mutant-minus-wt-residue"
CHECKPOINT_FILENAMES = (
    "config.json", "model.safetensors", "special_tokens_map.json",
    "tokenizer_config.json", "vocab.txt",
)
IDENTITY_FIELDS = (
    "candidate_id", "gene", "protein_variant", "position", "wt", "mut",
    "group_id", "split", "reference_sequence_sha256",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        while chunk := file.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_hash(value) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _inputs(candidates: list[dict], references_path: Path) -> tuple[str, list[dict]]:
    references = json.loads(Path(references_path).read_text())
    sequence = references["sequences"]["FGFR2"]
    if not sequence or any(residue not in protein_letters_1to3 for residue in sequence):
        raise ValueError("FGFR2 reference must contain only canonical amino acids")
    # Stay within the pretrained context; never silently crop or truncate.
    if len(sequence) > 1022:
        raise ValueError("FGFR2 reference exceeds the fixed full-sequence context")
    sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
    if references["sequence_sha256"]["FGFR2"] != sequence_hash:
        raise ValueError("reference sequence hash mismatch")
    if not candidates:
        raise ValueError("at least one eligible FGFR2 candidate is required")
    identities = []
    seen = set()
    for candidate in candidates:
        try:
            identity = {name: candidate[name] for name in IDENTITY_FIELDS}
        except KeyError as error:
            raise ValueError(f"missing candidate identity field: {error.args[0]}") from error
        candidate_id = identity["candidate_id"]
        if candidate_id in seen:
            raise ValueError("duplicate candidate identity")
        seen.add(candidate_id)
        position, wt, mutant = identity["position"], identity["wt"], identity["mut"]
        if identity["gene"] != "FGFR2" or identity["split"] not in {
            "development", "validation"
        }:
            raise ValueError("only evaluator-approved development/validation FGFR2 candidates")
        if type(position) is not int or position < 1 or position > len(sequence):
            raise ValueError("candidate position outside full reference sequence")
        if wt not in protein_letters_1to3 or mutant not in protein_letters_1to3 or wt == mutant:
            raise ValueError("candidate must specify one canonical amino-acid substitution")
        if sequence[position - 1] != wt:
            raise ValueError("candidate WT does not match exact assay sequence")
        expected_variant = f"p.{protein_letters_1to3[wt]}{position}{protein_letters_1to3[mutant]}"
        if identity["protein_variant"] != expected_variant or candidate_id != f"FGFR2:{expected_variant}":
            raise ValueError("candidate identity does not match its amino-acid substitution")
        if identity["reference_sequence_sha256"] != sequence_hash:
            raise ValueError("candidate reference sequence hash mismatch")
        if not isinstance(identity["group_id"], str) or not identity["group_id"]:
            raise ValueError("candidate requires its evaluator-owned group identity")
        identities.append(identity)
    return sequence, identities


def _mutant_sequence(sequence: str, identity: dict) -> str:
    index = identity["position"] - 1
    mutant = sequence[:index] + identity["mut"] + sequence[index + 1:]
    differences = [index for index, pair in enumerate(zip(sequence, mutant)) if pair[0] != pair[1]]
    if differences != [identity["position"] - 1]:
        raise ValueError("mutant sequence must differ at exactly the requested residue")
    return mutant


def _validate_hidden(hidden: np.ndarray, batch: int, sequence_length: int) -> np.ndarray:
    hidden = np.asarray(hidden)
    if hidden.shape != (batch, sequence_length + 2, EMBEDDING_DIM):
        raise ValueError("model representation shape does not match full sequence, BOS and EOS")
    if not np.isfinite(hidden).all():
        raise ValueError("nonfinite model representation")
    return hidden.astype(np.float32, copy=False)


class _FrozenEncoder:
    def __init__(self, tokenizer, model, torch, device: str):
        self.tokenizer, self.model, self.torch, self.device = tokenizer, model, torch, device

    def encode(self, sequences: list[str]) -> np.ndarray:
        encoded = self.tokenizer(
            sequences, return_tensors="pt", padding=True, truncation=False,
            add_special_tokens=True,
        )
        token_ids = encoded["input_ids"].cpu().numpy()
        for row, sequence in enumerate(sequences):
            expected = [self.tokenizer.cls_token_id]
            expected.extend(self.tokenizer.convert_tokens_to_ids(list(sequence)))
            expected.append(self.tokenizer.eos_token_id)
            if list(token_ids[row, :len(expected)]) != expected:
                raise ValueError("tokenizer failed exact amino-acid/BOS/EOS alignment")
            if encoded["attention_mask"][row].sum().item() != len(expected):
                raise ValueError("tokenizer truncated or changed the full assay sequence")
        encoded = {name: tensor.to(self.device) for name, tensor in encoded.items()}
        with self.torch.inference_mode():
            hidden = self.model(**encoded).last_hidden_state
        return hidden.detach().to(device="cpu", dtype=self.torch.float32).numpy()


def _load_encoder(*, device: str, cache_dir: Path | None):
    # Lazy dependencies keep ordinary project tests independent of a model runtime.
    import torch
    from huggingface_hub import snapshot_download
    from transformers import AutoModel, AutoTokenizer

    if device not in {"cpu", "mps"}:
        raise ValueError("device must be cpu or mps")
    if device == "mps" and not torch.backends.mps.is_available():
        raise ValueError("MPS was requested but is unavailable")
    snapshot = Path(snapshot_download(
        repo_id=MODEL_ID, revision=MODEL_REVISION, allow_patterns=list(CHECKPOINT_FILENAMES),
        cache_dir=str(cache_dir) if cache_dir else None,
    ))
    files = {}
    for name in CHECKPOINT_FILENAMES:
        path = snapshot / name
        if not path.is_file():
            raise ValueError(f"missing pinned checkpoint file: {name}")
        files[name] = {"path": str(path.absolute()), "sha256": _sha256(path)}
    tokenizer = AutoTokenizer.from_pretrained(
        str(snapshot), local_files_only=True, trust_remote_code=False,
    )
    model = AutoModel.from_pretrained(
        str(snapshot), local_files_only=True, trust_remote_code=False,
        use_safetensors=True, add_pooling_layer=False,
    )
    if model.config.hidden_size != EMBEDDING_DIM or model.config.num_hidden_layers != 12:
        raise ValueError("checkpoint architecture does not match pinned ESM2-35M")
    model.requires_grad_(False)
    model.to(device)
    model.eval()
    packages = {}
    for package in ("torch", "transformers", "huggingface-hub", "safetensors", "numpy"):
        packages[package] = importlib.metadata.version(package)
    provenance = {"files": files, "packages": packages, "device": device,
                  "python": platform.python_version(), "machine": platform.machine()}
    return _FrozenEncoder(tokenizer, model, torch, device), provenance


@dataclass
class SequenceFeatureCache:
    candidate_ids: tuple[str, ...]
    wt: np.ndarray
    delta: np.ndarray
    manifest: dict
    _row_by_id: dict[str, int] = field(init=False, repr=False)

    def __post_init__(self):
        self._row_by_id = {candidate_id: row for row, candidate_id in enumerate(self.candidate_ids)}

    def features(self, candidate_id: str) -> dict[str, np.ndarray]:
        row = self._row_by_id[candidate_id]
        return {"wt": self.wt[row].copy(), "delta": self.delta[row].copy()}

    def matrix(self, candidate_ids: list[str], representation: str = "delta") -> np.ndarray:
        rows = [self._row_by_id[candidate_id] for candidate_id in candidate_ids]
        if representation == "delta":
            return self.delta[rows].copy()
        if representation == "wt_delta":
            return np.concatenate((self.wt[rows], self.delta[rows]), axis=1)
        raise ValueError("representation must be delta or wt_delta")


def extract_features(
    candidates: list[dict], references_path: Path, output_dir: Path, *,
    device: str = "cpu", batch_size: int = 4, cache_dir: Path | None = None,
) -> SequenceFeatureCache:
    """Extract only the evaluator-supplied candidate list, or verify an existing cache."""
    sequence, identities = _inputs(candidates, references_path)
    output_dir = Path(output_dir)
    if type(batch_size) is not int or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    if device not in {"cpu", "mps"}:
        raise ValueError("device must be cpu or mps")
    if (output_dir / "features.npz").exists() or (output_dir / "manifest.json").exists():
        if not all((output_dir / name).is_file() for name in ("features.npz", "manifest.json")):
            raise ValueError("incomplete existing cache; refusing to overwrite it")
        return load_features(candidates, references_path, output_dir)
    started = time.perf_counter()
    encoder, provenance = _load_encoder(device=device, cache_dir=cache_dir)
    hidden_wt = _validate_hidden(encoder.encode([sequence]), 1, len(sequence))[0]
    positions = np.asarray([identity["position"] for identity in identities], dtype=np.int64)
    wt = hidden_wt[positions].copy()
    delta = np.empty_like(wt)
    for start in range(0, len(identities), batch_size):
        chunk = identities[start:start + batch_size]
        mutant_sequences = [_mutant_sequence(sequence, identity) for identity in chunk]
        hidden_mutant = _validate_hidden(encoder.encode(mutant_sequences), len(chunk), len(sequence))
        for row, identity in enumerate(chunk):
            delta[start + row] = hidden_mutant[row, identity["position"]] - wt[start + row]
        print(f"ESM2 FGFR2 features: {min(start + batch_size, len(identities))}/{len(identities)}",
              file=sys.stderr, flush=True)
    if not np.isfinite(delta).all():
        raise ValueError("nonfinite signed difference representation")
    candidate_ids = tuple(identity["candidate_id"] for identity in identities)
    output_dir.mkdir(parents=True, exist_ok=True)
    features_path = output_dir / "features.npz"
    # Open with exclusive creation so a concurrent extraction cannot overwrite a cache.
    with features_path.open("xb") as file:
        np.savez_compressed(file, candidate_ids=np.asarray(candidate_ids), wt=wt, delta=delta)
    manifest = {
        "schema_version": SCHEMA_VERSION, "model_id": MODEL_ID,
        "model_revision": MODEL_REVISION, "embedding_dim": EMBEDDING_DIM,
        "representation": REPRESENTATION, "outcomes_used": False,
        "model_frozen": True, "model_eval": True, "inference_mode": True,
        "trust_remote_code": False, "weight_format": "safetensors", "dtype": "float32",
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "sequence_length": len(sequence), "references_file_sha256": _sha256(Path(references_path)),
        "candidate_identities": identities, "candidate_identity_sha256": _canonical_hash(identities),
        "candidate_row_mapping": {candidate_id: row for row, candidate_id in enumerate(candidate_ids)},
        "features_sha256": _sha256(features_path), "checkpoint": provenance,
        "batch_size": batch_size, "created_at": datetime.now(UTC).isoformat(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    with (output_dir / "manifest.json").open("x") as file:
        json.dump(manifest, file, indent=2, sort_keys=True)
        file.write("\n")
    return load_features(candidates, references_path, output_dir)


def load_features(
    candidates: list[dict], references_path: Path, output_dir: Path, *,
    verify_checkpoint_files: bool = True,
) -> SequenceFeatureCache:
    """Fail closed on input, checkpoint, feature-array or provenance mismatches."""
    sequence, identities = _inputs(candidates, references_path)
    output_dir = Path(output_dir)
    manifest = json.loads((output_dir / "manifest.json").read_text())
    expected = {
        "schema_version": SCHEMA_VERSION, "model_id": MODEL_ID,
        "model_revision": MODEL_REVISION, "embedding_dim": EMBEDDING_DIM,
        "representation": REPRESENTATION, "outcomes_used": False, "model_frozen": True,
        "model_eval": True, "inference_mode": True, "trust_remote_code": False,
        "weight_format": "safetensors", "dtype": "float32",
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "sequence_length": len(sequence), "references_file_sha256": _sha256(Path(references_path)),
    }
    if any(manifest.get(name) != value for name, value in expected.items()):
        raise ValueError("cache provenance does not match the frozen representation")
    cached_identities = manifest.get("candidate_identities", [])
    current_map = {identity["candidate_id"]: identity for identity in identities}
    cached_map = {identity["candidate_id"]: identity for identity in cached_identities}
    if (len(cached_map) != len(cached_identities) or cached_map != current_map
            or manifest.get("candidate_identity_sha256") != _canonical_hash(cached_identities)):
        raise ValueError("cache candidate identities differ from the evaluator-approved list")
    candidate_ids = tuple(identity["candidate_id"] for identity in cached_identities)
    if manifest.get("candidate_row_mapping") != {
        candidate_id: row for row, candidate_id in enumerate(candidate_ids)
    }:
        raise ValueError("cache candidate row mapping is inconsistent")
    checkpoint = manifest.get("checkpoint", {})
    if (set(checkpoint.get("files", {})) != set(CHECKPOINT_FILENAMES)
            or not checkpoint.get("packages") or checkpoint.get("device") not in {"cpu", "mps"}):
        raise ValueError("missing checkpoint provenance")
    if verify_checkpoint_files:
        for entry in checkpoint["files"].values():
            path = Path(entry["path"])
            if not path.is_file() or _sha256(path) != entry["sha256"]:
                raise ValueError("checkpoint file hash mismatch")
    features_path = output_dir / "features.npz"
    if _sha256(features_path) != manifest.get("features_sha256"):
        raise ValueError("feature file hash mismatch")
    with np.load(features_path, allow_pickle=False) as contents:
        if set(contents.files) != {"candidate_ids", "wt", "delta"}:
            raise ValueError("unexpected cached feature arrays")
        if tuple(contents["candidate_ids"].tolist()) != candidate_ids:
            raise ValueError("cache candidate array mapping mismatch")
        wt, delta = contents["wt"].copy(), contents["delta"].copy()
    for values in (wt, delta):
        if values.shape != (len(candidate_ids), EMBEDDING_DIM) or values.dtype != np.float32:
            raise ValueError("cached feature shape or dtype mismatch")
        if not np.isfinite(values).all():
            raise ValueError("cached feature array is nonfinite")
    return SequenceFeatureCache(candidate_ids, wt, delta, manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True,
                        help="JSON list already filtered by the evaluator")
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", choices=("cpu", "mps"), default="cpu")
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--cache-dir", type=Path)
    args = parser.parse_args()
    cache = extract_features(
        json.loads(args.candidates.read_text()), args.references, args.output,
        device=args.device, batch_size=args.batch_size, cache_dir=args.cache_dir,
    )
    print(json.dumps({"candidates": len(cache.candidate_ids), "embedding_dim": EMBEDDING_DIM,
                      "cache": str(args.output.absolute()),
                      "elapsed_seconds": cache.manifest["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
