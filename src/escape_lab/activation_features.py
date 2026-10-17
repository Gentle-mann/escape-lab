"""Outcome-independent additions for the development-only activation pilot.

Context and motif proximity use exact assay sequences and their existing FGFR1
alignment. Proximity is along the sequence, not a simulated structural distance.
AlphaMissense and EVE are pretrained constraint/pathogenicity inputs, not labels
or signed activation estimates. No source assay outcomes or annotations are read.
"""

import csv
import hashlib
import json
import math
from pathlib import Path

from Bio.SeqUtils.ProtParamData import kd

from escape_lab.biological_features import AMINO_ACIDS, GENES, SCHEMA_VERSION
from escape_lab.data import aligned_position_map, parse_protein_variant

SEQUENCE_ANCHORS = {
    "brake_n546": (546, "N"),
    "brake_e562": (562, "E"),
    "brake_k638": (638, "K"),
    "activation_loop_y653": (653, "Y"),
    "activation_loop_y654": (654, "Y"),
}
SCORE_NAMES = ("alphamissense", "eve")
SCORE_STATUSES = frozenset({"available", "missing", "discordant_not_collapsed"})
IDENTITY_FIELDS = (
    "candidate_id", "gene", "protein_variant", "position", "aligned_position", "wt", "mut",
    "split",
)


class ActivationFeatureAugmenter:
    """Register permitted identities, then serve only outcome-independent additions.

    The caller owns the permitted development candidate list, including external
    group exclusions. CSV identity allowlisting occurs before score/status parsing.
    The base store is retained for interface composition; its features are not
    returned or mutated here. All feature preparation is outcome-independent.
    """

    def __init__(self, base_store, candidates, references_path, scores_path):
        self.base_store = base_store
        references = json.loads(Path(references_path).read_text())
        if references.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("Unsupported activation reference schema")
        sequences = references.get("sequences", {})
        hashes = references.get("sequence_sha256", {})
        if set(sequences) != set(GENES) or set(hashes) != set(GENES):
            raise ValueError("Activation references require all exact assay genes")
        for gene, sequence in sequences.items():
            if (not isinstance(sequence, str) or not sequence
                    or any(residue not in AMINO_ACIDS for residue in sequence)
                    or hashlib.sha256(sequence.encode()).hexdigest() != hashes[gene]):
                raise ValueError("Invalid activation reference sequence/hash")
        for name, (position, residue) in SEQUENCE_ANCHORS.items():
            reference = sequences["FGFR1"]
            if position > len(reference) or reference[position - 1] != residue:
                raise ValueError(f"Activation sequence anchor does not match reference: {name}")
        if base_store is not None and getattr(base_store, "_sequences", sequences) != sequences:
            raise ValueError("Base feature store and activation references disagree")
        self._sequences = dict(sequences)
        self._alignment = {
            gene: aligned_position_map(sequences["FGFR1"], sequence)
            for gene, sequence in sequences.items()
        }
        self._candidates = {}
        for candidate in candidates:
            identity = self._validated_identity(candidate)
            candidate_id = identity["candidate_id"]
            if candidate_id in self._candidates:
                raise ValueError("Duplicate activation candidate identity")
            self._candidates[candidate_id] = identity
        if not self._candidates:
            raise ValueError("Activation features require eligible development candidates")
        self._scores = self._read_scores(scores_path)
        self.score_coverage = {
            score: {
                status: sum(values[score][1] == status for values in self._scores.values())
                for status in sorted(SCORE_STATUSES)
            }
            for score in SCORE_NAMES
        }

    def _validated_identity(self, candidate):
        if any(field not in candidate for field in IDENTITY_FIELDS):
            raise ValueError("Incomplete activation candidate identity")
        identity = {field: candidate[field] for field in IDENTITY_FIELDS}
        gene, position = identity["gene"], identity["position"]
        wt, mut = identity["wt"], identity["mut"]
        if identity["split"] not in {"development", "validation"}:
            raise ValueError("Activation features only permit development/validation candidates")
        if (not isinstance(gene, str) or gene not in self._sequences
                or type(position) is not int or not 1 <= position <= len(self._sequences[gene])
                or not isinstance(wt, str) or len(wt) != 1 or wt not in AMINO_ACIDS
                or not isinstance(mut, str) or len(mut) != 1 or mut not in AMINO_ACIDS
                or type(identity["aligned_position"]) is not int):
            raise ValueError("Invalid activation candidate identity")
        if self._sequences[gene][position - 1] != wt:
            raise ValueError("Activation candidate WT differs from exact assay reference")
        if identity["aligned_position"] != self._alignment[gene][position]:
            raise ValueError("Activation candidate aligned position differs from assay mapping")
        if (not isinstance(identity["protein_variant"], str)
                or parse_protein_variant(identity["protein_variant"]) != (wt, position, mut)
                or identity["candidate_id"] != f"{gene}:{identity['protein_variant']}"):
            raise ValueError("Activation candidate protein name differs from reference identity")
        return identity

    def _read_scores(self, path):
        result = {
            candidate_id: {name: (math.nan, "missing") for name in SCORE_NAMES}
            for candidate_id in self._candidates
        }
        required = {"candidate_id_name_match", "gene", "protein_variant"}
        for name in SCORE_NAMES:
            required.update({name, f"{name}_status"})
        seen = set()
        with Path(path).open(newline="") as handle:
            reader = csv.DictReader(handle)
            if not required.issubset(reader.fieldnames or []):
                raise ValueError("Missing activation predictor score columns")
            for row in reader:
                candidate_id = row["candidate_id_name_match"]
                # This guard deliberately precedes every numeric or status parse.
                if candidate_id not in self._candidates:
                    continue
                if candidate_id in seen:
                    raise ValueError("Duplicate allowed predictor score row")
                seen.add(candidate_id)
                identity = self._candidates[candidate_id]
                protein = f"{identity['wt']}{identity['position']}{identity['mut']}"
                if row["gene"] != identity["gene"] or row["protein_variant"] != protein:
                    raise ValueError("Predictor score identity differs from exact assay candidate")
                for name in SCORE_NAMES:
                    status = row[f"{name}_status"]
                    if status not in SCORE_STATUSES:
                        raise ValueError(f"Unknown {name} score status")
                    value = math.nan
                    if status == "available":
                        value = float(row[name])
                        if not math.isfinite(value) or not 0 <= value <= 1:
                            raise ValueError(f"Invalid {name} probability score")
                    result[candidate_id][name] = value, status
        return result

    def features(self, candidate):
        identity = self._validated_identity(candidate)
        candidate_id = identity["candidate_id"]
        if self._candidates.get(candidate_id) != identity:
            raise ValueError("Query identity is not an approved activation candidate")
        gene, position = identity["gene"], identity["position"]
        aligned = identity["aligned_position"]
        wt, mut = identity["wt"], identity["mut"]
        sequence = self._sequences[gene]
        result = {"activation_aligned_position": aligned / 1000}
        for offset in range(-3, 4):
            name = f"minus{-offset}" if offset < 0 else f"plus{offset}"
            index = position - 1 + offset
            result[f"activation_context_{name}"] = sequence[index] if 0 <= index < len(sequence) else "PAD"
        for center in range(450, 811, 20):
            result[f"activation_position_rbf_{center}"] = math.exp(-0.5 * ((aligned - center) / 20) ** 2)
        deltas = {
            "charge_delta": (int(mut in "KR") - int(mut in "DE"))
            - (int(wt in "KR") - int(wt in "DE")),
            "hydropathy_delta": (kd[mut] - kd[wt]) / 10,
            "proline_delta": int(mut == "P") - int(wt == "P"),
        }
        for name, (anchor, _) in SEQUENCE_ANCHORS.items():
            proximity = math.exp(-abs(aligned - anchor) / 10)
            feature = f"activation_{name}_sequence_proximity"
            result[feature] = proximity
            for delta, value in deltas.items():
                result[f"{feature}_x_{delta}"] = proximity * value
        # Substitution heuristic only: no claim of phosphorylation or faithful mimicry.
        result["activation_st_to_de_substitution"] = int(wt in "ST" and mut in "DE")
        result["activation_tyrosine_loss"] = int(wt == "Y" and mut != "Y")
        for name, (value, status) in self._scores[candidate_id].items():
            result[f"activation_{name}"] = value
            result[f"activation_{name}_missing"] = int(status != "available")
            result[f"activation_{name}_conflict"] = int(status == "discordant_not_collapsed")
        return result
