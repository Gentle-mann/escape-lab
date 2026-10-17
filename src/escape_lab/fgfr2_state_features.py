"""Outcome-free paired FGFR2 conformation descriptors, not physical energies.

Both contexts use identical positions and heavy-atom names. Mutation chemistry
interacts with fixed reference geometry; no mutant relaxation is simulated.
Author-provided modeled structures are state proxies, not measured populations.
"""

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParamData import kd

from escape_lab.fgfr2_sequence_features import IDENTITY_FIELDS, _inputs

SCHEMA_VERSION = "fgfr2-paired-state-features-v1"
ARMS = ("inactive", "active", "contrast")
ANCHORS = {"brake": {549: "N", 565: "E", 641: "K"},
           "dfg": {644: "D", 645: "F", 646: "G"}, "alpha_c": {534: "E"}}
CANONICAL = {"PTR": "Y", "HID": "H", "HIE": "H", "HIP": "H",
             "HSD": "H", "HSE": "H", "HSP": "H", "CYX": "C"}
SIDECHAIN_HEAVY_ATOMS = dict(zip("ACDEFGHIKLMNPQRSTVWY",
                               (1, 2, 4, 5, 7, 0, 6, 4, 5, 4, 4, 4, 3, 5, 7, 2, 3, 3, 10, 8),
                               strict=True))
DONORS = frozenset("KRNHQSTWY")
ACCEPTORS = frozenset("DENQSTY")
CONTEXT_NAMES = ("degree", "hydropathy", "charge", "aromatic", "donor", "acceptor")
CONTEXT_KEYS = (*CONTEXT_NAMES, "heavy_degree", "heavy_hydropathy", "heavy_charge",
                "anchor_brake", "anchor_dfg", "anchor_alpha_c")
PHOSPHATE_ATOMS = frozenset({"P", "O1P", "O2P", "O3P", "OP1", "OP2", "OP3"})


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _write_json(path, value):
    with Path(path).open("x") as handle:
        json.dump(value, handle, sort_keys=True, separators=(",", ":"), allow_nan=False)
        handle.write("\n")


def _properties(residue):
    return {"degree": 1.0, "hydropathy": kd[residue] / 4.5,
            "charge": float(int(residue in "KR") - int(residue in "DE")),
            "aromatic": float(residue in "FWY"), "donor": float(residue in DONORS),
            "acceptor": float(residue in ACCEPTORS)}


def _load_state(specification, sequence):
    path = Path(specification["path"]).resolve()
    if file_hash(path) != specification["sha256"]:
        raise ValueError("state source hash mismatch")
    mapping = {}
    for author, position in specification["author_to_assay"].items():
        if (not isinstance(author, str) or str(int(author)) != author
                or type(position) is not int or not 1 <= position <= len(sequence)):
            raise ValueError("invalid author-to-assay mapping")
        mapping[int(author)] = position
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("author-to-assay mapping must be injective")
    excluded = set(specification.get("excluded_positions", []))
    if any(type(position) is not int or not 1 <= position <= len(sequence)
           for position in excluded):
        raise ValueError("invalid excluded reference position")
    file_format = specification.get("format")
    if file_format not in {"pdb", "cif"}:
        raise ValueError("state format must be pdb or cif")
    parser = PDBParser(QUIET=True) if file_format == "pdb" else MMCIFParser(QUIET=True)
    structure = parser.get_structure("state", str(path))
    if len(structure) != 1:
        raise ValueError("exactly one prepared model is required")
    chain_name = specification["chain"]
    if chain_name not in structure[0]:
        raise ValueError("declared state chain is unavailable")
    output = {}
    for residue in structure[0][chain_name]:
        author = residue.id[1]
        if author not in mapping or "CA" not in residue:
            continue
        if residue.id[2] != " ":
            raise ValueError("insertion codes require explicit separate mapping")
        position = mapping[author]
        canonical = seq1(residue.resname, custom_map=CANONICAL)
        if position in excluded:
            continue
        if canonical != sequence[position - 1]:
            raise ValueError("state residue identity differs from exact assay reference")
        if position in output:
            raise ValueError("duplicate mapped residue geometry")
        atoms = {}
        for atom in residue.get_atoms():
            if atom.element.upper() in {"H", "D", "P"} or atom.name in PHOSPHATE_ATOMS:
                continue
            coordinate = np.asarray(atom.coord, dtype=float)
            if coordinate.shape != (3,) or not np.isfinite(coordinate).all():
                raise ValueError("invalid state atomic coordinate")
            atoms[atom.name] = coordinate
        if "CA" in atoms:
            output[position] = atoms
    if not output:
        raise ValueError("state has no accepted reference geometry")
    return output, excluded


def paired_contexts(states, sequence, excluded_positions=()):
    """Return matched position contexts plus geometry coverage, with no outcomes.

    Input states map assay positions to atom-name/XYZ dictionaries. Contacts
    exclude sequence neighbors within two residues. All sums divide by a fixed
    10, never a data-derived neighborhood normalization. Distances are Angstroms.
    """
    if set(states) != {"inactive", "active"}:
        raise ValueError("exactly inactive and active states are required")
    common = sorted(set(states["inactive"]) & set(states["active"])
                    - set(excluded_positions))
    if not common or any(type(position) is not int or not 1 <= position <= len(sequence)
                         for position in common):
        raise ValueError("paired states require valid common assay positions")
    for anchor in ANCHORS.values():
        for position, residue in anchor.items():
            if position > len(sequence) or sequence[position - 1] != residue:
                raise ValueError("structural anchor differs from exact assay reference")
            if position not in common:
                raise ValueError("structural anchor missing from paired geometry")
    atom_names = {}
    for position in common:
        names = sorted(set(states["inactive"][position]) & set(states["active"][position])
                       - PHOSPHATE_ATOMS)
        if "CA" not in names:
            raise ValueError("every paired residue requires a common CA atom")
        atom_names[position] = names
    indices = {position: index for index, position in enumerate(common)}
    properties = {name: np.array([_properties(sequence[position - 1])[name]
                                 for position in common]) for name in CONTEXT_NAMES}
    separation = np.abs(np.subtract.outer(common, common)) > 2
    output = {}
    for state_name, state in states.items():
        xyz = np.asarray([state[position]["CA"] for position in common], dtype=float)
        if xyz.shape != (len(common), 3) or not np.isfinite(xyz).all():
            raise ValueError("invalid common CA coordinates")
        ca_distances = np.linalg.norm(xyz[:, None] - xyz[None, :], axis=-1)
        ca_contacts = 1 / (1 + (ca_distances / 8) ** 6) * separation
        heavy_contacts = np.zeros((len(common), len(common)))
        heavy_coordinates = []
        for position in common:
            atoms = np.asarray([state[position][name] for name in atom_names[position]], dtype=float)
            if atoms.shape != (len(atom_names[position]), 3) or not np.isfinite(atoms).all():
                raise ValueError("invalid common heavy coordinates")
            heavy_coordinates.append(atoms)
        for i, atoms in enumerate(heavy_coordinates):
            for j in range(i + 1, len(common)):
                if not separation[i, j]:
                    continue
                distance = float(np.linalg.norm(
                    atoms[:, None] - heavy_coordinates[j][None, :], axis=-1).min())
                contact = 1 / (1 + (distance / 5) ** 6)
                heavy_contacts[i, j] = heavy_contacts[j, i] = contact
        state_contexts = {}
        for i, position in enumerate(common):
            values = {name: float(ca_contacts[i] @ properties[name] / 10)
                      for name in CONTEXT_NAMES}
            for name in ("degree", "hydropathy", "charge"):
                values[f"heavy_{name}"] = float(heavy_contacts[i] @ properties[name] / 10)
            for name, anchor in ANCHORS.items():
                distance = min(ca_distances[i, indices[other]] for other in anchor)
                values[f"anchor_{name}"] = math.exp(-float(distance) / 8)
            state_contexts[position] = values
        output[state_name] = state_contexts
    coverage = {"common_positions": common, "common_position_count": len(common),
                "common_heavy_atom_count": sum(len(names) for names in atom_names.values()),
                "excluded_positions": sorted(set(excluded_positions))}
    return output, coverage


def _variant_values(context, wt, mutant):
    original, changed = _properties(wt), _properties(mutant)
    delta = {name: changed[name] - original[name] for name in CONTEXT_NAMES}
    delta["size_proxy"] = (SIDECHAIN_HEAVY_ATOMS[mutant] - SIDECHAIN_HEAVY_ATOMS[wt]) / 10
    delta["proline"] = float(mutant == "P") - float(wt == "P")
    delta["glycine"] = float(mutant == "G") - float(wt == "G")
    values = {f"state_context_{name}": context[name] for name in CONTEXT_NAMES}
    interactions = (("degree", "size_proxy"), ("degree", "proline"), ("degree", "glycine"),
                    ("hydropathy", "hydropathy"), ("charge", "charge"),
                    ("aromatic", "aromatic"), ("donor", "acceptor"), ("acceptor", "donor"))
    for name, chemical in interactions:
        values[f"state_{name}_x_{chemical}_delta"] = context[name] * delta[chemical]
    for name in ANCHORS:
        context_name = f"anchor_{name}"
        values[f"state_{context_name}"] = context[context_name]
        for chemical in ("hydropathy", "charge", "proline"):
            values[f"state_{context_name}_x_{chemical}_delta"] = (
                context[context_name] * delta[chemical])
    values["state_heavy_degree"] = context["heavy_degree"]
    for name, chemical in (("degree", "size_proxy"), ("hydropathy", "hydropathy"),
                           ("charge", "charge")):
        values[f"state_heavy_{name}_x_{chemical}_delta"] = (
            context[f"heavy_{name}"] * delta[chemical])
    return values


FEATURE_NAMES = tuple(sorted([*_variant_values(dict.fromkeys(CONTEXT_KEYS, 0.0), "A", "K"),
                              "state_paired_available"]))


def mutation_features(candidate, contexts, arm):
    """Use signed active minus inactive context with an identical availability flag."""
    if arm not in ARMS:
        raise ValueError("unknown paired-state feature arm")
    position = candidate["position"]
    present = position in contexts["inactive"] and position in contexts["active"]
    if not present:
        return dict.fromkeys(FEATURE_NAMES, 0.0)
    if arm == "contrast":
        context = {name: contexts["active"][position][name] - contexts["inactive"][position][name]
                   for name in CONTEXT_KEYS}
    else:
        context = contexts[arm][position]
    values = _variant_values(context, candidate["wt"], candidate["mut"])
    values["state_paired_available"] = 1.0
    return {name: float(values[name]) for name in FEATURE_NAMES}


def build_state_features(spec_path, references_path, candidates, output_dir):
    """Create a new immutable cache from approved identities and structural sources."""
    started = time.perf_counter()
    spec_path, references_path = Path(spec_path).resolve(), Path(references_path).resolve()
    sequence, identities = _inputs(list(candidates), references_path)
    specification = json.loads(spec_path.read_text())
    if set(specification["states"]) != {"inactive", "active"}:
        raise ValueError("state specification requires inactive and active")
    states = {}
    excluded = set(specification.get("excluded_positions", []))
    for name, state_spec in specification["states"].items():
        states[name], state_excluded = _load_state(state_spec, sequence)
        excluded.update(state_excluded)
    context, coverage = paired_contexts(states, sequence, excluded)
    cache = {identity["candidate_id"]: {
        "identity": identity,
        "features": {arm: mutation_features(identity, context, arm) for arm in ARMS}}
             for identity in identities}
    coverage["candidate_count"] = len(identities)
    coverage["paired_candidates"] = sum(
        int(row["features"]["contrast"]["state_paired_available"]) for row in cache.values())
    coverage["unpaired_candidates"] = len(identities) - coverage["paired_candidates"]
    contrasts = np.array([[row["features"]["contrast"][name] for name in FEATURE_NAMES]
                          for row in cache.values()])
    quantiles = {name: list(np.quantile(np.abs(contrasts[:, index]), [0, .25, .5, .75, .95, 1]))
                 for index, name in enumerate(FEATURE_NAMES) if name != "state_paired_available"}
    source_hashes = {str(spec_path): file_hash(spec_path),
                     str(references_path): file_hash(references_path),
                     str(Path(__file__).resolve()): file_hash(__file__)}
    for state_spec in specification["states"].values():
        source_hashes[str(Path(state_spec["path"]).resolve())] = state_spec["sha256"]
    manifest = {
        "schema_version": SCHEMA_VERSION, "outcomes_used": False,
        "arms": list(ARMS), "feature_names": list(FEATURE_NAMES), "feature_count": len(FEATURE_NAMES),
        "reference_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "source_sha256": source_hashes, "coverage": coverage,
        "absolute_raw_contrast_quantiles": quantiles,
        "quantile_levels": [0, .25, .5, .75, .95, 1],
        "conventions": {"ca_kernel_angstrom": 8, "heavy_kernel_angstrom": 5,
                        "kernel_power": 6, "sequence_neighbor_exclusion": 2,
                        "fixed_sum_divisor": 10, "anchor_exponential_scale_angstrom": 8,
                        "size_proxy": "sidechain-heavy-atom-count-change-divided-by-10",
                        "contrast": "active-minus-inactive",
                        "common_atom_names_only": True, "no_mutant_relaxation": True,
                        "donor_acceptor": "residue-propensity-flags-not-hydrogen-bond-counts"},
        "state_specification": specification,
        "elapsed_seconds": time.perf_counter() - started,
    }
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    _write_json(output_dir / "feature-cache.json", cache)
    manifest["feature_cache_sha256"] = file_hash(output_dir / "feature-cache.json")
    _write_json(output_dir / "manifest.json", manifest)
    return manifest


class FGFR2StateFeatureStore:
    """Hash-checked, identity-allowlisted geometry. No labels are accessible here."""

    def __init__(self, output_dir):
        output_dir = Path(output_dir)
        self.manifest = json.loads((output_dir / "manifest.json").read_text())
        manifest = self.manifest
        if (manifest.get("schema_version") != SCHEMA_VERSION
                or manifest.get("outcomes_used") is not False
                or manifest.get("arms") != list(ARMS)
                or manifest.get("feature_names") != list(FEATURE_NAMES)
                or manifest.get("feature_count") != len(FEATURE_NAMES)):
            raise ValueError("invalid paired-state cache provenance")
        for path, digest in manifest["source_sha256"].items():
            if file_hash(path) != digest:
                raise ValueError("paired-state source hash mismatch")
        path = output_dir / "feature-cache.json"
        if file_hash(path) != manifest["feature_cache_sha256"]:
            raise ValueError("paired-state feature cache hash mismatch")
        self._cache = json.loads(path.read_text())
        for candidate_id, row in self._cache.items():
            if row["identity"]["candidate_id"] != candidate_id or set(row["features"]) != set(ARMS):
                raise ValueError("malformed paired-state candidate identity")
            availability = []
            for values in row["features"].values():
                if set(values) != set(FEATURE_NAMES) or any(
                    type(value) not in (int, float) or not math.isfinite(value)
                    for value in values.values()
                ):
                    raise ValueError("invalid paired-state cached features")
                availability.append(values["state_paired_available"])
            if len(set(availability)) != 1 or availability[0] not in (0, 1):
                raise ValueError("paired-state availability differs between arms")

    def features(self, candidate, arm):
        if arm not in ARMS:
            raise ValueError("unknown paired-state feature arm")
        identity = {name: candidate.get(name) for name in IDENTITY_FIELDS}
        row = self._cache.get(identity["candidate_id"])
        if row is None or identity != row["identity"]:
            raise ValueError("candidate identity is not approved for paired-state cache")
        return dict(row["features"][arm])
