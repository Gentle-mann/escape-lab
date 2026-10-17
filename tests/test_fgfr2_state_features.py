"""Paired-state geometry must respect identities, shared atoms, and invariances."""

import copy
import hashlib
import json

import numpy as np
import pytest
from Bio.Data.IUPACData import protein_letters_1to3

from escape_lab import fgfr2_state_features as module


def sequence():
    residues = list("A" * 700)
    for position, residue in {491: "C", 534: "E", 549: "N", 565: "E",
                              641: "K", 644: "D", 645: "F", 646: "G"}.items():
        residues[position - 1] = residue
    return "".join(residues)


def identity(position=10):
    return {"candidate_id": f"FGFR2:p.Ala{position}Lys", "gene": "FGFR2",
            "protein_variant": f"p.Ala{position}Lys", "position": position,
            "wt": "A", "mut": "K", "group_id": f"block-{position}",
            "split": "development",
            "reference_sequence_sha256": hashlib.sha256(sequence().encode()).hexdigest()}


def atom_states():
    positions = (10, 15, 30, 491, 534, 549, 565, 641, 644, 645, 646)
    inactive = {}
    for index, position in enumerate(positions):
        ca = np.array([index * 1.5, index % 3, index % 2], dtype=float)
        inactive[position] = {"CA": ca, "N": ca + [0.2, 0.1, 0.4]}
    active = copy.deepcopy(inactive)
    active[15]["CA"] += [4.0, 0.2, 0.1]
    active[15]["N"] += [3.0, 0.2, 0.1]
    return {"inactive": inactive, "active": active}


def contexts(states):
    return module.paired_contexts(states, sequence(), excluded_positions={491})[0]


def test_identical_states_zero_contrast_and_all_arms_share_mask():
    states = atom_states()
    states["active"] = copy.deepcopy(states["inactive"])
    context = contexts(states)
    for arm in module.ARMS:
        values = module.mutation_features(identity(), context, arm)
        assert len(values) == 31
        assert values["state_paired_available"] == 1
        if arm == "contrast":
            assert all(value == 0 for key, value in values.items()
                       if key != "state_paired_available")
    missing = module.mutation_features(identity(400), context, "contrast")
    assert len(missing) == 31 and set(missing.values()) == {0.0}


def test_swapping_states_reverses_contrast_but_not_availability():
    states = atom_states()
    original = module.mutation_features(identity(), contexts(states), "contrast")
    swapped = {"inactive": states["active"], "active": states["inactive"]}
    reversed_values = module.mutation_features(identity(), contexts(swapped), "contrast")
    assert any(abs(value) > 1e-5 for key, value in original.items()
               if key != "state_paired_available")
    for key in original:
        assert reversed_values[key] == pytest.approx(
            original[key] if key == "state_paired_available" else -original[key])


def test_independent_rigid_rotations_leave_all_geometry_unchanged():
    states = atom_states()
    transformed = copy.deepcopy(states)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    for atoms in transformed["active"].values():
        for name, coordinate in atoms.items():
            atoms[name] = coordinate @ rotation + [100, -10, 8]
    for atoms in transformed["inactive"].values():
        for name, coordinate in atoms.items():
            atoms[name] = coordinate @ rotation.T + [-20, 80, 11]
    original, moved = contexts(states), contexts(transformed)
    for arm in module.ARMS:
        assert module.mutation_features(identity(), moved, arm) == pytest.approx(
            module.mutation_features(identity(), original, arm), abs=1e-12)


def test_only_common_residues_and_common_atom_names_contribute():
    states = atom_states()
    original = contexts(states)
    states["active"][400] = {"CA": np.zeros(3), "N": np.zeros(3)}
    states["active"][10]["P"] = np.zeros(3)  # No counterpart, never a shared atom.
    states["inactive"][491]["CA"] = np.zeros(3)  # Engineered site is excluded.
    changed = contexts(states)
    assert module.mutation_features(identity(), changed, "contrast") == pytest.approx(
        module.mutation_features(identity(), original, "contrast"))


def test_heavy_contacts_respond_to_sidechain_motion_with_same_ca():
    states = atom_states()
    states["active"] = copy.deepcopy(states["inactive"])
    states["inactive"][15]["CB"] = np.array([1.0, 0, 0])
    states["active"][15]["CB"] = np.array([15.0, 0, 0])
    values = module.mutation_features(identity(), contexts(states), "contrast")
    assert values["state_context_degree"] == 0
    assert abs(values["state_heavy_degree"]) > 1e-4


def _pdb_text(atoms_by_position, residue_overrides=None):
    text = []
    serial = 1
    overrides = residue_overrides or {}
    for position, atoms in sorted(atoms_by_position.items()):
        residue = overrides.get(position, protein_letters_1to3[sequence()[position - 1]].upper())
        for name, point in atoms.items():
            element = name[0]
            text.append(f"ATOM  {serial:5d} {name:^4s} {residue:>3s} A{position:4d}    "
                        f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}"
                        f"{1:6.2f}{20:6.2f}          {element:>2s}  ")
            serial += 1
    return "\n".join([*text, "TER", "END", ""])


@pytest.fixture
def prepared(tmp_path):
    reference = tmp_path / "references.json"
    reference.write_text(json.dumps({"sequences": {"FGFR2": sequence()},
                                    "sequence_sha256": {"FGFR2": identity()[
                                        "reference_sequence_sha256"]}}))
    states = {}
    for name, atoms in atom_states().items():
        path = tmp_path / f"{name}.pdb"
        path.write_text(_pdb_text(atoms, {491: "ALA"}))
        states[name] = {"path": str(path), "sha256": module.file_hash(path),
                        "format": "pdb", "chain": "A", "excluded_positions": [491],
                        "author_to_assay": {str(position): position for position in atoms}}
    spec = tmp_path / "state-spec.json"
    spec.write_text(json.dumps({"states": states}))
    return spec, reference, tmp_path / "cache"


def test_build_load_and_exact_identity_allowlist(prepared):
    spec, reference, output = prepared
    manifest = module.build_state_features(spec, reference, [identity(), identity(400)], output)
    assert manifest["outcomes_used"] is False
    assert manifest["coverage"]["paired_candidates"] == 1
    cache = module.FGFR2StateFeatureStore(output)
    for arm in module.ARMS:
        assert cache.features(identity(), arm)["state_paired_available"] == 1
        assert cache.features(identity(400), arm)["state_paired_available"] == 0
    altered = identity()
    altered["group_id"] = "other"
    with pytest.raises(ValueError, match="identity"):
        cache.features(altered, "active")


@pytest.mark.parametrize("field,value", [("wt", "C"), ("split", "final"),
                                         ("gene", "FGFR1"), ("position", 0)])
def test_invalid_candidates_rejected_before_cache(prepared, field, value):
    spec, reference, output = prepared
    candidate = identity()
    candidate[field] = value
    with pytest.raises(ValueError):
        module.build_state_features(spec, reference, [candidate], output)
    assert not output.exists()


def test_unexcluded_mismatch_rejected(prepared):
    spec, reference, output = prepared
    contents = json.loads(spec.read_text())
    contents["states"]["active"]["excluded_positions"] = []
    spec.write_text(json.dumps(contents))
    with pytest.raises(ValueError, match="residue identity"):
        module.build_state_features(spec, reference, [identity()], output)


def test_noninjective_mapping_rejected(prepared):
    spec, reference, output = prepared
    contents = json.loads(spec.read_text())
    contents["states"]["active"]["author_to_assay"]["15"] = 10
    spec.write_text(json.dumps(contents))
    with pytest.raises(ValueError, match="injective"):
        module.build_state_features(spec, reference, [identity()], output)


@pytest.mark.parametrize("which", ["cache", "source", "reference", "spec"])
def test_changed_provenance_or_cache_rejected(prepared, which):
    spec, reference, output = prepared
    module.build_state_features(spec, reference, [identity()], output)
    paths = {"cache": output / "feature-cache.json", "source": spec.parent / "active.pdb",
             "reference": reference, "spec": spec}
    with paths[which].open("a") as handle:
        handle.write("\n ")
    with pytest.raises(ValueError, match="hash"):
        module.FGFR2StateFeatureStore(output)


def test_reference_hash_mismatch_rejected(prepared):
    spec, reference, output = prepared
    contents = json.loads(reference.read_text())
    contents["sequence_sha256"]["FGFR2"] = "bad"
    reference.write_text(json.dumps(contents))
    with pytest.raises(ValueError, match="hash"):
        module.build_state_features(spec, reference, [identity()], output)


def test_source_hash_mismatch_rejected(prepared):
    spec, reference, output = prepared
    with (spec.parent / "active.pdb").open("a") as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="hash"):
        module.build_state_features(spec, reference, [identity()], output)


def test_engineered_site_is_neither_query_geometry_nor_context(prepared):
    spec, reference, output = prepared
    candidate = identity(491)
    candidate.update({"wt": "C", "protein_variant": "p.Cys491Lys",
                      "candidate_id": "FGFR2:p.Cys491Lys"})
    module.build_state_features(spec, reference, [candidate], output)
    cache = module.FGFR2StateFeatureStore(output)
    for arm in module.ARMS:
        assert set(cache.features(candidate, arm).values()) == {0.0}


def test_missing_critical_anchor_fails_instead_of_inventing_distance():
    states = atom_states()
    del states["active"][534]
    with pytest.raises(ValueError, match="anchor missing"):
        contexts(states)


def test_reference_anchor_identity_is_checked():
    bad_reference = sequence()[:533] + "A" + sequence()[534:]
    with pytest.raises(ValueError, match="anchor differs"):
        module.paired_contexts(atom_states(), bad_reference, excluded_positions={491})


def test_nonfinite_common_coordinates_rejected():
    states = atom_states()
    states["active"][10]["N"][0] = float("nan")
    with pytest.raises(ValueError, match="coordinates"):
        contexts(states)


def test_phosphate_atoms_never_enter_common_heavy_contacts():
    states = atom_states()
    original = contexts(states)
    states["inactive"][15]["O1P"] = np.zeros(3)
    states["active"][15]["O1P"] = np.array([99.0, 99.0, 99.0])
    changed = contexts(states)
    for arm in module.ARMS:
        assert module.mutation_features(identity(), changed, arm) == pytest.approx(
            module.mutation_features(identity(), original, arm))


def test_unknown_arm_rejected(prepared):
    spec, reference, output = prepared
    module.build_state_features(spec, reference, [identity()], output)
    with pytest.raises(ValueError, match="arm"):
        module.FGFR2StateFeatureStore(output).features(identity(), "chosen-after-results")
