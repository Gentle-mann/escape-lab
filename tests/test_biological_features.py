import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Residue

from escape_lab.biological_features import (
    MISSING_DISTANCE,
    BiologicalFeatureStore,
    _chain_geometry,
    _sequence_mapping,
    build_feature_store,
)

ROOT = Path(__file__).resolve().parents[1]


def atom(name, coordinate, element="C", serial=1):
    return Atom.Atom(name, np.asarray(coordinate, dtype=float), 0, 1, " ", name, serial, element)


def test_geometry_uses_heavy_atoms_and_protein_neighbors_not_nearby_hydrogens():
    chain = Chain.Chain("A")
    for number, x in [(1, 0), (2, 6), (3, 12)]:
        residue = Residue.Residue((" ", number, " "), "ALA", " ")
        residue.add(atom("CA", [x, 0, 0]))
        residue.add(atom("CB", [x, 1, 0], serial=2))
        chain.add(residue)
    ligand = Residue.Residue(("H_LIG", 99, " "), "LIG", " ")
    ligand.add(atom("C1", [3, 0, 0]))
    ligand.add(atom("H1", [0, 0, 0.1], element="H", serial=2))
    geometry = _chain_geometry(chain, ligand, {"test": [1]})
    first = geometry[1]
    assert first["ligand_min_heavy_distance"] == pytest.approx(3.0)
    assert first["ligand_ca_distance"] == pytest.approx(3.0)
    assert first["ligand_atom_pairs_4a"] == 2
    assert first["ca_neighbors_8a"] == 1
    assert first["sasa"] > 0
    assert first["motif_test_ca_distance"] == 0
    assert geometry[2]["motif_test_ca_distance"] == pytest.approx(6.0)


def test_sequence_mapping_handles_insertions_without_filling_gaps_from_previous_position():
    template = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
    reference = "QQQ" + template[:20] + "PPP" + template[20:] + "QQQ"
    mapping, ambiguous = _sequence_mapping(reference, template)
    assert not ambiguous
    assert mapping[4] == 1
    assert mapping[23] == 20
    assert all(position not in mapping for position in [24, 25, 26])
    assert mapping[27] == 21
    assert mapping[3 + len(template) + 3] == len(template)


def test_ambiguous_repeated_domain_mapping_is_not_silently_chosen():
    domain = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
    mapping, ambiguous = _sequence_mapping(domain + "GGGGGGGG" + domain, domain)
    assert ambiguous
    assert not mapping


@pytest.fixture(scope="module")
def built_store(tmp_path_factory):
    raw = ROOT / "data/raw"
    structures = ROOT / "data/features-research"
    required = [raw / "FGFR1-pem.json", structures / "7WCL.cif", structures / "6MZW.cif"]
    if not all(path.exists() for path in required):
        pytest.skip("Optional local public reference metadata/coordinates are not downloaded")
    output = tmp_path_factory.mktemp("biological-features")
    manifest = build_feature_store(raw, structures, output)
    return output, manifest, BiologicalFeatureStore(output)


def candidate(gene="FGFR1", position=488, wt="C", mut="S"):
    return {"gene": gene, "position": position, "wt": wt, "mut": mut}


def test_actual_templates_capture_the_covalent_drug_contact(built_store):
    _, manifest, store = built_store
    pem = store.features(candidate(), "pem", "structure")
    fut = store.features(candidate(), "fut", "structure")
    assert 9 < pem["structure_ligand_min_heavy_distance"] < 10
    assert 1.6 < fut["structure_ligand_min_heavy_distance"] < 2.0
    assert fut["structure_motif_covalent_cys_ca_distance"] == 0
    assert pem["structure_motif_covalent_cys_ca_distance"] == 0
    assert manifest["templates"]["pem"]["pdb_id"] == "7WCL"
    assert manifest["templates"]["fut"]["pdb_id"] == "6MZW"
    assert manifest["templates"]["fut"]["author_chain"] == "A"
    assert manifest["coverage"]["FGFR1"]["fut"]["geometry_available"] == 282


def test_alignment_uses_assay_numbering_for_every_paralog(built_store):
    _, _, store = built_store
    for gene, position in [("FGFR1", 488), ("FGFR2", 491), ("FGFR3", 482), ("FGFR4", 477)]:
        features = store.features(candidate(gene, position), "fut", "structure")
        assert features["structure_motif_covalent_cys_ca_distance"] == 0
        assert features["structure_ligand_min_heavy_distance"] < 2
        assert features["structure_homolog_template"] == (gene != "FGFR1")
    for gene, position in [("FGFR1", 641), ("FGFR2", 644), ("FGFR3", 635), ("FGFR4", 630)]:
        features = store.features(candidate(gene, position, "D", "A"), "pem", "structure")
        assert features["structure_motif_dfg_ca_distance"] == 0


def test_engineered_and_native_template_mismatches_mask_sidechain_geometry(built_store):
    _, _, store = built_store
    changed = store.features(candidate(position=584, wt="C"), "pem", "structure")
    assert changed["structure_geometry_available"] == 1
    assert changed["structure_engineered_site"] == 1
    assert changed["structure_sidechain_geometry_available"] == 0
    assert changed["structure_ligand_min_heavy_distance"] == MISSING_DISTANCE
    assert changed["structure_sasa"] == -1
    assert changed["structure_ligand_ca_distance"] != MISSING_DISTANCE
    # FGFR4's native gatekeeper-hinge region includes Cys where FGFR1 has Tyr.
    mismatch = store.features(candidate("FGFR4", 552, "C", "S"), "pem", "structure")
    assert mismatch["structure_template_residue_matches"] == 0
    assert mismatch["structure_sidechain_geometry_available"] == 0
    assert mismatch["structure_ligand_min_heavy_distance"] == MISSING_DISTANCE
    assert mismatch["structure_ligand_ca_distance"] < MISSING_DISTANCE


def test_outside_structure_remains_missing_and_interactions_do_not_treat_it_as_contact(built_store):
    output, _, store = built_store
    sequence = json.loads((output / "references.json").read_text())["sequences"]["FGFR1"]
    features = store.features(candidate(position=1, wt=sequence[0]), "pem", "interactions")
    assert features["structure_geometry_available"] == 0
    assert features["structure_alignment_available"] == 0
    assert features["structure_ligand_ca_distance"] == MISSING_DISTANCE
    assert all(value == 0 for name, value in features.items() if name.startswith("interaction_"))
    assert all(type(value) in (int, float) and math.isfinite(value) for value in features.values())


def test_features_are_cumulative_numeric_and_ignore_outcome_poison(built_store):
    _, _, store = built_store
    clean = candidate()
    poisoned = {**clean, "label": True, "score": 9000, "split": "final", "candidate_id": "positive",
                "aligned_position": -999, "context": "leaked", "drug": "fut"}
    chemistry = store.features(clean, "pem", "chemistry")
    structure = store.features(clean, "pem", "structure")
    interactions = store.features(clean, "pem", "interactions")
    assert chemistry.items() <= structure.items() <= interactions.items()
    assert len([name for name in interactions if name.startswith("interaction_")]) <= 12
    assert store.features(poisoned, "pem", "interactions") == interactions
    assert chemistry == store.features(clean, "fut", "chemistry")
    assert all(type(value) in (int, float) and math.isfinite(value) for value in interactions.values())


@pytest.mark.parametrize("patch", [
    {"wt": "A"}, {"position": True}, {"position": 0}, {"position": 99999},
    {"gene": "FGFR5"}, {"mut": "X"}, {"wt": None},
])
def test_invalid_sequence_identity_is_rejected(built_store, patch):
    _, _, store = built_store
    with pytest.raises(ValueError):
        store.features({**candidate(), **patch}, "pem")


def test_unknown_drug_or_family_is_rejected(built_store):
    _, _, store = built_store
    with pytest.raises(ValueError):
        store.features(candidate(), "unknown")
    with pytest.raises(ValueError):
        store.features(candidate(), "pem", "hidden_labels")


def test_cache_hash_tampering_is_rejected_before_feature_use(built_store, tmp_path):
    output, _, _ = built_store
    for name in ["manifest.json", "references.json", "feature-cache.json"]:
        (tmp_path / name).write_bytes((output / name).read_bytes())
    cache = json.loads((tmp_path / "feature-cache.json").read_text())
    cache["FGFR1"]["488"]["pem"]["geometry_available"] = 999
    (tmp_path / "feature-cache.json").write_text(json.dumps(cache))
    with pytest.raises(ValueError, match="hash|checksum"):
        BiologicalFeatureStore(tmp_path)


def test_build_reads_only_named_metadata_and_coordinates_never_measurements(built_store, tmp_path, monkeypatch):
    original_read = Path.read_text
    original_bytes = Path.read_bytes
    opened = []

    def checked_read(path, *args, **kwargs):
        opened.append(str(path))
        assert path.suffix != ".csv" and path.name not in {"outcomes.json", "candidates.json"}
        return original_read(path, *args, **kwargs)

    def checked_bytes(path, *args, **kwargs):
        opened.append(str(path))
        assert path.suffix != ".csv" and path.name not in {"outcomes.json", "candidates.json"}
        return original_bytes(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", checked_read)
    monkeypatch.setattr(Path, "read_bytes", checked_bytes)
    manifest = build_feature_store(ROOT / "data/raw", ROOT / "data/features-research", tmp_path)
    assert opened
    assert manifest["feature_cache_sha256"] == built_store[1]["feature_cache_sha256"]
    assert manifest["references_sha256"] == built_store[1]["references_sha256"]


def test_builder_refuses_different_assay_references_and_does_not_overwrite(built_store, tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    for gene in ["FGFR1", "FGFR2", "FGFR3", "FGFR4"]:
        for drug in ["pem", "fut"]:
            source = ROOT / "data/raw" / f"{gene}-{drug}.json"
            (raw / source.name).write_bytes(source.read_bytes())
    wrong = raw / "FGFR1-fut.json"
    data = copy.deepcopy(json.loads(wrong.read_text()))
    target = data["targetGenes"][0]["targetSequence"]
    target["sequence"] = "GCT" + target["sequence"][3:]
    wrong.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="reference|sequence"):
        build_feature_store(raw, ROOT / "data/features-research", tmp_path / "output")
    assert not (tmp_path / "output" / "manifest.json").exists()
    with pytest.raises(FileExistsError):
        build_feature_store(ROOT / "data/raw", ROOT / "data/features-research", built_store[0])
