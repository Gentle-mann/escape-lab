"""Outcome-independent features from exact assay sequences and two FGFR1 structures.

The structures are homolog templates, not simulated mutant or patient structures.
Only named MaveDB metadata JSON files supply sequences; no scores are opened.
Geometry is computed once, then served by a hash-validated public feature store.
"""

import copy
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.SASA import ShrakeRupley
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParamData import kd

GENES = ("FGFR1", "FGFR2", "FGFR3", "FGFR4")
DRUGS = ("pem", "fut")
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
BLOSUM = substitution_matrices.load("BLOSUM62")
MISSING_DISTANCE = 999.0
SCHEMA_VERSION = "biological-features-v1"
SCORE_SETS = {
    "FGFR1": {"pem": "10", "fut": "4"},
    "FGFR2": {"pem": "9", "fut": "7"},
    "FGFR3": {"pem": "5", "fut": "2"},
    "FGFR4": {"pem": "8", "fut": "6"},
}
TEMPLATES = {
    "pem": {"pdb_id": "7WCL", "ligand": "8ZF", "author_chain": "A",
            "first_protein_label_seq_id": 3, "last_protein_label_seq_id": 310},
    "fut": {"pdb_id": "6MZW", "ligand": "TZ0", "author_chain": "A",
            "first_protein_label_seq_id": 5, "last_protein_label_seq_id": 311},
}
# Author FGFR1 positions, validated against exact sequence and deposited atoms.
# These describe general kinase machinery, not measured resistance hotspots.
ANCHORS = {
    "covalent_cys": (488,), "catalytic_lys": (514,), "alpha_c_glu": (531,),
    "gatekeeper": (561,), "p_loop": tuple(range(484, 493)),
    "hrd": (621, 622, 623), "dfg": (641, 642, 643), "ape": (668, 669, 670),
}
ANCHOR_SEQUENCE = {488: "C", 514: "K", 531: "E", 561: "V",
                   621: "H", 622: "R", 623: "D", 641: "D", 642: "F", 643: "G",
                   668: "A", 669: "P", 670: "E"}
PUBLIC_FIELDS = frozenset({"gene", "position", "wt", "mut"})


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json_write(path, value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
    with Path(path).open("x") as handle:
        handle.write(encoded)


def _sequence_mapping(reference, template):
    """Map 1-based target positions to template positions, omitting ambiguous sites.

    Local alignment prevents terminal tags/tails being forced into the kinase.
    Consensus across all optimal alignments (bounded to 64) avoids an arbitrary
    first alignment. Insertions have no mapping, rather than inheriting a neighbor.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = BLOSUM
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(reference, template)
    try:
        count = len(alignments)
    except OverflowError:
        count = 65
    if not count:
        return {}, set()
    if count > 64:
        return {}, set(range(1, len(reference) + 1))
    maps = []
    for alignment in alignments:
        mapping = {}
        for ref_block, template_block in zip(*alignment.aligned, strict=True):
            for ref_index, template_index in zip(range(*ref_block), range(*template_block), strict=True):
                mapping[int(ref_index) + 1] = int(template_index) + 1
        maps.append(mapping)
    all_positions = set().union(*(mapping.keys() for mapping in maps))
    ambiguous = {position for position in all_positions
                 if len({mapping.get(position) for mapping in maps}) != 1}
    return {position: maps[0][position] for position in all_positions - ambiguous}, ambiguous


def _heavy_atoms(residue):
    # Bio.PDB selects the highest-occupancy alternative; ties retain file order.
    return [atom for atom in residue.get_atoms() if atom.element not in {"H", "D"}]


def _chain_geometry(chain, ligand, anchors):
    """Raw geometry for one protein chain and its ligand, keyed by author position.

    Distances are Angstroms, SASA is square Angstroms. Contacts count atom pairs
    with strict distance <4/<5 Angstroms; density counts other CAs within 8A.
    SASA excludes the ligand, waters, other chains and hydrogen atoms.
    """
    protein = copy.deepcopy(chain)
    for residue in list(protein):
        if residue.id[0] != " ":
            protein.detach_child(residue.id)
            continue
        for atom in list(residue):
            if atom.element in {"H", "D"}:
                residue.detach_child(atom.id)
    if not len(protein):
        raise ValueError("Template chain has no protein residues")
    ligand_atoms = _heavy_atoms(ligand)
    if not ligand_atoms:
        raise ValueError("Ligand has no heavy atoms")
    ligand_xyz = np.asarray([atom.coord for atom in ligand_atoms], dtype=float)
    ShrakeRupley(n_points=240, probe_radius=1.4).compute(protein, level="R")
    ca = {residue.id[1]: np.asarray(residue["CA"].coord, dtype=float)
          for residue in protein if "CA" in residue}
    output = {}
    for residue in protein:
        position = residue.id[1]
        if residue.id[2] != " ":
            raise ValueError("Template insertion codes require explicit mapping support")
        atoms = _heavy_atoms(residue)
        xyz = np.asarray([atom.coord for atom in atoms], dtype=float)
        distances = np.linalg.norm(xyz[:, None, :] - ligand_xyz[None, :, :], axis=2)
        backbone = [atom.coord for atom in atoms if atom.name in {"N", "CA", "C", "O", "OXT"}]
        result = {
            "ligand_min_heavy_distance": float(distances.min()),
            "ligand_atom_pairs_4a": int((distances < 4).sum()),
            "ligand_atom_pairs_5a": int((distances < 5).sum()),
            "ligand_backbone_distance": MISSING_DISTANCE,
            "ligand_ca_distance": MISSING_DISTANCE,
            "sasa": float(residue.sasa), "ca_available": int(position in ca),
            "ca_neighbors_8a": -1,
        }
        if backbone:
            result["ligand_backbone_distance"] = float(np.linalg.norm(
                np.asarray(backbone)[:, None, :] - ligand_xyz[None, :, :], axis=2).min())
        if position in ca:
            result["ligand_ca_distance"] = float(np.linalg.norm(ligand_xyz - ca[position], axis=1).min())
            result["ca_neighbors_8a"] = sum(
                1 for other, coord in ca.items()
                if other != position and np.linalg.norm(coord - ca[position]) < 8)
        for name, positions in anchors.items():
            points = [ca[value] for value in positions if value in ca]
            available = position in ca and bool(points)
            result[f"motif_{name}_available"] = int(available)
            result[f"motif_{name}_ca_distance"] = (
                float(np.linalg.norm(np.asarray(points) - ca[position], axis=1).min())
                if available else MISSING_DISTANCE)
        output[position] = result
    return output


def _empty_geometry():
    result = {
        "alignment_available": 0, "alignment_ambiguous": 0, "geometry_available": 0,
        "ca_available": 0, "sidechain_geometry_available": 0,
        "template_residue_matches": 0, "homolog_template": 0, "engineered_site": 0,
        "ligand_min_heavy_distance": MISSING_DISTANCE,
        "ligand_backbone_distance": MISSING_DISTANCE, "ligand_ca_distance": MISSING_DISTANCE,
        "ligand_atom_pairs_4a": -1, "ligand_atom_pairs_5a": -1,
        "sasa": -1, "ca_neighbors_8a": -1,
    }
    for name in ANCHORS:
        result[f"motif_{name}_available"] = 0
        result[f"motif_{name}_ca_distance"] = MISSING_DISTANCE
    return result


def _read_assay_sequences(raw_dir):
    sequences, hashes = {}, {}
    for gene in GENES:
        for drug in DRUGS:
            path = raw_dir / f"{gene}-{drug}.json"
            metadata = json.loads(path.read_text())
            expected = "urn:mavedb:00001256-a-" + SCORE_SETS[gene][drug]
            drug_name = {"pem": "pemigatinib", "fut": "futibatinib"}[drug]
            if (metadata.get("urn") != expected or len(metadata.get("targetGenes", [])) != 1
                    or metadata["targetGenes"][0].get("name") != gene
                    or metadata.get("title", "").lower() != f"{gene} {drug_name} resistance".lower()):
                raise ValueError("Assay reference source identity mismatch")
            target = metadata["targetGenes"][0]["targetSequence"]
            dna = target["sequence"].upper()
            if len(dna) % 3 or any(base not in "ACGT" for base in dna):
                raise ValueError("Assay reference must be unambiguous complete codons")
            sequence = str(Seq(dna).translate()).removesuffix("*")
            if not sequence or any(aa not in AMINO_ACIDS for aa in sequence):
                raise ValueError("Assay reference contains unsupported protein residues")
            if gene in sequences and sequences[gene] != sequence:
                raise ValueError("Drug assays have different reference sequences")
            sequences[gene] = sequence
            hashes[path.name] = _hash(path)
    return sequences, hashes


def _load_template(structure_dir, drug, fgfr1_reference):
    specification = TEMPLATES[drug]
    code, author_chain = specification["pdb_id"], specification["author_chain"]
    path = structure_dir / f"{code}.cif"
    mmcif = MMCIF2Dict(str(path))
    if mmcif.get("_entry.id") != [code]:
        raise ValueError("Structure identity mismatch")
    entity_ids = mmcif["_entity_poly.entity_id"]
    if "1" not in entity_ids:
        raise ValueError("Expected protein entity 1")
    entity_index = entity_ids.index("1")
    full_sequence = mmcif["_entity_poly.pdbx_seq_one_letter_code_can"][entity_index].replace("\n", "")
    first = specification["first_protein_label_seq_id"]
    last = specification["last_protein_label_seq_id"]
    sequence = full_sequence[first - 1:last]
    if len(full_sequence) != last or any(aa not in AMINO_ACIDS for aa in sequence):
        raise ValueError("Unexpected template construct or expression tag boundary")
    # Validate the previously verified SIFTS bounds independently against the
    # actual assay sequence. Do not use the stale 6MZW _struct_ref_seq offsets.
    mapping, ambiguous = _sequence_mapping(fgfr1_reference, sequence)
    if ambiguous or len(mapping) != len(sequence):
        raise ValueError("FGFR1 template cannot be aligned completely and unambiguously")
    differences = {position: (fgfr1_reference[position - 1], sequence[index - 1])
                   for position, index in mapping.items()
                   if fgfr1_reference[position - 1] != sequence[index - 1]}
    if differences != {584: ("C", "S")}:
        raise ValueError("Template substitutions differ from verified C584S construct")
    label_to_author = {}
    label_to_residue = {}
    for label, author, chain, entity, component, insertion in zip(
            mmcif["_atom_site.label_seq_id"], mmcif["_atom_site.auth_seq_id"],
            mmcif["_atom_site.auth_asym_id"], mmcif["_atom_site.label_entity_id"],
            mmcif["_atom_site.label_comp_id"], mmcif["_atom_site.pdbx_PDB_ins_code"], strict=True):
        if chain != author_chain or entity != "1" or label in {".", "?"}:
            continue
        if insertion not in {".", "?"}:
            raise ValueError("Unexpected insertion code in template")
        label, author = int(label), int(author)
        if label in label_to_author and label_to_author[label] != author:
            raise ValueError("Ambiguous template residue identity")
        residue = seq1(component)
        if not 1 <= label <= len(full_sequence) or residue != full_sequence[label - 1]:
            raise ValueError("Template atom identity disagrees with entity sequence")
        label_to_author[label] = author
        label_to_residue[label] = residue
    for position, residue in ANCHOR_SEQUENCE.items():
        if fgfr1_reference[position - 1] != residue or position not in mapping:
            raise ValueError("Kinase anchor does not match assay reference")
        label = mapping[position] + first - 1
        if label_to_author.get(label) != position or label_to_residue.get(label) != residue:
            raise ValueError("Kinase anchor does not match observed structure")
    structure = MMCIFParser(QUIET=True).get_structure(code, str(path))
    chain = structure[0][author_chain]
    ligands = [residue for residue in chain if residue.resname == specification["ligand"]]
    if len(ligands) != 1:
        raise ValueError("Expected exactly one target ligand in template author chain")
    if drug == "fut":
        bonds = zip(mmcif.get("_struct_conn.conn_type_id", []),
                    mmcif.get("_struct_conn.ptnr1_auth_seq_id", []),
                    mmcif.get("_struct_conn.ptnr1_label_atom_id", []),
                    mmcif.get("_struct_conn.ptnr2_label_comp_id", []), strict=True)
        if not any(kind == "covale" and position == "488" and atom == "SG" and ligand == "TZ0"
                   for kind, position, atom, ligand in bonds):
            raise ValueError("Futibatinib template lacks verified covalent cysteine connection")
    geometry = _chain_geometry(chain, ligands[0], ANCHORS)
    return {"sequence": sequence, "offset": first - 1, "label_to_author": label_to_author,
            "label_to_residue": label_to_residue, "geometry": geometry,
            "engineered_labels": {mapping[584] + first - 1}, "sha256": _hash(path)}


def build_feature_store(raw_dir, structure_dir, output_dir):
    """Build a deterministic public cache without reading any assay outcomes.

    Existing artifacts are never overwritten. Full sequence metadata is trusted
    input, but only explicit sequence/identity fields enter the feature products.
    """
    raw_dir, structure_dir, output_dir = Path(raw_dir), Path(structure_dir), Path(output_dir)
    names = ("feature-cache.json", "references.json", "manifest.json")
    if any((output_dir / name).exists() for name in names):
        raise FileExistsError("Feature artifacts already exist; use a new output directory")
    sequences, source_hashes = _read_assay_sequences(raw_dir)
    templates = {drug: _load_template(structure_dir, drug, sequences["FGFR1"]) for drug in DRUGS}
    cache, mappings, coverage = {}, {}, {}
    for gene, reference in sequences.items():
        cache[gene] = {str(position): {} for position in range(1, len(reference) + 1)}
        mappings[gene], coverage[gene] = {}, {}
        for drug, template in templates.items():
            mapping, ambiguous = _sequence_mapping(reference, template["sequence"])
            if len(mapping) < 250:
                raise ValueError("Insufficient unambiguous kinase template alignment")
            mapping_records = {}
            for position in range(1, len(reference) + 1):
                features = _empty_geometry()
                features["homolog_template"] = int(gene != "FGFR1")
                features["alignment_ambiguous"] = int(position in ambiguous)
                if position in mapping:
                    index = mapping[position]
                    label = index + template["offset"]
                    author = template["label_to_author"].get(label)
                    template_residue = template["sequence"][index - 1]
                    engineered = label in template["engineered_labels"]
                    same = reference[position - 1] == template_residue
                    features.update(alignment_available=1, template_residue_matches=int(same),
                                    engineered_site=int(engineered))
                    mapping_records[str(position)] = {
                        "template_label_seq_id": label, "template_auth_seq_id": author,
                        "template_residue": template_residue, "target_residue": reference[position - 1],
                        "observed": author in template["geometry"], "engineered": engineered,
                    }
                    if author in template["geometry"]:
                        features.update(template["geometry"][author])
                        features["geometry_available"] = 1
                        features["sidechain_geometry_available"] = int(same and not engineered)
                        if not same or engineered:
                            features["ligand_min_heavy_distance"] = MISSING_DISTANCE
                            features["ligand_atom_pairs_4a"] = -1
                            features["ligand_atom_pairs_5a"] = -1
                            features["sasa"] = -1
                cache[gene][str(position)][drug] = features
            mappings[gene][drug] = mapping_records
            coverage[gene][drug] = {"reference_positions": len(reference), **{
                field: sum(cache[gene][str(position)][drug][field] for position in range(1, len(reference) + 1))
                for field in ("alignment_available", "alignment_ambiguous", "geometry_available",
                              "sidechain_geometry_available", "engineered_site")}}
    references = {"schema_version": SCHEMA_VERSION, "sequences": sequences,
                  "sequence_sha256": {gene: hashlib.sha256(sequence.encode()).hexdigest()
                                      for gene, sequence in sequences.items()},
                  "mappings": mappings}
    output_dir.mkdir(parents=True, exist_ok=True)
    _json_write(output_dir / "feature-cache.json", cache)
    _json_write(output_dir / "references.json", references)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "feature_cache_sha256": _hash(output_dir / "feature-cache.json"),
        "references_sha256": _hash(output_dir / "references.json"),
        "implementation_sha256": _hash(__file__),
        "source_metadata_sha256": source_hashes,
        "structure_sha256": {TEMPLATES[drug]["pdb_id"] + ".cif": template["sha256"]
                             for drug, template in templates.items()},
        "templates": TEMPLATES, "coverage": coverage,
        "structure_feature_names": sorted(_empty_geometry()),
        "methods": {"alignment": "local BLOSUM62, gap open -10, extend -0.5; optimal consensus",
                    "sasa": "ShrakeRupley, 240 points, 1.4A probe, protein chain only, no hydrogens",
                    "contacts": "heavy atom pairs <4A/<5A; C-alpha neighbors <8A, self excluded",
                    "sidechain_mask": "mask unequal template residue or engineered C584S site",
                    "missing_distance": MISSING_DISTANCE, "missing_area_or_count": -1,
                    "alternate_atoms": "Bio.PDB highest occupancy; first in file on ties",
                    "interpretation": "FGFR1 bound-structure proxy, not modeled mutant/paralog binding affinity"},
        "sources": {
            "pem_structure": "https://files.rcsb.org/download/7WCL.cif",
            "fut_structure": "https://files.rcsb.org/download/6MZW.cif",
            "pem_mapping": "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/7wcl",
            "fut_mapping": "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/6mzw",
            "pem_paper": "https://doi.org/10.1038/s42004-022-00718-z",
            "fut_paper": "https://doi.org/10.1002/cmdc.201800719",
            "hydropathy": "https://doi.org/10.1016/0022-2836(82)90515-0",
            "blosum": "https://doi.org/10.1073/pnas.89.22.10915",
        },
        "licenses": {"pdb_coordinates": "CC0 1.0", "assay_source_metadata": "CC BY-SA 4.0"},
        "excluded_inputs": ["assay scores", "resistance annotations", "split assignments",
                            "candidate identifiers", "learned variant effect scores"],
    }
    _json_write(output_dir / "manifest.json", manifest)
    return manifest


def _chemistry(wt, mut):
    charge_wt = int(wt in "KR") - int(wt in "DE")
    charge_mut = int(mut in "KR") - int(mut in "DE")
    result = {
        "blosum": float(BLOSUM[wt, mut]) / 10,
        "wt_hydropathy": kd[wt] / 10, "mut_hydropathy": kd[mut] / 10,
        "hydropathy_delta": (kd[mut] - kd[wt]) / 10,
        "wt_charge": charge_wt, "mut_charge": charge_mut, "charge_delta": charge_mut - charge_wt,
        "aromatic_delta": int(mut in "FWY") - int(wt in "FWY"),
        "proline_delta": int(mut == "P") - int(wt == "P"),
        "glycine_delta": int(mut == "G") - int(wt == "G"),
        "cysteine_loss": int(wt == "C" and mut != "C"), "same_residue": int(wt == mut),
    }
    for residue in AMINO_ACIDS:
        result[f"wt_{residue}"] = int(wt == residue)
        result[f"mut_{residue}"] = int(mut == residue)
    return result


def _interactions(chemistry, structure):
    # Interaction masks are evaluated before the sentinel can enter arithmetic.
    proximity = math.exp(-structure["ligand_ca_distance"] / 5) if structure["ca_available"] else 0
    burial = 1 / (1 + structure["sasa"] / 50) if structure["sidechain_geometry_available"] else 0
    density = (structure["ca_neighbors_8a"] / (structure["ca_neighbors_8a"] + 5)
               if structure["ca_available"] else 0)
    result = {f"{name}_x_ligand_proximity": chemistry[name] * proximity
              for name in ("hydropathy_delta", "charge_delta", "aromatic_delta", "proline_delta",
                           "glycine_delta", "blosum")}
    for name in ("hydropathy_delta", "charge_delta"):
        result[f"{name}_x_burial"] = chemistry[name] * burial
    result["proline_delta_x_density"] = chemistry["proline_delta"] * density
    for name, motif in (("cysteine_loss", "covalent_cys"), ("glycine_delta", "p_loop"),
                        ("aromatic_delta", "gatekeeper")):
        distance = structure[f"motif_{motif}_ca_distance"]
        value = math.exp(-distance / 5) if structure[f"motif_{motif}_available"] else 0
        result[f"{name}_x_{motif}_proximity"] = chemistry[name] * value
    return result


class BiologicalFeatureStore:
    """Immutable-by-interface, pickleable cache. No raw data is accessed at inference."""

    def __init__(self, output_dir):
        output_dir = Path(output_dir)
        manifest = json.loads((output_dir / "manifest.json").read_text())
        if manifest.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("Unsupported biological feature schema")
        for name, key in (("feature-cache.json", "feature_cache_sha256"),
                          ("references.json", "references_sha256")):
            if _hash(output_dir / name) != manifest.get(key):
                raise ValueError(f"Feature artifact checksum/hash mismatch: {name}")
        if _hash(__file__) != manifest.get("implementation_sha256"):
            raise ValueError("Feature implementation checksum/hash mismatch")
        cache = json.loads((output_dir / "feature-cache.json").read_text())
        references = json.loads((output_dir / "references.json").read_text())
        if references.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("Unsupported reference schema")
        sequences = references.get("sequences", {})
        if set(cache) != set(GENES) or set(sequences) != set(GENES):
            raise ValueError("Feature cache genes differ from exact reference genes")
        expected_keys = set(_empty_geometry())
        for gene, sequence in sequences.items():
            if (not isinstance(sequence, str) or not sequence or any(aa not in AMINO_ACIDS for aa in sequence)
                    or hashlib.sha256(sequence.encode()).hexdigest() != references["sequence_sha256"][gene]):
                raise ValueError("Invalid reference sequence/hash")
            if set(cache[gene]) != {str(position) for position in range(1, len(sequence) + 1)}:
                raise ValueError("Incomplete feature cache positions")
            for by_drug in cache[gene].values():
                if set(by_drug) != set(DRUGS):
                    raise ValueError("Incomplete drug feature cache")
                for features in by_drug.values():
                    if set(features) != expected_keys or any(
                            type(value) not in (int, float) or not math.isfinite(value)
                            for value in features.values()):
                        raise ValueError("Unexpected or nonfinite feature values")
        self.manifest = manifest
        self._sequences = sequences
        self._cache = cache

    def features(self, candidate, drug, family="structure"):
        if drug not in DRUGS or family not in {"chemistry", "structure", "interactions"}:
            raise ValueError("Unknown drug or biological feature family")
        public = {name: candidate[name] for name in PUBLIC_FIELDS if name in candidate}
        if set(public) != PUBLIC_FIELDS:
            raise ValueError("Candidate must provide gene, position, wt and mut")
        gene, position, wt, mut = public["gene"], public["position"], public["wt"], public["mut"]
        if (not isinstance(gene, str) or gene not in self._sequences or type(position) is not int
                or not 1 <= position <= len(self._sequences[gene])
                or not isinstance(wt, str) or len(wt) != 1 or wt not in AMINO_ACIDS
                or not isinstance(mut, str) or len(mut) != 1 or mut not in AMINO_ACIDS):
            raise ValueError("Invalid candidate reference identity or substitution")
        if self._sequences[gene][position - 1] != wt:
            raise ValueError("Candidate WT does not match exact assay reference sequence")
        chemistry = _chemistry(wt, mut)
        result = {f"chemistry_{name}": value for name, value in chemistry.items()}
        if family == "chemistry":
            return result
        structure = self._cache[gene][str(position)][drug]
        result.update({f"structure_{name}": value for name, value in structure.items()})
        if family == "interactions":
            result.update({f"interaction_{name}": value
                           for name, value in _interactions(chemistry, structure).items()})
        return result
