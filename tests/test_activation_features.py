import csv
import hashlib
import json
import math

import pytest

from escape_lab.activation_features import ActivationFeatureAugmenter


def candidate(position=500, wt="A", mut="G", gene="FGFR1", split="development"):
    names = {"A": "Ala", "G": "Gly", "P": "Pro", "S": "Ser", "D": "Asp", "Y": "Tyr"}
    protein = f"p.{names[wt]}{position}{names[mut]}"
    return {
        "candidate_id": f"{gene}:{protein}", "gene": gene, "protein_variant": protein,
        "position": position, "aligned_position": position, "wt": wt, "mut": mut,
        "split": split,
    }


@pytest.fixture
def references(tmp_path):
    sequence = list("A" * 830)
    for position, residue in {520: "S", 530: "Y", 546: "N", 562: "E", 638: "K",
                              653: "Y", 654: "Y"}.items():
        sequence[position - 1] = residue
    sequence = "".join(sequence)
    content = {
        "schema_version": "biological-features-v1",
        "sequences": {gene: sequence for gene in ("FGFR1", "FGFR2", "FGFR3", "FGFR4")},
        "sequence_sha256": {
            gene: hashlib.sha256(sequence.encode()).hexdigest()
            for gene in ("FGFR1", "FGFR2", "FGFR3", "FGFR4")
        },
    }
    path = tmp_path / "references.json"
    path.write_text(json.dumps(content))
    return path


def scores(tmp_path, rows):
    fields = ["gene", "protein_variant", "candidate_id_name_match", "alphamissense",
              "alphamissense_status", "eve", "eve_status", "prediction_target", "revel"]
    path = tmp_path / "scores.csv"
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def score_row(value=None, **changes):
    item = candidate() if value is None else value
    row = {
        "gene": item["gene"], "protein_variant": f"{item['wt']}{item['position']}{item['mut']}",
        "candidate_id_name_match": item["candidate_id"], "alphamissense": "0.7",
        "alphamissense_status": "available", "eve": "0.4", "eve_status": "available",
        "prediction_target": "untrusted annotation", "revel": "DO NOT PARSE",
    }
    row.update(changes)
    return row


def make_store(tmp_path, references, items=None, rows=None):
    items = [candidate()] if items is None else items
    rows = [score_row()] if rows is None else rows
    return ActivationFeatureAugmenter(None, items, references, scores(tmp_path, rows))


def test_features_include_context_scores_and_sequence_proximity_only(tmp_path, references):
    store = make_store(tmp_path, references)
    values = store.features(candidate())
    assert values["activation_alphamissense"] == 0.7
    assert values["activation_eve"] == 0.4
    assert values["activation_context_minus3"] == "A"
    assert values["activation_aligned_position"] == 0.5
    assert values["activation_brake_n546_sequence_proximity"] == pytest.approx(math.exp(-4.6))
    assert not any("revel" in name or "target" in name or "label" in name for name in values)
    assert not any("3d" in name or "distance" in name for name in values)


def test_positions_change_features_despite_identical_substitution(tmp_path, references):
    items = [candidate(), candidate(510)]
    store = make_store(tmp_path, references, items, [score_row(item) for item in items])
    first, second = [store.features(item) for item in items]
    assert first["activation_position_rbf_510"] != second["activation_position_rbf_510"]
    assert first["activation_brake_n546_sequence_proximity"] != second[
        "activation_brake_n546_sequence_proximity"]


@pytest.mark.parametrize("status", ["missing", "discordant_not_collapsed"])
def test_missing_or_conflicted_scores_are_not_parsed(tmp_path, references, status):
    store = make_store(tmp_path, references, rows=[score_row(eve="POISON", eve_status=status)])
    values = store.features(candidate())
    assert math.isnan(values["activation_eve"])
    assert values["activation_eve_missing"] == 1
    assert values["activation_eve_conflict"] == int(status == "discordant_not_collapsed")


def test_absent_score_rows_are_missing_not_zero(tmp_path, references):
    values = make_store(tmp_path, references, rows=[]).features(candidate())
    assert math.isnan(values["activation_alphamissense"])
    assert values["activation_alphamissense_missing"] == 1


def test_excluded_rows_are_skipped_before_status_and_float_parsing(tmp_path, references):
    excluded = candidate(510, split="final")
    rows = [score_row(), score_row(excluded, alphamissense="POISON", eve="POISON",
                                   alphamissense_status="POISON", eve_status="POISON")]
    assert make_store(tmp_path, references, rows=rows).features(candidate())["activation_eve"] == 0.4


@pytest.mark.parametrize("changes", [
    {"wt": "S"}, {"position": 501}, {"aligned_position": 501},
    {"protein_variant": "p.Ala500Pro"}, {"candidate_id": "FGFR1:p.Ala500Pro"},
    {"split": "final"}, {"split": "external"},
])
def test_candidate_identity_is_validated(tmp_path, references, changes):
    bad = candidate() | changes
    with pytest.raises(ValueError):
        make_store(tmp_path, references, items=[bad])


@pytest.mark.parametrize("changes", [
    {"gene": "FGFR2"}, {"protein_variant": "S500G"}, {"protein_variant": "A501G"},
    {"alphamissense_status": "unknown"}, {"eve": "nan"}, {"eve": "1.1"},
])
def test_allowed_score_rows_require_exact_identity_and_valid_scores(tmp_path, references, changes):
    with pytest.raises(ValueError):
        make_store(tmp_path, references, rows=[score_row(**changes)])


def test_anchor_sequence_is_checked_even_if_hash_was_updated(tmp_path, references):
    content = json.loads(references.read_text())
    sequence = list(content["sequences"]["FGFR1"])
    sequence[545] = "A"
    content["sequences"]["FGFR1"] = "".join(sequence)
    content["sequence_sha256"]["FGFR1"] = hashlib.sha256("".join(sequence).encode()).hexdigest()
    references.write_text(json.dumps(content))
    with pytest.raises(ValueError, match="anchor"):
        make_store(tmp_path, references)


def test_annotations_and_outcomes_do_not_enter_features(tmp_path, references):
    store = make_store(tmp_path, references)
    clean = store.features(candidate())
    annotated = candidate() | {"activation": 999, "label": True, "prediction_target": "GOF",
                               "human_annotation": "resistant", "context": "POISON"}
    assert clean == store.features(annotated)


def test_query_cannot_change_registered_identity_or_add_other_candidates(tmp_path, references):
    store = make_store(tmp_path, references)
    for bad in [candidate(510), candidate() | {"aligned_position": 501},
                candidate() | {"split": "final"}]:
        with pytest.raises(ValueError):
            store.features(bad)


def test_st_to_de_substitution_flag_excludes_tyrosine(tmp_path, references):
    items = [candidate(520, "S", "D"), candidate(530, "Y", "D")]
    store = make_store(tmp_path, references, items, [score_row(item) for item in items])
    serine, tyrosine = [store.features(item) for item in items]
    assert serine["activation_st_to_de_substitution"] == 1
    assert tyrosine["activation_st_to_de_substitution"] == 0
    assert tyrosine["activation_tyrosine_loss"] == 1


def test_duplicate_allowed_score_rows_are_rejected(tmp_path, references):
    with pytest.raises(ValueError, match="Duplicate"):
        make_store(tmp_path, references, rows=[score_row(), score_row()])
