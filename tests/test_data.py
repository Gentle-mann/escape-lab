import pytest

from escape_lab.data import (
    aligned_position_map,
    assign_split,
    collapse_labels,
    parse_measurement,
    parse_protein_variant,
    resistance_label,
    validate_source_identity,
)


def test_resistance_requires_valid_replicated_enrichment():
    assert resistance_label([1.4, 1.7, 1.8]) is True
    assert resistance_label([1.0, 1.5, 1.6]) is False
    assert resistance_label([1.0, float('nan'), 2.0]) is None
    assert resistance_label([1.0, 2.0]) is None
    assert resistance_label([1.0, -1.0, 2.0]) is None
    assert resistance_label([1.0, None, 2.0]) is None


def test_missing_measurements_are_preserved_not_imputed_as_negative():
    assert parse_measurement('NA') is None
    assert parse_measurement('nan') is None
    assert parse_measurement('1.25') == 1.25


def test_conflicting_protein_representations_remain_inconclusive():
    assert collapse_labels([True, True]) is True
    assert collapse_labels([False, False]) is False
    assert collapse_labels([True, False]) is None
    assert collapse_labels([True, None]) is None


def test_protein_identity_has_no_outcome_fields():
    assert parse_protein_variant('p.Asn549Lys') == ('N', 549, 'K')
    with pytest.raises(ValueError):
        parse_protein_variant('p.Asn549Ter')


def test_alignment_groups_homologous_positions_even_with_insertion():
    mapping = aligned_position_map('ACDEFGHIK', 'ACDWEFGHIK')
    assert mapping[5] == 4
    assert mapping[10] == 9
    assert assign_split(mapping[5]) == assign_split(4)


def test_split_is_deterministic_and_block_grouped():
    assert assign_split(101) == assign_split(101)
    assert assign_split(96) == assign_split(95)
    assert assign_split(101) in {'development', 'validation', 'final'}


def test_source_identity_rejects_swapped_drug_or_gene():
    metadata = {'urn': 'urn:mavedb:00001256-a-9', 'targetGenes': [{'name': 'FGFR2'}],
                'title': 'FGFR2 pemigatinib resistance'}
    assert validate_source_identity(metadata, 'FGFR2', 'pem') == metadata['urn']
    with pytest.raises(ValueError):
        validate_source_identity(metadata, 'FGFR2', 'fut')
    with pytest.raises(ValueError):
        validate_source_identity(metadata, 'FGFR1', 'pem')
