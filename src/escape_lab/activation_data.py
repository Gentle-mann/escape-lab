"""Development-only activation ingestion, preserving replicate summaries.

The binary target is a protein-level adaptation of the study's nucleotide-level
definition: at least three of four per-protein replicate summaries are >= 1.5.
It must not be described as an exact reconstruction of published variant labels.
"""

import math
import re
import statistics

import numpy as np
import openpyxl

ACTIVATION_HEADER = (
    "FGFR", "Index", "Transcript ID", "NT MUT", "AA MUT",
    "Activation BR1", "Activation BR2", "Activation BR3", "Activation BR4",
)


def _candidate_lookup(candidates):
    if any(c["split"] not in ("development", "validation") for c in candidates):
        raise ValueError("Activation ingestion only permits development candidates")
    if len({c["candidate_id"] for c in candidates}) != len(candidates):
        raise ValueError("Duplicate candidate identity")
    lookup = {(c["gene"], c["wt"], c["position"], c["mut"]): c["candidate_id"]
              for c in candidates}
    if len(lookup) != len(candidates):
        raise ValueError("Duplicate candidate protein identity")
    return lookup


def activation_records_from_rows(rows, candidates):
    """Return complete protein-level activation records and aggregate coverage.

    Candidate eligibility must already exclude reserved groups. This importer
    rejects original final candidates and looks up gene/protein identity before
    accessing any nucleotide identity or activation value. It never reads source
    columns beyond I. Codons are collapsed separately within each replicate;
    missing measurements are not zero. The continuous target matches pilot v1.
    """
    candidates = list(candidates)
    lookup = _candidate_lookup(candidates)
    replicates, seen = {}, set()
    count = 0
    for row in rows:
        gene, aa = row[0], row[4]
        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", str(aa))
        if match is None:
            continue
        wt, position, mut = match.groups()
        cid = lookup.get((gene, wt, int(position), mut))
        if cid is None:
            continue
        identity = (gene, row[3])
        if identity in seen:
            raise ValueError("Duplicate nucleotide identity in activation source")
        seen.add(identity)
        count += 1
        bucket = replicates.setdefault(cid, [[], [], [], []])
        for i in range(4):
            value = row[5 + i]
            if value is None or str(value).strip().lower() in ("", "na", "nan", "n/a"):
                continue
            value = float(value)
            if not math.isfinite(value) or value < 0:
                raise ValueError("Activation requires finite nonnegative enrichment")
            bucket[i].append(value)

    records, spreads = {}, []
    for c in candidates:
        cid = c["candidate_id"]
        reps = replicates.get(cid, [[], [], [], []])
        if not all(reps):
            continue
        medians = [statistics.median(rep) for rep in reps]
        records[cid] = {
            "value": math.log2(1 + statistics.median(medians)),
            "replicates": medians,
            "active": sum(value >= 1.5 for value in medians) >= 3,
        }
        spreads.append(float(np.std(np.log2(1 + np.asarray(medians)))))

    audit = {
        "eligible_nucleotide_rows": count,
        "eligible_protein_variants": len(candidates),
        "complete_protein_variants": len(records),
        "incomplete_protein_variants": len(candidates) - len(records),
        "median_within_variant_replicate_sd_log2_1p": statistics.median(spreads)
        if spreads else None,
        "aggregation": "median across codons within replicate, median of four replicates",
        "transform": "log2(1 + enrichment)",
        "active_definition": "at least three of four protein replicate medians >= 1.5",
        "active_definition_scope": "protein-level adaptation, not exact published nucleotide labels",
    }
    return records, audit


def load_development_activation(workbook_path, candidates):
    """Read only S10 identity/activation columns A:I with an exact header check."""
    candidates = list(candidates)
    _candidate_lookup(candidates)
    workbook = openpyxl.load_workbook(workbook_path, read_only=True, data_only=True)
    try:
        sheet = workbook["S10_Screen data table"]
        header = next(sheet.iter_rows(min_row=10, max_row=10, max_col=9, values_only=True))
        if tuple(header) != ACTIVATION_HEADER:
            raise ValueError("Unexpected activation source columns")
        return activation_records_from_rows(
            sheet.iter_rows(min_row=11, max_col=9, values_only=True), candidates)
    finally:
        workbook.close()
