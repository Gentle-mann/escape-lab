import copy
import math

import numpy as np
import pytest

from escape_lab.mechanistic_pilot import (
    ActivationRegressor,
    ActivationStore,
    activation_from_rows,
    advancement_gate,
    crossfit_activation,
    eligible_candidates,
    evaluate_pilot,
)
from escape_lab.mechanistic_pilot_run import paired_intervals


def candidate(i, group=None):
    return {"candidate_id": f"FGFR1:p.Ala{i + 1}Gly", "gene": "FGFR1", "wt": "A",
            "mut": "G", "position": i + 1, "group_id": group or f"block-{i // 4}",
            "split": "development"}


class Store:
    def features(self, c, drug, family):
        return {"chemistry_x": float(c["position"] % 7),
                "structure_sasa": float(c["position"]),
                "structure_ca_available": 1,
                "structure_ligand_ca_distance": 4.0}


def test_exclusions_remove_entire_groups_and_final():
    rows = [candidate(0, "reserved"), candidate(1, "reserved"), candidate(2, "ok")]
    rows.append({**candidate(3, "ok"), "split": "final"})
    assert eligible_candidates(rows, {"reserved"}) == [rows[2]]


def test_import_aggregates_codons_within_replicates_and_excludes_before_values():
    class Forbidden:
        def __getitem__(self, index):
            if index == 0:
                return "FGFR1"
            if index == 4:
                return "A99G"
            raise AssertionError("Excluded value accessed")

    rows = [Forbidden(), ["FGFR1", 1, "transcript", "G1A", "A1G", 1, 2, 3, 4],
            ["FGFR1", 2, "transcript", "G1C", "A1G", 3, 4, 5, 6]]
    values, audit = activation_from_rows(rows, [candidate(0)])
    assert values[candidate(0)["candidate_id"]] == pytest.approx(math.log2(4.5))
    assert audit["eligible_nucleotide_rows"] == 2
    assert audit["complete_protein_variants"] == 1


def test_import_missing_is_not_zero_and_negative_is_rejected():
    row = ["FGFR1", 1, "tx", "G1A", "A1G", None, 1, 1, 1]
    values, audit = activation_from_rows([row], [candidate(0)])
    assert values == {}
    assert audit["incomplete_protein_variants"] == 1
    row[5] = -1
    with pytest.raises(ValueError, match="nonnegative"):
        activation_from_rows([row], [candidate(0)])


def test_duplicate_nucleotide_identity_rejected():
    row = ["FGFR1", 1, "tx", "G1A", "A1G", 1, 1, 1, 1]
    with pytest.raises(ValueError, match="Duplicate nucleotide"):
        activation_from_rows([row, row], [candidate(0)])


def test_activation_feature_contract_excludes_ligand_information():
    model = ActivationRegressor(Store())
    features = model.features([candidate(0)])[0]
    assert "structure_ligand_ca_distance" not in features
    assert features["structure_sasa"] == 1


@pytest.mark.parametrize("shuffled", [False, True])
def test_nested_crossfit_never_uses_target_group_labels(shuffled):
    rows = [candidate(i) for i in range(24)]
    truth = {c["candidate_id"]: float(c["position"] % 5) for c in rows}

    class MeanRegressor:
        def fit(self, training, values):
            self.mean = np.mean([values[c["candidate_id"]] for c in training])
            return self

        def predict(self, pool):
            return np.full(len(pool), self.mean)

    first, audit = crossfit_activation(rows, truth, MeanRegressor, "test", shuffled)
    altered = dict(truth)
    for c in rows:
        if c["group_id"] == rows[0]["group_id"]:
            altered[c["candidate_id"]] += 1000
    second, _ = crossfit_activation(rows, altered, MeanRegressor, "test", shuffled)
    assert first[rows[0]["candidate_id"]] == second[rows[0]["candidate_id"]]
    assert all(not set(a["training_groups"]) & set(a["prediction_groups"]) for a in audit)


def test_activation_store_requires_exact_available_identity():
    c = candidate(0)
    store = ActivationStore(Store(), {c["candidate_id"]: 2.0})
    feature = store.features(c, "pem", "structure")
    assert feature["activation_x_ligand_proximity"] == pytest.approx(2 * math.exp(-4 / 5))
    with pytest.raises(KeyError):
        store.features(candidate(1), "pem", "structure")


def test_outer_test_activation_cannot_change_sequence_only_predictions():
    rows = [candidate(i) for i in range(36)]
    activation = {c["candidate_id"]: float(c["position"] % 5) for c in rows}
    labels = {c["candidate_id"]: {"pem": i % 3 == 0, "fut": i % 4 == 0}
              for i, c in enumerate(rows)}
    first = evaluate_pilot(rows, labels, activation, Store())
    heldout = {r["candidate_id"] for r in first["records"]["predicted_activation"]
               if r["fold"] == "pilot-0"}
    altered = {cid: value + (1000 if cid in heldout else 0)
               for cid, value in activation.items()}
    second = evaluate_pilot(rows, labels, altered, Store())
    for arm in ("structure", "predicted_activation", "shuffled_activation"):
        a = [r for r in first["records"][arm] if r["fold"] == "pilot-0"]
        b = [r for r in second["records"][arm] if r["fold"] == "pilot-0"]
        assert a == b
    populations = [{(r["candidate_id"], r["drug"]) for r in records}
                   for records in first["records"].values()]
    assert all(p == populations[0] for p in populations)


def test_gate_handles_undefined_metrics_without_crashing():
    metric = {"macro_ap": None, "per_drug": {d: {"average_precision": None, "brier": None}
                                             for d in ("pem", "fut")}, "by_fold": {}}
    result = {"metrics": {arm: copy.deepcopy(metric) for arm in
                          ("structure", "predicted_activation", "shuffled_activation")}}
    assert advancement_gate(result, {})["status"] == "inconclusive"


def test_paired_bootstrap_identical_scores_have_zero_difference():
    rows = [{"candidate_id": str(i), "group_id": str(i // 2), "drug": drug,
             "label": i % 2 == 0, "probability": (i + 1) / 10}
            for i in range(6) for drug in ("pem", "fut")]
    records = {arm: rows for arm in ("structure", "predicted_activation", "measured_activation",
                                    "shuffled_activation")}
    intervals = paired_intervals(records, replicates=10)
    assert all(v == {"percentile_2_5": 0.0, "percentile_97_5": 0.0}
               for v in intervals["macro_ap_difference_vs_structure"].values())
