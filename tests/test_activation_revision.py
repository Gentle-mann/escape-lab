import numpy as np
import pytest

from escape_lab.activation_revision import (
    ActivationRevisionModel,
    DualActivationStore,
    crossfit_outputs,
    evaluate_revision,
    resistance_gate,
    summarize_activation,
)
from escape_lab.activation_revision_run import paired_uncertainty


def candidates(n=36):
    return [{"candidate_id": f"FGFR1:p.Ala{i+1}Gly", "gene": "FGFR1", "wt": "A", "mut": "G",
             "position": i+1, "group_id": f"block-{i//4}", "split": "development"}
            for i in range(n)]


class Store:
    def features(self, c, drug, family):
        return {"chemistry_x": float(c["position"] % 7), "structure_sasa": c["position"],
                "structure_ca_available": 1, "structure_ligand_ca_distance": 4.}


class Augmenter:
    def features(self, c):
        return {"context_test": str(c["position"] % 5), "novel_missing": np.nan}


def targets(rows):
    return {c["candidate_id"]: {"value": float(i % 5), "active": i % 3 == 0,
                               "replicates": [float(i % 5)] * 4}
            for i, c in enumerate(rows)}


@pytest.mark.parametrize("arm", ["original", "capacity", "biology", "biology_tail"])
def test_model_rejects_activation_labels_outside_training(arm):
    rows = candidates()
    model = ActivationRevisionModel(Store(), Augmenter(), arm)
    with pytest.raises(ValueError, match="exactly"):
        model.fit(rows[:-1], targets(rows))


def test_two_head_continuous_output_matches_enriched_regression():
    rows = candidates(80)
    a = ActivationRevisionModel(Store(), Augmenter(), "biology").fit(rows, targets(rows))
    b = ActivationRevisionModel(Store(), Augmenter(), "biology_tail").fit(rows, targets(rows))
    np.testing.assert_array_equal(a.predict(rows)["value"], b.predict(rows)["value"])
    assert np.all((b.predict(rows)["active_score"] >= 0) & (b.predict(rows)["active_score"] <= 1))


def test_tail_class_fallback_preserves_all_predictions():
    rows = candidates()
    truth = targets(rows)
    for row in truth.values():
        row["active"] = False
    model = ActivationRevisionModel(Store(), Augmenter(), "biology_tail").fit(rows, truth)
    assert model.tail_fallback == "single_class_training"
    assert len(model.predict(rows)["active_score"]) == len(rows)
    assert np.all(model.predict(rows)["active_score"] == 0)


def test_crossfit_two_heads_cannot_use_the_prediction_groups_targets():
    rows = candidates()
    truth = targets(rows)

    class MeanModel:
        tail_fallback = None

        def fit(self, training, outcomes):
            assert set(outcomes) == {c["candidate_id"] for c in training}
            self.value = np.mean([r["value"] for r in outcomes.values()])
            self.active = np.mean([r["active"] for r in outcomes.values()])
            return self

        def predict(self, pool):
            return {"value": np.full(len(pool), self.value),
                    "active_score": np.full(len(pool), self.active)}

    first, audit = crossfit_outputs(rows, truth, MeanModel, "test")
    for c in rows:
        if c["group_id"] == rows[0]["group_id"]:
            truth[c["candidate_id"]] = {"value": 1000., "active": True}
    second, _ = crossfit_outputs(rows, truth, MeanModel, "test")
    assert first[rows[0]["candidate_id"]] == second[rows[0]["candidate_id"]]
    assert all(not set(a["training_groups"]) & set(a["prediction_groups"]) for a in audit)


def test_extra_head_is_only_added_to_two_head_resistance_model():
    row = candidates(1)[0]
    values = {row["candidate_id"]: {"value": 2., "active_score": .8}}
    one = DualActivationStore(Store(), values, False).features(row, "pem", "structure")
    two = DualActivationStore(Store(), values, True).features(row, "pem", "structure")
    assert "activation_active_score" not in one
    assert two["activation_active_score"] == .8


def test_undefined_gene_class_metric_is_not_silently_removed_from_macro():
    rows = [{"candidate_id": str(i), "gene": gene, "group_id": str(i), "fold": "f",
             "truth": float(i), "prediction": float(i), "active": False,
             "active_score": float(i)} for gene in ("FGFR1", "FGFR2", "FGFR3", "FGFR4")
            for i in range(4)]
    metrics = summarize_activation(rows)
    assert metrics["macro_gene_spearman"] == 1.
    assert metrics["macro_gene_active_ap"] is None


@pytest.mark.parametrize("value", [None, float("nan")])
def test_resistance_gate_undefined_values_cannot_pass(value):
    metric = {"macro_ap": value, "per_drug": {d: {"average_precision": value, "brier": .1}
                                              for d in ("pem", "fut")}, "by_fold": {}}
    result = {"resistance_metrics": {"structure": metric, "biology_tail": metric}}
    assert resistance_gate(result)["status"] == "inconclusive"


def test_outer_labels_cannot_change_either_activation_head_or_resistance_predictions():
    rows = candidates(80)
    for i, c in enumerate(rows):
        c["gene"] = f"FGFR{i % 4 + 1}"
    truth = targets(rows)
    resistance = {c["candidate_id"]: {"pem": i % 3 == 0, "fut": i % 5 == 0}
                  for i, c in enumerate(rows)}
    first = evaluate_revision(rows, truth, resistance, Store(), Augmenter())
    heldout = {r["candidate_id"] for r in first["activation_records"]["original"]
               if r["fold"] == "pilot-0"}
    for cid in heldout:
        truth[cid]["value"] += 1000
        truth[cid]["active"] = not truth[cid]["active"]
    second = evaluate_revision(rows, truth, resistance, Store(), Augmenter())
    for arm in ("original", "capacity", "biology", "biology_tail"):
        for first_row, second_row in zip(first["activation_records"][arm],
                                         second["activation_records"][arm], strict=True):
            if first_row["fold"] == "pilot-0":
                assert first_row["prediction"] == second_row["prediction"]
                assert first_row["active_score"] == second_row["active_score"]
        a = [r for r in first["resistance_records"][arm] if r["fold"] == "pilot-0"]
        b = [r for r in second["resistance_records"][arm] if r["fold"] == "pilot-0"]
        assert a == b


def test_bootstrap_identical_predictors_produce_zero_deltas():
    activation = []
    resistance = []
    for gene in ("FGFR1", "FGFR2", "FGFR3", "FGFR4"):
        for i in range(8):
            cid = f"{gene}-{i}"
            activation.append({"candidate_id": cid, "gene": gene, "group_id": str(i // 2),
                               "truth": float(i), "prediction": float(i), "active": i % 2 == 0,
                               "active_score": float(i)})
            for drug in ("pem", "fut"):
                resistance.append({"candidate_id": cid, "group_id": str(i // 2), "drug": drug,
                                   "label": i % 2 == 0, "probability": (i + 1) / 10})
    raw = {"activation_records": {"original": activation, "biology_tail": activation},
           "resistance_records": {"structure": resistance, "biology_tail": resistance}}
    report = paired_uncertainty(raw, draws=10)
    assert all(m["percentile_2_5"] == 0 and m["percentile_97_5"] == 0
               for m in report["contrasts"].values())
