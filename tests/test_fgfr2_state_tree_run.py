import copy

import pytest

from escape_lab.fgfr2_state_tree_run import fixed_partitions, paired_uncertainty, screening_gate


def candidates():
    return [{"candidate_id": f"{gene}-{i}", "gene": gene, "split": "development",
             "group_id": f"g{i}"} for i in range(3) for gene in ("FGFR1", "FGFR2")]


def test_saved_partitions_keep_pooled_training_and_only_score_fgfr2():
    data = candidates()
    tasks = fixed_partitions(data, {f"pilot-{i}": [f"g{i}"] for i in range(3)})
    for _, train, test in tasks:
        assert len(train) == 4
        assert len(test) == 1 and test[0]["gene"] == "FGFR2"
        assert not {c["group_id"] for c in train} & {c["group_id"] for c in test}
        assert train == [c for c in data if c["group_id"] != test[0]["group_id"]]


def test_saved_partition_rejects_final_or_missing_mapping():
    data = candidates()
    mapping = {f"pilot-{i}": [f"g{i}"] for i in range(3)}
    with pytest.raises(ValueError):
        fixed_partitions(data, {**mapping, "pilot-0": ["g1"]})
    data[0]["split"] = "final"
    with pytest.raises(ValueError):
        fixed_partitions(data, mapping)


def metric(rho, mae=1., hits=3):
    return {"mean_fold_spearman": rho, "mean_fold_mae": mae, "high_hits": hits,
            "by_fold": {str(i): {"spearman": rho} for i in range(3)}}


def test_fixed_gate_requires_ranking_error_and_tail_guards():
    baseline = metric(.38)
    assert screening_gate({"baseline": baseline, "contrast": metric(.41)})["status"] == "pass"
    for contrast in (metric(.40), metric(.42, 1.06), metric(.42, hits=2), metric(None)):
        assert screening_gate({"baseline": baseline, "contrast": contrast})["status"] != "pass"


def test_paired_block_intervals_identity_and_alignment():
    rows = [{"candidate_id": str(i), "group_id": "g", "fold": "pilot-0",
             "truth": float(i), "score": float(i), "high": i == 2, "strong": False}
            for i in range(3)]
    arms = {"baseline": rows, "contrast": copy.deepcopy(list(reversed(rows)))}
    result = paired_uncertainty(arms, draws=10)
    for m in result["contrasts"].values():
        assert m["percentile_2_5"] == 0
        assert m["percentile_97_5"] == 0
        assert m["valid_draws"] == 10
    arms["contrast"][0]["group_id"] = "bad"
    with pytest.raises(ValueError):
        paired_uncertainty(arms, draws=10)
