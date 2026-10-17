import numpy as np

from escape_lab.baseline import NumericalSearch, candidate_features


def candidate(i, mutant='K'):
    return {'candidate_id': f'G:p.Ala{i}{mutant}', 'gene': 'FGFR1', 'wt': 'A',
            'mut': mutant, 'position': i, 'aligned_position': i, 'context': 'AAAGAAA',
            'group_id': f'block-{i // 12}'}


def test_features_cannot_learn_from_outcome_columns_or_ids():
    clean = candidate(500)
    poisoned = {**clean, 'dual_resistant': True, 'score': 999, 'split': 'final',
                'candidate_id': 'secret-positive', 'evidence_ids': ['answer']}
    assert candidate_features(clean) == candidate_features(poisoned)


def test_cold_start_and_single_class_are_explicit_and_finite():
    pool = [candidate(500), candidate(520, 'E')]
    model = NumericalSearch(seed=1)
    model.fit([], [])
    assert np.isfinite(model.predict(pool)[0]).all()
    model.fit(pool, [False, False])
    probability, uncertainty = model.predict(pool)
    assert np.isfinite(probability).all()
    assert np.isfinite(uncertainty).all()
    assert ((probability >= 0) & (probability <= 1)).all()


def test_rankings_exclude_measured_candidates_and_are_reproducible():
    pool = [candidate(i, 'K' if i % 2 else 'E') for i in range(500, 525)]
    model = NumericalSearch(seed=4, trees=20)
    model.fit(pool[:10], [bool(i % 2) for i in range(10)])
    a = model.shortlists(pool, {pool[12]['candidate_id']}, limit=5)
    b = model.shortlists(pool, {pool[12]['candidate_id']}, limit=5)
    assert a == b
    for values in a.values():
        assert len(values) == 5
        assert all(row['candidate_id'] != pool[12]['candidate_id'] for row in values)
        assert all(0 <= row['numerical_probability'] <= 1 for row in values)
