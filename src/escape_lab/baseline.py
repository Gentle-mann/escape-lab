"""Numerical active search using sequence-only features and disclosed labels."""

from collections import Counter
from operator import itemgetter

import numpy as np
from Bio.Align import substitution_matrices
from Bio.SeqUtils.ProtParamData import kd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression

BLOSUM = substitution_matrices.load('BLOSUM62')


def candidate_features(candidate):
    """Explicit allowlist: no scores, labels, split IDs or data-derived fields."""
    wt, mut = candidate['wt'], candidate['mut']
    position = candidate['aligned_position']
    features = {
        'gene': candidate['gene'], 'wt': wt, 'mut': mut, 'substitution': wt + mut,
        'position': position / 1000, 'blosum': float(BLOSUM[wt, mut]) / 10,
        'hydropathy_delta': (kd[mut] - kd[wt]) / 10,
        'charge_delta': (int(mut in 'KR') - int(mut in 'DE')
                         - int(wt in 'KR') + int(wt in 'DE')),
        'aromatic_delta': int(mut in 'FWY') - int(wt in 'FWY'),
        'proline_delta': int(mut == 'P') - int(wt == 'P'),
        'glycine_delta': int(mut == 'G') - int(wt == 'G'),
    }
    # Smooth position features transfer across held-out aligned blocks rather
    # than encode an exact residue identity learned from hidden assay values.
    for center in range(450, 811, 20):
        features[f'position_rbf_{center}'] = float(np.exp(-((position-center)/20)**2))
    for offset, residue in enumerate(candidate.get('context', '')):
        features[f'context_{offset}'] = residue
    return features


class NumericalSearch:
    def __init__(self, kind='extra_trees', seed=0, trees=100, leaf=3):
        if kind not in {'extra_trees', 'logistic'}:
            raise ValueError('Unknown numerical model')
        self.kind, self.seed, self.trees, self.leaf = kind, seed, trees, leaf
        self.vectorizer = DictVectorizer(sparse=False)
        self.model = None
        self.prior = 0.5
        self.gene_rates = {}
        self.fit_count = 0

    def fit(self, candidates, labels):
        if len(candidates) != len(labels) or any(type(x) is not bool for x in labels):
            raise ValueError('Training requires one disclosed boolean label per candidate')
        self.fit_count += 1
        self.prior = (sum(labels) + 1) / (len(labels) + 2)
        gene_n, gene_positive = Counter(), Counter()
        for candidate, label in zip(candidates, labels, strict=True):
            gene_n[candidate['gene']] += 1
            gene_positive[candidate['gene']] += int(label)
        self.gene_rates = {gene: (gene_positive[gene] + 10*self.prior)/(count+10)
                           for gene, count in gene_n.items()}
        self.model = None
        if len(set(labels)) < 2:
            return self
        x = self.vectorizer.fit_transform([candidate_features(c) for c in candidates])
        if self.kind == 'extra_trees':
            model = ExtraTreesClassifier(n_estimators=self.trees, min_samples_leaf=self.leaf,
                                         max_features=0.8, random_state=self.seed, n_jobs=1)
        else:
            model = LogisticRegression(C=1.0, max_iter=1000, random_state=self.seed)
        model.fit(x, np.asarray(labels, dtype=int))
        self.model = model
        return self

    def predict(self, candidates):
        if not candidates:
            return np.array([]), np.array([])
        if self.model is None:
            return np.full(len(candidates), self.prior), np.full(len(candidates), 0.5)
        x = self.vectorizer.transform([candidate_features(c) for c in candidates])
        probability = self.model.predict_proba(x)[:, 1]
        if self.kind == 'extra_trees':
            predictions = [tree.predict_proba(x)[:, 1] for tree in self.model.estimators_]
            uncertainty = np.std(predictions, axis=0)
        else:
            # This is a ranking diagnostic, not a calibrated epistemic interval.
            uncertainty = np.sqrt(probability * (1-probability))
        return probability, uncertainty

    def shortlists(self, pool, revealed_ids, limit=8, explored_groups=None):
        remaining = [c for c in pool if c['candidate_id'] not in revealed_ids]
        probability, uncertainty = self.predict(remaining)
        explored_groups = explored_groups or set()
        rows = []
        for candidate, p, u in zip(remaining, probability, uncertainty, strict=True):
            rows.append({**candidate, 'numerical_probability': float(p),
                         'model_disagreement': float(u),
                         'exploration_score': float(p + 0.5*u),
                         'gene_prior': float(self.gene_rates.get(candidate['gene'], self.prior))})
        exploit = sorted(rows, key=itemgetter('numerical_probability', 'candidate_id'), reverse=True)
        explore = sorted(rows, key=itemgetter('exploration_score', 'candidate_id'), reverse=True)
        diverse, selected_groups = [], set(explored_groups)
        for row in exploit:
            if row['group_id'] not in selected_groups:
                diverse.append(row)
                selected_groups.add(row['group_id'])
        chosen = {row['candidate_id'] for row in diverse}
        diverse.extend(row for row in exploit if row['candidate_id'] not in chosen)
        gene_prior = sorted(rows, key=itemgetter('gene_prior', 'candidate_id'), reverse=True)
        return {'exploit': exploit[:limit], 'explore': explore[:limit],
                'diverse': diverse[:limit], 'gene_prior': gene_prior[:limit]}
