"""Read-only, aggregate replay of published fixed-development activation experiments.

These are reused development results. Replay executes no fits, reads no row-level
outcomes, and cannot be presented as a new unseen evaluation or fresh discovery.
"""

import hashlib
import json
from pathlib import Path

import numpy as np

from escape_lab.activation_planner import _safe_result
from escape_lab.planner_diagnostics import aggregate_fit_diagnostics, sanitize_diagnostics


class DevelopmentReplay:
    scope = 'fixed-development-replay'
    metric = 'mean_fold_high_activation_ap'
    new_fits = 0

    def __init__(self, model_report, ranking_report, *, source_hashes):
        metrics = model_report['metrics']
        ranking = ranking_report['metrics']
        self.source_hashes = dict(source_hashes)
        # Reject inconsistent fold/count summaries; aggregate reports cannot prove variant identity.
        expected = ranking['margin']['by_fold']
        for name in ('legacy_both', 'legacy_wt', 'legacy_delta'):
            observed = metrics[name]['by_fold']
            if set(observed) != set(expected):
                raise ValueError('Replay folds differ across reports')
            for fold, row in observed.items():
                target = expected[fold]
                if (row['n'], row['budget'], row['high']['positives']) != (
                        target['n'], target['budget'], target['high']['positives']):
                    raise ValueError('Replay case counts or endpoint differ')
        if not np.isclose(metrics['legacy_both']['mean_fold_high_ap'],
                          ranking['margin']['mean_fold_high_ap'], rtol=0, atol=1e-12):
            raise ValueError('Baseline reproduction differs from ranking repair')
        comparisons = {'probability_ap': ranking['probability']['mean_fold_high_ap'],
                       'ranking_ap': ranking['margin']['mean_fold_high_ap']}
        comparisons['delta_ap'] = comparisons['ranking_ap'] - comparisons['probability_ap']
        ablation = {'combined_ap': metrics['legacy_both']['mean_fold_high_ap'],
                    'wt_only_ap': metrics['legacy_wt']['mean_fold_high_ap'],
                    'delta_only_ap': metrics['legacy_delta']['mean_fold_high_ap']}
        if len(_safe_result('compare_rankings', comparisons)) != 3:
            raise ValueError('Invalid replay comparison metrics')
        if len(_safe_result('request_sequence_ablation', ablation)) != 3:
            raise ValueError('Invalid replay ablation metrics')
        scores, scales = [], []
        for row in model_report.get('diagnostics', []):
            if row.get('arm') != 'legacy_both':
                continue
            n = row.get('n')
            if type(n) is not int or n <= 0:
                raise ValueError('Invalid aggregate diagnostic count')
            # Archives count exact ones but not exact zeros. Do not pretend these are all ties.
            scores.append(sanitize_diagnostics({
                'score_count': n,
                'probability_exact_one_fraction': row['probability_exact_one'] / n,
                'probability_unique_fraction': row['probability_unique'] / n,
                'ranking_score_abs_max': row['margin_abs_max']}))
            scale = row.get('model', {}).get('biology_scaling', {})
            scales.append(sanitize_diagnostics({
                'minimum_training_scale': scale.get('minimum_training_scale'),
                'below_floor_feature_count': scale.get('below_floor_count')}))
        fit = aggregate_fit_diagnostics(scores)
        feature = aggregate_fit_diagnostics(scales)
        self.results = {'inspect_score_saturation': fit, 'inspect_feature_shift': feature,
                        'compare_rankings': comparisons, 'request_sequence_ablation': ablation}
        self.executed = set()

    @property
    def available_actions(self):
        return set(self.results)

    def execute(self, action):
        if action not in self.results or action in self.executed:
            raise ValueError('Unavailable or duplicate replay tool')
        self.executed.add(action)
        return dict(self.results[action])


def load_development_replay(root):
    root = Path(root)
    paths = {
        'model_improvements': root / 'docs/archive/reports/model-improvements/results.json',
        'ranking_repair': root / 'docs/archive/reports/mechanistic-pilot/fgfr2-sequence-ap-v1/results.json',
    }
    payloads, hashes = {}, {}
    for name, path in paths.items():
        data = path.read_bytes()
        hashes[name] = hashlib.sha256(data).hexdigest()
        payloads[name] = json.loads(data)
    return DevelopmentReplay(payloads['model_improvements'], payloads['ranking_repair'],
                             source_hashes=hashes)
