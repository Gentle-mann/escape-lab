> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Autonomous FGFR prediction run

Run: `20261017T194348Z-5e66cbab`
Outcome: **prototype_margin_not_met**

## Development and validation

| Method | Development macro AP | Validation macro AP | Promoted |
|---|---:|---:|---|
| fixed | 0.0600 | 0.0548 | False |
| random | 0.0620 | 0.0481 | False |
| agent | 0.0705 | 0.0572 | True |

## Frozen final evaluation

| Method | Pemigatinib AP | Futibatinib AP | Macro AP |
|---|---:|---:|---:|
| fixed | 0.1419 | 0.0765 | 0.1092 |
| random | 0.1419 | 0.0765 | 0.1092 |
| agent | 0.1178 | 0.0718 | 0.0948 |

## Attempts

- fixed 0: completed; accepted=False.
- fixed 1: completed; accepted=False.
- fixed 2: completed; accepted=False.
- random 1: completed; accepted=False.
- random 2: completed; accepted=True.
- agent 1: completed; accepted=True. Incumbent logistic (full features, gene included) is the best tried head (macro AP 0.0600; fut 0.0502, pem 0.0697). Diagnosis: fut is the constraint-binding, weaker head and is weakest by gene on FGFR3 (AP 0.031) and FGFR2 (0.038) despite adequate n, while fold-2 collapses (fut AP 0.025, only 10 positives) from positive scarcity. Proposed single change: keep the incumbent logistic/full/gene setup but weaken the penalty (regularization C=10 vs 1) so gene/sequence coefficients for fut on FGFR2/FGFR3 can express more strongly. Expected tradeoff: higher variance on the ultra-sparse fold-2, but the per-drug floor gives fut headroom (0.0502 vs min allowed 0.0302) and pem is comfortably above its floor (0.0697 vs 0.0497). All other keys held at incumbent/inactive defaults to isolate the change; no labels, thresholds, splits, or comparator are altered.
- agent 2: completed; accepted=False. Incumbent (logistic, full features, gene included, C=10) is the best tried: macro AP 0.0705, fut 0.0662, pem 0.0748. Less regularization (C=10) beat C=1 and both tree variants on both heads, so I keep C=10 (the max allowed) and the full+gene representation. Fut is the limiting head; its constraint floor is 0.0462. The weakest cells are fold 2 (only 10 fut / 20 pem positives) and FGFR3 (fut 0.0416, pem 0.0435), i.e., sparse positional signal. The single change is position_scale 20->40: coarser position bins pool sparse mutation counts for more stable per-position estimates, aiming to lift fut/FGFR3 without removing gene or sequence information. Tradeoff: coarser bins sacrifice fine positional resolution (could slightly blunt sharp position-specific effects) in exchange for lower variance in sparse folds; expected fut drop is well under the 0.02 allowance, with pem roughly unchanged.
- fixed validation: completed; accepted=False.
- random validation: completed; accepted=False.
- agent validation: completed; accepted=False.

Elapsed: 100.7 seconds. Within budget: True.

The ordinary search uses two random configurations, not an exhaustive AutoML search. All methods share training data, grouped folds, and evaluation rules. Matched trial budgets are not identical hardware time; actual timings are recorded.

Sequence blocks are not independent tumors or laboratories. Results are retrospective assay predictions, not treatment recommendations. Public-data pretraining exposure cannot be ruled out. No claim is made against unavailable historical CAGI submissions.
