> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Autonomous FGFR prediction run

Run: `20261017T193042Z-7992cb52`
Outcome: **no_agent_validation_improvement**

## Development and validation

| Method | Development macro AP | Validation macro AP | Promoted |
|---|---:|---:|---|
| fixed | 0.0600 | 0.0548 | False |
| random | 0.0620 | 0.0481 | False |
| agent | 0.0600 | 0.0548 | False |

Final test remains sealed. No complete, validation-improving agent revision qualified, or final evaluation was disabled/budget-limited.

## Attempts

- fixed 0: completed; accepted=False.
- fixed 1: completed; accepted=False.
- fixed 2: completed; accepted=False.
- random 1: completed; accepted=False.
- random 2: completed; accepted=True.
- agent 1: proposal_failed; accepted=False. Crusoe request failed; no pipeline was proposed.
- agent 2: proposal_failed; accepted=False. Crusoe request failed; no pipeline was proposed.
- fixed validation: completed; accepted=False.
- random validation: completed; accepted=False.

Elapsed: 94.8 seconds. Within budget: True.

The ordinary search uses two random configurations, not an exhaustive AutoML search. All methods share training data, grouped folds, and evaluation rules. Matched trial budgets are not identical hardware time; actual timings are recorded.

Sequence blocks are not independent tumors or laboratories. Results are retrospective assay predictions, not treatment recommendations. Public-data pretraining exposure cannot be ruled out. No claim is made against unavailable historical CAGI submissions.
