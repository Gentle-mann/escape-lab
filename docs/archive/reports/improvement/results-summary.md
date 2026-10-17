> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# Escape Lab: measured autonomous improvement result

**The autonomous loop completed, but its selected revision performed worse on the
untouched final test. The predictive improvement goal was not achieved.**

Completed run: `20261017T194348Z-5e66cbab` on October 17, 2026.
Model: `nvidia/NVIDIA-Nemotron-3.5-Super-VL`, served by Crusoe.

## What was built and tested

The pipeline ingests the real FGFR assay dataset, fits separate pemigatinib and
futibatinib resistance predictors, and evaluates three fixed model families.
Nemotron then receives aggregate development errors and autonomously proposes
two bounded configurations. The numerical evaluator trains each proposal,
accepts or rejects it, and returns development feedback for the next revision.
A two-trial random search uses the same configuration space and folds.

Validation never enters the agent prompt. Models, source, data, and evaluation
rules are frozen before a one-shot final evaluation. Final predictions are
committed before the evaluator requests labels. This is an autonomous batch
pipeline; no researcher conversation is required during a run.

## Scores

Average precision (AP) measures ranking of resistant variants; higher is better.
It is not percent classification accuracy. Macro AP averages the two drug scores.

| Selection method | Development macro AP | Validation macro AP | Final macro AP |
|---|---:|---:|---:|
| Fixed baseline, selected from three families | 0.05999 | 0.05481 | **0.10923** |
| Two-trial random search | 0.06202 | 0.04805 | **0.10923*** |
| Two-proposal Nemotron loop | **0.07051** | **0.05722** | 0.09481 |

*Random search failed validation promotion and retained the original fixed
baseline, as specified in advance. Its final predictions are therefore identical
to the fixed baseline. They are not two independent confirmations. The development
and validation columns show its unsuccessful candidate before that fallback.

The agent's final macro AP was **0.01442 lower** than the fixed comparator. It did
not meet the preset success requirement of a gain of at least 0.05 over both
controls. It also breached the maximum allowed 0.02 per-drug AP deterioration.

| Final prediction target | Fixed AP | Agent AP | Agent minus fixed |
|---|---:|---:|---:|
| Pemigatinib | 0.14194 | 0.11783 | −0.02411 |
| Futibatinib | 0.07653 | 0.07180 | −0.00473 |

Probability quality worsened too. Brier score measures squared probability
error, so lower is better: pemigatinib rose from **0.08560 to 0.12490**;
futibatinib rose from **0.05943 to 0.06482**. A development AP gain does not imply
better calibrated probabilities, and Brier was a diagnostic rather than a
promotion criterion in this protocol.

The final set contained 1,819 protein variants across seven aligned sequence
blocks. After excluding inconclusive labels separately, pemigatinib had 1,802
scorable variants with 152 positives; futibatinib had 1,807 with 108 positives.
The constant-score AP references are approximately 0.08435 and 0.05977. These
contextual references do not turn the experiment into a clinically accurate test.

## What the agent actually changed

1. **Revision 1: accepted on development.** Keep logistic regression and the
   existing features; change inverse regularization `C` from 1 to 10. Development
   macro AP increased from 0.05999 to 0.07051. Validation also increased slightly,
   from 0.05481 to 0.05722, so the frozen promotion rule admitted it to final testing.
2. **Revision 2: rejected.** Keep `C=10` and widen the positional RBF feature scale
   from 20 to 40. Development macro AP fell to 0.06128, so revision 1 was retained.

Both proposals were generated and numerically scored, with no baseline substituted
for a failed request in the completed run. The agent's public rationales are
unverified explanations: revision 2 incorrectly described the RBF width as
position bins and referred to a floor relative to its incumbent. The evaluator
used the actual feature implementation and the unchanged fixed-reference floor.
Those explanatory mistakes did not change scoring or acceptance.

The observed pattern is consistent with a development improvement that did not
generalize. It does not establish a universal conclusion about LLMs or biology.
The fixed comparator is a development-selected baseline, not an exhaustive AutoML
optimization. A larger ordinary hyperparameter search could also find `C=10`.

## Runtime, operational failures, and verification

The completed scientific run took **100.74 seconds** and two inference calls,
using **19,277 provider-reported tokens**. Both calls completed with valid tool
arguments and required no format repair.

Three preceding real runs failed operationally before any agent configuration
could be scored. They took 94.79, 71.68, and 85.97 seconds, totaling **252.44 seconds**.
All four real attempts therefore consumed **353.18 seconds of run time**, excluding
setup, debugging, and synthetic diagnostics. Measured usage across the four was
68,920 tokens, plus unknown actual usage for the first two failed requests.
Their 21,059-token conservative budget reserve is an estimate, not billed usage.
Three separate synthetic schema probes took another 37.85 seconds and 14,954
reported tokens. Earlier smoke tests and engineering time are not included.

The preserved [retry record](inference-retry-note.md) documents the pinned model,
larger reasoning/output allowance, and inlined tool schema that resolved the
observed compatibility failure. Scientific data, folds, feature choices, seeds,
two-configuration budgets, and acceptance rules were not changed in response to
validation scores. The first failed run incidentally evaluated validation; its
inherited-baseline agent row is not an agent prediction result. Reuse during
operational recovery is not independent validation.

Verification of the completed implementation:

- **140 tests passed; Ruff passed.** Tests cover strict configuration parsing,
  group separation, failed proposals, scoring, and final-evaluation boundaries.
- Every required control evaluation and both agent evaluations completed.
- Frozen source/data/protocol hashes and four distinct fitted model artifacts
  verified after completion.
- The completed final receipt matches the frozen manifest and the committed
  prediction checksum.
- Independent code and aggregate-result reviews found no material scoring error.

## What this means for the next iteration

The engineering mechanism is working; the biological prediction result failed.
The final holdout is now consumed, and the CLI blocks fresh tuning against this
same dataset state. Do not delete the receipt or retry until a favorable final
score appears.

A new research iteration should investigate the weak biological representation,
probability quality, and cross-block generalization using development evidence,
then establish a fresh untouched evaluation before claiming improvement. The
current test cannot serve as that fresh evidence. Do not present this run as a
cancer discovery, a treatment recommendation, or a successful predictive advance.

Only seven final blocks from one retrospective assay dataset were available.
The recorded paired-block bootstrap is descriptive sensitivity analysis, not
independent laboratory or patient replication. Public-data pretraining exposure
cannot be ruled out.

## Evidence

- Machine report and complete metrics (`20261017T194348Z-5e66cbab/report.json`; local artifact)
- [Run report](20261017T194348Z-5e66cbab/report.md)
- Predeclared protocol (`20261017T194348Z-5e66cbab/protocol.json`; local artifact)
- First proposal (`20261017T194348Z-5e66cbab/proposal-1.json`; local artifact)
- Second proposal (`20261017T194348Z-5e66cbab/proposal-2.json`; local artifact)
- Frozen selection and checksums (`20261017T194348Z-5e66cbab/freeze.json`; local artifact)
- [Final aggregate scores](20261017T194348Z-5e66cbab/final-results.json)
- [Current protocol documentation](../../docs/revised-prediction-protocol.md)
