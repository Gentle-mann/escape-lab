> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# Experiments 1 and 2: completed 17 October 2026

**Neither experiment demonstrated improvement under the predeclared criteria.** Both formal verdicts are inconclusive. You.com evidence gave a very small positive difference versus no evidence and essentially tied unrelated evidence. The no-evidence agent performed worse than fixed and random selection on every repeat; its small positive mean versus the numerical policy was insufficient.

The batch completed 90 isolated searches: ten repeats across three grouped outer folds and three evidence arms. Every search had three proposal opportunities. All selected models and unlabeled predictions were committed together before evaluator scoring. Proposal agents received only their own inner-training feedback and assigned evidence cards. They never received evaluation outcomes, other agents' history, raw papers, or filesystem/web tools. An independent audit checked all 270 requests.

The cohort contains 3,730 variants in 15 aligned groups. These are reused development data. The ten repeats measure procedural variability on one fixed split, not ten independent biological validations. No final or external panel was scored.

The primary score is the mean of ten repeat-level pooled out-of-fold macro average precisions (AP), averaging pemigatinib and futibatinib. Higher AP is better; lower Brier error is better. AP is a ranking metric, not classification accuracy.

| Method | Macro AP | Pemigatinib AP | Futibatinib AP | Macro Brier |
|---|---:|---:|---:|---:|
| Agent, no evidence | 0.02897045 | 0.03962548 | 0.01831541 | 0.03503282 |
| Agent + You.com evidence | 0.02914259 | 0.03977374 | 0.01851144 | 0.03504838 |
| Agent + unrelated evidence | 0.02914202 | 0.03978038 | 0.01850365 | 0.03504858 |
| Fixed structural model | 0.03462765 | 0.04922321 | 0.02003209 | 0.03578062 |
| Random selection (3 tries) | 0.03226587 | 0.04389715 | 0.02063458 | 0.03461819 |
| Numerical selection (3 tries) | 0.02880030 | 0.03945396 | 0.01814664 | 0.03501657 |
| Full grid (11 tries) | 0.03493831 | 0.04756983 | 0.02230679 | 0.03459281 |

**Experiment 1: relevant evidence versus no evidence and placebo.** Relevant evidence improved macro AP by only +0.000172143 over control and +0.000000571 over placebo. Both fall below the predeclared +0.002 minimum, and both adjusted uncertainty intervals include zero. Relevant/control wins, ties and losses were 6, 1 and 3; relevant/placebo were 4, 4 and 2. This does not establish a useful benefit from the two reviewed You.com cards.

**Experiment 2: agent selection versus simpler selection.** The no-evidence agent was specified as the primary agent before execution. It scored 0.02897045, below fixed structural selection (0.03462765) and random selection (0.03226587) in all ten repeats. The agent-minus-random adjusted interval is entirely negative, supporting worse performance in this conditional development comparison. The agent-minus-fixed adjusted interval includes zero, so the repeated losses are descriptive rather than proof of general inferiority. Its +0.000170146 mean advantage over the numerical policy comes largely from one better repeat and does not meet the gate.

| Predeclared contrast | Absolute macro AP difference | Adjusted 99% interval | Improvement gate passed |
|---|---:|---:|---|
| Agent + You.com evidence minus Agent, no evidence | +0.000172143 | [-0.000803, +0.001742] | No |
| Agent + You.com evidence minus Agent + unrelated evidence | +0.000000571 | [-0.001373, +0.001485] | No |
| Agent, no evidence minus Fixed structural model | -0.005657204 | [-0.014721, +0.001738] | No |
| Agent, no evidence minus Random selection (3 tries) | -0.003295420 | [-0.010566, -0.000564] | No |
| Agent, no evidence minus Numerical selection (3 tries) | +0.000170146 | [-0.000175, +0.001166] | No |

Intervals use 10,000 paired bootstrap draws, resampling aligned groups and repeat IDs with shared draws across methods. The five primary contrasts use 99% percentile intervals for the predeclared multiplicity adjustment. All per-drug AP and Brier safeguards passed; the failed requirements were minimum gain and a positive lower interval bound. The intervals condition on the fitted predictions and do not include retraining or independent-study uncertainty.

The full grid is a larger-budget reference, with eleven new evaluations. Agent, random and numerical selection each had three opportunities beyond the common reference. Fixed structural selection had no adaptive search cost. The numerical arm tests the specified coarse-to-fine policy, not all numerical optimizers. Numerical model fits were cached across identical configurations, so this batch's wall time is not a fair standalone comparison of optimizer cost.

**What changed in the decisions.** Relevant evidence changed 9 of 30 fold-level choices compared with control and 6 of 30 compared with placebo. All differences occurred in outer fold 1. In folds 0 and 2, every agent retained the common chemistry model in every repeat. Evidence changed some choices but did not produce a reliable prediction gain. Exact choices are preserved in choices.json and summarized in the audit JSON.

**Execution and integrity.** All 90 searches entered scoring, including failures. Of 270 proposal opportunities, 264 produced valid proposals; five responses were truncated and one control request failed with a connection error. No scientific result triggered a rerun. Relevant produced 90/90 valid proposals, control 88/90 and placebo 86/90. All 314 run-artifact hashes, 388 joint-commitment entries, 77 protocol-input hashes and 292 preserved prior-artifact hashes verified. The commitment contains all search logs, three prediction banks and 28 fitted-model files. The pre-run suite passed 583 tests and Ruff; the final focused check passed 41 tests. Frozen source hashes remained unchanged throughout execution.

The trusted evaluator reads a monolithic local outcome archive internally and returns only the requested partition. Thus the guarantee is isolation of proposal-agent inputs and delayed evaluation scoring, not an assertion that the trusted process never parsed other archive bytes.

| Arm | Searches | Valid proposals | Failed opportunities | Tokens with reported usage | Estimated inference cost |
|---|---:|---:|---:|---:|---:|
| Agent, no evidence | 30 | 88/90 | 2 | 940,922 | $2.67732 |
| Agent + You.com evidence | 30 | 90/90 | 0 | 1,049,996 | $3.10683 |
| Agent + unrelated evidence | 30 | 86/90 | 4 | 1,020,087 | $2.97578 |

Runtime was 21.58 minutes. Reported usage covers 3,011,005 tokens and implies $8.75993 inference cost at the frozen estimate of $1.40 per million input tokens and $4.40 per million output tokens. One connection-error request has unknown provider usage; its budget ledger conservatively charged 29,264 estimated tokens without an input/output split. The cost table excludes that unknown usage and cache discounts. The new placebo retrieval used four You.com searches and four requested Contents pages over two extraction calls, estimated at $0.024. Existing relevant evidence was reused. Combined known-usage inference plus new retrieval is approximately $8.78; this is an estimate, not a complete invoice.

This experiment tests two fixed evidence cards within a twelve-configuration catalog. It does not test unrestricted literature synthesis, new biological features, or transfer to a fresh study. The results justify neither claiming research-assisted improvement nor claiming that research can never help. They provide no evidence to replace the simpler baselines with this agent-selection policy.

- [Machine results](20261017-experiments-1-2/results.json)
- Frozen protocol (`reports/repeated-research/20261017-experiments-1-2/protocol.json`; local artifact)
- [Human-readable protocol](../../docs/repeated-research-protocol.md)
- Joint prediction commitment (`reports/repeated-research/20261017-experiments-1-2/joint-freeze.json`; local artifact)
- Selections (`reports/repeated-research/20261017-experiments-1-2/choices.json`; local artifact)
- [Verification and choice summary](20261017-experiments-1-2-audit.json)
- Relevant evidence dossier (`research/you-evidence/you-pilot-20261017/evidence.json`; local artifact)
- Placebo evidence dossier (`research/repeated-research-placebo-evidence.json`; local artifact)
