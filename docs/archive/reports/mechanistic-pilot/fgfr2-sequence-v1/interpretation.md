> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Interpretation of the FGFR2 sequence pilot

The primary hypothesis did not pass. Frozen ESM2 features improved average
precision relative to the matched biology-only classifier, but the primary
model retrieved fewer high-activation variants than the previous best ranking
within the experimental nomination budget. It should not replace that reference.

This was a fixed development comparison on 668 FGFR2 variants, ten occupied
aligned-position blocks and the same three global folds used previously. The
literature-derived primary endpoint was replicate consistency plus median
enrichment >5, adapted to protein-level replicate summaries. There were 26
positives. All arms nominated 28, 18 and 22 variants in the three folds, totaling
68. Random selection at those budgets would recover about 2.65 positives on
average; this is an expectation, not a fitted or simulated competing model.

| Arm | Mean fold AP | High positives recovered /26 | Strong positives recovered /7 |
|---|---:|---:|---:|
| Previous enriched regression | 0.0719 | 3 | 0 |
| FGFR2 biology-only classifier | 0.0471 | 0 | 0 |
| ESM + biology classifier, primary | 0.0973 | 1 | 0 |
| ESM + biology regression, secondary | 0.1171 | 4 | 0 |

Average precision evaluates the ranking across thresholds. It can rise while
the number recovered at a particular budget falls. The primary model improved
AP in only one of three folds; its increase over the previous ranking was
0.02546, below the predeclared 0.03 gate, while retrieval fell from three to one.
The descriptive paired block-bootstrap AP difference interval versus the old
ranking spans -0.0687 to +0.2677. It also spans zero versus the matched
biology-only classifier. These are uncertain signals on reused development
blocks, not independent confirmation.

The secondary regression arm is a weak follow-up lead: its AP is higher and it
recovers one additional positive. It still misses 22 of 26 high activators,
misses all seven strongest activators, and would fail the three-extra-recoveries
criterion even if selected after the fact. There is no matched biology-only
ridge arm, so the regression result cannot isolate a causal benefit from ESM
features. It must not be promoted as a newly selected winning primary model.

The clearest unresolved failure is pilot-2: all four arms recover zero of its
12 higher-activation variants within the top 22 nominations. Thus the new
sequence representation has not solved transfer to difficult positions. This
run does not identify a unique cause. Limited high-activation supervision,
feature compression, the linear prediction heads, and a missing representation
of activation mechanisms remain plausible explanations.

The research motivated mutation-position WT context and signed mutant-minus-WT
context rather than treating evolutionary plausibility as an activation score.
The implementation uses frozen ESM2-35M, independent train-only PCA16 for the
WT and delta branches, existing biological features, and fixed regularized
classification/regression heads. No large-model fine-tuning or LLM calls were
used. See [research and protocol](../../../research/fgfr2-sequence-pilot/research-and-protocol.md)
and [interpretation boundaries](../../../research/fgfr2-sequence-pilot/interpretation-boundaries.md).

Do not spend the untouched external datasets on this result. A scientifically
distinct next hypothesis would retain a continuous activation target and add
carefully audited active-versus-inactive structural features, with a matched
head ablation. The existing sequence cache can be reused. Such structural
features would remain proxies, not a physical simulation or a validated
prediction of cancer growth or drug resistance.

Operationally, extraction took 49.62 seconds after first-time dependency and
checkpoint setup. Fitting/scoring all nine new fold-model combinations took
0.64 seconds; the evaluator including ingestion, checks and uncertainty took
4.60 seconds. First-time downloads were slow and are separately documented in
the runtime benchmark. All 399 scoped tests passed; all 518 protected artifacts
remained unchanged. Independent arithmetic reproduced all twelve arm/fold AP
and retrieval values. No post-score model changes, external outcome evaluation,
original final evaluation, or resistance-label use occurred.
