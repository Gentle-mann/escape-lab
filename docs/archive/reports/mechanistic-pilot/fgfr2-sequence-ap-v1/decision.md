> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Decision: preserve raw classifier margins for AP ranking

The single targeted test improved mean-fold AP from 0.09732 to 0.11135 and high-activation retrieval from 1 to 4 of 26 at the same 68-pick budget. The original probabilities reproduced exactly for all 668 variants. The sole change was ranking by the unchanged classifier decision score before sigmoid rounding. Retain this ranking method for the development AP evaluator; do not describe margins as calibrated probabilities.

The original matched sequence ablation improved mean AP from 0.04709 to 0.09732, mostly through pilot-1 (0.05792 to 0.21570). Thus sequence representations helped mainly in one existing development fold. Comparison with the pooled tree also changes objective and training cohort, so it does not isolate a sequence effect.

Pilot-0 had 136 predictions rounded to 1.0. Raw margins resolved this tie, improved its AP from 0.02463 to 0.06671 and recovered three additional high activators. The other folds were unchanged. No fitting settings, input features, folds or labels changed.

This is a technical ranking correction, not biological model improvement or independent validation. Unsafe feature extrapolation remains: some raw decision margins exceed one million. The revised model still finds only 4 of 26 high activators and no strongest activators. Its descriptive reused-block AP-gain interval is approximately 0.000 to 0.033; the boundary at zero and one-fold concentration limit the conclusion.

Three fits/scoring took 0.31 seconds, evaluator total 2.92 seconds. All 29 focused tests passed, all 594 protected files were unchanged, and separate manual AP/retrieval and weighted-bootstrap arithmetic verified the saved results. The same agent performed that arithmetic; no independent-agent review is claimed. Reserved external assay outcomes and original final evaluation were not used. Stop after this one sensitivity; no tuning sweep was performed.

See report.md and verification.json for full results and receipts.
