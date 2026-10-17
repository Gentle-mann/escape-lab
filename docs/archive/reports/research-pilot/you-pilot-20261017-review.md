> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# You.com pilot result

The integration ran successfully, but this small paired pilot showed **no improvement**.
Both GLM-5.3 arms selected `interactions-logistic-c0.1` and produced identical predictions.

| Measure | Control | You.com evidence |
|---|---:|---:|
| Mean external Spearman correlation | 0.29695 | 0.29695 |
| Pem correlation | 0.50856 | 0.50856 |
| Fut correlation | 0.08535 | 0.08535 |
| Inner selection AP | 0.04217 | 0.04217 |
| Valid proposals / opportunities | 2 / 3 | 2 / 3 |
| Inference tokens charged | 33,899 | 33,960 |

The research-minus-control difference is exactly zero for these committed predictions.
Its paired bootstrap interval is also [0, 0] because the prediction vectors are identical;
this is not evidence that research can never help. The fixed structure baseline scored
0.40568, above both agents. Inner AP is a selection diagnostic, not a separate test.
Both arms' third responses reached the 8,192-token completion limit. No scientific
restart or outcome-based retry occurred.

## What You.com contributed

Four frozen searches returned 20 results. One Contents request selected four source
pages; three were returned. Two reviewed general biological cards passed the gate:

- Solvent environment: [Barelier et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23874896).
- Limits of static structures for kinetics: [Wong, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10926948/).

The conformational-plasticity page was not returned. The atomic-interaction page
supported its proposed claim but lacked its own publication date in the retrieved
text, so both cards were omitted. Missing evidence was not fabricated or replaced.
The research hypotheses referenced the reviewed solvent/burial context, but the
sequence of evaluated configurations and final model matched the control.

## Evaluation boundary and verification

Both fresh clients had the same prompt, catalog, numerical seeds, budgets and scoring
rules. Control received an empty evidence list. Research received only the two fixed
claims, canonical citations and years. The six saved request payloads were audited:
no raw source text, supporting excerpts, test outcomes or evaluator artifacts entered
the proposal context. Only aggregate inner-CV feedback was available during selection.

Both arms' model files and predictions were jointly hash-committed before S13 scoring.
All 50 agent prediction entries per arm are identical. Verified 75 run-artifact hashes,
56 joint-commitment hashes, and 394 preserved prior artifacts. All four numerical
baselines reproduced the original results to tolerance 1e-12. No credential appeared
in the new research or run artifacts. The integration source remains byte-identical
to the version that passed 408 tests and Ruff; the nine pilot tests also passed in
an independent execution audit.

This was one stochastic pair using an already spent retrospective biochemical panel.
It tests two curated cards within a fixed model catalog, not unrestricted literature
synthesis or new biological features. It cannot establish an average benefit or
absence of benefit. No fresh holdout was consumed.

The pilot took 179.8 seconds. Estimated You.com cost was
$0.024; inference was $0.19627; total **$0.22027**.
These are usage-based estimates, not invoices, and inference excludes cache discounts.

Artifacts:

- [Scored run summary](you-pilot-20261017/results-summary.md)
- Frozen protocol (`reports/research-pilot/you-pilot-20261017/protocol.json`; local artifact)
- Joint prediction commitment (`reports/research-pilot/you-pilot-20261017/joint-freeze.json`; local artifact)
- Released evidence cards (`research/you-evidence/you-pilot-20261017/evidence.json`; local artifact)
- Curation decisions (`research/you-evidence/you-pilot-20261017/curation-decisions.json`; local artifact)
