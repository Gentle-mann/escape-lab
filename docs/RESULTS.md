# Escape Lab: six hours of building and experiments

In approximately six hours, we built a research pipeline around real FGFR data: autonomous experiment proposals, numerical models, independent verification, and an auditable record of every revision. We integrated Nemotron and GLM-5.3 on Crusoe, added reviewed evidence from You.com, and tested sequence and structural representations.

## Measured ranking improvement

On **668 FGFR2 variants** across three fixed development folds, we improved high-activation ranking while keeping the same features, fitted coefficients and evaluation budget.

| Ranking | Mean-fold AP ↑ | Highly activating variants found in 68 picks | Strongly activating variants found |
|---|---:|---:|---:|
| Original sequence classifier probabilities | 0.09732 | 1 / 26 | 0 / 7 |
| **Same classifier, raw decision margin** | **0.11135** | **4 / 26** | 0 / 7 |
| Historical biology tree reference | 0.07186 | 3 / 26 | 0 / 7 |

**AP increased by 0.01402, a 14.4% relative gain. Three additional highly activating variants were recovered within the same budget.** AP measures ranking quality.

The improvement came from a numerical repair we implemented: 136 probabilities in the first fold had rounded to exactly 1.0. Ranking by the classifier's raw decision margin restored their ordering. The other two folds were unchanged.

The verification run took **2.92 seconds with cached features**, reproduced the original probabilities exactly, and preserved all **594 protected artifacts**. [Implementation](../src/escape_lab/fgfr2_sequence_ap.py) · [Measured results](archive/reports/mechanistic-pilot/fgfr2-sequence-ap-v1/report.md).

## Agent experiments

We also tested whether changing the agent or adding literature evidence improved **drug-resistance prediction**.

| Experiment | Comparator | Agent result |
|---|---:|---:|
| Original final resistance test, macro AP | Fixed model 0.10923 | Nemotron 0.09481 |
| Revised internal evaluation, macro AP | Fixed structure 0.03463; random/grid 0.03494 | Nemotron and GLM-5.3: 0.03052 |
| S13 biochemical panel, mean Spearman | Fixed structure 0.40568 | Nemotron 0.20892; GLM-5.3 0.29695 |
| You.com evidence pilot, reused S13 | GLM without evidence 0.29695 | GLM with evidence 0.29695 |
| Experiment 1: repeated evidence, mean macro AP | No evidence 0.02897045; placebo 0.02914202 | You.com evidence 0.02914259 |
| Experiment 2: repeated selection, mean macro AP | Fixed 0.03462765; random 0.03226587; numerical 0.02880030 | No-evidence agent 0.02897045 |

**Experiments 1 and 2 are complete; both are inconclusive under their predeclared improvement gates.** Relevant evidence improved AP by only +0.000172143 over no evidence and essentially tied placebo. The no-evidence agent scored below fixed and random selection in all ten repeats; its small gain over the numerical policy was insufficient. These comparisons did not establish an agent advantage over the numerical controls.

All 90 searches entered scoring, with 264 valid proposals and six failed opportunities across 270 requests. Proposal inputs were isolated from evaluation outcomes, and models and unlabeled predictions were committed before scoring. The study reused 3,730 development variants across fixed folds; its ten repetitions are not independent biological validations. No final or external panel was scored. See the [completed review](archive/reports/repeated-research/20261017-experiments-1-2-review.md) and [complete experiment ledger](EXPERIMENTS.md#3-youcom-evidence-and-repeated-selector-studies) for uncertainty, failures, costs and evaluation limits.

## Reported external-score references

Drug-resistance ranking scores for pretrained variant-effect methods:

| Method | Pemigatinib AP | Futibatinib AP | Average AP |
|---|---:|---:|---:|
| [EVE](https://www.nature.com/articles/s41586-021-04043-8) | 0.0894 | 0.0648 | 0.0771 |
| [REVEL](https://pmc.ncbi.nlm.nih.gov/articles/PMC5065685/) | 0.0671 | 0.0499 | 0.0585 |
| [AlphaMissense](https://deepmind.google/research/publications/21083/) | 0.0627 | 0.0500 | 0.0564 |

These are resistance reference scores; our latest 0.11135 measures activation. The reference scoring protocol is not archived here, so a direct comparison requires matching the endpoint, variants and evaluation procedure.

## Evaluation notes and next steps

The activation gain is a development result from a scoring repair, with a descriptive paired interval of approximately **[0, 0.03261]**. Feature extrapolation remains unresolved; earlier ridge models had higher AP but unstable continuous predictions. Details are in the [experiment ledger](EXPERIMENTS.md#5-fgfr2-sequence-and-high-activation-ranking).

The original resistance final set and S13 panel have been evaluated. RLY and MANO remain reserved for an independent, assay-compatible test. Next steps are stronger biological representations and external validation. Private specialized models may perform better; they have not been evaluated here.

Preclinical research prototype. No patient-level or treatment recommendations.
