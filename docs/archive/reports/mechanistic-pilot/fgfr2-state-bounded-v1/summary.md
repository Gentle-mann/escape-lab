> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# FGFR2 active versus inactive structures: pilot outcome

The static active-minus-inactive representation gives a small ranking improvement over an equally corrected ridge baseline. It does not beat our stronger historical tree model and does not establish generalization. External assay outcomes and the original final evaluation were not opened for this work.

## Structures actually used

The literature makes the naive crystal comparison unreliable: the 2017 mechanistic study describes 2PSQ as partially active. The 2023 PNAS workflow instead built inactive FGFR2 from FGFR1 3KY2 and prepared active FGFR2 from 2PVF. We used the authors’ prepared starting coordinates, extracted without reading trajectories or mutant archives. Both have the same C491A background, ATP, and unphosphorylated activation-loop tyrosines. These are structural proxies; the inactive state is a homology model.

Only common residue and heavy-atom identities were used. Position 491 and uncertain modeled loops were masked in both states. This leaves 279 positions and 2,236 common heavy atoms, covering 528 of the 668 evaluated variants. All 140 unsupported variants remain in the cohort with zero geometry and an explicit availability flag. Masking the engineered site does not remove possible longer-range effects of that preparation.

Features encode contact environments, side-chain packing, mutation chemistry and proximity to the molecular brake/DFG/alphaC motifs. The contrast subtracts inactive features from active features. This is not a mutant energy calculation, relaxed structure, catalytic simulation or cellular growth model.

Primary sources: [2017 eLife](https://elifesciences.org/articles/21137), [2023 PNAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC9974468/), [authors’ structural archive index](https://yzhang.hpc.nyu.edu/MD_data_FGFR/). [2024 paired-state kinase work](https://www.nature.com/articles/s41467-024-50812-0) motivates the hypothesis but uses substantially richer computations and does not validate these FGFR2 descriptors.

## Why two runs were necessary

The first fixed run faithfully reproduced the old ridge predictions, then failed both gates. Continuous errors exposed a pre-existing extrapolation problem hidden by ranking-only evaluation: one localized position feature had a training SD of 1.17e-8, so an ordinary held-out value became 78.8 million after standardization. Baseline mean-fold MAE was 3,933 and contrast MAE was 21,495 on the log activation scale. Those outputs are unusable for continuous activation estimation. The complete failed run is preserved in ../fgfr2-state-v1/.

We then performed exactly one explicitly post-result numerical-stability sensitivity. Every biology feature was limited to its own training minimum/maximum at inference, identically for all four arms. Training matrices, fitted ridge coefficients, sequence PCA, targets and folds remained unchanged. Every arm’s unbounded scores reproduced all 668 first-run predictions exactly before bounded results were accepted. No clipping threshold, model hyperparameter or feature set was selected by a search.

This correction is substantial: at least one feature was limited for 76–98% of baseline fold rows and 79–98% of contrast fold rows. It prevents unsafe extrapolation but can suppress real novel signal. Corrected versus uncorrected gains cannot be attributed to the structural hypothesis.

## Matched corrected comparison

All trained arms use the same 668 FGFR2 development variants, ten occupied position blocks, three original folds, continuous target log2(1 + median activation), and ridge model. The historical tree row is an unchanged cached reference.

| Model | Mean-fold Spearman ↑ | Mean-fold MAE ↓ | High activators found /26 | Strong activators found /7 |
|---|---:|---:|---:|---:|
| Matched ridge baseline | 0.23373 | 0.53482 | 2 | 0 |
| Baseline + inactive state | 0.23298 | 0.53251 | 3 | 0 |
| Baseline + active state | 0.22697 | 0.54118 | 2 | 0 |
| Baseline + active − inactive | 0.26784 | 0.53265 | 2 | 1 |
| Historical biology tree model | 0.38350 | 0.49053 | 3 | 0 |

The retrieval budget is 68 total nominations. Spearman is a ranking correlation, not percentage accuracy.

The contrast gains 0.03411 Spearman over the matched corrected baseline and slightly improves mean MAE (0.53482 to 0.53265). It passes both original numerical gates when applied to this sensitivity run. However, the gain is not robust: the paired block-resampling descriptive interval is −0.070 to +0.092, and fold three gets worse (0.18961 to 0.08580). The interval is conditional on reused development predictions, not an independent confidence claim. High-activation retrieval remains 2 of 26; one strong activator is a secondary observation, not a new success criterion.

The stronger historical tree remains better overall: correlation 0.38350 and MAE 0.49053. The active-minus-inactive idea has limited developmental support, but the current implementation is not ready for external validation as our best predictor.

## Verification and timing

- Original ridge baseline reproduced exactly; all four sensitivity models reproduced original unbounded scores exactly.
- 63 focused structural/model/evaluation tests passed; the earlier 449-case regression suite also passed. Scoped Ruff checks passed.
- All 558 first-run protected files and 570 sensitivity-run protected files were unchanged at completion. See the independent verification artifacts for subsequent checks.
- First run: 1.13 seconds fitting/scoring, 3.89 seconds evaluator total. Sensitivity: 1.57 seconds fitting/scoring, 4.54 seconds total. Structural feature construction took 0.30 seconds. Research and implementation time are excluded.
- RLY2023 S5, MANO2021 Data1 and the original final evaluation remain reserved. Research structural files are separate from those assay outcomes.

## Next decision

A defensible next development experiment would test these same frozen contrast descriptors with the stronger existing tree baseline, without changing folds or opening external outcomes. First-run failure and the numerical correction must remain visible. A more ambitious state-population or free-energy approach would require additional mechanistic modeling; these results neither validate nor disprove that broader hypothesis.
