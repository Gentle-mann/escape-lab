# Experiment ledger

This ledger records completed experiments, failed attempts, and pending work through 17 October 2026, including the completed repeated evidence and selector experiments. It includes the parallel GLM-5.3 and You.com research work. The [results summary](RESULTS.md) is the concise reading path.

AP means average precision, a ranking metric. Higher AP and Spearman correlation are better; lower Brier error and mean absolute error (MAE) are better. Pooled macro AP across the two drugs and mean-fold FGFR2 high-activation AP are different endpoints. No table below measures percentage accuracy, patient outcomes, cancer cures, or mutation emergence.

## Data, task, and evaluation boundaries

The primary source is the [FGFR mutational scanning study, DOI 10.1038/s41588-025-02431-8](https://doi.org/10.1038/s41588-025-02431-8), with [open article and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC12807871/). Eight public MAVEdb resistance score sets yield 8,407 nucleotide-variant records normalized to 7,522 protein variants. Drug resistance is measured separately for pemigatinib and futibatinib in engineered cell assays. Activation pilots use the study's four-replicate activation measurements. The endpoints and cell backgrounds differ.

The evaluator keeps representations of the same protein variant and homologous aligned-position blocks together. Missing or ambiguous measurements are not converted to negative labels. Training transforms are fitted within training partitions. Agents receive permitted aggregate development feedback, never outcome tables or authority to change labels, folds, or metrics. This is a logical access boundary around the agent; the trusted evaluator may parse a monolithic archive. It is not a claim of OS-level isolation or secret public literature.

**The original final holdout was consumed by v1. S13 was consumed by v2.** Subsequent GLM and You.com S13 scores are retrospective. Activation, sequence, structural, and AP repair pilots reuse development folds and exclude original final and prospective RLY/MANO blocks. Reserved RLY-4008 S5 and MANO Data 1 numerical outcomes were not evaluated by those pilots. Public-data pretraining exposure cannot be ruled out.

Source run paths below identify preserved local evidence. The publication excludes raw outcomes, private predictions, fitted models, credentials, downloaded third-party archives, and local caches. Source code and protocols describe the experiments; a repeated execution on the same data is reproduction, not new independent evidence.

## Experiment index

| Family / run | Question or change | Outcome |
|---|---|---|
| Initial ingestion, experiment vault, numerical models | Can permitted archival reveals and bounded numerical comparisons run reproducibly? | Infrastructure implemented and tested; not itself evidence of scientific improvement |
| v1 attempt `20261017T193042Z-7992cb52` | First Nemotron configuration proposals | Request failures; no scored agent proposal |
| v1 attempt `20261017T193621Z-d6e64c8c` | Pin model and increase timeout | Completion truncation; no scored agent proposal |
| v1 attempt `20261017T193915Z-ee008058` | Increase completion allowance | Tool-schema incompatibility; no scored agent proposal |
| v1 completed `20261017T194348Z-5e66cbab` | Two autonomous revisions versus fixed/random | Development improved; untouched final AP worsened |
| Revision v2 `20261017T201321Z-a29b17a8` | Biological/structural features, nested selection, stronger guards | Fixed structural external point estimate improved; agent advantage unproven |
| Model swap `20261017T210043Z-a4800896` | Same v2 protocol, GLM-5.3 in place of Nemotron | Internal predictions identical; external point estimate higher but uncertain; six proposal failures |
| `you-pilot-20261017` | Reviewed You.com biological evidence versus empty evidence | Identical selected model and predictions |
| `20261017-experiments-1-2` | Repeated relevant/control/placebo evidence; budget-matched selectors | Completed 90 searches; both experiments inconclusive under predeclared improvement gates |
| `activation-transfer-v1` | Measured versus predicted activation for resistance | Measured activation informative; predicted activation failed advancement |
| `activation-revision-v2` | Capacity, biological features, activation-tail head | Better broad activation correlation, no demonstrated resistance gain |
| `fgfr2-sequence-v1` | Frozen ESM2 plus supervised FGFR2 heads | Higher classifier AP, insufficient retrieval; primary gate failed |
| `fgfr2-state-v1` | Active/inactive/contrast geometry with matched ridge | Failed; severe extrapolation exposed |
| `fgfr2-state-bounded-v1` | Training-range bounding sensitivity | Nominal exploratory gate passed; still below historical tree, uncertain |
| `fgfr2-state-tree-v1` | Same contrast features added to stronger matched tree | Worse all three folds; failed |
| `fgfr2-sequence-ap-v1` | Raw classifier margin instead of saturated probability | Technical ranking screen passed; AP 0.09732 to 0.11135 |

## 1. Original autonomous resistance experiment

Protocol: [revised prediction protocol](archive/docs/revised-prediction-protocol.md). Three fixed families (logistic regression, Extra Trees, histogram gradient boosting) select one development reference. Two random configurations and two Nemotron proposals share the allowed space and folds. A separate validation promotion rule precedes a one-shot final evaluation. Final predictions and source/data/protocol hashes are committed before final scoring.

| Method | Development macro AP | Validation macro AP | Final macro AP |
|---|---:|---:|---:|
| Fixed family selection | 0.05999 | 0.05481 | 0.10923 |
| Two-trial random search | 0.06202 | 0.04805 | 0.10923 |
| Two-proposal Nemotron loop | 0.07051 | 0.05722 | 0.09481 |

Random failed validation promotion and reverted to the fixed model, so its identical final score is not an independent confirmation. Nemotron increased logistic `C` from 1 to 10; widening positional RBF scale from 20 to 40 was subsequently rejected. Final AP dropped by 0.01442, missing the preset +0.05 improvement criterion and the per-drug guard. Pemigatinib AP fell 0.14194 to 0.11783; futibatinib fell 0.07653 to 0.07180. Brier error also worsened (0.08560 to 0.12490 and 0.05943 to 0.06482).

The final cohort contained 1,819 protein variants in seven aligned blocks. Defined labels numbered 1,802/152 positives for pemigatinib and 1,807/108 for futibatinib. Sequence blocks are not independent tumors or laboratories. This result supports a development-overfitting failure, not a universal conclusion about language models.

Operational attempts are retained:

| Run suffix | Model / operational change | Scored proposals | Seconds | Provider-reported tokens |
|---|---|---:|---:|---:|
| `193042Z-7992cb52` | Automatically selected Nemotron-3-Super-120B-A12B; requests failed | 0 | 94.79 | Actual usage unavailable for two failed requests |
| `193621Z-d6e64c8c` | Pinned Nemotron-3.5-Super-VL, 60-second timeout; 2,048-token ceiling truncated replies | 0 | 71.68 | 22,548 |
| `193915Z-ee008058` | 8,192-token allowance; referenced schema returned invalid string-valued config | 0 | 85.97 | 27,095 |
| `194348Z-5e66cbab` | Inlined equivalent strict schema; complete scientific run | 2/2 | 100.74 | 19,277 |

The first failed attempt incidentally evaluated validation using an inherited baseline, so its agent row is not an agent prediction. The next two skipped validation/final. Three synthetic schema probes used fabricated scores, 37.85 seconds, and 14,954 reported tokens; they are engineering diagnostics, not biological experiments. Across four real attempts, runtime was 353.18 seconds and known reported usage was 68,920 tokens, with the two earliest failed requests' actual usage unknown. No unsuccessful run was erased.

Evidence: [v1 complete summary](archive/reports/improvement/results-summary.md), [operational retry record](archive/reports/improvement/inference-retry-note.md), and all four timestamped run directories.

## 2. Biological revision v2 and GLM-5.3 model swap

Revision v2 adds substitution chemistry and drug-contact/packing/motif geometry from [FGFR1 7WCL](https://www.rcsb.org/structure/7WCL) and [FGFR1 6MZW](https://www.rcsb.org/structure/6MZW). FGFR2–4 use aligned FGFR1 proxies, not experimentally observed mutant poses. Training-only preprocessing and selection guards cover AP, probability error, and fold consistency. Twelve fixed configurations define the catalog; agents receive three opportunities per search.

The development cohort excludes original-final and all nine aligned blocks containing the external variants, leaving 3,730 variants in 15 groups. Three outer grouped folds each contain a fresh three-fold inner search; an additional full-data search selects the external model. The independent v2 endpoint was residual biochemical kinase activity from [Saifudeen et al., Supplementary Table S13](https://doi.org/10.1038/s41587-026-03090-8): 25 eligible single variants per drug across nine blocks. This is a different assay from cellular resistance. No external AP threshold was invented; its metric is rank correlation.

| Method | Internal pooled macro AP | Internal macro Brier | S13 mean Spearman |
|---|---:|---:|---:|
| v1 fixed, retrained on common cohort | 0.034413 | 0.035374 | -0.018063 |
| Fixed structural logistic C0.1 | 0.034628 | 0.035781 | 0.405681 |
| Three-trial random selector | 0.034938 | 0.034593 | 0.283740 |
| Twelve-candidate grid | 0.034938 | 0.034593 | 0.182598 |
| Nemotron-3.5-Super-VL selector | 0.030516 | 0.035174 | 0.208915 |
| GLM-5.3 selector, retrospective swap | 0.030516 | 0.035174 | 0.296952 |

All original v2 predeclared comparison intervals included zero. Fixed structure's external gain over retrained v1 was +0.42374 with descriptive interval [-0.08953, +0.69897]. The original internal reference AP 0.10923 is from another split and is not a comparator here. Constant-score internal macro AP was approximately 0.03389, so internal rankings remained weak.

The swap preserves data, folds, seeds, numerical catalog, agent instructions, and opportunity/repair budgets. It changes the model ID to **GLM-5.3**, not GLM-5.2. Both agents selected identical models in all three outer folds, producing identical internal predictions. Full-data selection differed: Nemotron selected structure-logistic-C1; GLM selected interactions-logistic-C0.1. GLM's external mean-correlation increase over Nemotron was +0.088037 with descriptive interval [-0.358545, +0.368979], below fixed structure. This is one stochastic retrospective comparison, not evidence of a general stronger-model benefit.

| Agent | Valid proposals / opportunities | Requests / repairs | Reported tokens | Run seconds | Estimated inference USD |
|---|---:|---:|---:|---:|---:|
| Nemotron | 12/12 | 14/2 | 148,958 | 243.55 | 0.03090 |
| GLM-5.3 | 6/12 | 17/5 | 165,500 | 393.63 | 0.47414 |

GLM had four truncated and two invalid opportunities; each consumed its budget and retained its current model. Both used an 8,192-token completion ceiling. Model-specific reasoning/sampling defaults remained unspecified. Timings include local/concurrent work and are not controlled serving-speed measurements. Cost estimates used uncached provider list prices recorded on the experiment date, not invoices. Bootstrap intervals condition on fixed predictions, omit refitting and LLM-selection variability, and are not biological replication.

Evidence: [v2 summary](archive/reports/revision-v2/results-summary.md), [GLM model-swap summary](archive/reports/model-swap-v2/20261017T210043Z-a4800896/results-summary.md). Protocols: [v2 research](archive/docs/revision-v2-research.md), [model swap](archive/docs/model-swap-v2-protocol.md).

## 3. You.com evidence and repeated selector studies

One paired GLM-5.3 pilot compared empty evidence with two reviewed general biological cards. Four frozen searches produced 20 results; one Contents request selected four pages and returned three. Only supported, dated cards entered the agent context. The accepted cards concerned solvent environment and limitations of static structures for kinetics. Raw retrieved text, excerpts, test outcomes, and evaluator artifacts were excluded from the proposal prompt.

| Arm | Selected model | Inner selection AP | S13 mean Spearman | Valid proposals | Seconds |
|---|---|---:|---:|---:|---:|
| No evidence | interactions-logistic-C0.1 | 0.04217 | 0.29695 | 2/3 | 86.26 |
| Reviewed You.com evidence | interactions-logistic-C0.1 | 0.04217 | 0.29695 | 2/3 | 80.12 |

Both selected predictions were identical; correlation difference and its conditional bootstrap interval were exactly zero. This does not establish that literature research never helps. Each third opportunity truncated at the token limit. Total pilot runtime was 179.8 seconds; estimated retrieval cost $0.024 plus inference $0.19627, total $0.22027. The already-spent S13 panel was reused; no fresh external holdout was consumed. All four numerical controls reproduced original selections and predictions to 1e-12.

Evidence: [You.com paired summary](archive/reports/research-pilot/you-pilot-20261017/results-summary.md), [execution review](archive/reports/research-pilot/you-pilot-20261017-review.md). Protocol: [You.com pilot](archive/docs/you-research-pilot-protocol.md).

### Completed experiments 1 and 2

The follow-on [repeated research protocol](archive/docs/repeated-research-protocol.md) completed all **90 isolated searches**: ten planner repetitions across three grouped outer folds and three evidence arms (empty, relevant, truthful unrelated placebo). The cohort contains 3,730 variants in 15 aligned groups. The primary score averages ten repeat-level pooled out-of-fold macro AP values, each averaging the two drugs.

| Method | Mean macro AP | Mean macro Brier |
|---|---:|---:|
| Agent, no evidence | 0.02897045 | 0.03503282 |
| Agent + You.com evidence | 0.02914259 | 0.03504838 |
| Agent + unrelated evidence | 0.02914202 | 0.03504858 |
| Fixed structural model | 0.03462765 | 0.03578062 |
| Random selection, three tries | 0.03226587 | 0.03461819 |
| Numerical selection, three tries | 0.02880030 | 0.03501657 |
| Full grid, eleven tries | 0.03493831 | 0.03459281 |

**Experiment 1, evidence benefit: inconclusive.** Relevant evidence increased macro AP by only +0.000172143 over no evidence and +0.000000571 over placebo. Both gains missed the predeclared +0.002 minimum, and both adjusted 99% intervals included zero. The two reviewed You.com cards did not demonstrate a useful prediction benefit.

**Experiment 2, agent selection: inconclusive.** The predeclared primary agent was the no-evidence arm. It performed below fixed and random selection in all ten repeats, with mean differences of -0.005657204 and -0.003295420. The agent-minus-random adjusted 99% interval was entirely negative, [-0.010566, -0.000564]; the agent-minus-fixed interval included zero. Its +0.000170146 mean gain over the specified numerical policy also missed the improvement gate. These results do not establish an advantage over simpler selection.

All five primary contrasts failed both the minimum-gain and positive-lower-bound requirements; their per-drug AP and Brier safeguards passed. Intervals use 10,000 paired bootstrap draws over aligned groups and repeat IDs and condition on fitted predictions. Agent, random and numerical selection each had three opportunities beyond a common reference. Fixed selection had no adaptive search cost; full grid is a larger-budget reference. The numerical arm tests one specified policy, not every numerical optimizer.

All 90 searches entered scoring, including failures. **270 requests produced 264 valid proposals and six failed opportunities**: five truncated responses and one connection error. No scientific result triggered a rerun. Models and unlabeled predictions were committed together before evaluator scoring. The independent request audit found no evaluation outcomes, cross-arm history or raw source text in proposal inputs; agents had no filesystem or web tools. The trusted evaluator can parse a monolithic local outcome archive, so this is an input-isolation and delayed-scoring guarantee, not OS-level isolation.

Runtime was 21.58 minutes, with 3,011,005 reported tokens and an estimated $8.75993 in known-usage inference plus $0.024 for new retrieval. One failed request has unknown provider usage; this is not a complete invoice. Cached numerical fits make wall time unsuitable as a standalone optimizer-cost comparison. Incomplete placebo-source acquisition attempts remain recorded locally.

These are reused development data. Ten repetitions measure procedural variability on the same biological folds, not ten independent biological validations. No original final, S13, RLY or MANO outcome was scored. The study tests two fixed evidence cards within a twelve-configuration catalog; it does not establish fresh-study generalization or that literature research can never help.

Evidence: [completed review](archive/reports/repeated-research/20261017-experiments-1-2-review.md), [aggregate results](archive/reports/repeated-research/20261017-experiments-1-2/results.json), and [verification and choice summary](archive/reports/repeated-research/20261017-experiments-1-2-audit.json). Raw outcomes, predictions and request logs remain local.

## 4. Activation as an intermediate resistance feature

`activation-transfer-v1` used 2,719 complete variants from 2,749 eligible examples, across 11 aligned blocks and three fixed folds. The same complete cohort was used for all comparisons. Predicted activation was generated with cross-fitting; measured activation assumes the activation assay is already available for the held-out variant and is therefore a different deployment task.

| Resistance input | Macro AP | Interpretation |
|---|---:|---|
| Existing structural features | 0.06016 | Matched reference |
| Structure + predicted activation | 0.05961 | Failed advancement |
| Structure + measured activation | 0.11718 | Assay-assisted signal; not sequence-only prediction |
| Structure + predicted shuffled activation | 0.05763 | Single shuffle capacity control |

Predicted activation's AP difference interval was [-0.00513, +0.00320]; measured activation's was [+0.01308, +0.13811]. These are descriptive resamples of eleven related blocks. Predicted versus measured activation correlated 0.376 overall but only 0.195 in FGFR2. The measured arm's Brier error improved 0.04949 to 0.04164, yet a post-hoc training-prevalence check scored 0.03083, so accurate probabilities were not established. Runtime: 9.60 seconds fitting, 14.26 seconds total. No LLM calls.

`activation-revision-v2` used the identical cohort. It separated extra tree capacity from added flanking sequence, smooth position, motif proximity, mutation interactions, and AlphaMissense/EVE features. The evolutionary/pathogenicity scores are auxiliary predictions, not drug-resistance or activation labels.

| Arm | Macro within-gene activation rho | Activation MAE | Activation-class macro AP | Resistance macro AP |
|---|---:|---:|---:|---:|
| Original | 0.33902 | 0.24706 | 0.04632 | 0.05961 |
| More capacity, original features | 0.36730 | 0.26181 | 0.04817 | 0.05817 |
| Capacity + biology | 0.38101 | 0.26329 | 0.05527 | 0.06012 |
| Biology + activation classifier | 0.38101 | 0.26329 | 0.05428 | 0.06067 |

FGFR2 correlation improved 0.19510 to 0.34645, but the tail head did not reliably retrieve activation extremes. The broad activation definition had 122 positives: FGFR1 5/690, FGFR2 96/668, FGFR3 2/686, FGFR4 19/675. One FGFR3 fold lacked positives, making the component gate inconclusive. The resistance gate failed: gain over structural AP was only 0.00052, with one fold losing 0.02729 and interval [-0.00553, +0.00684]. Runtime: 69.08 seconds fitting, 75.00 total; no LLM calls. These broader activation labels differ from the later 26-variant high-activation target.

Evidence: [activation-transfer report](archive/reports/mechanistic-pilot/activation-transfer-v1/report.md), [activation-revision report](archive/reports/mechanistic-pilot/activation-revision-v2/report.md).

## 5. FGFR2 sequence and high-activation ranking

The remaining pilots share 668 FGFR2 examples, ten occupied aligned blocks, and three fixed folds of 275/180/213 rows. “High” means at least three of four replicate summaries >=1.5 and median >5; there are 26 positives (7/7/12 by fold). “Strong” means median >10 with the same consistency rule; seven positives, including a zero-positive fold. Ranking budgets are 28/18/22, totaling 68. Mean-fold AP is primary; pooled scoring is a separate diagnostic.

Frozen ESM2-35M supplies wild-type residue and mutant-minus-wild-type vectors. Separate training-only PCA16 branches combine with biological features. Logistic uses C0.1 with balanced classes; ridge uses alpha100. No protein-model fine-tuning or live LLM is involved.

| Arm | Mean-fold high AP | High hits /26 | Strong hits /7 |
|---|---:|---:|---:|
| Cached pooled biology HGB | 0.071864 | 3 | 0 |
| FGFR2 biology logistic | 0.047092 | 0 | 0 |
| ESM + biology logistic, probability | 0.097321 | 1 | 0 |
| ESM + biology ridge | 0.117085 | 4 | 0 |

The sequence classifier's original primary gate failed: fewer high hits than the historical tree, AP improvement under the preset 0.03 threshold, and improvement in only one fold. The matched biology-only classifier isolates sequence representation; the historical tree also differs in training cohort and objective. ESM extraction took 49.62 seconds; head fitting/scoring 0.64 seconds. Later diagnostics revealed severe continuous ridge extrapolation, so its AP must not be advertised as evidence of accurate activity estimates.

Evidence: [FGFR2 sequence report](archive/reports/mechanistic-pilot/fgfr2-sequence-v1/report.md).

## 6. Active versus inactive structure experiments

Prepared active coordinates derive from FGFR2 2PVF; inactive coordinates are a published FGFR1-3KY2-templated FGFR2 homology model. Both retain engineered C491A and lack activation-loop phosphotyrosine. The common mask excludes position 491 and uncertain loops, leaving paired geometry for 528/668 variants; the remainder retain explicit missingness. Thirty-one contact, packing, motif-distance and substitution-interaction descriptors represent each structural arm. These are static proxies, not free energies, state populations, relaxed mutant structures, or molecular-dynamics predictions.

Background sources: [eLife molecular-brake study](https://elifesciences.org/articles/21137), [PNAS structural-state study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9974468/), and [authors' prepared structures](https://yzhang.hpc.nyu.edu/MD_data_FGFR/).

The original matched ridge experiment used continuous log2(1 + median activation):

| Unbounded ridge arm | Mean-fold rho | Mean-fold MAE | High AP | High hits /26 |
|---|---:|---:|---:|---:|
| Sequence/biology baseline | 0.19585 | 3933.44035 | 0.117085 | 4 |
| + Inactive geometry | 0.20696 | 1102.51852 | 0.122270 | 5 |
| + Active geometry | 0.19078 | 1777.89393 | 0.119663 | 5 |
| + Active-minus-inactive geometry | 0.22328 | 21494.56594 | 0.083132 | 5 |
| Historical biology HGB | 0.38350 | 0.49053 | 0.071864 | 3 |

Both structural gates failed. A near-constant training positional feature became enormous after standardization on a held-out region. Thus attractive ranking numbers coexisted with unusable numerical predictions. Baseline ridge reproduction was exact. Runtime: 1.13 seconds fitting, 3.89 total.

A declared **post-result sensitivity** limited standardized biology features to their training ranges at inference, retaining the same fitted coefficients, PCA, and training matrices:

| Bounded ridge arm | Mean-fold rho | Mean-fold MAE | High AP | High hits /26 | Strong hits /7 |
|---|---:|---:|---:|---:|---:|
| Baseline | 0.23373 | 0.53482 | 0.067906 | 2 | 0 |
| + Inactive | 0.23298 | 0.53251 | 0.067103 | 3 | 0 |
| + Active | 0.22697 | 0.54118 | 0.063924 | 2 | 0 |
| + Contrast | 0.26784 | 0.53265 | 0.071182 | 2 | 1 |

Nominal incremental/state-contrast gates passed, but contrast's rho gain interval was [-0.07004, +0.09244], the third fold worsened, and the historical tree remained stronger on correlation and error. This is sensitivity analysis, not fresh confirmation. Runtime: 1.57 seconds fitting, 4.54 total.

The final tiny structural test retrained the historical HGB with exactly the same pooled FGFR1–4 cohort and hyperparameters in both arms, adding the cached contrast features only to the challenger. Six total fits, no tuning:

| Matched tree arm | Mean-fold rho | Mean-fold MAE | High AP | High hits /26 |
|---|---:|---:|---:|---:|
| Biology HGB | 0.38350 | 0.49053 | 0.071864 | 3 |
| Biology HGB + contrast | 0.32900 | 0.54613 | 0.054112 | 1 |

Every fold's correlation and MAE worsened; all four screening checks failed. Baseline reproduction was exact. Runtime: 5.28 seconds fitting, 6.67 total. This result does not support promoting these static contrast features.

Evidence: [unbounded structure report](archive/reports/mechanistic-pilot/fgfr2-state-v1/report.md), [bounded sensitivity report](archive/reports/mechanistic-pilot/fgfr2-state-bounded-v1/report.md), [matched tree report](archive/reports/mechanistic-pilot/fgfr2-state-tree-v1/report.md).

## 7. Targeted AP numerical repair

The most recent change keeps the original fitted ESM/biology logistic classifier and uses raw decision margins for ranking instead of sigmoid probabilities. In pilot-0, 136 probabilities rounded to exactly 1.0; margins retained 275 unique scores. The other folds had no exact-one saturation.

| Fold | Probability AP | Margin AP | High hits before / after |
|---|---:|---:|---:|
| pilot-0 | 0.024633 | 0.066708 | 0 / 3 |
| pilot-1 | 0.215701 | 0.215701 | 1 / 1 |
| pilot-2 | 0.051628 | 0.051628 | 0 / 0 |
| **Mean AP / total hits** | **0.097321** | **0.111346** | **1 / 4** |

The technical screen passed: gain >=0.01, no fold AP loss, and no retrieval loss. Both rankings retrieved zero of seven strong variants. AP gain was 0.014025 with descriptive paired-block percentile interval approximately [0, 0.032615]. Three fits took 0.31 seconds; total 2.92 seconds. The original probabilities reproduced exactly and 594 protected files stayed unchanged.

This recovers ranking information; it does not add biological information, improve probability calibration, fix unsafe extrapolation, or establish generalization. It was implemented during development, not discovered autonomously by the runtime LLM loop. It is not the highest historical AP across all arms, and the raw margin's largest value still exceeded one million in the first fold. A robust next model must address that extrapolation without choosing changes against the reserved external outcomes.

Evidence: [AP repair report](archive/reports/mechanistic-pilot/fgfr2-sequence-ap-v1/report.md). Implementation: [fgfr2_sequence_ap.py](../src/escape_lab/fgfr2_sequence_ap.py).

## Research audits that were not predictive experiments

- External-prediction inventory reconciled 7,522 protein names and paper-provided AlphaMissense/EVE scores. These are general pathogenicity scores, not historical drug-specific predictions. CIViC records and paper-curated resistance identities are previously observed evidence, not blind forecasts. No CAGI participant prediction files were obtained, and no superiority over them is claimed.
- The mechanistic [data-to-model map](archive/docs/mechanistic-data-to-model-map.md) distinguishes measured activation, binding proxies, missing kinetic/cellular parameters, and reserved external assays. It is a feasibility map, not a whole-cell simulator.
- Structural audits rejected treating 2PSQ as a clean inactive comparator, documented prepared-template limitations, and fixed masks before structural scoring.
- v3 development-readiness, extraction/runtime probes, schema smoke tests, and artifact-preservation checks are preparation or verification, not additional scientific wins.

The completed repeated-study results are included above. Any new external test requires a separately frozen protocol and should remain visibly separate from these reused development results.
