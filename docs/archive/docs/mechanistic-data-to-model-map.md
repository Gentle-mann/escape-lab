> Historical research record. See [current results](../../RESULTS.md) for the release interpretation.

# Escape Lab: mechanistic data-to-model map

Research snapshot: 17 October 2026. This is a proposed next experiment, separate from the exact GLM-versus-Nemotron v2 rerun. No model was trained or evaluated for this audit. Original runtime files, outcomes, receipts and frozen experiments were not changed.

## Decision

We can build a mechanism-informed predictor now. We cannot yet identify or validate a complete molecular-to-cell growth simulator from the available measurements. The important missing evidence is a matched chain of baseline function, target inhibition, downstream signaling and growth in the same variant, cell background, concentration and exposure schedule.

Use the unused activation measurements as an auxiliary phenotype, verified molecular structures as drug-specific features, and separate models for each measured endpoint. Test cellular target inhibition and cell growth separately. A favorable result would establish predictive value for the added biology, not that we have recovered the entire causal mechanism.

## 1. What data we have

| Layer | Verified resource | Appropriate use | Status / limitation |
|---|---|---|---|
| Mutation identity | Original paired manifest: 8,407 nucleotide records, 7,522 protein variants | Canonical sequence/variant identity and grouped splits | Multiple nucleotide records are not independent protein examples. |
| Baseline functional phenotype | Four activation score sets, four replicate columns, aggregate and SDV | Predict activation-assay enrichment; optionally supply measured activation in a separately labeled assay-assisted task | Not ingested by the current resistance pipeline. MCF10A cellular selection, not purified catalytic activity. |
| Drug chemistry | PubChem pemigatinib CID 86705695; futibatinib CID 71621331 | Free-molecule identity, stereochemistry and descriptors | Structure is an input, not measured binding potency. Two drugs do not support learning broad chemical generalization. |
| Protein–drug geometry | Existing PDB 7WCL and 6MZW; available FGFR2–futibatinib structure 8W2X | Contact geometry, mutant-versus-WT hypotheses, covalent-site features | Existing templates are FGFR1 proxies; both have engineered C584S. 8W2X is a potential FGFR2-specific improvement, not yet integrated. |
| Biochemical inhibition | Already evaluated S13: 25 eligible mutants × pemigatinib/futibatinib | Development-only auxiliary endpoint or diagnostic, excluding any newly held-out positions | Spent. Residual enzyme activity is not cellular survival or an automatically measured IC50. |
| Cellular response | Existing three-replicate drug-enrichment endpoints | Learn continuous assay response and/or the declared binary target | NCI-H1581; one concentration per drug. Original final test is spent. |
| Small cross-context follow-ups | Original S4: seven FGFR2 mutants + WT + parental, serum-response ED50. S5: seven mutants + WT across six inhibitors | Development calibration and endpoint checks | Same-study selected validations; not an independent final test. Numeric cells were not inspected in this audit. |
| Existing pretrained scores | Local AlphaMissense 7,522; EVE 5,554; consistent REVEL 7,274 | Ordinary auxiliary-feature baselines with missingness flags | General variant scores, not drug-resistance truth. |

Local count/status sources: `data/processed/manifest.json`, `src/escape_lab/data.py`, `research/external-predictions/2026-10-17-aggregate-sources/README.md`, and `reports/revision-v2/results-summary.md`. The creation manifest still says the final set is sealed; the completed receipts are authoritative. `data/processed/final-evaluation.json` and `data/revision-v2/external-receipt.json` both report completed evaluations.

Activation metadata verified through the live MAVEdb API:

| Gene | Score-set accession suffix under `urn:mavedb:00001256-a-` | Nucleotide records |
|---|---|---:|
| FGFR1 | 12 | 2,103 |
| FGFR2 | 1 | 2,122 |
| FGFR3 | 3 | 2,102 |
| FGFR4 | 11 | 2,080 |
| Total | | 8,407 |

Identifier-only reconciliation of original S10 gives 7,522 protein identities. These are additional measurements of existing variants, not 7,522 new independent final-test examples. The working experiment endpoint is [MAVEdb 00001256-a](https://api.mavedb.org/api/v1/experiments/urn:mavedb:00001256-a); the older accession in the paper did not resolve during this audit.

The [original study](https://www.nature.com/articles/s41588-025-02431-8) uses different cells for activation and resistance. Its detailed drug protocol specifies 20 nM pemigatinib or 10 nM futibatinib, 72 hours followed by another 96 hours after refreshment. Preserve these explicit intervals rather than interpreting an inclusive day count as an exact duration. S10 contains endpoint enrichment, not a matched untreated NCI-H1581 growth trajectory. Raw sequencing is available under [PRJNA1136966](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1136966), but reprocessing the same experiment does not create an independent test.

## 2. Model components and assumptions

### Activation branch

Input: sequence, mutation and permitted structural descriptors. Output: predicted activation-assay phenotype, with uncertainty.

Do not rename this output catalytic activity. Enrichment can combine receptor expression, signaling and proliferation. Transfer between cell contexts is a hypothesis to test.

For sequence-only prediction, train this branch within the training partition and generate grouped out-of-fold predictions for downstream training. Keep all nucleotide representations, aligned positions and auxiliary labels together. Feeding the true activation measurement for a held-out variant changes the task to assay-assisted prediction; report that experiment separately.

### Inhibition branch

Use checked free-drug chemistry, gene-specific structure where available, mutant-versus-WT geometry and explicit reversible/covalent features. Sources: [pemigatinib](https://pubchem.ncbi.nlm.nih.gov/compound/86705695), [futibatinib](https://pubchem.ncbi.nlm.nih.gov/compound/71621331), [7WCL](https://www.rcsb.org/structure/7WCL), [6MZW](https://www.rcsb.org/structure/6MZW), [8W2X](https://www.rcsb.org/structure/8W2X).

Audit sequence numbering, construct differences, missing coordinates and ligand identity before feature generation. In a covalent structure, the bound chemical component may represent an adduct; do not mistake its chemical representation for the unreacted drug. Docking or structural relaxation can generate features and hypotheses, but needs experimental calibration before it becomes a potency prediction.

For reversible competition, a simplified enzyme relationship is:

`IC50 = Ki × (1 + ATP / Km_ATP)`.

This is conditional on kinetic assumptions, not a universal cellular formula. For a two-step covalent inhibitor, a common simplified relationship is:

`k_obs(D) = k_inact × D / (K_I + D)`.

One endpoint at one dose does not separately determine binding affinity, ATP affinity, reaction rate, Hill slope, receptor turnover and intracellular drug exposure. Relevant primary references: [Cheng–Prusoff](https://pubmed.ncbi.nlm.nih.gov/4202581/) and [time-dependent inhibition analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC11089659/).

[Saifudeen / KiRHub](https://www.nature.com/articles/s41587-026-03090-8) provides a large biochemical resource, but assay products must remain distinct: the mutant residual-activity screen and WT dose-response measurements are not the same evidence. The full panel uses a radiometric assay with controlled ATP conditions; single-concentration mutant values cannot be relabeled as measured dose-response constants. Confirm each table's ATP, exposure and concentration metadata before combining studies.

The [RLY study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10481131/) supplies limited kinetic anchors through time-dependent intact-mass measurements. These do not constitute a full mutant-specific kinetic panel. If WT parameters are used to calibrate an inhibitor, call the result WT-calibrated transfer to withheld mutants, not completely zero-shot prediction. Do not confuse drug-metabolism enzyme kinetics with FGFR kinetics.

### Cellular-response branch

Use a small regularized model to relate predicted activation, inhibition-related features and explicit assay context to the actual observed endpoint. Initially predict assay-specific response rankings or enrichment. Do not output absolute cell counts, division/death rates or calibrated cellular IC50 unless those targets are independently supported by training data.

If concentration and time are constant in training, their effects cannot be learned by merely adding them as input fields. A mechanistic dose/time model would need experimental variation to constrain those parameters.

The relevant missing measurements include receptor abundance, untreated phosphorylation, pathway dependence, intracellular exposure, initial cell counts and untreated growth rate. DepMap offers [genomic and expression resources](https://depmap.org/portal/data_page/), but measurements from a reference cell line are not automatically matched to the engineered cells in our assays. Such annotations are context priors, not substitutes for missing paired measurements.

Growth-rate normalization requires the appropriate measurements. [Hafner et al.](https://pubmed.ncbi.nlm.nih.gov/27135972/) explains why conventional response metrics can change with baseline division rate. Our pipeline must not invent that normalization from endpoint enrichment alone.

## 3. External test candidates

### A. Ready intermediate challenge: reserved RLY study S5

- Existing metadata identifies 11 eligible FGFR2 single mutants in six aligned blocks.
- The header explicitly includes futibatinib and pemigatinib, alongside RLY-4008, infigratinib and erdafitinib.
- Endpoint: fold change in IC50 relative to FGFR2 WT, with SD, for cellular pFGFR2 inhibition after two hours in engineered 293T cells.
- Primary question: does the model rank mutation-dependent loss of target inhibition for our existing drugs?
- It does not test cell survival or the complete growth pathway. New-drug RLY-4008 extrapolation should be a separate, optional experiment.

[S5 asset](https://aacr.figshare.com/articles/journal_contribution/24091984), [direct DOCX](https://ndownloader.figshare.com/files/42276655), [paper](https://doi.org/10.1158/2159-8290.CD-23-0475). Header-only inspection was performed; numeric outcome cells were not extracted. The local metadata-only exclusion audit leaves 4,546 original development/validation variants if these six blocks alone are excluded. A final joint holdout design may leave fewer.

### B. Broader cellular candidate: Nakamura MANO study

[Nakamura et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8285406/) studies FGFR variants in NIH-3T3 cells and includes both target drugs. Its pooled drug assay uses five-day exposure, a concentration range and triplicates. The headline describes 160 mutants and ten fusions, but that is not our eligible kinase-domain single-mutant test count.

The verified quantitative source is **Supplementary Data 1**, `41698_2021_204_MOESM3_ESM.xlsx`, supplemented in the paper by Figure 4. **Supplementary Data 2 is a pairwise-comparison table**, not the response matrix. Raw barcode sequencing is author-request-only. The archive was successfully retrieved through [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8285406/supplementaryFiles). The header-only audit confirmed a numeric IC50 table, with futibatinib in column G and pemigatinib in H.

Data 1 has 70 distinct single-mutant labels. Fourteen match our current WT-position/variant identifiers: five FGFR1 and nine FGFR2, spanning 12 existing aligned blocks. These are provisional matches, not a completed sequence/construct or finite-outcome audit. The other 56 labels must not be assumed eligible or automatically classified as outside-domain without checking. Figure 4 contains common-variant results separately; it was not digitized. Use a Data-1-only test scope initially rather than promising the headline study size.

This is an independently performed cellular experiment, but transformation-based selection and pooled competition limit transfer. It cannot be assumed to represent all mutations or all cell backgrounds. There are at most 28 matched mutant–drug entries before checking missingness/censoring. For this test use within-drug IC50 ranking; WT-normalized potency shifts require verified matched WT controls. Include a gene-only control, since gene-level potency differences can otherwise masquerade as mutation-effect prediction.

### C. Smaller cellular replication panels

| Study | Verified design | Role / unresolved issue |
|---|---|---|
| [Facchinetti 2024](https://doi.org/10.1158/1078-0432.CCR-24-1834) | 17 FGFR2 mutants, Ba/F3 FGFR2::BICC1, 48-hour viability, both drugs | Controlled fusion-context replication. The public S5 asset is a figure/PPTX; numerical recoverability needs a format audit. |
| [LC-F2-1 2026](https://doi.org/10.1038/s41467-026-76339-0) | 12 FGFR2 mutants, Ba/F3 TEL-FGFR2, 72-hour proliferation, nine-point threefold dilution, duplicates, both drugs | Independent-lab replication with a named Source Data workbook, `41467_2026_76339_MOESM4_ESM.xlsx`. Confirm eligible rows and exact data availability before freezing. |
| [Goyal 2025](https://doi.org/10.1016/j.annonc.2024.12.011) | Biochemical, structural and cellular phosphorylation studies | Valuable mechanistic calibration candidate; supplement eligibility and potential reuse across studies still need an audit. |

Choose primary and replication roles based on metadata and accessibility before viewing model performance. Do not promote whichever study gives the best result. Study independence also requires checking shared constructs, copied measurements and overlapping investigators; a different DOI alone is insufficient.

## 4. What “untouched” can honestly mean here

The original 1,819-variant final split and the v2 S13 panel have already been evaluated. They can support development diagnostics, not a fresh final claim.

Unused activation endpoints are new training information about existing variants. A new external assay can supply independent measurements even when its mutations overlap old experiments. For a claim about unseen positions, exclude those positions from every training modality, including activation and already-spent biochemical data.

The reserved RLY numeric table and newly identified cellular tables have not been used to select this proposed model. However, published abstracts and search snippets exposed some qualitative associations and isolated values. Therefore call a future run a preregistered external numerical evaluation with declared literature exposure, not perfectly literature-blind discovery. Pretraining contamination cannot be excluded.

Maintain an exposure ledger with source, asset hash, metadata read, qualitative findings exposed, numerical fields accessed, and whether any outcomes affected development. Mechanistic literature may inform the design; the final proposing agent must not retrieve the test papers or view test labels.

## 5. Evaluation that distinguishes better data from a better mechanism

First freeze an endpoint contract. Cellular enrichment, biochemical residual activity, phosphorylation IC50 and viability IC50 are separate measurements. Do not concatenate them into one resistance probability target.

Compare these models under the same permitted data, grouped splits and tuning budget:

| Model | Purpose |
|---|---|
| Existing v2 numerical predictor | Historical reference |
| Drug/gene-only reference | Check whether broad gene potency differences explain the result |
| Regularized empirical model with the proposed new features | Fair test of whether extra data alone explains improvement |
| Activation-only addition | Contribution of the functional phenotype |
| Structure/drug-feature addition | Contribution of molecular features |
| Combined mechanism-informed model | Value of the complete representation |
| Same inputs with versus without explicit mechanistic constraints | Whether the constraints themselves help |

A hybrid beating only the old model does not establish that simulation helped. A flexible model with the same inputs may be stronger. Keep the simpler result if the mechanistic assumptions reduce predictive performance.

Suggested endpoints, to lock after the metadata audit but before outcomes:

- Activation: rank correlation and error on the declared transformed enrichment scale; optional GoF classification assessed separately.
- RLY S5: per-drug Spearman correlation with reported potency shift; pairwise comparison against the strongest preregistered empirical control.
- Cellular panel: per-drug response ranking; use potency-shift error only if the model was trained to predict that quantity. Treat greater-than/less-than values as censored observations, not exact numbers.
- Report coverage, omitted rows, confidence intervals and drug-specific degradation. Bootstrap independent mutation/position groups, not duplicate rows. Small selected panels remain weak evidence even with many assay replicates.

Freeze all model selections and predictions for both external stages before opening either outcome table. Never use the intermediate test result to tune the growth prediction and still call the growth test a single untouched evaluation of the original design.

The metadata-only audit finds 3,328 original development/validation variants remaining after excluding Data 1's mapped blocks, or **2,749 after jointly excluding those blocks and the reserved RLY blocks**. This retains the existing 12-residue aligned-block scheme. These counts concern the original training pool only; final sequence mapping or additional development sources can change the eligible pool. Apply the same exclusions to auxiliary activation, S4/S5 and biochemical training observations. Do not silently relax grouping or recover excluded examples from another modality.

Success means a reproducible gain over the matched empirical control on independent measurements, with uncertainty supporting the claim and a separately reported cellular result. Improved target-inhibition ranking alone supports only that intermediate task. Positive point estimates with wide intervals are inconclusive. Better development scores without external improvement are a predictive failure, even if the agent loop executes correctly.

## 6. Build sequence and data contract

1. Finalize an outcome-blind inventory of the candidate external assets and eligible variants. Fix test scope before importing auxiliary training labels.
2. Ingest activation replicates into a separate trusted table; preserve original nucleotide-to-protein relationships, units, replicate identities and uncertainty.
3. Add checked free-drug structures and gene-specific geometry; record every mismatch and missing feature. Avoid expensive molecular dynamics until a specific missing observable and a way to validate it are defined.
4. Fit activation and response models with grouped cross-fitting. Add only low-dimensional inhibition parameters that the development data can constrain.
5. Fit the matched empirical controls and component ablations. Freeze preprocessing, models, predictions, scripts and source hashes.
6. Score the reserved target-inhibition and cellular-response assays once, then report all prespecified results.

Every observation should retain: study/asset ID and hash; canonical sequence and variant; cell line/species/construct/fusion; drug identity and chemical form; measured endpoint and unit; concentration, ATP and exposure metadata where available; replicate and censoring flags; measured/predicted/assumed origin; and train/development/test/exposure status. Missing information stays missing. Biological parameters supplied by assumption must have ranges and a sensitivity analysis, not fabricated precision.

For eventual validation of a true mechanistic simulator, obtain matched variant-specific baseline receptor/signaling measurements, dose/time inhibition curves, downstream signaling, untreated and treated cell counts, and preferably intracellular target engagement or exposure. That is the experiment needed to connect the full chain.

Redistribution requires the license of each exact asset. Local MAVE-derived measurements are marked CC BY-SA 4.0. MANO is CC BY 4.0. AACR article and Figshare supplement licenses may differ; the RLY S5 asset metadata reports CC BY 4.0. Do not infer a supplement's license from the journal alone.

## Final metadata audit note

Completed and saved to `research/mechanistic-map/mano-metadata-audit.json`. The workbook SHA-256 is `db674c47d458f7eec5aff88a8c0f3def33fca3b1a62d91803fffbed29fdcef93`. The audit retained headers, identifier matches, cell-type counts and exclusion counts. It did not extract numerical response values, inspect predictive performance, or establish response completeness/censoring. Downloading and mechanically inventorying the workbook is not a claim that it was never accessed.

Recommended next test pairing: RLY S5 for cellular target inhibition, and the eligible MANO Data 1 subset for cellular response, with both predictions committed before either score is revealed. The smaller independent panels above remain possible prespecified replications. We have not automatically reserved or spent any additional evaluation receipt in this research task.

After writing this report, all 28 source/data artifacts bound by the original v2 protocol were checked and remained hash-identical. This audit did not change the active GLM comparison.
