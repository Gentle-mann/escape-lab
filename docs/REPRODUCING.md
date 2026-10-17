# Reproduction guide

## What a clean clone can do

`uv sync --frozen --python 3.12`, `uv run --frozen pytest -q`, and `uv run --frozen ruff check src tests scripts` validate the implementation using synthetic fixtures. The local original workspace had 583 passing tests at publication preparation. A clean clone can skip tests that explicitly need optional local reference caches. No credentials or scientific evaluation are needed for these checks.

The repository includes aggregate results and source code, not the private artifact chain required to rerun every historical evaluation. Do not represent the public clone as a one-command reproduction of the headline AP. Historical runners intentionally require prior protocol hashes, frozen cohorts, features and receipts.

## Original FGFR ingestion and replay

```sh
uv run --frozen escape-lab prepare
uv run --frozen escape-lab status
uv run --frozen escape-lab baseline --stage development --seeds 0 1 2
uv run --frozen escape-lab serve --port 8017
```

`prepare` downloads eight public MAVEdb score sets and source metadata, normalizing 8,407 paired nucleotide variants into 7,522 protein variants. The data manifest records exact URLs, licenses and checksums. No row-level outcomes are shipped in Git. The source study is https://doi.org/10.1038/s41588-025-02431-8 and an example accession is `urn:mavedb:00001256-a-9`.

## Agent experiments

Copy `.env.example` to `.env` and configure keys locally. Crusoe model access and price can change; discover available Nemotron IDs using `uv run --frozen escape-lab crusoe-models`. The separate `model_swap` adapter discovers GLM across the full provider catalog. The historical IDs were `nvidia/NVIDIA-Nemotron-3.5-Super-VL` and `zai-org/GLM-5.3`.

The source entry points are `escape_lab.improvement_agent`, `escape_lab.revision_cli`, `escape_lab.model_swap`, `escape_lab.research_pilot`, and `escape_lab.research_repeat_run`. Their protocols are preserved in `docs/archive/docs/`. They depend on frozen local feature and evaluation artifacts. Do not delete receipts, substitute different inputs silently, or rerun a final-scoring command to make a demo work. Charged inference and ordinary deterministic code tests are separate operations.

## Activation, sequence and structure pilots

The activation target comes from Supplementary Table S10 of the FGFR study, using codon-collapsed replicate measurements and the cohort exclusions documented in the protocols. Its workbook is stored locally at `research/external-predictions/2026-10-17-aggregate-sources/raw/fgfr-supplementary-tables.xlsx`; it is not distributed here. Obtain the original supplement from the study, check its recorded checksum in the original experiment record, and preserve the same eligibility and grouping logic before any new run. External S13, RLY and MANO panels are different measurements and must not be substituted for S10.

The sequence experiment requires exact FGFR2 references, eligible candidate metadata, earlier cohort/fold artifacts, and a frozen ESM cache. The pinned checkpoint is `facebook/esm2_t12_35M_UR50D`, revision `6fbf070e65b0b7291e7bbcd451118c216cff79d8`. The original optional environment was:

```sh
HF_HUB_DISABLE_XET=1 TOKENIZERS_PARALLELISM=false uv run --frozen \
  --with 'torch==2.8.0' --with 'transformers==4.57.1' \
  python -m escape_lab.fgfr2_sequence_features --help
```

Torch and Transformers are optional overlay dependencies; their transitive dependencies include Hugging Face Hub and Safetensors. They are not part of the core `uv.lock`. The historical run used Apple MPS, batch size 4; CPU inference is supported but timing will differ. No fine-tuning of ESM occurred. See [original sequence run instructions](archive/research/fgfr2-sequence-pilot/run-instructions.md).

The structure comparison uses prepared active/inactive models, a common residue mask, and engineered-background caveats. It does not calculate molecular dynamics or free energies. See [structure audit](archive/research/fgfr2-state-pilot/structure-audit.md).

The latest AP entry point is `python -m escape_lab.fgfr2_sequence_ap --help`. Actual evaluation requires all original protected artifacts. It trains the same three heads and compares raw margins with sigmoid probabilities; the original wall time was 2.92 seconds after features were available. Do not expect that timing to include data preparation or ESM extraction.

## Evaluation boundaries

Original v1 final: consumed. S13 biochemical panel: consumed and retrospectively reused. RLY2023 S5 and MANO2021 Data1: reserved. The FGFR2 pilots reuse development folds and are exploratory. Public data may have been in pretrained model data; runtime outcome isolation does not prove absence of memorization.

For scientific confirmation, freeze a selected procedure and a baseline, predeclare the appropriate external assay endpoint, then run a separate evaluation. The publication and frontend work do not authorize consuming these reserved assays.
