# Reproduce Escape Lab

Our ambition is an autonomous laboratory that learns how to make better scientific predictions. Reproduction is how another researcher can check each step. The public runner now rebuilds the FGFR2 development experiment from public source data and fits the same three classifier heads without private historical predictions, credentials, or LLM calls.

## Start in a clean clone

```sh
git clone https://github.com/Gentle-mann/escape-lab.git
cd escape-lab
uv sync --frozen --python 3.12
uv run --frozen escape-lab-reproduce doctor
uv run --frozen escape-lab-reproduce saturation-demo
uv run --frozen pytest -q
uv run --frozen ruff check src tests scripts
```

`doctor` checks the packaged frozen cohort and installed dependencies without network access or loading outcomes. `saturation-demo` takes four synthetic margins, shows how sigmoid rounding collapses them to one probability, and compares their AP. Its 0.5-to-1.0 AP example explains the numerical mechanism; the scientific result comes from the real-data command below.

## Rebuild the headline development experiment

```sh
uv sync --frozen --python 3.12 --extra sequence
uv run --frozen --extra sequence escape-lab-reproduce prepare \
  --output reports/reproduction/inputs-001
HF_HUB_DISABLE_XET=1 TOKENIZERS_PARALLELISM=false \
  uv run --frozen --extra sequence escape-lab-reproduce headline \
  --inputs reports/reproduction/inputs-001 \
  --output reports/reproduction/run-001 \
  --device mps
```

Use `--device cpu` on machines without Apple MPS. MPS, batch size 4, and Python 3.12 reproduce the original extraction configuration. The optional `sequence` extra pins Torch 2.8.0, Transformers 4.57.1, Hugging Face Hub 0.36.2, Safetensors 0.8.0, and their full dependency resolution in `uv.lock`. The standard Torch distribution can install substantial platform-specific dependencies, especially on Linux. Model weights are downloaded from the pinned public checkpoint, separately from Python packages.

The commands perform this fixed recipe:

1. Download eight MAVEdb metadata records, two RCSB structures (7WCL and 6MZW), and the study supplement from Europe PMC. Verify exact coordinate and workbook checksums; rebuild and verify the original biological feature and reference hashes.
2. Use the packaged list of **668 FGFR2 development variants**, ten occupied residue blocks, and the original three folds. The list contains identities and grouping metadata only. Reserved groups are explicitly excluded.
3. Extract AlphaMissense and EVE predictor scores from S10 using an identity allowlist before parsing scores. Collapse codon representations only when their scores agree. Derive activation labels with the existing development-only S10 importer.
4. Extract frozen ESM2-35M wild-type context and mutant-minus-wild-type features, using checkpoint revision `6fbf070e65b0b7291e7bbcd451118c216cff79d8` and verified checkpoint file hashes. Fit three fixed legacy classifiers. Scaling and PCA fit inside each training fold.
5. Rank each held-out fold with the same classifier's probabilities and decision margins. Write aggregate results, fitting time, package versions, checksums, diagnostics, and a comparison with the published values.

The published aggregate targets are:

| Ranking | Mean-fold activation AP | High-activation hits | Total picks |
|---|---:|---:|---:|
| Probability | 0.0973209165 | 1 | 68 |
| Decision margin | 0.1113457484 | 4 | 68 |

**Verified from a source-only temporary workspace:** public input preparation took 40.82 seconds, followed by fresh ESM extraction and evaluation in 59.73 seconds, including 0.31 seconds for the three heads. Both AP values and the 1-to-4 hit improvement reproduced exactly. These timings exclude first-time dependency installation. [The verification receipt](archive/reports/reproduction/verification.json) records input and checkpoint provenance, package versions and the retained failed preparation attempt. No historical embeddings or predictions were supplied.

`results.json` reports `comparison.matches_published_aggregate`; it must be true before calling a run an aggregate reproduction. The `headline` command exits with code 1 if its completed evaluation does not match the published aggregate. CPU/GPU backend differences can change floating-point embeddings and downstream rankings, especially in the historical near-constant feature branch. Every run reports its measured result and whether it matches. No hyperparameter search or automatic retry changes the result.

Both preparation and evaluation require a new output directory. Failed attempts retain their status. All generated data remain under ignored `reports/`; individual predictions are written with owner-only permissions. No previous experiment receipt, private prediction file, external assay, or original final split is needed. This reproduces the fixed development comparison; independent biological confirmation remains a separate experiment.

For an already downloaded original workbook, pass `--workbook /absolute/path/to/workbook.xlsx` to `prepare`; it must match the frozen checksum. For a previously extracted and verified ESM cache, add `--feature-cache /absolute/path/to/cache` to `headline`. The result explicitly records cached features. Without that option, features are freshly extracted.

## Source and license provenance

The recipe is packaged at [`src/escape_lab/reproduction_assets/fgfr2-headline-v1.json`](../src/escape_lab/reproduction_assets/fgfr2-headline-v1.json), including cohort identities, folds, model revision, and input hashes. [The FGFR study](https://doi.org/10.1038/s41588-025-02431-8) supplies the measurements. Covered article material is CC BY 4.0; source MAVEdb metadata and derived cohort metadata are CC BY-SA 4.0; RCSB coordinates are CC0. Third-party predictor score terms remain applicable. Source workbooks, scores and model weights are retrieved locally rather than redistributed in this repository.

## Agent and earlier experiment paths

The public-source runner above is independent of paid inference. For Nemotron on Crusoe, copy `.env.example` to `.env` and set credentials locally. Discover available models with `uv run --frozen escape-lab crusoe-models`; historical IDs were `nvidia/NVIDIA-Nemotron-3.5-Super-VL` and `zai-org/GLM-5.3`.

Original resistance replay remains available:

```sh
uv run --frozen escape-lab prepare
uv run --frozen escape-lab status
uv run --frozen escape-lab baseline --stage development --seeds 0 1 2
uv run --frozen escape-lab serve --port 8017
```

That ingestion downloads the eight MAVEdb score sets and normalizes 8,407 paired nucleotide variants into 7,522 protein variants. It is a separate experiment from the 668-variant activation reproduction above.

Historical agent and structure runners retain their original artifact checks. Their protocols and aggregate results are in [the experiment ledger](EXPERIMENTS.md) and [the archive](archive/manifest.json). They do not replace the new public-source reproduction command.

Original v1 final and the S13 biochemical panel were previously consumed. RLY2023 S5 and MANO2021 Data1 remain reserved. Select and freeze a procedure, baseline and assay-compatible endpoint before a separate external evaluation.
