# kira-proteoqc

Deterministic, explainable proteostasis QC for single-cell expression data.

## Build requirements

- Rust >= 1.95

## Install

Install from crates.io:

```bash
cargo install kira-proteoqc
```


## Usage examples

Standalone run (single input, per-cell outputs):

```bash
kira-proteoqc run \
  --input ./data/pbmc3k \
  --out ./out/pbmc3k \
  --mode cell \
  --json \
  --tsv
```

Standalone run (sample mode):

```bash
kira-proteoqc run \
  --input ./data/inf \
  --out ./out/inf \
  --mode sample \
  --json \
  --tsv
```

Pipeline run (shared cache lookup + pipeline artifacts):

```bash
kira-proteoqc run \
  --input ./data/inf \
  --out ./out/inf \
  --mode sample \
  --run-mode pipeline
```

Validation command:

```bash
kira-proteoqc validate \
  --input ./data/inf \
  --run-mode pipeline
```

Geneset listing:

```bash
kira-proteoqc geneset show
```

Geneset coverage on input:

```bash
kira-proteoqc geneset show --input ./data/inf
```

## Modes

- `--run-mode standalone` (default): existing standalone behavior and outputs.
- `--run-mode pipeline`: pipeline contract mode for `kira-organelle`.

## Pipeline cache lookup

In pipeline mode, `kira-proteoqc` first searches the input directory for shared cache as specified in [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md):

- no prefix: `kira-organelle.bin`
- prefixed dataset: `<PREFIX>.kira-organelle.bin`

Prefix is detected non-recursively from names like `<PREFIX>_matrix.mtx(.gz)`, `<PREFIX>_features.tsv(.gz)`, `<PREFIX>_barcodes.tsv(.gz)`.

Behavior:

- cache exists and valid: use shared cache path.
- cache missing: warn once and fall back to MTX input.
- cache exists but invalid: hard error (no silent fallback).

## Pipeline output contract

In pipeline mode, outputs are written to:

- `--out <DIR>` -> `<DIR>/kira-proteoqc/`

Required artifacts:

- `proteoqc.tsv` (per-cell contract table)
- `summary.json` (run-level aggregates)
- `panels_report.tsv` (panel audit)
- `pipeline_step.json` (ingestion manifest for `kira-organelle`)

All TSV float values are fixed `%.6f`.

## Shared cache specification

- Cache format specification: [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md)
- Reader validates header/magic/version/endian/header-size/file-bytes, header CRC64-ECMA, section bounds, string tables, and CSC invariants.
