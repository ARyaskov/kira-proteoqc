# kira-proteoqc

Deterministic, explainable proteostasis QC for single-cell expression data.

## Build requirements

- Rust >= 1.95

## Modes

- `--run-mode standalone` (default): existing standalone behavior and outputs.
- `--run-mode pipeline`: pipeline contract mode for `kira-organelle`.

## Pipeline cache lookup

In pipeline mode, `kira-proteoqc` first searches the input directory for shared cache as specified in `CACHE_FILE.md`:

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

- Canonical format: `CACHE_FILE.md`
- Reader validates header/magic/version/endian/header-size/file-bytes, header CRC64-ECMA, section bounds, string tables, and CSC invariants.
