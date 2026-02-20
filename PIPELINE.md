# kira-proteoqc Pipeline Contract

## Run mode

- `standalone`: default local workflow.
- `pipeline`: strict pipeline IO contract for `kira-organelle`.

## Shared cache input (pipeline mode)

Lookup is non-recursive in `--input` directory:

1. detect optional prefix from `*_matrix.mtx(.gz)`, `*_features.tsv(.gz)`, `*_barcodes.tsv(.gz)`
2. expected cache path:
   - no prefix: `kira-organelle.bin`
   - with prefix: `<PREFIX>.kira-organelle.bin`

Semantics:

- valid cache: used directly
- missing cache: warning emitted, MTX fallback used
- invalid cache: hard error

Cache format specification: [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md).

## Pipeline output layout

Given `--out <OUT>`, pipeline mode writes:

- `<OUT>/kira-proteoqc/proteoqc.tsv`
- `<OUT>/kira-proteoqc/summary.json`
- `<OUT>/kira-proteoqc/panels_report.tsv`
- `<OUT>/kira-proteoqc/pipeline_step.json`

## File contracts

### `proteoqc.tsv`

Per-cell table with exact deterministic column order:

`barcode sample condition species libsize nnz expressed_genes proteostasis_load misfolded_protein_burden chaperone_capacity proteasome_activity_proxy protein_quality_balance stress_proteostasis_index regime flags confidence`

### `panels_report.tsv`

Panel-level audit table:

`panel_id panel_name panel_group panel_size_defined panel_size_mappable missing_genes coverage_median coverage_p10 sum_median sum_p90 sum_p99`

### `summary.json`

Run-level machine-readable aggregate:

- `tool {name, version, simd}`
- `input {n_cells, species}`
- `distributions {proteostasis_load, misfolded_protein_burden, stress_proteostasis_index}`
- `regimes {counts, fractions}`
- `qc {low_confidence_fraction, low_chaperone_signal_fraction}`

### `pipeline_step.json`

Aggregator ingest manifest:

- fixed `tool` identity with `stage = "proteostasis"`
- artifact paths: `proteoqc.tsv`, `summary.json`, `panels_report.tsv`
- exported cell metrics list
- exported regimes list

## Determinism

- fixed file naming and directory layout in pipeline mode
- fixed TSV column order
- stable regime labels
- sorted-map JSON sections where map ordering matters
- no silent fallback on malformed shared cache
