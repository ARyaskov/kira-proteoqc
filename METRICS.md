# kira-proteoqc Metrics Specification

This document defines the canonical metric conventions, formulas, units, and JSON output contracts implemented in `kira-proteoqc`.

## Scope

- Covers metric computation currently implemented in code.
- Covers JSON contracts produced by:
  - standalone mode (`proteoqc.json`, schema v1)
  - pipeline mode (`summary.json`, `pipeline_step.json`)
- Does not define external aggregator behavior beyond emitted fields.

## Canonical Conventions

- Numeric base type in compute path: `f32`.
- JSON export type for summary metrics: `f64`.
- Missing/undefined numeric values:
  - standalone JSON uses optional fields where applicable (`Option` -> `null` when absent).
  - pipeline JSON emits concrete numeric values (no `NaN`/`Inf`).
- `NaN` handling:
  - any `NaN` in core score vectors is a hard error.
- Determinism:
  - fixed formulas and fixed operation order.
  - deterministic reduction over canonical CSC expression cache.
  - stable sorting for percentile-like stats (`total_cmp`).

## Units and Scale

- Axis raw and integrated raw scores:
  - Unit: arbitrary expression-derived unit (mean raw count per geneset, then linear combinations).
  - Range: unbounded real values.
- Robust z-scores:
  - Unit: standard deviations in robust space (dimensionless).
  - Range: unbounded real values.
- Pipeline proxy metrics (`proteostasis_load`, etc.):
  - Unit: dimensionless.
  - Range: `[0, 1]` by sigmoid transform.
- Confidence:
  - Unit: dimensionless.
  - Range: `[0, 1]`.

## Statistical Primitives

### Trimmed Mean

For vector `v`, length `n`, trim fraction `p`:

1. sort ascending
2. `k = floor(p * n)`
3. if `n < 2*k + 1`: use arithmetic mean over all values
4. else use mean over `v[k..n-k)`

### Median

- odd `n`: middle value
- even `n`: average of two middle values

### MAD

`MAD = median(|v_i - median(v)|)`

### Robust z-score

`z = (x - median) / (1.4826 * MAD)`

Rule: if `MAD == 0`, return `0`.

## Geneset Reduction Convention

Per geneset, per cell:

1. sum expression values for all genes in geneset
2. divide by `gene_count`

This produces geneset mean expression per cell (raw-count based in current implementation).

## Axis Raw Metrics

All axis metrics are computed per cell (or as single sample vector in sample mode).

- `P_core = mean(proteasome_core geneset)`
- `P_reg = mean(proteasome_regulator geneset)`
- `U_axis = mean(ubiquitin_axis geneset)`
- `E3_axis = mean(e3_ligases geneset)`
- `DUB_axis = mean(dubs geneset)`
- `HSP70 = mean(chaperone_hsp70 geneset)`
- `HSP90 = mean(chaperone_hsp90 geneset)`
- `HSP40 = mean(chaperone_hsp40 geneset)`
- `ERAD = mean(erad geneset)`
- `Ribo = mean(ribosome_load geneset)`

Formulas:

- `PCS_raw = 0.6*P_core + 0.4*P_reg`
- `UTP_raw = 0.5*U_axis + 0.35*E3_axis - 0.15*DUB_axis`
- `CLS_raw = 0.45*HSP70 + 0.35*HSP90 + 0.20*HSP40`
- `ERAD_raw = ERAD`
- `Ribo_raw = Ribo`

## Integrated Metrics

- `Capacity_raw = 0.55*PCS_raw + 0.25*CLS_raw + 0.20*ERAD_raw`
- `PII_raw = Ribo_raw - Capacity_raw`
- `PFS_raw = 0.40*PII_raw + 0.25*UTP_raw + 0.20*Ribo_raw - 0.15*PCS_raw`

PFS contribution decomposition:

- `c_pii = 0.40*PII_raw`
- `c_utp = 0.25*UTP_raw`
- `c_ribo = 0.20*Ribo_raw`
- `c_pcs = -0.15*PCS_raw`

Z-score fields:

- In cell mode: `capacity_z`, `pii_z`, `pfs_z` are computed.
- In sample mode: z-score fields are `null` / absent (`None`).

## Risk Flag Metrics

Flags are deterministic rule-based predicates:

- `fragile_high`
- `proteasome_addiction`
- `proteotoxic_stress`
- `er_degradation_overdrive`

Cell mode evaluates on per-cell vectors and reports fraction in details.
Sample mode evaluates sample-level predicates; `fragile_high` uses derived per-cell PFS estimation path and top-10% condition.

## Pipeline Per-cell Proxy Metrics (`proteoqc.tsv`)

Pipeline table columns include:

- `proteostasis_load = sigmoid(0.5*PII_raw + 0.5*Ribo_raw)`
- `misfolded_protein_burden = sigmoid(PII_raw)`
- `chaperone_capacity = sigmoid(CLS_raw)`
- `proteasome_activity_proxy = sigmoid(PCS_raw)`
- `protein_quality_balance = sigmoid(Capacity_raw - PII_raw)`
- `stress_proteostasis_index = sigmoid(PFS_raw)`

Where `sigmoid(x) = 1 / (1 + exp(-x))`.

Additional pipeline fields:

- `regime`: categorical class from threshold rules over stress/misfolded/proteasome proxies.
- `flags`: comma-separated warnings (`LOW_CONFIDENCE`, `LOW_CHAPERONE_SIGNAL`).
- `confidence`: deterministic score in `[0,1]`, based on warnings, chaperone proxy, and geneset coverage.

## JSON Contract: `proteoqc.json` (Standalone, schema v1)

Top-level required fields:

- `tool: string`
- `version: string`
- `schema_version: "v1"`
- `input_meta`
- `scores`
- `risk_flags`
- `explainability`
- `timecourse` (nullable)

`input_meta`:

- `genes: u64|null`
- `cells: u64|null`
- `nnz: u64|null`
- `normalization: {log1p: bool, cp10k: bool, raw: bool}`
- `mode: "cell"|"sample"`
- `timecourse: bool`

`scores`:

- `per_sample: [ { id, PCS_raw, UTP_raw, CLS_raw, ERAD_raw, Ribo_raw, Capacity_raw, PII_raw, PFS_raw } ] | null`
- `per_cell_tsv_path: string|null`
- `distributions: object|null` (placeholder in current schema)

`risk_flags` entries:

- `name: string`
- `fired: bool`
- `threshold: string|null`
- `details: string|null`

`explainability`:

- `component_contributions: array|null`
- `geneset_coverage: [ { geneset, found, total, fraction } ]`
- `pfs_contributions: { pii, utp, ribo, pcs } | null`

`timecourse` (if present):

- `timepoints: [ { label, pfs, pii, pcs, cls, utp } ]`
- `deltas: [ { from, to, delta_pfs, delta_pii, delta_pcs, delta_cls } ]`
- `trajectory: string`

## JSON Contract: `summary.json` (Pipeline mode)

Top-level required fields:

- `tool: { name, version, simd }`
- `input: { n_cells, species }`
- `distributions: { proteostasis_load, misfolded_protein_burden, stress_proteostasis_index }`
- `regimes: { counts, fractions }`
- `qc: { low_confidence_fraction, low_chaperone_signal_fraction }`

Distribution object type:

- `{ median: f64, p90: f64, p99: f64 }`

Rounding:

- Pipeline summary numeric outputs are rounded to 6 decimals where implemented.

## JSON Contract: `pipeline_step.json` (Pipeline mode)

Top-level required fields:

- `tool: { name, stage, version }`
- `artifacts: { summary, primary_metrics, panels }`
- `cell_metrics: { file, id_column, regime_column, confidence_column, flag_column }`
- `regimes: [string]`

Current canonical artifact names:

- `summary = "summary.json"`
- `primary_metrics = "proteoqc.tsv"`
- `panels = "panels_report.tsv"`

## Field Naming Rules

- JSON fields use `snake_case` except explicit legacy names in standalone score payload (`PCS_raw`, etc.).
- Metric IDs are stable API surface; changing names requires schema/version bump.
