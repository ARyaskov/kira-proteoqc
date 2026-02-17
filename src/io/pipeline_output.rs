use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use serde::Serialize;

use crate::ctx::Ctx;
use crate::math::reduce::GeneSetReducer;

const PIPELINE_DIR: &str = "kira-proteoqc";
const IO_BUF_CAPACITY: usize = 1 << 20; // 1 MiB

#[derive(Debug, Clone, Serialize)]
struct ToolMeta {
    name: String,
    version: String,
    simd: String,
}

#[derive(Debug, Clone, Serialize)]
struct SummaryInput {
    n_cells: usize,
    species: String,
}

#[derive(Debug, Clone, Serialize)]
struct DistStats {
    median: f64,
    p90: f64,
    p99: f64,
}

#[derive(Debug, Clone, Serialize)]
struct Regimes {
    counts: BTreeMap<String, usize>,
    fractions: BTreeMap<String, f64>,
}

#[derive(Debug, Clone, Serialize)]
struct SummaryQc {
    low_confidence_fraction: f64,
    low_chaperone_signal_fraction: f64,
}

#[derive(Debug, Clone, Serialize)]
struct PipelineSummary {
    tool: ToolMeta,
    input: SummaryInput,
    distributions: SummaryDistributions,
    regimes: Regimes,
    qc: SummaryQc,
}

#[derive(Debug, Clone, Serialize)]
struct SummaryDistributions {
    proteostasis_load: DistStats,
    misfolded_protein_burden: DistStats,
    stress_proteostasis_index: DistStats,
}

#[derive(Debug, Clone, Serialize)]
struct PipelineStepTool {
    name: String,
    stage: String,
    version: String,
}

#[derive(Debug, Clone, Serialize)]
struct PipelineArtifacts {
    summary: String,
    primary_metrics: String,
    panels: String,
}

#[derive(Debug, Clone, Serialize)]
struct PipelineCellMetrics {
    file: String,
    id_column: String,
    regime_column: String,
    confidence_column: String,
    flag_column: String,
}

#[derive(Debug, Clone, Serialize)]
struct PipelineStep {
    tool: PipelineStepTool,
    artifacts: PipelineArtifacts,
    cell_metrics: PipelineCellMetrics,
    regimes: Vec<String>,
}

pub fn ensure_pipeline_out_dir(base_out: &Path) -> Result<std::path::PathBuf> {
    let out = base_out.join(PIPELINE_DIR);
    std::fs::create_dir_all(&out).with_context(|| format!("failed to create {}", out.display()))?;
    Ok(out)
}

pub fn write_pipeline_outputs(ctx: &Ctx, out_dir: &Path) -> Result<()> {
    write_pipeline_cell_tsv(ctx, &out_dir.join("proteoqc.tsv"))?;
    write_panels_report(ctx, &out_dir.join("panels_report.tsv"))?;
    write_summary_json(ctx, &out_dir.join("summary.json"))?;
    write_pipeline_step_json(&out_dir.join("pipeline_step.json"))?;
    Ok(())
}

fn write_pipeline_cell_tsv(ctx: &Ctx, path: &Path) -> Result<()> {
    let axis = ctx.axis_raw.as_ref().context("axis raw scores missing")?;
    let integrated = ctx
        .integrated_scores
        .as_ref()
        .context("integrated scores missing")?;
    let n_cells = ctx.cells.len();
    let row_mode = detect_row_mode(
        &[
            ("axis pcs", axis.pcs.len()),
            ("axis cls", axis.cls.len()),
            ("axis ribo", axis.ribo.len()),
            ("integrated capacity", integrated.capacity_raw.len()),
            ("integrated pii", integrated.pii_raw.len()),
            ("integrated pfs", integrated.pfs_raw.len()),
        ],
        n_cells,
    )?;
    let n_rows = if row_mode == RowMode::PerCell {
        n_cells
    } else {
        1
    };

    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::with_capacity(IO_BUF_CAPACITY, file);

    writeln!(
        w,
        "barcode\tsample\tcondition\tspecies\tlibsize\tnnz\texpressed_genes\tproteostasis_load\tmisfolded_protein_burden\tchaperone_capacity\tproteasome_activity_proxy\tprotein_quality_balance\tstress_proteostasis_index\tregime\tflags\tconfidence"
    )?;

    let row_indices = if row_mode == RowMode::PerCell {
        let mut ordered = (0..n_rows).collect::<Vec<usize>>();
        ordered.sort_unstable_by(|a, b| ctx.cells[*a].cmp(&ctx.cells[*b]));
        ordered
    } else {
        vec![0usize]
    };

    for i in row_indices {
        let idx = if row_mode == RowMode::PerCell { i } else { 0 };
        let libsize = 0u64;
        let nnz = 0u64;
        let expressed_genes = 0u64;

        let misfolded = sigmoid(integrated.pii_raw[idx]) as f64;
        let chaperone = sigmoid(axis.cls[idx]) as f64;
        let proteasome = sigmoid(axis.pcs[idx]) as f64;
        let load = sigmoid(0.5 * integrated.pii_raw[idx] + 0.5 * axis.ribo[idx]) as f64;
        let balance = sigmoid(integrated.capacity_raw[idx] - integrated.pii_raw[idx]) as f64;
        let stress = sigmoid(integrated.pfs_raw[idx]) as f64;

        let confidence = compute_confidence(ctx, chaperone);
        let regime = classify_regime(stress, misfolded, proteasome);
        let flags = build_flags(confidence, chaperone);
        let barcode = if row_mode == RowMode::PerCell {
            ctx.cells[i].as_str()
        } else {
            "sample"
        };

        writeln!(
            w,
            "{}\tsample\tunknown\tunknown\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}",
            barcode,
            libsize,
            nnz,
            expressed_genes,
            load,
            misfolded,
            chaperone,
            proteasome,
            balance,
            stress,
            regime,
            flags,
            confidence
        )?;
    }
    Ok(())
}

fn write_panels_report(ctx: &Ctx, path: &Path) -> Result<()> {
    let collection = ctx.genesets.as_ref().context("genesets missing")?;
    let expr = ctx.expr_reader()?;
    let mut reducer = GeneSetReducer::new(&expr, ctx.threads, ctx.cache_block, ctx.prefetch);

    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::with_capacity(IO_BUF_CAPACITY, file);
    writeln!(
        w,
        "panel_id\tpanel_name\tpanel_group\tpanel_size_defined\tpanel_size_mappable\tmissing_genes\tcoverage_median\tcoverage_p10\tsum_median\tsum_p90\tsum_p99"
    )?;

    for gs in &collection.resolved {
        let mut raw = vec![0.0f32; expr.n_cells()];
        if !gs.gene_ids.is_empty() {
            reducer.per_cell_raw(&gs.gene_ids, &mut raw)?;
        }
        let mut sums = raw
            .iter()
            .map(|v| *v as f64 * gs.gene_ids.len() as f64)
            .collect::<Vec<_>>();
        let mut cov = vec![
            if gs.total == 0 {
                0.0
            } else {
                gs.gene_ids.len() as f64 / gs.total as f64
            };
            expr.n_cells()
        ];

        sort_f64(&mut cov);
        sort_f64(&mut sums);
        let cov_median = percentile_sorted(&cov, 0.50);
        let cov_p10 = percentile_sorted(&cov, 0.10);
        let sum_median = percentile_sorted(&sums, 0.50);
        let sum_p90 = percentile_sorted(&sums, 0.90);
        let sum_p99 = percentile_sorted(&sums, 0.99);

        let missing = if gs.missing.is_empty() {
            String::new()
        } else {
            gs.missing.join(",")
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            gs.id,
            gs.id,
            gs.axis,
            gs.total,
            gs.gene_ids.len(),
            missing,
            cov_median,
            cov_p10,
            sum_median,
            sum_p90,
            sum_p99
        )?;
    }
    Ok(())
}

fn write_summary_json(ctx: &Ctx, path: &Path) -> Result<()> {
    let axis = ctx.axis_raw.as_ref().context("axis raw scores missing")?;
    let integrated = ctx
        .integrated_scores
        .as_ref()
        .context("integrated scores missing")?;
    let n_cells = ctx.cells.len();
    let row_mode = detect_row_mode(
        &[
            ("axis pcs", axis.pcs.len()),
            ("axis cls", axis.cls.len()),
            ("axis ribo", axis.ribo.len()),
            ("integrated capacity", integrated.capacity_raw.len()),
            ("integrated pii", integrated.pii_raw.len()),
            ("integrated pfs", integrated.pfs_raw.len()),
        ],
        n_cells,
    )?;
    let n = if row_mode == RowMode::PerCell {
        n_cells
    } else {
        1
    };
    let mut regimes_count: BTreeMap<String, usize> = BTreeMap::new();
    let mut low_conf = 0usize;
    let mut low_chaperone = 0usize;

    let mut load_vals = Vec::with_capacity(n);
    let mut misfolded_vals = Vec::with_capacity(n);
    let mut stress_vals = Vec::with_capacity(n);
    for i in 0..n {
        let idx = if row_mode == RowMode::PerCell { i } else { 0 };
        let misfolded = sigmoid(integrated.pii_raw[idx]) as f64;
        let chaperone = sigmoid(axis.cls[idx]) as f64;
        let proteasome = sigmoid(axis.pcs[idx]) as f64;
        let load = sigmoid(0.5 * integrated.pii_raw[idx] + 0.5 * axis.ribo[idx]) as f64;
        let stress = sigmoid(integrated.pfs_raw[idx]) as f64;
        let confidence = compute_confidence(ctx, chaperone);
        let regime = classify_regime(stress, misfolded, proteasome).to_string();

        *regimes_count.entry(regime).or_insert(0) += 1;
        if confidence < 0.5 {
            low_conf += 1;
        }
        if chaperone < 0.25 {
            low_chaperone += 1;
        }
        load_vals.push(load);
        misfolded_vals.push(misfolded);
        stress_vals.push(stress);
    }

    let dist = SummaryDistributions {
        proteostasis_load: stats_from_values(&mut load_vals),
        misfolded_protein_burden: stats_from_values(&mut misfolded_vals),
        stress_proteostasis_index: stats_from_values(&mut stress_vals),
    };

    let mut regimes_fraction = BTreeMap::new();
    for (k, v) in &regimes_count {
        regimes_fraction.insert(k.clone(), round6(*v as f64 / n.max(1) as f64));
    }
    for name in regime_names() {
        regimes_count.entry(name.to_string()).or_insert(0);
        regimes_fraction.entry(name.to_string()).or_insert(0.0);
    }

    let summary = PipelineSummary {
        tool: ToolMeta {
            name: "kira-proteoqc".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            simd: crate::simd::backend_name().to_string(),
        },
        input: SummaryInput {
            n_cells,
            species: "unknown".to_string(),
        },
        distributions: dist,
        regimes: Regimes {
            counts: regimes_count,
            fractions: regimes_fraction,
        },
        qc: SummaryQc {
            low_confidence_fraction: round6(low_conf as f64 / n.max(1) as f64),
            low_chaperone_signal_fraction: round6(low_chaperone as f64 / n.max(1) as f64),
        },
    };

    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let writer = BufWriter::with_capacity(IO_BUF_CAPACITY, file);
    serde_json::to_writer_pretty(writer, &summary)?;
    Ok(())
}

fn write_pipeline_step_json(path: &Path) -> Result<()> {
    let step = PipelineStep {
        tool: PipelineStepTool {
            name: "kira-proteoqc".to_string(),
            stage: "proteostasis".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
        },
        artifacts: PipelineArtifacts {
            summary: "summary.json".to_string(),
            primary_metrics: "proteoqc.tsv".to_string(),
            panels: "panels_report.tsv".to_string(),
        },
        cell_metrics: PipelineCellMetrics {
            file: "proteoqc.tsv".to_string(),
            id_column: "barcode".to_string(),
            regime_column: "regime".to_string(),
            confidence_column: "confidence".to_string(),
            flag_column: "flags".to_string(),
        },
        regimes: regime_names().iter().map(|s| s.to_string()).collect(),
    };
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let writer = BufWriter::with_capacity(IO_BUF_CAPACITY, file);
    serde_json::to_writer_pretty(writer, &step)?;
    Ok(())
}

fn classify_regime(stress: f64, misfolded: f64, proteasome: f64) -> &'static str {
    if stress < 0.25 && misfolded < 0.30 {
        "BalancedProteostasis"
    } else if stress < 0.50 && proteasome >= 0.50 {
        "CompensatedStress"
    } else if stress >= 0.50 && misfolded >= 0.55 && proteasome >= 0.45 {
        "ProteotoxicStress"
    } else if stress >= 0.60 && proteasome >= 0.70 {
        "ProteasomeOverload"
    } else if stress >= 0.75 && misfolded >= 0.70 && proteasome < 0.45 {
        "ProteostasisCollapse"
    } else {
        "Unclassified"
    }
}

fn regime_names() -> [&'static str; 6] {
    [
        "BalancedProteostasis",
        "CompensatedStress",
        "ProteotoxicStress",
        "ProteasomeOverload",
        "ProteostasisCollapse",
        "Unclassified",
    ]
}

fn build_flags(confidence: f64, chaperone_capacity: f64) -> String {
    let mut flags = Vec::new();
    if confidence < 0.5 {
        flags.push("LOW_CONFIDENCE");
    }
    if chaperone_capacity < 0.25 {
        flags.push("LOW_CHAPERONE_SIGNAL");
    }
    flags.join(",")
}

fn compute_confidence(ctx: &Ctx, chaperone_capacity: f64) -> f64 {
    let mut score = 1.0f64;
    if !ctx.warnings.is_empty() {
        score -= 0.15;
    }
    if chaperone_capacity < 0.25 {
        score -= 0.20;
    }
    if let Some(gs) = &ctx.genesets {
        let mut total = 0.0f64;
        for r in &gs.resolved {
            let frac = if r.total == 0 {
                0.0
            } else {
                r.gene_ids.len() as f64 / r.total as f64
            };
            total += frac;
        }
        if !gs.resolved.is_empty() {
            score *= (total / gs.resolved.len() as f64).clamp(0.0, 1.0);
        }
    }
    score.clamp(0.0, 1.0)
}

fn stats_from_values(values: &mut Vec<f64>) -> DistStats {
    sort_f64(values);
    DistStats {
        median: round6(percentile_sorted(values, 0.50)),
        p90: round6(percentile_sorted(values, 0.90)),
        p99: round6(percentile_sorted(values, 0.99)),
    }
}

fn round6(v: f64) -> f64 {
    (v * 1_000_000.0).round() / 1_000_000.0
}

fn sort_f64(values: &mut [f64]) {
    values.sort_unstable_by(f64::total_cmp);
}

fn percentile_sorted(values: &[f64], p: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let idx = ((values.len() - 1) as f64 * p).floor() as usize;
    values[idx]
}

fn sigmoid(x: f32) -> f32 {
    1.0 / (1.0 + (-x).exp())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RowMode {
    PerCell,
    PerSample,
}

fn detect_row_mode(lengths: &[(&str, usize)], n_cells: usize) -> Result<RowMode> {
    let all_cells = lengths.iter().all(|(_, len)| *len == n_cells);
    if all_cells {
        return Ok(RowMode::PerCell);
    }
    let all_sample = lengths.iter().all(|(_, len)| *len == 1);
    if all_sample {
        return Ok(RowMode::PerSample);
    }
    for (name, len) in lengths {
        if *len != n_cells && *len != 1 {
            bail!(
                "{} length mismatch: {} (expected {} or 1)",
                name,
                len,
                n_cells
            );
        }
    }
    bail!("inconsistent score vector lengths for pipeline output")
}
