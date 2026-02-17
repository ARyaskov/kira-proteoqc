use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::ctx::Ctx;
use crate::schema::v1::{
    DeltaSummary, Explainability, GenesetCoverage, InputMeta, Mode, Normalization, PerSampleScore,
    PfsContributions, ProteoQcV1, RiskFlag, Scores, TimecourseResult, TimepointSummary,
};

pub fn build_report(ctx: &Ctx) -> Result<ProteoQcV1> {
    let input_meta = InputMeta {
        genes: ctx.input_meta.genes,
        cells: ctx.input_meta.cells,
        nnz: ctx.input_meta.nnz,
        normalization: Normalization {
            log1p: ctx.log1p,
            cp10k: true,
            raw: true,
        },
        mode: ctx.mode.clone(),
        timecourse: ctx.timecourse,
    };

    let axis = ctx.axis_raw.as_ref().context("axis raw scores missing")?;
    let integrated = ctx
        .integrated_scores
        .as_ref()
        .context("integrated scores missing")?;

    let per_sample = vec![PerSampleScore {
        id: "sample".to_string(),
        pcs_raw: Some(mean_vec(&axis.pcs)? as f64),
        utp_raw: Some(mean_vec(&axis.utp)? as f64),
        cls_raw: Some(mean_vec(&axis.cls)? as f64),
        erad_raw: Some(mean_vec(&axis.erad)? as f64),
        ribo_raw: Some(mean_vec(&axis.ribo)? as f64),
        capacity_raw: Some(mean_vec(&integrated.capacity_raw)? as f64),
        pii_raw: Some(mean_vec(&integrated.pii_raw)? as f64),
        pfs_raw: Some(mean_vec(&integrated.pfs_raw)? as f64),
    }];

    let per_cell_tsv_path = if matches!(ctx.mode, Mode::Cell) && ctx.write_tsv {
        Some("proteoqc.tsv".to_string())
    } else {
        None
    };

    let scores = Scores {
        per_sample: Some(per_sample),
        per_cell_tsv_path,
        distributions: None,
    };

    let risk_flags = ctx
        .risk_flags
        .iter()
        .map(|f| RiskFlag {
            name: f.name.clone(),
            fired: f.fired,
            threshold: Some(f.threshold.clone()),
            details: f.details.clone(),
        })
        .collect::<Vec<_>>();

    let geneset_coverage = if let Some(gs) = &ctx.genesets {
        gs.resolved
            .iter()
            .map(|g| GenesetCoverage {
                geneset: g.id.clone(),
                found: g.gene_ids.len() as u64,
                total: g.total as u64,
                fraction: if g.total == 0 {
                    0.0
                } else {
                    g.gene_ids.len() as f64 / g.total as f64
                },
            })
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };

    let pfs_contrib = if let Some(c) = &ctx.pfs_contributions {
        Some(PfsContributions {
            pii: mean_vec(&c.c_pii)? as f64,
            utp: mean_vec(&c.c_utp)? as f64,
            ribo: mean_vec(&c.c_ribo)? as f64,
            pcs: mean_vec(&c.c_pcs)? as f64,
        })
    } else {
        None
    };

    let explainability = Explainability {
        component_contributions: None,
        geneset_coverage,
        pfs_contributions: pfs_contrib,
    };

    let timecourse = ctx.timecourse_result.as_ref().map(|tc| TimecourseResult {
        timepoints: tc
            .timepoints
            .iter()
            .map(|t| TimepointSummary {
                label: t.label.clone(),
                pfs: t.pfs,
                pii: t.pii,
                pcs: t.pcs,
                cls: t.cls,
                utp: t.utp,
            })
            .collect(),
        deltas: tc
            .deltas
            .iter()
            .map(|d| DeltaSummary {
                from: d.from.clone(),
                to: d.to.clone(),
                delta_pfs: d.delta_pfs,
                delta_pii: d.delta_pii,
                delta_pcs: d.delta_pcs,
                delta_cls: d.delta_cls,
            })
            .collect(),
        trajectory: tc.trajectory.clone(),
    });

    Ok(ProteoQcV1 {
        tool: "kira-proteoqc".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        schema_version: "v1".to_string(),
        input_meta,
        scores,
        risk_flags,
        explainability,
        timecourse,
    })
}

pub fn write_json(path: &Path, ctx: &Ctx) -> Result<()> {
    let report = build_report(ctx)?;
    let file = std::fs::File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let writer = std::io::BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &report)?;
    Ok(())
}

fn mean_vec(values: &[f32]) -> Result<f32> {
    if values.is_empty() {
        return Ok(0.0);
    }
    let mut sum = 0.0;
    for v in values {
        if v.is_nan() {
            bail!("NaN encountered in score vector");
        }
        sum += *v;
    }
    Ok(sum / values.len() as f32)
}
