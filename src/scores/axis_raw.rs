use anyhow::{Result, bail};
use tracing::warn;

use crate::ctx::Ctx;
use crate::geneset::ResolvedGeneset;
use crate::math::reduce::GeneSetReducer;
use crate::math::stats::trimmed_mean;
use crate::schema::v1::Mode;
use crate::scores::AxisRawScores;

#[cfg(feature = "fusion")]
use crate::fusion;

pub fn compute_axis_raw(ctx: &mut Ctx) -> Result<AxisRawScores> {
    compute_axis_raw_with_mode(ctx, ctx.mode.clone())
}

pub fn compute_axis_raw_with_mode(ctx: &mut Ctx, mode: Mode) -> Result<AxisRawScores> {
    let mut warnings = std::mem::take(&mut ctx.warnings);
    let mut scratch = std::mem::take(&mut ctx.scratch_cell_buf);

    let result = (|| -> Result<AxisRawScores> {
        let reader = ctx.expr_reader()?;
        let genesets = ctx
            .genesets
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("genesets not resolved"))?;

        #[cfg(feature = "fusion")]
        {
            if ctx.fusion != "off" {
                return compute_axis_raw_fusion(
                    &reader,
                    mode,
                    genesets.resolved.as_slice(),
                    &mut warnings,
                );
            }
        }

        let mut reducer = GeneSetReducer::new(&reader, ctx.threads, ctx.cache_block, ctx.prefetch);

        let n_cells = reader.n_cells();
        let cell_mode = matches!(mode, Mode::Cell);

        let mut pcs = vec![0.0f32; if cell_mode { n_cells } else { 1 }];
        let mut utp = vec![0.0f32; if cell_mode { n_cells } else { 1 }];
        let mut cls = vec![0.0f32; if cell_mode { n_cells } else { 1 }];
        let mut erad = vec![0.0f32; if cell_mode { n_cells } else { 1 }];
        let mut ribo = vec![0.0f32; if cell_mode { n_cells } else { 1 }];

        compute_weighted(
            &mut reducer,
            genesets.resolved.as_slice(),
            "proteasome_core",
            0.6,
            &mut pcs,
            &mut scratch,
            &mut warnings,
            mode.clone(),
        )?;
        compute_weighted(
            &mut reducer,
            genesets.resolved.as_slice(),
            "proteasome_regulator",
            0.4,
            &mut pcs,
            &mut scratch,
            &mut warnings,
            mode.clone(),
        )?;

        compute_weighted(
            &mut reducer,
            genesets.resolved.as_slice(),
            "ubiquitin_axis",
            0.5,
            &mut utp,
            &mut scratch,
            &mut warnings,
            mode.clone(),
        )?;
        compute_weighted(
            &mut reducer,
            genesets.resolved.as_slice(),
            "e3_ligases",
            0.35,
            &mut utp,
            &mut scratch,
            &mut warnings,
            mode.clone(),
        )?;
        compute_weighted(
            &mut reducer,
            genesets.resolved.as_slice(),
            "dubs",
            -0.15,
            &mut utp,
            &mut scratch,
            &mut warnings,
            mode.clone(),
        )?;

        compute_weighted(
            &mut reducer,
            genesets.resolved.as_slice(),
            "chaperone_hsp70",
            0.45,
            &mut cls,
            &mut scratch,
            &mut warnings,
            mode.clone(),
        )?;
        compute_weighted(
            &mut reducer,
            genesets.resolved.as_slice(),
            "chaperone_hsp90",
            0.35,
            &mut cls,
            &mut scratch,
            &mut warnings,
            mode.clone(),
        )?;
        compute_weighted(
            &mut reducer,
            genesets.resolved.as_slice(),
            "chaperone_hsp40",
            0.20,
            &mut cls,
            &mut scratch,
            &mut warnings,
            mode.clone(),
        )?;

        compute_weighted(
            &mut reducer,
            genesets.resolved.as_slice(),
            "erad",
            1.0,
            &mut erad,
            &mut scratch,
            &mut warnings,
            mode.clone(),
        )?;

        compute_weighted(
            &mut reducer,
            genesets.resolved.as_slice(),
            "ribosome_load",
            1.0,
            &mut ribo,
            &mut scratch,
            &mut warnings,
            mode,
        )?;

        check_nan(&pcs)?;
        check_nan(&utp)?;
        check_nan(&cls)?;
        check_nan(&erad)?;
        check_nan(&ribo)?;

        Ok(AxisRawScores {
            pcs,
            utp,
            cls,
            erad,
            ribo,
        })
    })();

    ctx.warnings = warnings;
    ctx.scratch_cell_buf = scratch;
    result
}

#[cfg(feature = "fusion")]
fn compute_axis_raw_fusion(
    reader: &crate::expr::reader::ExprReader<'_>,
    mode: Mode,
    resolved: &[ResolvedGeneset],
    warnings: &mut Vec<String>,
) -> Result<AxisRawScores> {
    let n_cells = reader.n_cells();
    let cell_mode = matches!(mode, Mode::Cell);

    let target_ids = [
        "proteasome_core",
        "proteasome_regulator",
        "ubiquitin_axis",
        "e3_ligases",
        "dubs",
        "chaperone_hsp70",
        "chaperone_hsp90",
        "chaperone_hsp40",
        "erad",
        "ribosome_load",
    ];

    let plan = fusion::build_plan(reader.n_genes(), n_cells, resolved, &target_ids);
    let mut out = vec![0.0f32; plan.targets.len() * n_cells];
    fusion::fused_reduce(reader, &plan, &mut out)?;

    let mut per_target: Vec<Vec<f32>> = Vec::new();
    for t in 0..plan.targets.len() {
        let start = t * n_cells;
        let mut v = out[start..start + n_cells].to_vec();
        let denom = plan.gene_counts[t].max(1) as f32;
        for x in v.iter_mut() {
            *x /= denom;
        }
        per_target.push(v);
    }

    let mut pcs = vec![0.0f32; if cell_mode { n_cells } else { 1 }];
    let mut utp = vec![0.0f32; if cell_mode { n_cells } else { 1 }];
    let mut cls = vec![0.0f32; if cell_mode { n_cells } else { 1 }];
    let mut erad = vec![0.0f32; if cell_mode { n_cells } else { 1 }];
    let mut ribo = vec![0.0f32; if cell_mode { n_cells } else { 1 }];

    let mut map = std::collections::HashMap::new();
    for (i, id) in plan.targets.iter().enumerate() {
        map.insert(id.as_str(), i);
    }

    add_weighted_from_target(
        "proteasome_core",
        0.6,
        &per_target,
        &map,
        &mut pcs,
        warnings,
        mode.clone(),
    )?;
    add_weighted_from_target(
        "proteasome_regulator",
        0.4,
        &per_target,
        &map,
        &mut pcs,
        warnings,
        mode.clone(),
    )?;
    add_weighted_from_target(
        "ubiquitin_axis",
        0.5,
        &per_target,
        &map,
        &mut utp,
        warnings,
        mode.clone(),
    )?;
    add_weighted_from_target(
        "e3_ligases",
        0.35,
        &per_target,
        &map,
        &mut utp,
        warnings,
        mode.clone(),
    )?;
    add_weighted_from_target(
        "dubs",
        -0.15,
        &per_target,
        &map,
        &mut utp,
        warnings,
        mode.clone(),
    )?;
    add_weighted_from_target(
        "chaperone_hsp70",
        0.45,
        &per_target,
        &map,
        &mut cls,
        warnings,
        mode.clone(),
    )?;
    add_weighted_from_target(
        "chaperone_hsp90",
        0.35,
        &per_target,
        &map,
        &mut cls,
        warnings,
        mode.clone(),
    )?;
    add_weighted_from_target(
        "chaperone_hsp40",
        0.20,
        &per_target,
        &map,
        &mut cls,
        warnings,
        mode.clone(),
    )?;
    add_weighted_from_target(
        "erad",
        1.0,
        &per_target,
        &map,
        &mut erad,
        warnings,
        mode.clone(),
    )?;
    add_weighted_from_target(
        "ribosome_load",
        1.0,
        &per_target,
        &map,
        &mut ribo,
        warnings,
        mode,
    )?;

    Ok(AxisRawScores {
        pcs,
        utp,
        cls,
        erad,
        ribo,
    })
}

#[cfg(feature = "fusion")]
fn add_weighted_from_target(
    id: &str,
    weight: f32,
    per_target: &[Vec<f32>],
    map: &std::collections::HashMap<&str, usize>,
    out: &mut [f32],
    warnings: &mut Vec<String>,
    mode: Mode,
) -> Result<()> {
    if let Some(&idx) = map.get(id) {
        if matches!(mode, Mode::Cell) {
            for (o, v) in out.iter_mut().zip(per_target[idx].iter()) {
                *o += weight * *v;
            }
        } else {
            let mean = trimmed_mean(&mut per_target[idx].clone(), SAMPLE_TRIM_P);
            out[0] += weight * mean;
        }
        Ok(())
    } else {
        warnings.push(format!("missing geneset '{}'", id));
        Ok(())
    }
}

const SAMPLE_TRIM_P: f32 = 0.0;

fn compute_weighted(
    reducer: &mut GeneSetReducer<'_>,
    resolved: &[ResolvedGeneset],
    id: &str,
    weight: f32,
    out: &mut [f32],
    scratch: &mut Vec<f32>,
    warnings: &mut Vec<String>,
    mode: Mode,
) -> Result<()> {
    let geneset = resolved.iter().find(|g| g.id == id);
    if let Some(gs) = geneset {
        if gs.gene_ids.is_empty() {
            warn!("geneset '{}' resolved to 0 genes", id);
            warnings.push(format!("geneset '{}' resolved to 0 genes", id));
            return Ok(());
        }
        let cell_mode = matches!(mode, Mode::Cell);
        if cell_mode {
            scratch.resize(out.len(), 0.0);
            reducer.per_cell_raw(&gs.gene_ids, scratch)?;
            for (o, v) in out.iter_mut().zip(scratch.iter()) {
                *o += weight * *v;
            }
        } else {
            scratch.resize(reducer.expr.n_cells(), 0.0);
            reducer.per_cell_raw(&gs.gene_ids, scratch)?;
            let mean = trimmed_mean(scratch, SAMPLE_TRIM_P);
            out[0] += weight * mean;
        }
        Ok(())
    } else {
        warn!("missing geneset '{}'", id);
        warnings.push(format!("missing geneset '{}'", id));
        Ok(())
    }
}

fn check_nan(values: &[f32]) -> Result<()> {
    for v in values {
        if v.is_nan() {
            bail!("NaN encountered in axis raw score");
        }
    }
    Ok(())
}
