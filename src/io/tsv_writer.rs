use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::ctx::Ctx;
use crate::schema::v1::Mode;

pub fn write_tsv(path: &Path, ctx: &Ctx) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::new(file);

    let axis = ctx.axis_raw.as_ref().context("axis raw scores missing")?;
    let integrated = ctx
        .integrated_scores
        .as_ref()
        .context("integrated scores missing")?;
    let extension = ctx.proteostasis_extension.as_ref().map(|r| &r.scores);

    match ctx.mode {
        Mode::Cell => {
            let n = ctx.cells.len();
            ensure_len(axis.pcs.len(), n, "pcs")?;
            ensure_len(axis.utp.len(), n, "utp")?;
            ensure_len(axis.cls.len(), n, "cls")?;
            ensure_len(axis.erad.len(), n, "erad")?;
            ensure_len(axis.ribo.len(), n, "ribo")?;
            ensure_len(integrated.capacity_raw.len(), n, "capacity")?;
            ensure_len(integrated.pii_raw.len(), n, "pii")?;
            ensure_len(integrated.pfs_raw.len(), n, "pfs")?;

            let pfs_z = integrated
                .pfs_z
                .as_ref()
                .context("PFS_z missing in per-cell mode")?;
            ensure_len(pfs_z.len(), n, "pfs_z")?;
            if let Some(ext) = extension {
                ensure_len(ext.chaperone_core.len(), n, "chaperone_core")?;
                ensure_len(ext.proteasome_core.len(), n, "proteasome_core")?;
                ensure_len(ext.upr_core.len(), n, "upr_core")?;
                ensure_len(ext.agg_core.len(), n, "agg_core")?;
                ensure_len(ext.cci.len(), n, "cci")?;
                ensure_len(ext.pci.len(), n, "pci")?;
                ensure_len(ext.upr_a.len(), n, "upr_a")?;
                ensure_len(ext.pls.len(), n, "pls")?;
                ensure_len(ext.sci.len(), n, "sci")?;
                ensure_len(ext.pcp.len(), n, "pcp")?;
                ensure_len(ext.chaperone_high.len(), n, "chaperone_high")?;
                ensure_len(ext.proteasome_high.len(), n, "proteasome_high")?;
                ensure_len(ext.upr_active.len(), n, "upr_active")?;
                ensure_len(ext.proteotoxic_high.len(), n, "proteotoxic_high")?;
                ensure_len(ext.imbalance_high.len(), n, "imbalance_high")?;
                ensure_len(ext.collapse_risk.len(), n, "collapse_risk")?;
            }

            writeln!(
                w,
                "cell_id\tPCS_raw\tUTP_raw\tCLS_raw\tERAD_raw\tRibo_raw\tCapacity_raw\tPII_raw\tPFS_raw\tPFS_z\tchaperone_core\tproteasome_core\tupr_core\tagg_core\tCCI\tPCI\tUPR_A\tPLS\tSCI\tPCP\tchaperone_high\tproteasome_high\tupr_active\tproteotoxic_high\timbalance_high\tcollapse_risk"
            )?;
            for i in 0..n {
                let (
                    chaperone_core,
                    proteasome_core,
                    upr_core,
                    agg_core,
                    cci,
                    pci,
                    upr_a,
                    pls,
                    sci,
                    pcp,
                    chaperone_high,
                    proteasome_high,
                    upr_active,
                    proteotoxic_high,
                    imbalance_high,
                    collapse_risk,
                ) = if let Some(ext) = extension {
                    (
                        ext.chaperone_core[i],
                        ext.proteasome_core[i],
                        ext.upr_core[i],
                        ext.agg_core[i],
                        ext.cci[i],
                        ext.pci[i],
                        ext.upr_a[i],
                        ext.pls[i],
                        ext.sci[i],
                        ext.pcp[i],
                        ext.chaperone_high[i],
                        ext.proteasome_high[i],
                        ext.upr_active[i],
                        ext.proteotoxic_high[i],
                        ext.imbalance_high[i],
                        ext.collapse_risk[i],
                    )
                } else {
                    (
                        f32::NAN,
                        f32::NAN,
                        f32::NAN,
                        f32::NAN,
                        f32::NAN,
                        f32::NAN,
                        f32::NAN,
                        f32::NAN,
                        f32::NAN,
                        f32::NAN,
                        false,
                        false,
                        false,
                        false,
                        false,
                        false,
                    )
                };
                writeln!(
                    w,
                    "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}",
                    ctx.cells[i],
                    axis.pcs[i],
                    axis.utp[i],
                    axis.cls[i],
                    axis.erad[i],
                    axis.ribo[i],
                    integrated.capacity_raw[i],
                    integrated.pii_raw[i],
                    integrated.pfs_raw[i],
                    pfs_z[i],
                    chaperone_core,
                    proteasome_core,
                    upr_core,
                    agg_core,
                    cci,
                    pci,
                    upr_a,
                    pls,
                    sci,
                    pcp,
                    chaperone_high,
                    proteasome_high,
                    upr_active,
                    proteotoxic_high,
                    imbalance_high,
                    collapse_risk
                )?;
            }
        }
        Mode::Sample => {
            writeln!(
                w,
                "sample\tPCS_raw\tUTP_raw\tCLS_raw\tERAD_raw\tRibo_raw\tCapacity_raw\tPII_raw\tPFS_raw\tchaperone_core\tproteasome_core\tupr_core\tagg_core\tCCI\tPCI\tUPR_A\tPLS\tSCI\tPCP\tchaperone_high\tproteasome_high\tupr_active\tproteotoxic_high\timbalance_high\tcollapse_risk"
            )?;
            let pcs = mean(&axis.pcs)?;
            let utp = mean(&axis.utp)?;
            let cls = mean(&axis.cls)?;
            let erad = mean(&axis.erad)?;
            let ribo = mean(&axis.ribo)?;
            let capacity = mean(&integrated.capacity_raw)?;
            let pii = mean(&integrated.pii_raw)?;
            let pfs = mean(&integrated.pfs_raw)?;
            let chaperone_core = mean_or_nan(extension.map(|e| e.chaperone_core.as_slice()))?;
            let proteasome_core = mean_or_nan(extension.map(|e| e.proteasome_core.as_slice()))?;
            let upr_core = mean_or_nan(extension.map(|e| e.upr_core.as_slice()))?;
            let agg_core = mean_or_nan(extension.map(|e| e.agg_core.as_slice()))?;
            let cci = mean_or_nan(extension.map(|e| e.cci.as_slice()))?;
            let pci = mean_or_nan(extension.map(|e| e.pci.as_slice()))?;
            let upr_a = mean_or_nan(extension.map(|e| e.upr_a.as_slice()))?;
            let pls = mean_or_nan(extension.map(|e| e.pls.as_slice()))?;
            let sci = mean_or_nan(extension.map(|e| e.sci.as_slice()))?;
            let pcp = mean_or_nan(extension.map(|e| e.pcp.as_slice()))?;
            let chaperone_high = frac_true(extension.map(|e| e.chaperone_high.as_slice())) >= 0.5;
            let proteasome_high = frac_true(extension.map(|e| e.proteasome_high.as_slice())) >= 0.5;
            let upr_active = frac_true(extension.map(|e| e.upr_active.as_slice())) >= 0.5;
            let proteotoxic_high =
                frac_true(extension.map(|e| e.proteotoxic_high.as_slice())) >= 0.5;
            let imbalance_high = frac_true(extension.map(|e| e.imbalance_high.as_slice())) >= 0.5;
            let collapse_risk = frac_true(extension.map(|e| e.collapse_risk.as_slice())) >= 0.5;
            writeln!(
                w,
                "sample\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}",
                pcs,
                utp,
                cls,
                erad,
                ribo,
                capacity,
                pii,
                pfs,
                chaperone_core,
                proteasome_core,
                upr_core,
                agg_core,
                cci,
                pci,
                upr_a,
                pls,
                sci,
                pcp,
                chaperone_high,
                proteasome_high,
                upr_active,
                proteotoxic_high,
                imbalance_high,
                collapse_risk
            )?;
        }
    }

    Ok(())
}

fn ensure_len(got: usize, expected: usize, name: &str) -> Result<()> {
    if got != expected {
        bail!("{} length mismatch: {} != {}", name, got, expected);
    }
    Ok(())
}

fn mean(values: &[f32]) -> Result<f32> {
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

fn mean_or_nan(values: Option<&[f32]>) -> Result<f32> {
    match values {
        Some(v) if !v.is_empty() => mean(v),
        _ => Ok(f32::NAN),
    }
}

fn frac_true(values: Option<&[bool]>) -> f32 {
    let Some(v) = values else {
        return 0.0;
    };
    if v.is_empty() {
        return 0.0;
    }
    let mut n = 0usize;
    for &flag in v {
        if flag {
            n += 1;
        }
    }
    n as f32 / v.len() as f32
}
