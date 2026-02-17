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

            writeln!(
                w,
                "cell_id\tPCS_raw\tUTP_raw\tCLS_raw\tERAD_raw\tRibo_raw\tCapacity_raw\tPII_raw\tPFS_raw\tPFS_z"
            )?;
            for i in 0..n {
                writeln!(
                    w,
                    "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                    ctx.cells[i],
                    axis.pcs[i],
                    axis.utp[i],
                    axis.cls[i],
                    axis.erad[i],
                    axis.ribo[i],
                    integrated.capacity_raw[i],
                    integrated.pii_raw[i],
                    integrated.pfs_raw[i],
                    pfs_z[i]
                )?;
            }
        }
        Mode::Sample => {
            writeln!(
                w,
                "sample\tPCS_raw\tUTP_raw\tCLS_raw\tERAD_raw\tRibo_raw\tCapacity_raw\tPII_raw\tPFS_raw"
            )?;
            let pcs = mean(&axis.pcs)?;
            let utp = mean(&axis.utp)?;
            let cls = mean(&axis.cls)?;
            let erad = mean(&axis.erad)?;
            let ribo = mean(&axis.ribo)?;
            let capacity = mean(&integrated.capacity_raw)?;
            let pii = mean(&integrated.pii_raw)?;
            let pfs = mean(&integrated.pfs_raw)?;
            writeln!(
                w,
                "sample\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                pcs, utp, cls, erad, ribo, capacity, pii, pfs
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
