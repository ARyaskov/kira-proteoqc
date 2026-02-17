use anyhow::{Result, bail};

use crate::ctx::Ctx;
use crate::schema::v1::Mode;

pub fn format_summary(ctx: &Ctx) -> Result<String> {
    let version = env!("CARGO_PKG_VERSION");
    let genes = ctx.input_meta.genes.unwrap_or(0);
    let cells = ctx.input_meta.cells.unwrap_or(0);
    let mode = match ctx.mode {
        Mode::Cell => "cell",
        Mode::Sample => "sample",
    };
    let integrated = ctx
        .integrated_scores
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("integrated scores missing"))?;
    let pfs = mean_vec(&integrated.pfs_raw)?;

    let mut out = String::new();
    out.push_str(&format!("kira-proteoqc v{}\n", version));
    out.push_str(&format!(
        "Input: {} genes, {} cells, mode={}\n",
        genes, cells, mode
    ));
    out.push_str(&format!("PFS: {:+.2}\n", pfs));

    if let Some(tc) = &ctx.timecourse_result {
        out.push_str(&format!("Trajectory: {}\n", tc.trajectory));
    }

    let fired: Vec<String> = ctx
        .risk_flags
        .iter()
        .filter(|f| f.fired)
        .map(|f| f.name.clone())
        .collect();
    if fired.is_empty() {
        out.push_str("Flags: none\n");
    } else {
        out.push_str(&format!("Flags: {}\n", fired.join(", ")));
    }

    Ok(out)
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
