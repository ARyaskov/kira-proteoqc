use anyhow::{Result, bail};

use crate::ctx::Ctx;
use crate::math::stats::{mad, median, robust_z};
use crate::schema::v1::Mode;
use crate::scores::axis_raw::compute_axis_raw_with_mode;
use crate::scores::{IntegratedScores, RiskFlag};

const FRAGILE_THRESHOLD: f32 = 1.5;

pub fn compute_risk_flags(ctx: &mut Ctx) -> Result<Vec<RiskFlag>> {
    let mut flags = Vec::new();
    match ctx.mode {
        Mode::Cell => {
            let axis = ctx
                .axis_raw
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("axis raw scores missing"))?;
            let integrated = ctx
                .integrated_scores
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("integrated scores missing"))?;

            let pfs_z = integrated
                .pfs_z
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("PFS_z missing in per-cell mode"))?;
            let pii_z = integrated
                .pii_z
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("PII_z missing in per-cell mode"))?;

            let pcs_z = zscore_vec(&axis.pcs)?;
            let utp_z = zscore_vec(&axis.utp)?;
            let cls_z = zscore_vec(&axis.cls)?;
            let erad_z = zscore_vec(&axis.erad)?;

            flags.push(flag_fragile_high_cell(pfs_z));
            flags.push(flag_proteasome_addiction_cell(&pcs_z, &utp_z, pfs_z));
            flags.push(flag_proteotoxic_stress_cell(&cls_z, pii_z));
            flags.push(flag_er_degradation_overdrive_cell(&erad_z, pii_z));
        }
        Mode::Sample => {
            let sample_pfs = ctx
                .integrated_scores
                .as_ref()
                .and_then(|s| s.pfs_raw.get(0).copied())
                .unwrap_or(0.0);
            flags.push(flag_fragile_high_sample(ctx, sample_pfs)?);

            let axis = ctx
                .axis_raw
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("axis raw scores missing"))?;
            let integrated = ctx
                .integrated_scores
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("integrated scores missing"))?;
            flags.push(flag_proteasome_addiction_sample(integrated, axis));
            flags.push(flag_proteotoxic_stress_sample(integrated, axis));
            flags.push(flag_er_degradation_overdrive_sample(integrated, axis));
        }
    }

    Ok(flags)
}

fn zscore_vec(values: &[f32]) -> Result<Vec<f32>> {
    let mut scratch = values.to_vec();
    let med = median(&mut scratch);
    let mut scratch2 = values.to_vec();
    let mad_val = mad(&mut scratch2, med);
    let mut out = Vec::with_capacity(values.len());
    for &v in values {
        if v.is_nan() {
            bail!("NaN encountered in axis raw values");
        }
        let z = robust_z(v, med, mad_val);
        out.push(z);
    }
    Ok(out)
}

fn flag_fragile_high_cell(pfs_z: &[f32]) -> RiskFlag {
    let mut count = 0usize;
    for &v in pfs_z {
        if v > FRAGILE_THRESHOLD {
            count += 1;
        }
    }
    let frac = fraction(count, pfs_z.len());
    RiskFlag {
        name: "fragile_high".to_string(),
        fired: count > 0,
        threshold: format!("PFS_z > {}", FRAGILE_THRESHOLD),
        details: Some(format!("fraction={:.4}", frac)),
    }
}

fn flag_fragile_high_sample(ctx: &mut Ctx, sample_pfs: f32) -> Result<RiskFlag> {
    // Compute per-cell PFS_raw on-demand for sample mode.
    let axis_cell = compute_axis_raw_with_mode(ctx, Mode::Cell)?;
    let (integrated_cell, _) =
        crate::scores::integrated::compute_integrated(&axis_cell, Mode::Cell)?;
    let pfs_raw = integrated_cell.pfs_raw;

    let mut count = 0usize;
    for &v in &pfs_raw {
        if v > FRAGILE_THRESHOLD {
            count += 1;
        }
    }
    let frac = fraction(count, pfs_raw.len());
    let fired = frac >= 0.10 && !pfs_raw.is_empty();
    Ok(RiskFlag {
        name: "fragile_high".to_string(),
        fired,
        threshold: format!("top10% PFS_raw > {}", FRAGILE_THRESHOLD),
        details: Some(format!(
            "fraction={:.4}, sample_pfs_raw={:.4}",
            frac, sample_pfs
        )),
    })
}

fn flag_proteasome_addiction_cell(pcs_z: &[f32], utp_z: &[f32], pfs_z: &[f32]) -> RiskFlag {
    let mut count = 0usize;
    for i in 0..pcs_z.len() {
        let pcs = pcs_z[i];
        let utp = utp_z[i];
        let pfs = pfs_z[i];
        if pcs > 1.0 && utp > 1.0 && pfs > -0.5 && pfs < 1.0 {
            count += 1;
        }
    }
    let frac = fraction(count, pcs_z.len());
    RiskFlag {
        name: "proteasome_addiction".to_string(),
        fired: count > 0,
        threshold: "PCS_z>1.0, UTP_z>1.0, PFS_z in (-0.5,1.0)".to_string(),
        details: Some(format!("fraction={:.4}", frac)),
    }
}

fn flag_proteasome_addiction_sample(
    integrated: &IntegratedScores,
    axis: &crate::scores::AxisRawScores,
) -> RiskFlag {
    let pcs = axis.pcs.get(0).copied().unwrap_or(0.0);
    let utp = axis.utp.get(0).copied().unwrap_or(0.0);
    let pfs = integrated.pfs_raw.get(0).copied().unwrap_or(0.0);
    let fired = pcs > 1.0 && utp > 1.0 && pfs > -0.5 && pfs < 1.0;
    RiskFlag {
        name: "proteasome_addiction".to_string(),
        fired,
        threshold: "PCS_raw>1.0, UTP_raw>1.0, PFS_raw in (-0.5,1.0)".to_string(),
        details: Some(format!("pcs={:.4}, utp={:.4}, pfs={:.4}", pcs, utp, pfs)),
    }
}

fn flag_proteotoxic_stress_cell(cls_z: &[f32], pii_z: &[f32]) -> RiskFlag {
    let mut count = 0usize;
    for i in 0..cls_z.len() {
        if cls_z[i] > 1.0 && pii_z[i] > 1.0 {
            count += 1;
        }
    }
    let frac = fraction(count, cls_z.len());
    RiskFlag {
        name: "proteotoxic_stress".to_string(),
        fired: count > 0,
        threshold: "CLS_z>1.0, PII_z>1.0".to_string(),
        details: Some(format!("fraction={:.4}", frac)),
    }
}

fn flag_proteotoxic_stress_sample(
    integrated: &IntegratedScores,
    axis: &crate::scores::AxisRawScores,
) -> RiskFlag {
    let cls = axis.cls.get(0).copied().unwrap_or(0.0);
    let pii = integrated.pii_raw.get(0).copied().unwrap_or(0.0);
    let fired = cls > 1.0 && pii > 1.0;
    RiskFlag {
        name: "proteotoxic_stress".to_string(),
        fired,
        threshold: "CLS_raw>1.0, PII_raw>1.0".to_string(),
        details: Some(format!("cls={:.4}, pii={:.4}", cls, pii)),
    }
}

fn flag_er_degradation_overdrive_cell(erad_z: &[f32], pii_z: &[f32]) -> RiskFlag {
    let mut count = 0usize;
    for i in 0..erad_z.len() {
        if erad_z[i] > 1.0 && pii_z[i] > 1.0 {
            count += 1;
        }
    }
    let frac = fraction(count, erad_z.len());
    RiskFlag {
        name: "er_degradation_overdrive".to_string(),
        fired: count > 0,
        threshold: "ERAD_z>1.0, PII_z>1.0".to_string(),
        details: Some(format!("fraction={:.4}", frac)),
    }
}

fn flag_er_degradation_overdrive_sample(
    integrated: &IntegratedScores,
    axis: &crate::scores::AxisRawScores,
) -> RiskFlag {
    let erad = axis.erad.get(0).copied().unwrap_or(0.0);
    let pii = integrated.pii_raw.get(0).copied().unwrap_or(0.0);
    let fired = erad > 1.0 && pii > 1.0;
    RiskFlag {
        name: "er_degradation_overdrive".to_string(),
        fired,
        threshold: "ERAD_raw>1.0, PII_raw>1.0".to_string(),
        details: Some(format!("erad={:.4}, pii={:.4}", erad, pii)),
    }
}

fn fraction(count: usize, total: usize) -> f32 {
    if total == 0 {
        0.0
    } else {
        count as f32 / total as f32
    }
}
