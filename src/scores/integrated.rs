use anyhow::{Result, bail};

use crate::math::stats::{mad, median, robust_z};
use crate::schema::v1::Mode;
use crate::scores::{AxisRawScores, IntegratedScores, PfsContributions};

pub fn compute_integrated(
    axis: &AxisRawScores,
    mode: Mode,
) -> Result<(IntegratedScores, PfsContributions)> {
    let n = axis.pcs.len();
    if axis.utp.len() != n || axis.cls.len() != n || axis.erad.len() != n || axis.ribo.len() != n {
        bail!("axis raw vectors length mismatch");
    }

    let mut capacity_raw = vec![0.0f32; n];
    let mut pii_raw = vec![0.0f32; n];
    let mut pfs_raw = vec![0.0f32; n];

    let mut c_pii = vec![0.0f32; n];
    let mut c_utp = vec![0.0f32; n];
    let mut c_ribo = vec![0.0f32; n];
    let mut c_pcs = vec![0.0f32; n];

    for i in 0..n {
        let pcs = axis.pcs[i];
        let cls = axis.cls[i];
        let erad = axis.erad[i];
        let ribo = axis.ribo[i];
        let utp = axis.utp[i];

        if pcs.is_nan() || cls.is_nan() || erad.is_nan() || ribo.is_nan() || utp.is_nan() {
            bail!("NaN encountered in axis raw inputs");
        }

        let capacity = 0.55 * pcs + 0.25 * cls + 0.20 * erad;
        let pii = ribo - capacity;
        let pfs = 0.40 * pii + 0.25 * utp + 0.20 * ribo - 0.15 * pcs;

        capacity_raw[i] = capacity;
        pii_raw[i] = pii;
        pfs_raw[i] = pfs;

        c_pii[i] = 0.40 * pii;
        c_utp[i] = 0.25 * utp;
        c_ribo[i] = 0.20 * ribo;
        c_pcs[i] = -0.15 * pcs;
    }

    let (capacity_z, pii_z, pfs_z) = if matches!(mode, Mode::Cell) {
        (
            Some(zscore(&capacity_raw)?),
            Some(zscore(&pii_raw)?),
            Some(zscore(&pfs_raw)?),
        )
    } else {
        (None, None, None)
    };

    Ok((
        IntegratedScores {
            capacity_raw,
            pii_raw,
            pfs_raw,
            capacity_z,
            pii_z,
            pfs_z,
        },
        PfsContributions {
            c_pii,
            c_utp,
            c_ribo,
            c_pcs,
        },
    ))
}

fn zscore(values: &[f32]) -> Result<Vec<f32>> {
    let mut scratch = values.to_vec();
    let med = median(&mut scratch);
    let mut scratch2 = values.to_vec();
    let mad_val = mad(&mut scratch2, med);

    let mut out = Vec::with_capacity(values.len());
    for &v in values {
        if v.is_nan() {
            bail!("NaN encountered in integrated raw scores");
        }
        let z = robust_z(v, med, mad_val);
        out.push(z);
    }
    Ok(out)
}
