use anyhow::{Result, bail};

use crate::scores::{DeltaSummary, TimecourseResult, TimepointSummary};

pub fn compute_timecourse(timepoints: Vec<TimepointSummary>) -> Result<TimecourseResult> {
    if timepoints.len() < 2 {
        bail!("timecourse requires at least 2 timepoints");
    }

    for tp in &timepoints {
        if tp.pfs.is_nan()
            || tp.pii.is_nan()
            || tp.pcs.is_nan()
            || tp.cls.is_nan()
            || tp.utp.is_nan()
        {
            bail!("NaN encountered in timepoint summary");
        }
    }

    let mut deltas = Vec::new();
    for win in timepoints.windows(2) {
        let a = &win[0];
        let b = &win[1];
        deltas.push(DeltaSummary {
            from: a.label.clone(),
            to: b.label.clone(),
            delta_pfs: b.pfs - a.pfs,
            delta_pii: b.pii - a.pii,
            delta_pcs: b.pcs - a.pcs,
            delta_cls: b.cls - a.cls,
        });
    }

    let first = &timepoints[0];
    let last = &timepoints[timepoints.len() - 1];

    let adaptive = last.pfs < first.pfs;
    let collapse = last.pfs > first.pfs && last.pcs < first.pcs;
    let addicted =
        last.pcs > first.pcs && last.utp > first.utp && (last.pfs - first.pfs).abs() < 0.5;

    let trajectory = if collapse {
        "collapse_trajectory"
    } else if addicted {
        "addicted_survival"
    } else if adaptive {
        "adaptive_recovery"
    } else {
        "undefined"
    };

    Ok(TimecourseResult {
        timepoints,
        deltas,
        trajectory: trajectory.to_string(),
    })
}
