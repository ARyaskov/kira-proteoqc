use anyhow::Result;
use tracing::info;

use crate::ctx::Ctx;
use crate::pipeline::Stage;
use crate::scores::timecourse::compute_timecourse;

pub struct Stage9Timecourse;

impl Stage9Timecourse {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage9Timecourse {
    fn name(&self) -> &'static str {
        "stage9_timecourse"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        if !ctx.timecourse {
            return Ok(());
        }
        if ctx.timecourse_points.len() < 2 {
            anyhow::bail!("timecourse requires at least 2 timepoints");
        }
        let result = compute_timecourse(ctx.timecourse_points.clone())?;
        ctx.timecourse_result = Some(result);
        info!("timecourse_ready");
        Ok(())
    }
}
