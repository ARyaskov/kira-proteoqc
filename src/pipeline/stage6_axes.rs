use anyhow::Result;
use tracing::info;

use crate::ctx::Ctx;
use crate::pipeline::Stage;
use crate::scores::axis_raw::compute_axis_raw;

pub struct Stage6Axes;

impl Stage6Axes {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage6Axes {
    fn name(&self) -> &'static str {
        "stage6_axes"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        let axis_raw = compute_axis_raw(ctx)?;
        ctx.axis_raw = Some(axis_raw);
        info!("axis_raw_ready");
        Ok(())
    }
}
