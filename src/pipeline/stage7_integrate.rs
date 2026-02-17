use anyhow::Result;
use tracing::info;

use crate::ctx::Ctx;
use crate::pipeline::Stage;
use crate::scores::integrated::compute_integrated;

pub struct Stage7Integrate;

impl Stage7Integrate {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage7Integrate {
    fn name(&self) -> &'static str {
        "stage7_integrate"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        let axis = ctx
            .axis_raw
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("axis raw scores missing"))?;
        let (integrated, contrib) = compute_integrated(axis, ctx.mode.clone())?;
        ctx.integrated_scores = Some(integrated);
        ctx.pfs_contributions = Some(contrib);
        info!("integrated_scores_ready");
        Ok(())
    }
}
