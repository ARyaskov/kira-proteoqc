use anyhow::Result;
use tracing::info;

use crate::ctx::Ctx;
use crate::pipeline::Stage;
use crate::scores::risk::compute_risk_flags;

pub struct Stage8Risk;

impl Stage8Risk {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage8Risk {
    fn name(&self) -> &'static str {
        "stage8_risk"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        let flags = compute_risk_flags(ctx)?;
        ctx.risk_flags = flags;
        info!("risk_flags_ready");
        Ok(())
    }
}
