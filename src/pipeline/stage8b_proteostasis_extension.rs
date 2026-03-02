use anyhow::Result;
use tracing::info;

use crate::ctx::Ctx;
use crate::metrics::proteostasis_extension::compute_extension;
use crate::pipeline::Stage;

pub struct Stage8bProteostasisExtension;

impl Stage8bProteostasisExtension {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage8bProteostasisExtension {
    fn name(&self) -> &'static str {
        "stage8b_proteostasis_extension"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        let result = compute_extension(ctx)?;
        ctx.proteostasis_extension = Some(result);
        info!("proteostasis_extension_ready");
        Ok(())
    }
}
