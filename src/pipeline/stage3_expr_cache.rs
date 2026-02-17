use anyhow::Result;
use tracing::info;

use crate::ctx::Ctx;
use crate::expr::writer;
use crate::pipeline::Stage;

pub struct Stage3ExprCache;

impl Stage3ExprCache {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage3ExprCache {
    fn name(&self) -> &'static str {
        "stage3_expr_cache"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        writer::ensure_expr_cache(ctx)?;
        info!(expr = %ctx.expr_path.display(), "expr_cache_ready");
        Ok(())
    }
}
