use anyhow::{Result, bail};
use tracing::info;

use crate::ctx::Ctx;
use crate::pipeline::Stage;

pub struct Stage5Math;

impl Stage5Math {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage5Math {
    fn name(&self) -> &'static str {
        "stage5_math"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        let header = ctx
            .expr_header
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("expr header missing"))?;
        if header.n_cells as usize != ctx.cells.len() {
            bail!("expr header n_cells does not match ctx");
        }
        if header.n_genes as usize != ctx.genes.len() {
            bail!("expr header n_genes does not match ctx");
        }
        if header.nnz as usize != ctx.nnz {
            bail!("expr header nnz does not match ctx");
        }
        if ctx.genesets.is_none() {
            bail!("genesets not resolved before Stage 5");
        }

        ctx.scratch_cell_buf.resize(ctx.cells.len(), 0.0);

        info!("math_primitives_ready");
        Ok(())
    }
}
