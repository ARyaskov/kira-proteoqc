use anyhow::Result;
use std::fs;
use tracing::info;

use crate::ctx::Ctx;
use crate::pipeline::Stage;

pub struct Stage0Scaffold;

impl Stage0Scaffold {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage0Scaffold {
    fn name(&self) -> &'static str {
        "stage0_scaffold"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        fs::create_dir_all(&ctx.output.out_dir)?;
        info!(
            out_dir = %ctx.output.out_dir.display(),
            "output_dir_ready"
        );

        ctx.report.input_meta.genes = ctx.input_meta.genes;
        ctx.report.input_meta.cells = ctx.input_meta.cells;
        ctx.report.input_meta.nnz = ctx.input_meta.nnz;

        if ctx.write_tsv {
            ctx.report.scores.per_cell_tsv_path = Some(
                ctx.output
                    .tsv_path
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .to_string(),
            );
        }

        Ok(())
    }
}
