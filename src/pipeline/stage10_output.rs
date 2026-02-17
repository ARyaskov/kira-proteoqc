use anyhow::Result;
use tracing::info;

use crate::ctx::{Ctx, RunMode};
use crate::io::{json_writer, pipeline_output, tsv_writer};
use crate::pipeline::Stage;

pub struct Stage10Output;

impl Stage10Output {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage10Output {
    fn name(&self) -> &'static str {
        "stage10_output"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        if matches!(ctx.run_mode, RunMode::Pipeline) {
            let out_dir = pipeline_output::ensure_pipeline_out_dir(&ctx.output.out_dir)?;
            pipeline_output::write_pipeline_outputs(ctx, &out_dir)?;
            info!(out_dir = %out_dir.display(), "stage10_pipeline_output_ready");
            return Ok(());
        }

        let report = json_writer::build_report(ctx)?;
        ctx.report = report;

        if ctx.write_json {
            json_writer::write_json(&ctx.output.json_path, ctx)?;
        }
        if ctx.write_tsv {
            tsv_writer::write_tsv(&ctx.output.tsv_path, ctx)?;
        }

        info!("stage10_output_ready");
        Ok(())
    }
}
