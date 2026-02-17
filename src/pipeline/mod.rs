use anyhow::Result;
use std::time::Instant;
use tracing::{info, warn};

use crate::ctx::Ctx;

pub mod stage0_scaffold;
pub mod stage10_output;
pub mod stage1_input;
pub mod stage2_h5ad;
pub mod stage3_expr_cache;
pub mod stage4_geneset;
pub mod stage5_math;
pub mod stage6_axes;
pub mod stage7_integrate;
pub mod stage8_risk;
pub mod stage9_timecourse;

pub trait Stage {
    fn name(&self) -> &'static str;
    fn run(&self, ctx: &mut Ctx) -> Result<()>;
}

pub struct Pipeline {
    stages: Vec<Box<dyn Stage>>,
}

impl Pipeline {
    pub fn new(stages: Vec<Box<dyn Stage>>) -> Self {
        Self { stages }
    }

    pub fn run(&self, ctx: &mut Ctx) -> Result<()> {
        info!(
            simd_backend = %crate::simd::backend_name(),
            "compute backend selected"
        );
        for stage in &self.stages {
            let start = Instant::now();
            info!(stage = stage.name(), "stage started");
            if let Err(err) = stage.run(ctx) {
                let elapsed_ms = start.elapsed().as_millis();
                warn!(
                    stage = stage.name(),
                    elapsed_ms = elapsed_ms as u64,
                    "stage failed"
                );
                return Err(err);
            }
            let elapsed_ms = start.elapsed().as_millis();
            info!(
                stage = stage.name(),
                elapsed_ms = elapsed_ms as u64,
                "stage finished"
            );
        }
        Ok(())
    }
}
