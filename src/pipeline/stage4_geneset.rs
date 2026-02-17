use anyhow::Result;
use tracing::info;

use crate::ctx::Ctx;
use crate::geneset::{load_builtin, load_user, merge_defs, resolve_collection};
use crate::pipeline::Stage;
use crate::schema::v1::GenesetCoverage;

pub struct Stage4Geneset;

impl Stage4Geneset {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage4Geneset {
    fn name(&self) -> &'static str {
        "stage4_geneset"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        let mut collection = load_builtin()?;
        if let Some(path) = &ctx.geneset_path {
            let user_defs = load_user(path)?;
            let merged = merge_defs(collection.defs, user_defs);
            collection.defs = merged;
        }

        let resolved = resolve_collection(collection, &ctx.gene_index);
        let mut coverage = Vec::with_capacity(resolved.resolved.len());
        let mut warnings = Vec::new();

        for gs in &resolved.resolved {
            let found = gs.gene_ids.len();
            let total = gs.total;
            let fraction = if total == 0 {
                0.0
            } else {
                found as f64 / total as f64
            };
            if found == 0 {
                warnings.push(format!("geneset '{}' resolved to 0 genes", gs.id));
            }
            coverage.push(GenesetCoverage {
                geneset: gs.id.clone(),
                found: found as u64,
                total: total as u64,
                fraction,
            });
        }

        ctx.warnings.extend(warnings);
        ctx.report.explainability.geneset_coverage = coverage.clone();
        ctx.genesets = Some(resolved);

        info!("geneset_resolved");
        Ok(())
    }
}
