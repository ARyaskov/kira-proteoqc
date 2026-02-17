use std::collections::HashMap;

use anyhow::Result;
use tracing::info;

use crate::ctx::{Ctx, InputFormat};
use crate::io::h5ad;
use crate::pipeline::Stage;

pub struct Stage2H5ad;

impl Stage2H5ad {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage2H5ad {
    fn name(&self) -> &'static str {
        "stage2_h5ad"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        if ctx.input_format != InputFormat::H5ad {
            return Ok(());
        }

        let summary = h5ad::read_h5ad_summary(&ctx.input)?;
        info!(
            nrows = summary.nrows,
            ncols = summary.ncols,
            nnz = summary.nnz,
            "h5ad_summary"
        );

        let (gene_index, mut warnings) = build_gene_index(&summary.genes);
        warnings.extend(summary.warnings.into_iter());

        ctx.genes = summary.genes;
        ctx.cells = summary.cells;
        ctx.nnz = summary.nnz;
        ctx.gene_index = gene_index;
        ctx.warnings = warnings;

        ctx.input_meta.genes = Some(summary.nrows as u64);
        ctx.input_meta.cells = Some(summary.ncols as u64);
        ctx.input_meta.nnz = Some(summary.nnz as u64);

        ctx.report.input_meta.genes = ctx.input_meta.genes;
        ctx.report.input_meta.cells = ctx.input_meta.cells;
        ctx.report.input_meta.nnz = ctx.input_meta.nnz;
        ctx.report.input_meta.normalization.log1p = ctx.log1p;
        ctx.report.input_meta.normalization.cp10k = true;
        ctx.report.input_meta.normalization.raw = true;
        ctx.report.input_meta.mode = ctx.mode.clone();
        ctx.report.input_meta.timecourse = ctx.timecourse;

        Ok(())
    }
}

fn build_gene_index(genes: &[String]) -> (HashMap<String, usize>, Vec<String>) {
    let mut index = HashMap::new();
    let mut warnings = Vec::new();

    for (i, symbol) in genes.iter().enumerate() {
        if let Some(first) = index.get(symbol) {
            warnings.push(format!(
                "duplicate gene symbol '{}' at row {} (kept first at row {})",
                symbol,
                i + 1,
                first + 1
            ));
        } else {
            index.insert(symbol.clone(), i);
        }
    }

    (index, warnings)
}
