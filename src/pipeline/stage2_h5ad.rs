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
    let mut duplicate_order: Vec<String> = Vec::new();
    let mut duplicate_seen: HashMap<String, usize> = HashMap::new();

    for (i, symbol) in genes.iter().enumerate() {
        if index.get(symbol).is_some() {
            let entry = duplicate_seen.entry(symbol.clone()).or_insert(0usize);
            if *entry == 0 {
                duplicate_order.push(symbol.clone());
            }
            *entry += 1;
        } else {
            index.insert(symbol.clone(), i);
        }
    }

    if !duplicate_order.is_empty() {
        let mut parts = Vec::with_capacity(duplicate_order.len());
        for symbol in duplicate_order {
            let first_row = index
                .get(&symbol)
                .copied()
                .unwrap_or_default()
                .saturating_add(1);
            let dup_count = duplicate_seen.get(&symbol).copied().unwrap_or(0);
            let total_count = dup_count + 1;
            parts.push(format!(
                "{}(x{}, first_row={})",
                symbol, total_count, first_row
            ));
        }
        warnings.push(format!(
            "duplicate gene symbols detected (kept first occurrence): {}",
            parts.join(", ")
        ));
    }

    (index, warnings)
}
