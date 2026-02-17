use std::collections::HashMap;

use anyhow::{Context, Result, bail};
use tracing::{info, warn};

use crate::ctx::{Ctx, RunMode};
use crate::input;
use crate::io::{barcodes, features, mtx, shared_cache};
use crate::pipeline::Stage;

pub struct Stage1Input;

impl Stage1Input {
    pub fn new() -> Self {
        Self
    }
}

impl Stage for Stage1Input {
    fn name(&self) -> &'static str {
        "stage1_input"
    }

    fn run(&self, ctx: &mut Ctx) -> Result<()> {
        if ctx.input_format == crate::ctx::InputFormat::H5ad {
            return Ok(());
        }
        let (genes, cells, nrows, ncols, nnz) = match ctx.run_mode {
            RunMode::Standalone => load_from_mtx(ctx, None)?,
            RunMode::Pipeline => load_pipeline_input(ctx)?,
        };

        if nrows != genes.len() {
            bail!(
                "MTX rows ({}) do not match features.tsv lines ({})",
                nrows,
                genes.len()
            );
        }
        if ncols != cells.len() {
            bail!(
                "MTX cols ({}) do not match barcodes.tsv lines ({})",
                ncols,
                cells.len()
            );
        }

        // Gene index is first-win for duplicates to preserve deterministic ordering.
        let (gene_index, warnings) = build_gene_index(&genes);

        ctx.genes = genes;
        ctx.cells = cells;
        ctx.nnz = nnz;
        ctx.gene_index = gene_index;
        ctx.warnings.extend(warnings);

        ctx.input_meta.genes = Some(nrows as u64);
        ctx.input_meta.cells = Some(ncols as u64);
        ctx.input_meta.nnz = Some(nnz as u64);

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

fn load_pipeline_input(ctx: &mut Ctx) -> Result<(Vec<String>, Vec<String>, usize, usize, usize)> {
    if let Some(cache_override) = &ctx.cache_override {
        let cache = shared_cache::SharedCache::open(cache_override).with_context(|| {
            format!(
                "failed to validate shared cache {}",
                cache_override.display()
            )
        })?;
        ctx.shared_cache_path = Some(cache_override.clone());
        ctx.shared_cache_used = true;
        info!(cache = %cache_override.display(), "shared_cache_loaded_override");
        return Ok((
            cache.genes,
            cache.barcodes,
            cache.header.n_genes as usize,
            cache.header.n_cells as usize,
            cache.header.nnz as usize,
        ));
    }

    let prefix = input::detect_prefix(&ctx.input)?;
    ctx.input_prefix = prefix.clone();
    let expected_cache = input::resolve_shared_cache_path(&ctx.input, prefix.as_deref());

    if expected_cache.exists() {
        let cache = shared_cache::SharedCache::open(&expected_cache).with_context(|| {
            format!(
                "failed to validate shared cache {}",
                expected_cache.display()
            )
        })?;
        ctx.shared_cache_path = Some(expected_cache.clone());
        ctx.shared_cache_used = true;
        info!(cache = %expected_cache.display(), "shared_cache_loaded");
        return Ok((
            cache.genes,
            cache.barcodes,
            cache.header.n_genes as usize,
            cache.header.n_cells as usize,
            cache.header.nnz as usize,
        ));
    }

    warn!(
        expected_cache = %expected_cache.display(),
        "shared cache not found; falling back to MTX input"
    );
    ctx.warnings.push(format!(
        "shared cache not found at {}; using MTX fallback",
        expected_cache.display()
    ));

    load_from_mtx(ctx, prefix.as_deref())
}

fn load_from_mtx(
    ctx: &mut Ctx,
    prefix: Option<&str>,
) -> Result<(Vec<String>, Vec<String>, usize, usize, usize)> {
    ctx.shared_cache_used = false;
    ctx.shared_cache_path = None;
    let input_dir = &ctx.input;
    let (matrix_path, features_path, barcodes_path) =
        input::resolve_mtx_input_files(input_dir, prefix);
    let matrix_path = matrix_path.context("missing matrix.mtx or matrix.mtx.gz")?;
    let features_path = features_path.context("missing features.tsv/genes.tsv (or .gz)")?;
    let barcodes_path = barcodes_path.context("missing barcodes.tsv or barcodes.tsv.gz")?;

    info!(
        matrix = %matrix_path.display(),
        features = %features_path.display(),
        barcodes = %barcodes_path.display(),
        "input_files"
    );

    ctx.mtx_matrix_path = Some(matrix_path.clone());
    ctx.mtx_features_path = Some(features_path.clone());
    ctx.mtx_barcodes_path = Some(barcodes_path.clone());

    let genes = features::read_features(&features_path)?;
    let cells = barcodes::read_barcodes(&barcodes_path)?;
    let mtx_summary = mtx::read_mtx_summary(&matrix_path)?;
    Ok((
        genes,
        cells,
        mtx_summary.nrows,
        mtx_summary.ncols,
        mtx_summary.nnz,
    ))
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
