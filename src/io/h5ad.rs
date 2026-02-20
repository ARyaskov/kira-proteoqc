use std::path::Path;

use anyhow::Result;
use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

#[derive(Debug)]
pub struct H5adSummary {
    pub genes: Vec<String>,
    pub cells: Vec<String>,
    pub nnz: usize,
    pub nrows: usize,
    pub ncols: usize,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Copy)]
pub struct H5adSparseMeta {
    pub nrows: usize,
    pub ncols: usize,
    pub nnz: usize,
}

pub fn read_h5ad_summary(path: &Path) -> Result<H5adSummary> {
    let reader = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::H5ad),
            strict: true,
        },
    );
    let canonical = reader.read_all().map_err(|e| anyhow::anyhow!(e.message))?;

    Ok(H5adSummary {
        genes: canonical.metadata.gene_symbols,
        cells: canonical.metadata.barcodes,
        nnz: canonical.matrix.values.len(),
        nrows: canonical.matrix.n_genes,
        ncols: canonical.matrix.n_cells,
        warnings: Vec::new(),
    })
}

pub fn read_h5ad_sparse(
    path: &Path,
) -> Result<(H5adSparseMeta, String, usize, Vec<u32>, Vec<f32>, Vec<u64>)> {
    let matrix = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::H5ad),
            strict: true,
        },
    )
    .read_matrix()
    .map_err(|e| anyhow::anyhow!(e.message))?;

    let indices = matrix
        .row_idx
        .iter()
        .map(|&x| u32::try_from(x).map_err(|_| anyhow::anyhow!("row index exceeds u32")))
        .collect::<Result<Vec<_>>>()?;
    let indptr = matrix
        .col_ptr
        .iter()
        .map(|&x| u64::try_from(x).map_err(|_| anyhow::anyhow!("indptr exceeds u64")))
        .collect::<Result<Vec<_>>>()?;

    let meta = H5adSparseMeta {
        nrows: matrix.n_genes,
        ncols: matrix.n_cells,
        nnz: matrix.values.len(),
    };

    Ok((
        meta,
        "csc_matrix".to_string(),
        indptr.len(),
        indices,
        matrix.values,
        indptr,
    ))
}
