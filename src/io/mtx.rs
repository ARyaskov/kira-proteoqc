use std::path::Path;

use anyhow::{Result, bail};
use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

#[derive(Debug, Clone, Copy)]
pub struct MtxSummary {
    pub nrows: usize,
    pub ncols: usize,
    pub nnz: usize,
}

pub fn read_mtx_summary(path: &Path) -> Result<MtxSummary> {
    let matrix = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_matrix()
    .map_err(|e| anyhow::anyhow!(e.message))?;

    Ok(MtxSummary {
        nrows: matrix.n_genes,
        ncols: matrix.n_cells,
        nnz: matrix.values.len(),
    })
}

pub fn for_each_entry<F>(path: &Path, mut f: F) -> Result<MtxSummary>
where
    F: FnMut(usize, usize, f32) -> Result<()>,
{
    let matrix = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_matrix()
    .map_err(|e| anyhow::anyhow!(e.message))?;

    for (col, window) in matrix.col_ptr.windows(2).enumerate() {
        for idx in window[0]..window[1] {
            let row = matrix.row_idx[idx];
            let value = matrix.values[idx];
            if row >= matrix.n_genes || col >= matrix.n_cells {
                bail!("MTX index out of bounds: ({}, {})", row + 1, col + 1);
            }
            f(row, col, value)?;
        }
    }

    Ok(MtxSummary {
        nrows: matrix.n_genes,
        ncols: matrix.n_cells,
        nnz: matrix.values.len(),
    })
}
