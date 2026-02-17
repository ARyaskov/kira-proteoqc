use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::io::open_maybe_gz;

#[derive(Debug, Clone, Copy)]
pub struct MtxSummary {
    pub nrows: usize,
    pub ncols: usize,
    pub nnz: usize,
}

pub fn read_mtx_summary(path: &Path) -> Result<MtxSummary> {
    let summary = for_each_entry(path, |_, _, _| Ok(()))?;
    Ok(summary)
}

pub fn for_each_entry<F>(path: &Path, mut f: F) -> Result<MtxSummary>
where
    F: FnMut(usize, usize, f32) -> Result<()>,
{
    let reader = open_maybe_gz(path)?;
    let mut reader = BufReader::new(reader);

    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            line.clear();
            continue;
        }
        if trimmed.starts_with('%') {
            line.clear();
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() != 3 {
            bail!("MTX header must have 3 fields: nrows ncols nnz");
        }
        let nrows: usize = parts[0].parse().context("invalid MTX nrows")?;
        let ncols: usize = parts[1].parse().context("invalid MTX ncols")?;
        let expected_nnz: usize = parts[2].parse().context("invalid MTX nnz")?;

        let mut nnz = 0usize;
        line.clear();
        while reader.read_line(&mut line)? > 0 {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                line.clear();
                continue;
            }
            if trimmed.starts_with('%') {
                line.clear();
                continue;
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() < 3 {
                bail!("MTX entry must have at least 3 fields: row col value");
            }
            // MTX indices are 1-based; convert to 0-based for bounds checks.
            let row: usize = parts[0].parse().context("invalid MTX row index")?;
            let col: usize = parts[1].parse().context("invalid MTX col index")?;
            if row == 0 || col == 0 {
                bail!("MTX indices are 1-based; found 0");
            }
            let row0 = row - 1;
            let col0 = col - 1;
            if row0 >= nrows || col0 >= ncols {
                bail!("MTX index out of bounds: ({}, {})", row, col);
            }
            let value: f32 = parts[2].parse().context("invalid MTX value")?;
            f(row0, col0, value)?;
            nnz += 1;
            line.clear();
        }

        if nnz != expected_nnz {
            bail!(
                "MTX nnz mismatch: header {} vs observed {}",
                expected_nnz,
                nnz
            );
        }

        return Ok(MtxSummary { nrows, ncols, nnz });
    }

    bail!("MTX file missing header")
}
