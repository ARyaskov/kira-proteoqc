use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result, bail};
use memmap2::Mmap;

use crate::expr::layout::{ExprHeaderV1, HEADER_SIZE, read_header};

pub struct ExprReader<'a> {
    header: &'a ExprHeaderV1,
    mmap: &'a Mmap,
}

impl<'a> ExprReader<'a> {
    pub fn new(header: &'a ExprHeaderV1, mmap: &'a Mmap) -> Self {
        Self { header, mmap }
    }

    pub fn n_cells(&self) -> usize {
        self.header.n_cells as usize
    }

    pub fn n_genes(&self) -> usize {
        self.header.n_genes as usize
    }

    pub fn gene_slice(&self, gene_id: usize) -> Result<(&[u32], &[f32])> {
        if gene_id >= self.n_genes() {
            bail!("gene_id out of range");
        }
        let gene_ptr = gene_ptr_slice(self.mmap, self.header);
        let cell_idx = cell_idx_slice(self.mmap, self.header);
        let values = values_slice(self.mmap, self.header);
        let start = gene_ptr[gene_id] as usize;
        let end = gene_ptr[gene_id + 1] as usize;
        Ok((&cell_idx[start..end], &values[start..end]))
    }
}

pub fn open_mmap(path: &Path) -> Result<(ExprHeaderV1, Mmap)> {
    let file = File::open(path).context("failed to open expr.bin")?;
    let mmap = unsafe { Mmap::map(&file).context("failed to mmap expr.bin")? };
    if mmap.len() < HEADER_SIZE {
        bail!("expr.bin too small");
    }
    let header = read_header(&mmap[..HEADER_SIZE])?;
    let expected = header.expected_len();
    if mmap.len() != expected {
        bail!(
            "expr.bin size mismatch: expected {}, got {}",
            expected,
            mmap.len()
        );
    }
    Ok((header, mmap))
}

pub fn gene_ptr_slice<'a>(mmap: &'a Mmap, header: &ExprHeaderV1) -> &'a [u64] {
    let start = header.gene_ptr_offset();
    let end = header.cell_idx_offset();
    let bytes = &mmap[start..end];
    unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u64, bytes.len() / 8) }
}

pub fn cell_idx_slice<'a>(mmap: &'a Mmap, header: &ExprHeaderV1) -> &'a [u32] {
    let start = header.cell_idx_offset();
    let end = header.values_offset();
    let bytes = &mmap[start..end];
    unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u32, bytes.len() / 4) }
}

pub fn values_slice<'a>(mmap: &'a Mmap, header: &ExprHeaderV1) -> &'a [f32] {
    let start = header.values_offset();
    let bytes = &mmap[start..];
    unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const f32, bytes.len() / 4) }
}
