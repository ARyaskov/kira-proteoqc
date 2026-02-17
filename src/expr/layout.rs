use std::io::{Read, Write};

use anyhow::{Context, Result, bail};

pub const MAGIC: [u8; 8] = *b"KIRAEXPR";
pub const VERSION: u32 = 1;
pub const LAYOUT_CSC: u32 = 1;
pub const HEADER_SIZE: usize = 32;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExprHeaderV1 {
    pub version: u32,
    pub n_genes: u32,
    pub n_cells: u32,
    pub nnz: u64,
    pub layout: u32,
}

impl ExprHeaderV1 {
    pub fn expected_len(&self) -> usize {
        let gene_ptr_bytes = (self.n_genes as usize + 1) * 8;
        let cell_idx_bytes = self.nnz as usize * 4;
        let values_bytes = self.nnz as usize * 4;
        HEADER_SIZE + gene_ptr_bytes + cell_idx_bytes + values_bytes
    }

    pub fn gene_ptr_offset(&self) -> usize {
        HEADER_SIZE
    }

    pub fn cell_idx_offset(&self) -> usize {
        self.gene_ptr_offset() + (self.n_genes as usize + 1) * 8
    }

    pub fn values_offset(&self) -> usize {
        self.cell_idx_offset() + self.nnz as usize * 4
    }
}

pub fn write_header<W: Write>(mut w: W, header: &ExprHeaderV1) -> Result<()> {
    w.write_all(&MAGIC)?;
    w.write_all(&header.version.to_le_bytes())?;
    w.write_all(&header.n_genes.to_le_bytes())?;
    w.write_all(&header.n_cells.to_le_bytes())?;
    w.write_all(&header.nnz.to_le_bytes())?;
    w.write_all(&header.layout.to_le_bytes())?;
    Ok(())
}

pub fn read_header<R: Read>(mut r: R) -> Result<ExprHeaderV1> {
    let mut magic = [0u8; 8];
    r.read_exact(&mut magic)?;
    if magic != MAGIC {
        bail!("expr.bin magic mismatch");
    }
    let version = read_u32(&mut r)?;
    if version != VERSION {
        bail!("unsupported expr.bin version {}", version);
    }
    let n_genes = read_u32(&mut r)?;
    let n_cells = read_u32(&mut r)?;
    let nnz = read_u64(&mut r)?;
    let layout = read_u32(&mut r)?;
    if layout != LAYOUT_CSC {
        bail!("unsupported expr.bin layout {}", layout);
    }
    Ok(ExprHeaderV1 {
        version,
        n_genes,
        n_cells,
        nnz,
        layout,
    })
}

fn read_u32<R: Read>(mut r: R) -> Result<u32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf).context("failed to read u32")?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64<R: Read>(mut r: R) -> Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf).context("failed to read u64")?;
    Ok(u64::from_le_bytes(buf))
}
