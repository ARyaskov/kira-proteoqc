use std::path::Path;

use anyhow::{Result, bail};

#[derive(Debug, Clone, Copy)]
pub struct SharedCacheHeaderV1 {
    pub version_major: u16,
    pub version_minor: u16,
    pub n_genes: u64,
    pub n_cells: u64,
    pub nnz: u64,
    pub genes_table_offset: u64,
    pub genes_table_bytes: u64,
    pub barcodes_table_offset: u64,
    pub barcodes_table_bytes: u64,
    pub col_ptr_offset: u64,
    pub row_idx_offset: u64,
    pub values_u32_offset: u64,
    pub n_blocks: u64,
    pub blocks_offset: u64,
    pub file_bytes: u64,
    pub header_crc64: u64,
    pub data_crc64: u64,
}

#[derive(Debug)]
pub struct SharedCache {
    pub header: SharedCacheHeaderV1,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
    col_ptr: Vec<u64>,
    row_idx: Vec<u32>,
    values_u32: Vec<u32>,
}

impl SharedCache {
    pub fn open(path: &Path) -> Result<Self> {
        let owned = kira_shared_sc_cache::read_shared_cache_owned(path)
            .map_err(|e| anyhow::anyhow!(e.to_string()))?;

        let header = SharedCacheHeaderV1 {
            version_major: 1,
            version_minor: 0,
            n_genes: owned.n_genes,
            n_cells: owned.n_cells,
            nnz: owned.nnz,
            genes_table_offset: 0,
            genes_table_bytes: 0,
            barcodes_table_offset: 0,
            barcodes_table_bytes: 0,
            col_ptr_offset: 0,
            row_idx_offset: 0,
            values_u32_offset: 0,
            n_blocks: 0,
            blocks_offset: 0,
            file_bytes: 0,
            header_crc64: 0,
            data_crc64: 0,
        };

        Ok(Self {
            header,
            genes: owned.genes,
            barcodes: owned.barcodes,
            col_ptr: owned.col_ptr,
            row_idx: owned.row_idx,
            values_u32: owned.values_u32,
        })
    }

    pub fn col_ptr(&self) -> Result<Vec<u64>> {
        Ok(self.col_ptr.clone())
    }

    pub fn row_idx_at(&self, idx: usize) -> Result<u32> {
        let value = self
            .row_idx
            .get(idx)
            .copied()
            .ok_or_else(|| anyhow::anyhow!("row_idx index out of bounds"))?;
        Ok(value)
    }

    pub fn value_u32_at(&self, idx: usize) -> Result<u32> {
        let value = self
            .values_u32
            .get(idx)
            .copied()
            .ok_or_else(|| anyhow::anyhow!("values_u32 index out of bounds"))?;
        Ok(value)
    }
}

pub fn write_shared_cache(
    path: &Path,
    genes: &[String],
    barcodes: &[String],
    col_ptr: &[u64],
    row_idx: &[u32],
    values_u32: &[u32],
) -> Result<()> {
    if row_idx.len() != values_u32.len() {
        bail!("row_idx and values_u32 length mismatch");
    }
    let input = kira_shared_sc_cache::SharedCacheWriteInput {
        genes,
        barcodes,
        col_ptr,
        row_idx,
        values_u32,
    };
    kira_shared_sc_cache::write_shared_cache(path, &input)
        .map_err(|e| anyhow::anyhow!(e.to_string()))
}
