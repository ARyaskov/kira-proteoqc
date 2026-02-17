use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result, bail};
use crc::{CRC_64_ECMA_182, Crc};
use memmap2::Mmap;

const HEADER_SIZE_V1: usize = 256;
const HEADER_CRC_OFFSET: usize = 120;
const MAGIC: &[u8; 4] = b"KORG";
const ENDIAN_TAG_LE: u32 = 0x1234_5678;

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
    mmap: Mmap,
    pub header: SharedCacheHeaderV1,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
}

impl SharedCache {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("failed to open shared cache at {}", path.display()))?;
        let mmap = unsafe { Mmap::map(&file).context("failed to mmap shared cache")? };
        validate_header(&mmap)?;
        let header = parse_header(&mmap[..HEADER_SIZE_V1])?;
        validate_crc(&mmap, &header)?;
        validate_section_bounds(&mmap, &header)?;

        let genes = parse_string_table(
            &mmap,
            header.genes_table_offset as usize,
            header.genes_table_bytes as usize,
            header.n_genes as usize,
            "genes",
        )?;
        let barcodes = parse_string_table(
            &mmap,
            header.barcodes_table_offset as usize,
            header.barcodes_table_bytes as usize,
            header.n_cells as usize,
            "barcodes",
        )?;

        validate_csc(&mmap, &header)?;

        Ok(Self {
            mmap,
            header,
            genes,
            barcodes,
        })
    }

    pub fn col_ptr(&self) -> Result<Vec<u64>> {
        let n = self.header.n_cells as usize + 1;
        let start = self.header.col_ptr_offset as usize;
        let end = start + n * 8;
        let bytes = &self.mmap[start..end];
        let mut out = Vec::with_capacity(n);
        for chunk in bytes.chunks_exact(8) {
            out.push(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
        Ok(out)
    }

    pub fn row_idx_at(&self, idx: usize) -> Result<u32> {
        if idx >= self.header.nnz as usize {
            bail!("row_idx index out of bounds");
        }
        let offset = self.header.row_idx_offset as usize + idx * 4;
        Ok(u32::from_le_bytes(
            self.mmap[offset..offset + 4].try_into().unwrap(),
        ))
    }

    pub fn value_u32_at(&self, idx: usize) -> Result<u32> {
        if idx >= self.header.nnz as usize {
            bail!("values_u32 index out of bounds");
        }
        let offset = self.header.values_u32_offset as usize + idx * 4;
        Ok(u32::from_le_bytes(
            self.mmap[offset..offset + 4].try_into().unwrap(),
        ))
    }
}

fn validate_header(mmap: &Mmap) -> Result<()> {
    if mmap.len() < HEADER_SIZE_V1 {
        bail!("shared cache file too small for v1 header");
    }
    Ok(())
}

fn parse_header(header: &[u8]) -> Result<SharedCacheHeaderV1> {
    if &header[0..4] != MAGIC {
        bail!("invalid shared cache magic");
    }
    let version_major = le_u16(header, 4)?;
    let version_minor = le_u16(header, 6)?;
    let endian_tag = le_u32(header, 8)?;
    let header_size = le_u32(header, 12)?;

    if version_major != 1 {
        bail!("unsupported shared cache major version: {}", version_major);
    }
    if version_minor != 0 {
        bail!("unsupported shared cache minor version: {}", version_minor);
    }
    if endian_tag != ENDIAN_TAG_LE {
        bail!("unsupported shared cache endian tag: 0x{endian_tag:08x}");
    }
    if header_size as usize != HEADER_SIZE_V1 {
        bail!("unexpected shared cache header size: {}", header_size);
    }

    Ok(SharedCacheHeaderV1 {
        version_major,
        version_minor,
        n_genes: le_u64(header, 16)?,
        n_cells: le_u64(header, 24)?,
        nnz: le_u64(header, 32)?,
        genes_table_offset: le_u64(header, 40)?,
        genes_table_bytes: le_u64(header, 48)?,
        barcodes_table_offset: le_u64(header, 56)?,
        barcodes_table_bytes: le_u64(header, 64)?,
        col_ptr_offset: le_u64(header, 72)?,
        row_idx_offset: le_u64(header, 80)?,
        values_u32_offset: le_u64(header, 88)?,
        n_blocks: le_u64(header, 96)?,
        blocks_offset: le_u64(header, 104)?,
        file_bytes: le_u64(header, 112)?,
        header_crc64: le_u64(header, 120)?,
        data_crc64: le_u64(header, 128)?,
    })
}

fn validate_crc(mmap: &Mmap, header: &SharedCacheHeaderV1) -> Result<()> {
    let mut tmp = mmap[..HEADER_SIZE_V1].to_vec();
    for b in &mut tmp[HEADER_CRC_OFFSET..HEADER_CRC_OFFSET + 8] {
        *b = 0;
    }
    let crc = Crc::<u64>::new(&CRC_64_ECMA_182);
    let computed = crc.checksum(&tmp);
    if computed != header.header_crc64 {
        bail!(
            "shared cache header CRC mismatch: expected {:016x}, got {:016x}",
            header.header_crc64,
            computed
        );
    }
    Ok(())
}

fn validate_section_bounds(mmap: &Mmap, header: &SharedCacheHeaderV1) -> Result<()> {
    if header.file_bytes as usize != mmap.len() {
        bail!(
            "shared cache file_bytes mismatch: header {} actual {}",
            header.file_bytes,
            mmap.len()
        );
    }
    if header.n_blocks != 0 || header.blocks_offset != 0 || header.data_crc64 != 0 {
        bail!("shared cache contains unsupported optional blocks");
    }

    let col_ptr_bytes = (header.n_cells as usize + 1) * 8;
    let row_idx_bytes = header.nnz as usize * 4;
    let values_bytes = header.nnz as usize * 4;

    check_range(
        mmap.len(),
        header.genes_table_offset as usize,
        header.genes_table_bytes as usize,
        "genes table",
    )?;
    check_range(
        mmap.len(),
        header.barcodes_table_offset as usize,
        header.barcodes_table_bytes as usize,
        "barcodes table",
    )?;
    check_range(
        mmap.len(),
        header.col_ptr_offset as usize,
        col_ptr_bytes,
        "col_ptr",
    )?;
    check_range(
        mmap.len(),
        header.row_idx_offset as usize,
        row_idx_bytes,
        "row_idx",
    )?;
    check_range(
        mmap.len(),
        header.values_u32_offset as usize,
        values_bytes,
        "values_u32",
    )?;

    Ok(())
}

fn validate_csc(mmap: &Mmap, header: &SharedCacheHeaderV1) -> Result<()> {
    let n_cells = header.n_cells as usize;
    let n_genes = header.n_genes as usize;
    let nnz = header.nnz as usize;

    let mut prev_ptr = 0u64;
    for c in 0..=n_cells {
        let ptr = read_u64_at(mmap, header.col_ptr_offset as usize + c * 8)?;
        if c == 0 && ptr != 0 {
            bail!("col_ptr[0] must be 0");
        }
        if ptr < prev_ptr {
            bail!("col_ptr must be monotonic");
        }
        prev_ptr = ptr;
    }
    if prev_ptr != header.nnz {
        bail!(
            "col_ptr[n_cells] must equal nnz: {} vs {}",
            prev_ptr,
            header.nnz
        );
    }

    for c in 0..n_cells {
        let start = read_u64_at(mmap, header.col_ptr_offset as usize + c * 8)? as usize;
        let end = read_u64_at(mmap, header.col_ptr_offset as usize + (c + 1) * 8)? as usize;
        if start > end || end > nnz {
            bail!("invalid col_ptr range for column {}", c);
        }
        let mut prev_row: Option<u32> = None;
        for k in start..end {
            let row = read_u32_at(mmap, header.row_idx_offset as usize + k * 4)?;
            if row as usize >= n_genes {
                bail!("row_idx out of bounds at nnz index {}", k);
            }
            if let Some(prev) = prev_row {
                if row <= prev {
                    bail!("row_idx must be strictly increasing within each column");
                }
            }
            prev_row = Some(row);
        }
    }

    Ok(())
}

fn parse_string_table(
    mmap: &Mmap,
    table_offset: usize,
    table_bytes: usize,
    expected_count: usize,
    label: &str,
) -> Result<Vec<String>> {
    check_range(mmap.len(), table_offset, table_bytes, label)?;
    let table = &mmap[table_offset..table_offset + table_bytes];
    if table.len() < 8 {
        bail!("{} table too small", label);
    }
    let count = le_u32(table, 0)? as usize;
    if count != expected_count {
        bail!(
            "{} table count mismatch: expected {} got {}",
            label,
            expected_count,
            count
        );
    }
    let offsets_bytes = (count + 1) * 4;
    let blob_offset = 4 + offsets_bytes;
    if table.len() < blob_offset {
        bail!("{} table offsets truncated", label);
    }

    let mut offsets = Vec::with_capacity(count + 1);
    for i in 0..=count {
        let off = le_u32(table, 4 + i * 4)? as usize;
        offsets.push(off);
    }
    for i in 1..offsets.len() {
        if offsets[i] < offsets[i - 1] {
            bail!("{} table offsets must be monotonic", label);
        }
    }
    let blob = &table[blob_offset..];
    if *offsets.last().unwrap_or(&0) != blob.len() {
        bail!(
            "{} table terminal offset mismatch: expected {} got {}",
            label,
            blob.len(),
            offsets.last().copied().unwrap_or(0)
        );
    }

    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let start = offsets[i];
        let end = offsets[i + 1];
        let s = std::str::from_utf8(&blob[start..end])
            .with_context(|| format!("{} table contains invalid UTF-8 at entry {}", label, i))?;
        out.push(s.to_string());
    }
    Ok(out)
}

fn check_range(file_len: usize, offset: usize, len: usize, label: &str) -> Result<()> {
    let Some(end) = offset.checked_add(len) else {
        bail!("{} range overflow", label);
    };
    if end > file_len {
        bail!("{} out of file bounds", label);
    }
    Ok(())
}

fn le_u16(buf: &[u8], offset: usize) -> Result<u16> {
    let end = offset + 2;
    if end > buf.len() {
        bail!("header truncated while reading u16");
    }
    Ok(u16::from_le_bytes(buf[offset..end].try_into().unwrap()))
}

fn le_u32(buf: &[u8], offset: usize) -> Result<u32> {
    let end = offset + 4;
    if end > buf.len() {
        bail!("header truncated while reading u32");
    }
    Ok(u32::from_le_bytes(buf[offset..end].try_into().unwrap()))
}

fn le_u64(buf: &[u8], offset: usize) -> Result<u64> {
    let end = offset + 8;
    if end > buf.len() {
        bail!("header truncated while reading u64");
    }
    Ok(u64::from_le_bytes(buf[offset..end].try_into().unwrap()))
}

fn read_u64_at(mmap: &Mmap, offset: usize) -> Result<u64> {
    check_range(mmap.len(), offset, 8, "u64 read")?;
    Ok(u64::from_le_bytes(
        mmap[offset..offset + 8].try_into().unwrap(),
    ))
}

fn read_u32_at(mmap: &Mmap, offset: usize) -> Result<u32> {
    check_range(mmap.len(), offset, 4, "u32 read")?;
    Ok(u32::from_le_bytes(
        mmap[offset..offset + 4].try_into().unwrap(),
    ))
}
