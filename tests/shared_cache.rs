use std::fs;
use std::path::Path;

use crc::{CRC_64_ECMA_182, Crc};
use kira_proteoqc::io::shared_cache::SharedCache;
use tempfile::TempDir;

#[test]
fn reads_valid_shared_cache_and_traverses_arrays() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("kira-organelle.bin");
    write_tiny_cache(&path);

    let cache = SharedCache::open(&path).unwrap();
    assert_eq!(cache.header.n_genes, 3);
    assert_eq!(cache.header.n_cells, 2);
    assert_eq!(cache.header.nnz, 3);
    assert_eq!(cache.genes, vec!["G1", "G2", "G3"]);
    assert_eq!(cache.barcodes, vec!["C1", "C2"]);
    assert_eq!(cache.col_ptr().unwrap(), vec![0, 2, 3]);
    assert_eq!(cache.row_idx_at(0).unwrap(), 0);
    assert_eq!(cache.row_idx_at(1).unwrap(), 2);
    assert_eq!(cache.row_idx_at(2).unwrap(), 1);
    assert_eq!(cache.value_u32_at(0).unwrap(), 5);
    assert_eq!(cache.value_u32_at(1).unwrap(), 1);
    assert_eq!(cache.value_u32_at(2).unwrap(), 7);
}

#[test]
fn rejects_tampered_header_crc() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("kira-organelle.bin");
    write_tiny_cache(&path);
    let mut bytes = fs::read(&path).unwrap();
    bytes[16] ^= 0x01;
    fs::write(&path, bytes).unwrap();
    let err = SharedCache::open(&path).unwrap_err().to_string();
    assert!(err.contains("CRC"), "unexpected error: {err}");
}

fn write_tiny_cache(path: &Path) {
    let genes = vec!["G1", "G2", "G3"];
    let barcodes = vec!["C1", "C2"];
    let col_ptr = vec![0u64, 2, 3];
    let row_idx = vec![0u32, 2, 1];
    let values = vec![5u32, 1, 7];

    let genes_table = make_table(&genes);
    let barcodes_table = make_table(&barcodes);

    let mut offset = 256usize;
    let genes_off = offset;
    offset += genes_table.len();
    offset = align64(offset);
    let barcodes_off = offset;
    offset += barcodes_table.len();
    offset = align64(offset);
    let col_ptr_off = offset;
    offset += col_ptr.len() * 8;
    offset = align64(offset);
    let row_idx_off = offset;
    offset += row_idx.len() * 4;
    offset = align64(offset);
    let values_off = offset;
    offset += values.len() * 4;
    let file_bytes = offset as u64;

    let mut bytes = vec![0u8; offset];
    bytes[0..4].copy_from_slice(b"KORG");
    bytes[4..6].copy_from_slice(&1u16.to_le_bytes());
    bytes[6..8].copy_from_slice(&0u16.to_le_bytes());
    bytes[8..12].copy_from_slice(&0x1234_5678u32.to_le_bytes());
    bytes[12..16].copy_from_slice(&256u32.to_le_bytes());
    bytes[16..24].copy_from_slice(&(genes.len() as u64).to_le_bytes());
    bytes[24..32].copy_from_slice(&(barcodes.len() as u64).to_le_bytes());
    bytes[32..40].copy_from_slice(&(values.len() as u64).to_le_bytes());
    bytes[40..48].copy_from_slice(&(genes_off as u64).to_le_bytes());
    bytes[48..56].copy_from_slice(&(genes_table.len() as u64).to_le_bytes());
    bytes[56..64].copy_from_slice(&(barcodes_off as u64).to_le_bytes());
    bytes[64..72].copy_from_slice(&(barcodes_table.len() as u64).to_le_bytes());
    bytes[72..80].copy_from_slice(&(col_ptr_off as u64).to_le_bytes());
    bytes[80..88].copy_from_slice(&(row_idx_off as u64).to_le_bytes());
    bytes[88..96].copy_from_slice(&(values_off as u64).to_le_bytes());
    bytes[96..104].copy_from_slice(&0u64.to_le_bytes());
    bytes[104..112].copy_from_slice(&0u64.to_le_bytes());
    bytes[112..120].copy_from_slice(&file_bytes.to_le_bytes());
    bytes[128..136].copy_from_slice(&0u64.to_le_bytes());

    bytes[genes_off..genes_off + genes_table.len()].copy_from_slice(&genes_table);
    bytes[barcodes_off..barcodes_off + barcodes_table.len()].copy_from_slice(&barcodes_table);

    for (i, v) in col_ptr.iter().enumerate() {
        let start = col_ptr_off + i * 8;
        bytes[start..start + 8].copy_from_slice(&v.to_le_bytes());
    }
    for (i, v) in row_idx.iter().enumerate() {
        let start = row_idx_off + i * 4;
        bytes[start..start + 4].copy_from_slice(&v.to_le_bytes());
    }
    for (i, v) in values.iter().enumerate() {
        let start = values_off + i * 4;
        bytes[start..start + 4].copy_from_slice(&v.to_le_bytes());
    }

    let crc = Crc::<u64>::new(&CRC_64_ECMA_182);
    let mut hdr = bytes[..256].to_vec();
    for b in &mut hdr[120..128] {
        *b = 0;
    }
    let crc_val = crc.checksum(&hdr);
    bytes[120..128].copy_from_slice(&crc_val.to_le_bytes());

    fs::write(path, bytes).unwrap();
}

fn make_table(strings: &[&str]) -> Vec<u8> {
    let mut blob = Vec::new();
    let mut offsets = Vec::with_capacity(strings.len() + 1);
    offsets.push(0u32);
    for s in strings {
        blob.extend_from_slice(s.as_bytes());
        offsets.push(blob.len() as u32);
    }
    let mut out = Vec::with_capacity(4 + offsets.len() * 4 + blob.len());
    out.extend_from_slice(&(strings.len() as u32).to_le_bytes());
    for off in offsets {
        out.extend_from_slice(&off.to_le_bytes());
    }
    out.extend_from_slice(&blob);
    out
}

fn align64(x: usize) -> usize {
    (x + 63) & !63
}
