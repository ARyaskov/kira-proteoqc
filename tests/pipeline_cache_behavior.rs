use std::fs;
use std::path::Path;

use crc::{CRC_64_ECMA_182, Crc};
use kira_proteoqc::ctx::{Ctx, RunMode};
use kira_proteoqc::pipeline::Pipeline;
use kira_proteoqc::pipeline::stage1_input::Stage1Input;
use kira_proteoqc::schema::v1::Mode;
use tempfile::TempDir;

#[test]
fn pipeline_mode_uses_cache_when_present() {
    let tmp = TempDir::new().unwrap();
    let cache = tmp.path().join("kira-organelle.bin");
    write_tiny_cache(&cache);
    let mut ctx = make_ctx(tmp.path());
    Pipeline::new(vec![Box::new(Stage1Input::new())])
        .run(&mut ctx)
        .unwrap();
    assert!(ctx.shared_cache_used);
    assert_eq!(ctx.genes.len(), 3);
    assert_eq!(ctx.cells.len(), 2);
    assert_eq!(ctx.nnz, 3);
}

#[test]
fn pipeline_mode_falls_back_to_mtx_when_cache_missing() {
    let tmp = TempDir::new().unwrap();
    write_10x(
        tmp.path(),
        "g1\tG1\ng2\tG2\n",
        "c1\nc2\n",
        "%%MatrixMarket matrix coordinate integer general\n2 2 1\n1 1 1\n",
    );
    let mut ctx = make_ctx(tmp.path());
    Pipeline::new(vec![Box::new(Stage1Input::new())])
        .run(&mut ctx)
        .unwrap();
    assert!(!ctx.shared_cache_used);
    assert!(!ctx.warnings.is_empty());
}

#[test]
fn pipeline_mode_errors_on_invalid_cache() {
    let tmp = TempDir::new().unwrap();
    fs::write(tmp.path().join("kira-organelle.bin"), b"bad").unwrap();
    let mut ctx = make_ctx(tmp.path());
    let err = Pipeline::new(vec![Box::new(Stage1Input::new())])
        .run(&mut ctx)
        .unwrap_err()
        .to_string();
    assert!(err.contains("shared cache"));
}

fn make_ctx(input: &Path) -> Ctx {
    let mut ctx = Ctx::new(
        input.to_path_buf(),
        input.join("out"),
        Mode::Cell,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );
    ctx.run_mode = RunMode::Pipeline;
    ctx
}

fn write_10x(dir: &Path, features: &str, barcodes: &str, mtx: &str) {
    fs::write(dir.join("features.tsv"), features).unwrap();
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
    fs::write(dir.join("matrix.mtx"), mtx).unwrap();
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
    let file_bytes = bytes.len() as u64;
    bytes[112..120].copy_from_slice(&file_bytes.to_le_bytes());
    bytes[128..136].copy_from_slice(&0u64.to_le_bytes());

    bytes[genes_off..genes_off + genes_table.len()].copy_from_slice(&genes_table);
    bytes[barcodes_off..barcodes_off + barcodes_table.len()].copy_from_slice(&barcodes_table);

    for (i, v) in col_ptr.iter().enumerate() {
        let s = col_ptr_off + i * 8;
        bytes[s..s + 8].copy_from_slice(&v.to_le_bytes());
    }
    for (i, v) in row_idx.iter().enumerate() {
        let s = row_idx_off + i * 4;
        bytes[s..s + 4].copy_from_slice(&v.to_le_bytes());
    }
    for (i, v) in values.iter().enumerate() {
        let s = values_off + i * 4;
        bytes[s..s + 4].copy_from_slice(&v.to_le_bytes());
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
    let mut offs = vec![0u32];
    for s in strings {
        blob.extend_from_slice(s.as_bytes());
        offs.push(blob.len() as u32);
    }
    let mut out = Vec::new();
    out.extend_from_slice(&(strings.len() as u32).to_le_bytes());
    for o in offs {
        out.extend_from_slice(&o.to_le_bytes());
    }
    out.extend_from_slice(&blob);
    out
}

fn align64(x: usize) -> usize {
    (x + 63) & !63
}
