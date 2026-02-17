use std::fs::File;
use std::io::{BufWriter, Write};

use kira_proteoqc::expr::layout::{ExprHeaderV1, LAYOUT_CSC, VERSION, write_header};
use kira_proteoqc::expr::reader::{ExprReader, open_mmap};
use kira_proteoqc::math::reduce::GeneSetReducer;
use tempfile::TempDir;

fn write_expr(path: &std::path::Path) {
    let header = ExprHeaderV1 {
        version: VERSION,
        n_genes: 3,
        n_cells: 4,
        nnz: 5,
        layout: LAYOUT_CSC,
    };
    let gene_ptr: Vec<u64> = vec![0, 2, 4, 5];
    let cell_idx: Vec<u32> = vec![0, 2, 1, 3, 2];
    let values: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let file = File::create(path).unwrap();
    let mut w = BufWriter::new(file);
    write_header(&mut w, &header).unwrap();
    for v in &gene_ptr {
        w.write_all(&v.to_le_bytes()).unwrap();
    }
    for v in &cell_idx {
        w.write_all(&v.to_le_bytes()).unwrap();
    }
    for v in &values {
        w.write_all(&v.to_le_bytes()).unwrap();
    }
    w.flush().unwrap();
}

#[test]
fn reducer_single_gene() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("expr.bin");
    write_expr(&path);

    let (header, mmap) = open_mmap(&path).unwrap();
    let reader = ExprReader::new(&header, &mmap);
    let mut reducer = GeneSetReducer::new(&reader, 1, 0, false);

    let mut out = vec![0.0f32; reader.n_cells()];
    reducer.per_cell_raw(&[1], &mut out).unwrap();
    assert_eq!(out, vec![0.0, 3.0, 0.0, 4.0]);
}

#[test]
fn reducer_multi_gene_mean() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("expr.bin");
    write_expr(&path);

    let (header, mmap) = open_mmap(&path).unwrap();
    let reader = ExprReader::new(&header, &mmap);
    let mut reducer = GeneSetReducer::new(&reader, 1, 0, false);

    let mut out = vec![0.0f32; reader.n_cells()];
    reducer.per_cell_raw(&[0, 2], &mut out).unwrap();
    assert_eq!(out, vec![0.5, 0.0, 3.5, 0.0]);
}

#[test]
fn reducer_determinism() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("expr.bin");
    write_expr(&path);

    let (header, mmap) = open_mmap(&path).unwrap();
    let reader = ExprReader::new(&header, &mmap);
    let mut reducer = GeneSetReducer::new(&reader, 1, 0, false);

    let mut out1 = vec![0.0f32; reader.n_cells()];
    let mut out2 = vec![0.0f32; reader.n_cells()];
    reducer.per_cell_raw(&[0, 2], &mut out1).unwrap();
    reducer.per_cell_raw(&[0, 2], &mut out2).unwrap();
    assert_eq!(out1, out2);
}
