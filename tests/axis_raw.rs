use std::fs::File;
use std::io::{BufWriter, Write};

use kira_proteoqc::ctx::Ctx;
use kira_proteoqc::expr::layout::{ExprHeaderV1, LAYOUT_CSC, VERSION, write_header};
use kira_proteoqc::expr::reader::open_mmap;
use kira_proteoqc::geneset::{GenesetCollection, ResolvedGeneset};
use kira_proteoqc::schema::v1::Mode;
use kira_proteoqc::scores::axis_raw::compute_axis_raw;
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

fn make_genesets() -> GenesetCollection {
    let resolved = vec![
        ResolvedGeneset {
            id: "proteasome_core".into(),
            axis: 'A',
            gene_ids: vec![0],
            missing: vec![],
            total: 1,
        },
        ResolvedGeneset {
            id: "proteasome_regulator".into(),
            axis: 'B',
            gene_ids: vec![1],
            missing: vec![],
            total: 1,
        },
        ResolvedGeneset {
            id: "ubiquitin_axis".into(),
            axis: 'D',
            gene_ids: vec![2],
            missing: vec![],
            total: 1,
        },
        ResolvedGeneset {
            id: "e3_ligases".into(),
            axis: 'D',
            gene_ids: vec![1],
            missing: vec![],
            total: 1,
        },
        ResolvedGeneset {
            id: "dubs".into(),
            axis: 'D',
            gene_ids: vec![0],
            missing: vec![],
            total: 1,
        },
        ResolvedGeneset {
            id: "chaperone_hsp70".into(),
            axis: 'E',
            gene_ids: vec![0],
            missing: vec![],
            total: 1,
        },
        ResolvedGeneset {
            id: "chaperone_hsp90".into(),
            axis: 'E',
            gene_ids: vec![1],
            missing: vec![],
            total: 1,
        },
        ResolvedGeneset {
            id: "chaperone_hsp40".into(),
            axis: 'E',
            gene_ids: vec![2],
            missing: vec![],
            total: 1,
        },
        ResolvedGeneset {
            id: "erad".into(),
            axis: 'F',
            gene_ids: vec![1],
            missing: vec![],
            total: 1,
        },
        ResolvedGeneset {
            id: "ribosome_load".into(),
            axis: 'F',
            gene_ids: vec![2],
            missing: vec![],
            total: 1,
        },
    ];
    GenesetCollection {
        version: "v1".into(),
        defs: vec![],
        resolved,
    }
}

fn assert_vec_close(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b.iter()) {
        assert!((x - y).abs() < 1e-6, "{} vs {}", x, y);
    }
}

#[test]
fn axis_raw_per_cell() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("expr.bin");
    write_expr(&path);
    let (header, mmap) = open_mmap(&path).unwrap();

    let mut ctx = Ctx::new(
        tmp.path().to_path_buf(),
        tmp.path().join("out"),
        Mode::Cell,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );
    ctx.expr_header = Some(header);
    ctx.expr_mmap = Some(mmap);
    ctx.genes = vec!["G1".into(), "G2".into(), "G3".into()];
    ctx.cells = vec!["C1".into(), "C2".into(), "C3".into(), "C4".into()];
    ctx.nnz = 5;
    ctx.genesets = Some(make_genesets());

    let axis = compute_axis_raw(&mut ctx).unwrap();
    assert_vec_close(&axis.pcs, &[0.6, 1.2, 1.2, 1.6]);
    assert_vec_close(&axis.utp, &[-0.15, 1.05, 2.2, 1.4]);
    assert_vec_close(&axis.cls, &[0.45, 1.05, 1.9, 1.4]);
    assert_vec_close(&axis.erad, &[0.0, 3.0, 0.0, 4.0]);
    assert_vec_close(&axis.ribo, &[0.0, 0.0, 5.0, 0.0]);
}

#[test]
fn axis_raw_per_sample_mean() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("expr.bin");
    write_expr(&path);
    let (header, mmap) = open_mmap(&path).unwrap();

    let mut ctx = Ctx::new(
        tmp.path().to_path_buf(),
        tmp.path().join("out"),
        Mode::Sample,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );
    ctx.expr_header = Some(header);
    ctx.expr_mmap = Some(mmap);
    ctx.genes = vec!["G1".into(), "G2".into(), "G3".into()];
    ctx.cells = vec!["C1".into(), "C2".into(), "C3".into(), "C4".into()];
    ctx.nnz = 5;
    ctx.genesets = Some(make_genesets());

    let axis = compute_axis_raw(&mut ctx).unwrap();
    let pcs_mean = (0.6 + 1.2 + 1.2 + 1.6) / 4.0;
    assert!((axis.pcs[0] - pcs_mean).abs() < 1e-6);
}

#[test]
fn axis_missing_geneset_zero() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("expr.bin");
    write_expr(&path);
    let (header, mmap) = open_mmap(&path).unwrap();

    let mut ctx = Ctx::new(
        tmp.path().to_path_buf(),
        tmp.path().join("out"),
        Mode::Cell,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );
    ctx.expr_header = Some(header);
    ctx.expr_mmap = Some(mmap);
    ctx.genes = vec!["G1".into(), "G2".into(), "G3".into()];
    ctx.cells = vec!["C1".into(), "C2".into(), "C3".into(), "C4".into()];
    ctx.nnz = 5;
    let mut gs = make_genesets();
    gs.resolved.retain(|g| g.id != "e3_ligases");
    ctx.genesets = Some(gs);

    let axis = compute_axis_raw(&mut ctx).unwrap();
    assert!((axis.utp[1] - 0.0).abs() < 1e-6);
    assert!(!ctx.warnings.is_empty());
}
