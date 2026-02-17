use std::fs::File;
use std::io::{BufWriter, Write};

use kira_proteoqc::ctx::Ctx;
use kira_proteoqc::expr::layout::{ExprHeaderV1, LAYOUT_CSC, VERSION, write_header};
use kira_proteoqc::expr::reader::open_mmap;
use kira_proteoqc::geneset::{GenesetCollection, ResolvedGeneset};
use kira_proteoqc::schema::v1::Mode;
use kira_proteoqc::scores::risk::compute_risk_flags;
use kira_proteoqc::scores::{AxisRawScores, IntegratedScores};
use tempfile::TempDir;

fn write_empty_expr(path: &std::path::Path, n_genes: u32, n_cells: u32) {
    let header = ExprHeaderV1 {
        version: VERSION,
        n_genes,
        n_cells,
        nnz: 0,
        layout: LAYOUT_CSC,
    };
    let gene_ptr: Vec<u64> = vec![0, 0];

    let file = File::create(path).unwrap();
    let mut w = BufWriter::new(file);
    write_header(&mut w, &header).unwrap();
    for v in &gene_ptr {
        w.write_all(&v.to_le_bytes()).unwrap();
    }
    w.flush().unwrap();
}

fn make_genesets() -> GenesetCollection {
    let mut resolved = Vec::new();
    let ids = vec![
        "proteasome_core",
        "proteasome_regulator",
        "ubiquitin_axis",
        "e3_ligases",
        "dubs",
        "chaperone_hsp70",
        "chaperone_hsp90",
        "chaperone_hsp40",
        "erad",
        "ribosome_load",
    ];
    for id in ids {
        resolved.push(ResolvedGeneset {
            id: id.to_string(),
            axis: 'A',
            gene_ids: vec![0],
            missing: vec![],
            total: 1,
        });
    }
    GenesetCollection {
        version: "v1".to_string(),
        defs: vec![],
        resolved,
    }
}

#[test]
fn risk_sample_mode_semantics() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("expr.bin");
    write_empty_expr(&path, 1, 10);
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
    ctx.genes = vec!["G1".into()];
    ctx.cells = (0..10).map(|i| format!("C{}", i)).collect();
    ctx.nnz = 0;
    ctx.genesets = Some(make_genesets());

    ctx.axis_raw = Some(AxisRawScores {
        pcs: vec![0.0],
        utp: vec![0.0],
        cls: vec![0.0],
        erad: vec![0.0],
        ribo: vec![0.0],
    });
    ctx.integrated_scores = Some(IntegratedScores {
        capacity_raw: vec![0.0],
        pii_raw: vec![0.0],
        pfs_raw: vec![0.0],
        capacity_z: None,
        pii_z: None,
        pfs_z: None,
    });

    let flags = compute_risk_flags(&mut ctx).unwrap();
    let fragile = flags.iter().find(|f| f.name == "fragile_high").unwrap();
    assert!(!fragile.fired);
}
