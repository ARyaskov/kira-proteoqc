use std::fs;

use kira_proteoqc::ctx::Ctx;
use kira_proteoqc::io::tsv_writer::write_tsv;
use kira_proteoqc::schema::v1::Mode;
use kira_proteoqc::scores::{AxisRawScores, IntegratedScores};
use tempfile::TempDir;

#[test]
fn tsv_per_cell_format() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("proteoqc.tsv");

    let mut ctx = Ctx::new(
        std::path::PathBuf::from("input"),
        tmp.path().to_path_buf(),
        Mode::Cell,
        false,
        None,
        true,
        true,
        true,
        "0.0.0-test",
    );
    ctx.cells = vec!["C1".into(), "C2".into()];
    ctx.axis_raw = Some(AxisRawScores {
        pcs: vec![1.0, 2.0],
        utp: vec![0.1, 0.2],
        cls: vec![0.3, 0.4],
        erad: vec![0.5, 0.6],
        ribo: vec![0.7, 0.8],
    });
    ctx.integrated_scores = Some(IntegratedScores {
        capacity_raw: vec![1.1, 1.2],
        pii_raw: vec![1.3, 1.4],
        pfs_raw: vec![1.5, 1.6],
        capacity_z: Some(vec![0.0, 0.0]),
        pii_z: Some(vec![0.0, 0.0]),
        pfs_z: Some(vec![0.1, 0.2]),
    });

    write_tsv(&path, &ctx).unwrap();
    let content = fs::read_to_string(&path).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 3);
    assert!(lines[0].starts_with("cell_id\tPCS_raw"));
}
