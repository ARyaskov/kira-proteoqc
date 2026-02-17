use kira_proteoqc::ctx::Ctx;
use kira_proteoqc::io::json_writer::build_report;
use kira_proteoqc::schema::v1::Mode;
use kira_proteoqc::scores::{AxisRawScores, IntegratedScores, PfsContributions, RiskFlag};
use serde_json::Value;

#[test]
fn json_report_populated() {
    let mut ctx = Ctx::new(
        std::path::PathBuf::from("input"),
        std::path::PathBuf::from("out"),
        Mode::Cell,
        false,
        None,
        true,
        true,
        true,
        "0.0.0-test",
    );
    ctx.input_meta.genes = Some(3);
    ctx.input_meta.cells = Some(2);
    ctx.input_meta.nnz = Some(5);

    ctx.axis_raw = Some(AxisRawScores {
        pcs: vec![1.0, 2.0],
        utp: vec![0.5, 0.5],
        cls: vec![1.0, 1.0],
        erad: vec![0.2, 0.2],
        ribo: vec![2.0, 2.0],
    });
    ctx.integrated_scores = Some(IntegratedScores {
        capacity_raw: vec![1.0, 1.0],
        pii_raw: vec![1.0, 1.0],
        pfs_raw: vec![1.0, 1.0],
        capacity_z: Some(vec![0.0, 0.0]),
        pii_z: Some(vec![0.0, 0.0]),
        pfs_z: Some(vec![0.0, 0.0]),
    });
    ctx.pfs_contributions = Some(PfsContributions {
        c_pii: vec![0.4, 0.4],
        c_utp: vec![0.1, 0.1],
        c_ribo: vec![0.2, 0.2],
        c_pcs: vec![-0.15, -0.15],
    });
    ctx.risk_flags = vec![RiskFlag {
        name: "fragile_high".to_string(),
        fired: false,
        threshold: ">1.5".to_string(),
        details: None,
    }];

    let report = build_report(&ctx).unwrap();
    let json = serde_json::to_value(report).unwrap();

    assert_eq!(json["tool"], "kira-proteoqc");
    assert_eq!(json["schema_version"], "v1");
    assert_eq!(json["input_meta"]["genes"], 3);
    assert_eq!(json["scores"]["per_cell_tsv_path"], "proteoqc.tsv");
    assert!(json["risk_flags"].is_array());
    assert!(json["explainability"]["pfs_contributions"].is_object());
}
