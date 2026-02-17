use kira_proteoqc::ctx::Ctx;
use kira_proteoqc::io::summary::format_summary;
use kira_proteoqc::schema::v1::Mode;
use kira_proteoqc::scores::{AxisRawScores, IntegratedScores, RiskFlag};

#[test]
fn summary_format() {
    let mut ctx = Ctx::new(
        std::path::PathBuf::from("input"),
        std::path::PathBuf::from("out"),
        Mode::Sample,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );
    ctx.input_meta.genes = Some(100);
    ctx.input_meta.cells = Some(10);
    ctx.axis_raw = Some(AxisRawScores {
        pcs: vec![1.0],
        utp: vec![1.0],
        cls: vec![1.0],
        erad: vec![1.0],
        ribo: vec![1.0],
    });
    ctx.integrated_scores = Some(IntegratedScores {
        capacity_raw: vec![1.0],
        pii_raw: vec![1.0],
        pfs_raw: vec![0.87],
        capacity_z: None,
        pii_z: None,
        pfs_z: None,
    });
    ctx.risk_flags = vec![RiskFlag {
        name: "fragile_high".to_string(),
        fired: true,
        threshold: "x".to_string(),
        details: None,
    }];

    let s = format_summary(&ctx).unwrap();
    assert!(s.contains("kira-proteoqc v"));
    assert!(s.contains("Input: 100 genes, 10 cells, mode=sample"));
    assert!(s.contains("PFS:"));
    assert!(s.contains("Flags: fragile_high"));
}
