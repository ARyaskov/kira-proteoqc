use std::fs;
use std::path::Path;

use assert_cmd::Command;
use serde_json::Value;
use tempfile::TempDir;

#[test]
fn pipeline_tsv_header_order_is_exact() {
    let tmp = TempDir::new().unwrap();
    let out = TempDir::new().unwrap();
    write_10x(tmp.path());
    run_pipeline(tmp.path(), out.path());

    let tsv = out.path().join("kira-proteoqc").join("proteoqc.tsv");
    let header = fs::read_to_string(tsv)
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_string();
    assert_eq!(
        header,
        "barcode\tsample\tcondition\tspecies\tlibsize\tnnz\texpressed_genes\tproteostasis_load\tmisfolded_protein_burden\tchaperone_capacity\tproteasome_activity_proxy\tprotein_quality_balance\tstress_proteostasis_index\tregime\tflags\tconfidence"
    );
}

#[test]
fn pipeline_summary_json_schema_fields_exist() {
    let tmp = TempDir::new().unwrap();
    let out = TempDir::new().unwrap();
    write_10x(tmp.path());
    run_pipeline(tmp.path(), out.path());

    let summary = out.path().join("kira-proteoqc").join("summary.json");
    let v: Value = serde_json::from_slice(&fs::read(summary).unwrap()).unwrap();
    assert_eq!(v["tool"]["name"], "kira-proteoqc");
    assert!(v["tool"]["version"].is_string());
    assert!(v["tool"]["simd"].is_string());
    assert!(v["input"]["n_cells"].is_number());
    assert!(v["input"]["species"].is_string());
    assert!(v["distributions"]["proteostasis_load"]["median"].is_number());
    assert!(v["regimes"]["counts"].is_object());
    assert!(v["regimes"]["fractions"].is_object());
    assert!(v["qc"]["low_confidence_fraction"].is_number());
    assert!(v["qc"]["low_chaperone_signal_fraction"].is_number());
}

#[test]
fn pipeline_step_json_matches_schema() {
    let tmp = TempDir::new().unwrap();
    let out = TempDir::new().unwrap();
    write_10x(tmp.path());
    run_pipeline(tmp.path(), out.path());

    let step = out.path().join("kira-proteoqc").join("pipeline_step.json");
    let v: Value = serde_json::from_slice(&fs::read(step).unwrap()).unwrap();
    assert_eq!(v["tool"]["name"], "kira-proteoqc");
    assert_eq!(v["tool"]["stage"], "proteostasis");
    assert!(v["tool"]["version"].is_string());
    assert_eq!(v["artifacts"]["summary"], "summary.json");
    assert_eq!(v["artifacts"]["primary_metrics"], "proteoqc.tsv");
    assert_eq!(v["artifacts"]["panels"], "panels_report.tsv");
    assert_eq!(v["cell_metrics"]["file"], "proteoqc.tsv");
    assert_eq!(v["cell_metrics"]["id_column"], "barcode");
    assert_eq!(v["cell_metrics"]["regime_column"], "regime");
    assert_eq!(v["cell_metrics"]["confidence_column"], "confidence");
    assert_eq!(v["cell_metrics"]["flag_column"], "flags");
    assert!(v["regimes"].is_array());
}

#[test]
fn pipeline_outputs_are_deterministic() {
    let in_dir = TempDir::new().unwrap();
    write_10x(in_dir.path());
    let out1 = TempDir::new().unwrap();
    let out2 = TempDir::new().unwrap();

    run_pipeline(in_dir.path(), out1.path());
    run_pipeline(in_dir.path(), out2.path());

    for name in [
        "proteoqc.tsv",
        "summary.json",
        "panels_report.tsv",
        "pipeline_step.json",
    ] {
        let a = fs::read(out1.path().join("kira-proteoqc").join(name)).unwrap();
        let b = fs::read(out2.path().join("kira-proteoqc").join(name)).unwrap();
        assert_eq!(a, b, "mismatch in {}", name);
    }
}

fn run_pipeline(input: &Path, out: &Path) {
    let mut cmd = Command::cargo_bin("kira-proteoqc").unwrap();
    cmd.args([
        "run",
        "--input",
        input.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--mode",
        "cell",
        "--run-mode",
        "pipeline",
    ]);
    cmd.assert().success();
}

fn write_10x(dir: &Path) {
    fs::write(
        dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 2 3\n1 1 5\n3 1 1\n2 2 7\n",
    )
    .unwrap();
    fs::write(dir.join("features.tsv"), "g1\tG1\ng2\tG2\ng3\tG3\n").unwrap();
    fs::write(dir.join("barcodes.tsv"), "C1\nC2\n").unwrap();
}
