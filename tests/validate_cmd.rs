use std::fs;
use std::path::Path;

use assert_cmd::Command;
use tempfile::TempDir;

fn write_10x(dir: &Path, features: &str, barcodes: &str, mtx: &str) {
    fs::write(dir.join("features.tsv"), features).unwrap();
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
    fs::write(dir.join("matrix.mtx"), mtx).unwrap();
}

#[test]
fn validate_command_ok() {
    let tmp = TempDir::new().unwrap();
    let mtx = "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 1\n2 2 1\n";
    let features = "g1\tGeneA\ng2\tGeneB\n";
    let barcodes = "cell1\ncell2\n";
    write_10x(tmp.path(), features, barcodes, mtx);

    let mut cmd = Command::cargo_bin("kira-proteoqc").unwrap();
    cmd.arg("validate").arg("--input").arg(tmp.path());
    cmd.assert().success();
}
