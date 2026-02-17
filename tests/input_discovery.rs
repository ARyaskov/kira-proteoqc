use std::fs;

use kira_proteoqc::input;
use tempfile::TempDir;

#[test]
fn detects_prefix_when_prefixed_files_exist() {
    let tmp = TempDir::new().unwrap();
    fs::write(tmp.path().join("ABC_matrix.mtx"), "x").unwrap();
    fs::write(tmp.path().join("ABC_features.tsv"), "x").unwrap();
    fs::write(tmp.path().join("ABC_barcodes.tsv"), "x").unwrap();
    let prefix = input::detect_prefix(tmp.path()).unwrap();
    assert_eq!(prefix.as_deref(), Some("ABC"));
}

#[test]
fn resolves_shared_cache_filename() {
    let tmp = TempDir::new().unwrap();
    let p1 = input::resolve_shared_cache_path(tmp.path(), None);
    let p2 = input::resolve_shared_cache_path(tmp.path(), Some("XYZ"));
    assert_eq!(
        p1.file_name().and_then(|s| s.to_str()),
        Some("kira-organelle.bin")
    );
    assert_eq!(
        p2.file_name().and_then(|s| s.to_str()),
        Some("XYZ.kira-organelle.bin")
    );
}
