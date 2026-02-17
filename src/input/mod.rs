use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;

const MATRIX_SUFFIXES: [&str; 2] = ["_matrix.mtx", "_matrix.mtx.gz"];
const FEATURES_SUFFIXES: [&str; 4] = [
    "_features.tsv",
    "_features.tsv.gz",
    "_genes.tsv",
    "_genes.tsv.gz",
];
const BARCODES_SUFFIXES: [&str; 2] = ["_barcodes.tsv", "_barcodes.tsv.gz"];

pub fn detect_prefix(input_dir: &Path) -> Result<Option<String>> {
    let mut prefixes = BTreeSet::new();
    for entry in fs::read_dir(input_dir)? {
        let entry = entry?;
        let Some(name) = entry.file_name().to_str().map(|s| s.to_string()) else {
            continue;
        };
        for suffix in MATRIX_SUFFIXES
            .iter()
            .chain(FEATURES_SUFFIXES.iter())
            .chain(BARCODES_SUFFIXES.iter())
        {
            if let Some(prefix) = name.strip_suffix(suffix) {
                if !prefix.is_empty() {
                    prefixes.insert(prefix.to_string());
                }
            }
        }
    }
    Ok(prefixes.into_iter().next())
}

pub fn resolve_shared_cache_path(input_dir: &Path, prefix: Option<&str>) -> PathBuf {
    match prefix {
        Some(prefix) if !prefix.is_empty() => {
            input_dir.join(format!("{prefix}.kira-organelle.bin"))
        }
        _ => input_dir.join("kira-organelle.bin"),
    }
}

pub fn resolve_mtx_input_files(
    input_dir: &Path,
    prefix: Option<&str>,
) -> (Option<PathBuf>, Option<PathBuf>, Option<PathBuf>) {
    (
        resolve_candidate(input_dir, prefix, "matrix.mtx", true),
        resolve_features_candidate(input_dir, prefix),
        resolve_candidate(input_dir, prefix, "barcodes.tsv", true),
    )
}

fn resolve_features_candidate(input_dir: &Path, prefix: Option<&str>) -> Option<PathBuf> {
    for base in ["features.tsv", "genes.tsv"] {
        if let Some(path) = resolve_candidate(input_dir, prefix, base, true) {
            return Some(path);
        }
    }
    None
}

fn resolve_candidate(
    input_dir: &Path,
    prefix: Option<&str>,
    base: &str,
    allow_gz: bool,
) -> Option<PathBuf> {
    if let Some(prefix) = prefix {
        let prefixed = input_dir.join(format!("{prefix}_{base}"));
        if prefixed.exists() {
            return Some(prefixed);
        }
        if allow_gz {
            let prefixed_gz = input_dir.join(format!("{prefix}_{base}.gz"));
            if prefixed_gz.exists() {
                return Some(prefixed_gz);
            }
        }
    }

    let plain = input_dir.join(base);
    if plain.exists() {
        return Some(plain);
    }
    if allow_gz {
        let gz = input_dir.join(format!("{base}.gz"));
        if gz.exists() {
            return Some(gz);
        }
    }
    None
}
