use std::path::{Path, PathBuf};

use anyhow::Result;

pub fn detect_prefix(input_dir: &Path) -> Result<Option<String>> {
    kira_scio::detect_prefix(input_dir).map_err(|e| anyhow::anyhow!(e.to_string()))
}

pub fn resolve_shared_cache_path(input_dir: &Path, prefix: Option<&str>) -> PathBuf {
    input_dir.join(kira_scio::resolve_shared_cache_filename(prefix))
}

pub fn resolve_mtx_input_files(
    input_dir: &Path,
    _prefix: Option<&str>,
) -> (Option<PathBuf>, Option<PathBuf>, Option<PathBuf>) {
    match kira_scio::discover(input_dir) {
        Ok(ds) => (Some(ds.matrix), ds.features.or(ds.genes), ds.barcodes),
        Err(_) => (None, None, None),
    }
}
