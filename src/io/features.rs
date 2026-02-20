use std::path::Path;

use anyhow::{Result, bail};
use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

pub fn read_features(path: &Path) -> Result<Vec<String>> {
    let md = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_metadata()
    .map_err(|e| anyhow::anyhow!(e.message))?;

    if md.gene_symbols.is_empty() {
        bail!("features.tsv is empty");
    }

    Ok(md.gene_symbols)
}
