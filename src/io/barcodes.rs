use std::path::Path;

use anyhow::{Result, bail};
use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

pub fn read_barcodes(path: &Path) -> Result<Vec<String>> {
    let md = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_metadata()
    .map_err(|e| anyhow::anyhow!(e.message))?;

    if md.barcodes.is_empty() {
        bail!("barcodes.tsv is empty");
    }

    Ok(md.barcodes)
}
