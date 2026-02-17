use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Result, bail};

use crate::io::open_maybe_gz;

pub fn read_barcodes(path: &Path) -> Result<Vec<String>> {
    let reader = open_maybe_gz(path)?;
    let mut reader = BufReader::new(reader);

    let mut barcodes = Vec::new();
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            line.clear();
            continue;
        }
        barcodes.push(trimmed.to_string());
        line.clear();
    }

    if barcodes.is_empty() {
        bail!("barcodes.tsv is empty");
    }

    Ok(barcodes)
}
