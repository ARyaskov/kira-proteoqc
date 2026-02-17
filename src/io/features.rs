use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Result, bail};

use crate::io::open_maybe_gz;

pub fn read_features(path: &Path) -> Result<Vec<String>> {
    let reader = open_maybe_gz(path)?;
    let mut reader = BufReader::new(reader);

    let mut features = Vec::new();
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            line.clear();
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.is_empty() {
            bail!("features.tsv line has no columns");
        }
        let symbol = if parts.len() > 1 && !parts[1].is_empty() {
            parts[1]
        } else {
            parts[0]
        };
        features.push(symbol.to_string());
        line.clear();
    }

    if features.is_empty() {
        bail!("features.tsv is empty");
    }

    Ok(features)
}
