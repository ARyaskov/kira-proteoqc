use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::geneset::GenesetDef;

pub fn load_builtin_v1() -> Result<Vec<GenesetDef>> {
    let content = include_str!("../../assets/genesets/proteoqc_v1.tsv");
    parse_geneset_tsv(content, "built-in v1")
}

pub fn load_geneset_tsv(path: &Path) -> Result<Vec<GenesetDef>> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read geneset TSV {}", path.display()))?;
    parse_geneset_tsv(&content, &path.display().to_string())
}

pub fn merge_defs(mut builtin: Vec<GenesetDef>, user: Vec<GenesetDef>) -> Vec<GenesetDef> {
    if user.is_empty() {
        return builtin;
    }
    let mut user_map: HashMap<String, GenesetDef> = HashMap::new();
    for def in user {
        user_map.insert(def.id.clone(), def);
    }

    let mut merged = Vec::with_capacity(builtin.len() + user_map.len());
    for def in builtin.drain(..) {
        if let Some(user_def) = user_map.remove(&def.id) {
            merged.push(user_def);
        } else {
            merged.push(def);
        }
    }
    for (_id, def) in user_map.into_iter() {
        merged.push(def);
    }
    merged
}

fn parse_geneset_tsv(content: &str, source: &str) -> Result<Vec<GenesetDef>> {
    let mut defs: HashMap<String, GenesetDef> = HashMap::new();
    let mut order: Vec<String> = Vec::new();

    for (idx, line) in content.lines().enumerate() {
        let line_no = idx + 1;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.len() != 3 {
            bail!("{}:{} malformed TSV (expected 3 columns)", source, line_no);
        }
        let id = parts[0].trim();
        let axis_str = parts[1].trim();
        let gene = parts[2].trim();
        if id.is_empty() || axis_str.is_empty() || gene.is_empty() {
            bail!("{}:{} empty field in TSV", source, line_no);
        }
        let axis = axis_str
            .chars()
            .next()
            .ok_or_else(|| anyhow::anyhow!("{}:{} invalid axis", source, line_no))?;
        if !('A'..='F').contains(&axis) {
            bail!("{}:{} axis must be A-F", source, line_no);
        }

        if !defs.contains_key(id) {
            order.push(id.to_string());
            defs.insert(
                id.to_string(),
                GenesetDef {
                    id: id.to_string(),
                    axis,
                    genes: Vec::new(),
                },
            );
        }

        let def = defs.get_mut(id).unwrap();
        if def.axis != axis {
            bail!("{}:{} axis mismatch for geneset '{}'", source, line_no, id);
        }
        def.genes.push(gene.to_string());
    }

    let mut out = Vec::with_capacity(order.len());
    for id in order {
        if let Some(def) = defs.remove(&id) {
            out.push(def);
        }
    }

    Ok(out)
}
