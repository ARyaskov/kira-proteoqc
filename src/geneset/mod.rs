mod loader;
mod resolve;

use std::path::Path;

use anyhow::Result;

pub use loader::{load_builtin_v1, load_geneset_tsv, merge_defs};
pub use resolve::{ResolvedGeneset, resolve_collection};

#[derive(Debug, Clone)]
pub struct GenesetDef {
    pub id: String,
    pub axis: char,
    pub genes: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct GenesetCollection {
    pub version: String,
    pub defs: Vec<GenesetDef>,
    pub resolved: Vec<ResolvedGeneset>,
}

pub fn load_builtin() -> Result<GenesetCollection> {
    let defs = load_builtin_v1()?;
    Ok(GenesetCollection {
        version: "v1".to_string(),
        defs,
        resolved: Vec::new(),
    })
}

pub fn load_user(path: &Path) -> Result<Vec<GenesetDef>> {
    load_geneset_tsv(path)
}
