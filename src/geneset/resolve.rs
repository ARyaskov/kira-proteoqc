use std::collections::HashSet;

use crate::geneset::{GenesetCollection, GenesetDef};

#[derive(Debug, Clone)]
pub struct ResolvedGeneset {
    pub id: String,
    pub axis: char,
    pub gene_ids: Vec<usize>,
    pub missing: Vec<String>,
    pub total: usize,
}

pub fn resolve_collection(
    mut collection: GenesetCollection,
    gene_index: &std::collections::HashMap<String, usize>,
) -> GenesetCollection {
    let mut resolved = Vec::with_capacity(collection.defs.len());
    for def in &collection.defs {
        resolved.push(resolve_def(def, gene_index));
    }
    collection.resolved = resolved;
    collection
}

fn resolve_def(
    def: &GenesetDef,
    gene_index: &std::collections::HashMap<String, usize>,
) -> ResolvedGeneset {
    let mut gene_ids = Vec::new();
    let mut missing = Vec::new();
    let mut seen = HashSet::new();

    for symbol in &def.genes {
        if let Some(&gid) = gene_index.get(symbol) {
            if seen.insert(gid) {
                gene_ids.push(gid);
            }
        } else {
            missing.push(symbol.clone());
        }
    }

    gene_ids.sort_unstable();

    ResolvedGeneset {
        id: def.id.clone(),
        axis: def.axis,
        gene_ids,
        missing,
        total: def.genes.len(),
    }
}
