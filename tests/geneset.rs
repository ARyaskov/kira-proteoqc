use std::collections::HashMap;
use std::fs;

use kira_proteoqc::geneset::{
    GenesetCollection, GenesetDef, load_geneset_tsv, merge_defs, resolve_collection,
};
use tempfile::TempDir;

#[test]
fn geneset_tsv_parse_order() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("gs.tsv");
    let content = "#geneset_id\taxis\tgene_symbol\nset_a\tA\tG1\nset_a\tA\tG2\nset_b\tB\tX1\n";
    fs::write(&path, content).unwrap();

    let defs = load_geneset_tsv(&path).unwrap();
    assert_eq!(defs.len(), 2);
    assert_eq!(defs[0].id, "set_a");
    assert_eq!(defs[0].genes, vec!["G1", "G2"]);
    assert_eq!(defs[1].id, "set_b");
}

#[test]
fn geneset_resolution() {
    let defs = vec![GenesetDef {
        id: "set_a".to_string(),
        axis: 'A',
        genes: vec!["G1".to_string(), "G2".to_string(), "G3".to_string()],
    }];
    let mut gene_index = HashMap::new();
    gene_index.insert("G1".to_string(), 2usize);
    gene_index.insert("G3".to_string(), 0usize);

    let collection = GenesetCollection {
        version: "v1".to_string(),
        defs,
        resolved: Vec::new(),
    };

    let resolved = resolve_collection(collection, &gene_index);
    let gs = &resolved.resolved[0];
    assert_eq!(gs.gene_ids, vec![0, 2]);
    assert_eq!(gs.missing, vec!["G2".to_string()]);
    assert_eq!(gs.total, 3);
}

#[test]
fn geneset_duplicate_override() {
    let builtin = vec![
        GenesetDef {
            id: "set_a".to_string(),
            axis: 'A',
            genes: vec!["G1".to_string()],
        },
        GenesetDef {
            id: "set_b".to_string(),
            axis: 'B',
            genes: vec!["G2".to_string()],
        },
    ];
    let user = vec![GenesetDef {
        id: "set_a".to_string(),
        axis: 'A',
        genes: vec!["GX".to_string()],
    }];

    let merged = merge_defs(builtin, user);
    assert_eq!(merged[0].genes, vec!["GX".to_string()]);
    assert_eq!(merged[1].id, "set_b");
}
