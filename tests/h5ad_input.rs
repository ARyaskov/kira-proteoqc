#![cfg(feature = "hdf5")]

use std::path::Path;

use hdf5::File;
use kira_proteoqc::ctx::Ctx;
use kira_proteoqc::pipeline::Pipeline;
use kira_proteoqc::pipeline::stage2_h5ad::Stage2H5ad;
use kira_proteoqc::schema::v1::Mode;
use tempfile::TempDir;

fn write_h5ad(path: &Path, genes: &[&str], cells: &[&str], encoding: &str) {
    let file = File::create(path).unwrap();

    let x = file.create_group("X").unwrap();
    let data: Vec<f32> = vec![1.0, 2.0];
    let indices: Vec<u32> = vec![0, 2];
    let indptr: Vec<u32> = vec![0, 1, 2];

    let ds_data = x
        .new_dataset::<f32>()
        .shape(data.len())
        .create("data")
        .unwrap();
    ds_data.write(&data).unwrap();

    let ds_indices = x
        .new_dataset::<u32>()
        .shape(indices.len())
        .create("indices")
        .unwrap();
    ds_indices.write(&indices).unwrap();

    let ds_indptr = x
        .new_dataset::<u32>()
        .shape(indptr.len())
        .create("indptr")
        .unwrap();
    ds_indptr.write(&indptr).unwrap();

    let shape = vec![genes.len() as u64, cells.len() as u64];
    let ds_shape = x
        .new_dataset::<u64>()
        .shape(shape.len())
        .create("shape")
        .unwrap();
    ds_shape.write(&shape).unwrap();

    let attr = x.new_attr::<String>().create("encoding-type").unwrap();
    attr.write_scalar(encoding).unwrap();

    let var = file.create_group("var").unwrap();
    let var_ds = var
        .new_dataset::<String>()
        .shape(genes.len())
        .create("gene_symbols")
        .unwrap();
    let genes_vec: Vec<String> = genes.iter().map(|s| s.to_string()).collect();
    var_ds.write(&genes_vec).unwrap();

    let obs = file.create_group("obs").unwrap();
    let obs_ds = obs
        .new_dataset::<String>()
        .shape(cells.len())
        .create("_index")
        .unwrap();
    let cells_vec: Vec<String> = cells.iter().map(|s| s.to_string()).collect();
    obs_ds.write(&cells_vec).unwrap();
}

#[test]
fn stage2_reads_h5ad_summary() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("tiny.h5ad");
    write_h5ad(&path, &["G1", "G2"], &["C1", "C2", "C3"], "csr_matrix");

    let mut ctx = Ctx::new(
        path.clone(),
        tmp.path().join("out"),
        Mode::Cell,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );

    let pipeline = Pipeline::new(vec![Box::new(Stage2H5ad::new())]);
    pipeline.run(&mut ctx).unwrap();

    assert_eq!(ctx.genes.len(), 2);
    assert_eq!(ctx.cells.len(), 3);
    assert_eq!(ctx.nnz, 2);
}

#[test]
fn stage2_duplicate_gene_symbol_warns() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("dup.h5ad");
    write_h5ad(&path, &["G1", "G1"], &["C1"], "csr_matrix");

    let mut ctx = Ctx::new(
        path.clone(),
        tmp.path().join("out"),
        Mode::Cell,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );

    let pipeline = Pipeline::new(vec![Box::new(Stage2H5ad::new())]);
    pipeline.run(&mut ctx).unwrap();

    assert_eq!(ctx.gene_index.get("G1"), Some(&0));
    assert_eq!(ctx.warnings.len(), 1);
}
