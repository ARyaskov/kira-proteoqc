#![cfg(feature = "hdf5")]

use std::path::Path;

use assert_cmd::Command;
use hdf5::File;
use tempfile::TempDir;

fn write_h5ad(path: &Path) {
    let file = File::create(path).unwrap();

    let x = file.create_group("X").unwrap();
    let data: Vec<f32> = vec![1.0];
    let indices: Vec<u32> = vec![0];
    let indptr: Vec<u32> = vec![0, 1];

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

    let shape = vec![1u64, 1u64];
    let ds_shape = x
        .new_dataset::<u64>()
        .shape(shape.len())
        .create("shape")
        .unwrap();
    ds_shape.write(&shape).unwrap();

    let attr = x.new_attr::<String>().create("encoding-type").unwrap();
    attr.write_scalar("csr_matrix").unwrap();

    let var = file.create_group("var").unwrap();
    let var_ds = var
        .new_dataset::<String>()
        .shape(1)
        .create("gene_symbols")
        .unwrap();
    var_ds.write(&vec!["G1".to_string()]).unwrap();

    let obs = file.create_group("obs").unwrap();
    let obs_ds = obs
        .new_dataset::<String>()
        .shape(1)
        .create("_index")
        .unwrap();
    obs_ds.write(&vec!["C1".to_string()]).unwrap();
}

#[test]
fn validate_command_h5ad_ok() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("tiny.h5ad");
    write_h5ad(&path);

    let mut cmd = Command::cargo_bin("kira-proteoqc").unwrap();
    cmd.arg("validate").arg("--input").arg(&path);
    cmd.assert().success();
}
