use std::fs;
use std::path::Path;

use kira_proteoqc::ctx::Ctx;
use kira_proteoqc::expr::reader;
use kira_proteoqc::pipeline::Pipeline;
use kira_proteoqc::pipeline::stage1_input::Stage1Input;
use kira_proteoqc::pipeline::stage3_expr_cache::Stage3ExprCache;
use kira_proteoqc::schema::v1::Mode;
use tempfile::TempDir;

#[cfg(feature = "hdf5")]
use hdf5::File;
#[cfg(feature = "hdf5")]
use kira_proteoqc::pipeline::stage2_h5ad::Stage2H5ad;

fn write_10x(dir: &Path, features: &str, barcodes: &str, mtx: &str) {
    fs::write(dir.join("features.tsv"), features).unwrap();
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
    fs::write(dir.join("matrix.mtx"), mtx).unwrap();
}

#[test]
fn mtx_expr_cache_roundtrip() {
    let tmp = TempDir::new().unwrap();
    let mtx =
        "%%MatrixMarket matrix coordinate integer general\n3 2 4\n1 1 1\n1 2 2\n3 1 3\n2 2 4\n";
    let features = "g1\tG1\ng2\tG2\ng3\tG3\n";
    let barcodes = "c1\nc2\n";
    write_10x(tmp.path(), features, barcodes, mtx);

    let mut ctx = Ctx::new(
        tmp.path().to_path_buf(),
        tmp.path().join("out"),
        Mode::Cell,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );

    let pipeline = Pipeline::new(vec![
        Box::new(Stage1Input::new()),
        Box::new(Stage3ExprCache::new()),
    ]);
    pipeline.run(&mut ctx).unwrap();

    let (header, mmap) = reader::open_mmap(&ctx.expr_path).unwrap();
    assert_eq!(header.n_genes, 3);
    assert_eq!(header.n_cells, 2);
    assert_eq!(header.nnz, 4);

    let gene_ptr = reader::gene_ptr_slice(&mmap, &header);
    let cell_idx = reader::cell_idx_slice(&mmap, &header);
    let values = reader::values_slice(&mmap, &header);

    assert_eq!(gene_ptr, &[0, 2, 3, 4]);
    assert_eq!(cell_idx, &[0, 1, 1, 0]);
    assert_eq!(values, &[1.0, 2.0, 4.0, 3.0]);
}

#[cfg(feature = "hdf5")]
fn write_h5ad_csc(path: &Path, genes: &[&str], cells: &[&str]) {
    let file = File::create(path).unwrap();

    let x = file.create_group("X").unwrap();
    let data: Vec<f32> = vec![1.0, 2.0, 3.0];
    let indices: Vec<u32> = vec![0, 1, 0];
    let indptr: Vec<u32> = vec![0, 1, 3];

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
    attr.write_scalar("csc_matrix").unwrap();

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

#[cfg(feature = "hdf5")]
#[test]
fn h5ad_csc_passthrough() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("csc.h5ad");
    write_h5ad_csc(&path, &["G1", "G2"], &["C1", "C2", "C3"]);

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

    let pipeline = Pipeline::new(vec![
        Box::new(Stage2H5ad::new()),
        Box::new(Stage3ExprCache::new()),
    ]);
    pipeline.run(&mut ctx).unwrap();

    let (header, mmap) = reader::open_mmap(&ctx.expr_path).unwrap();
    let gene_ptr = reader::gene_ptr_slice(&mmap, &header);
    let cell_idx = reader::cell_idx_slice(&mmap, &header);
    let values = reader::values_slice(&mmap, &header);

    assert_eq!(gene_ptr, &[0, 1, 3]);
    assert_eq!(cell_idx, &[0, 0, 1]);
    assert_eq!(values, &[1.0, 3.0, 2.0]);
}

#[cfg(feature = "hdf5")]
fn write_h5ad_csr_cells(path: &Path, genes: &[&str], cells: &[&str]) {
    let file = File::create(path).unwrap();

    let x = file.create_group("X").unwrap();
    let data: Vec<f32> = vec![5.0, 6.0];
    let indices: Vec<u32> = vec![0, 1];
    let indptr: Vec<u32> = vec![0, 1, 2, 2];

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
    attr.write_scalar("csr_matrix").unwrap();

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

#[cfg(feature = "hdf5")]
#[test]
fn h5ad_csr_transpose() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("csr.h5ad");
    write_h5ad_csr_cells(&path, &["G1", "G2"], &["C1", "C2", "C3"]);

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

    let pipeline = Pipeline::new(vec![
        Box::new(Stage2H5ad::new()),
        Box::new(Stage3ExprCache::new()),
    ]);
    pipeline.run(&mut ctx).unwrap();

    let (header, mmap) = reader::open_mmap(&ctx.expr_path).unwrap();
    let gene_ptr = reader::gene_ptr_slice(&mmap, &header);
    let cell_idx = reader::cell_idx_slice(&mmap, &header);
    let values = reader::values_slice(&mmap, &header);

    assert_eq!(gene_ptr, &[0, 1, 2]);
    assert_eq!(cell_idx, &[0, 1]);
    assert_eq!(values, &[5.0, 6.0]);
}

#[test]
fn cache_reuse_readonly() {
    let tmp = TempDir::new().unwrap();
    let mtx = "%%MatrixMarket matrix coordinate integer general\n2 2 1\n1 1 1\n";
    let features = "g1\tG1\ng2\tG2\n";
    let barcodes = "c1\nc2\n";
    write_10x(tmp.path(), features, barcodes, mtx);

    let mut ctx1 = Ctx::new(
        tmp.path().to_path_buf(),
        tmp.path().join("out"),
        Mode::Cell,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );
    let pipeline = Pipeline::new(vec![
        Box::new(Stage1Input::new()),
        Box::new(Stage3ExprCache::new()),
    ]);
    pipeline.run(&mut ctx1).unwrap();

    let expr_path = ctx1.expr_path.clone();
    let mut perms = fs::metadata(&expr_path).unwrap().permissions();
    perms.set_readonly(true);
    fs::set_permissions(&expr_path, perms).unwrap();

    let mut ctx2 = Ctx::new(
        tmp.path().to_path_buf(),
        tmp.path().join("out"),
        Mode::Cell,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );
    let pipeline = Pipeline::new(vec![
        Box::new(Stage1Input::new()),
        Box::new(Stage3ExprCache::new()),
    ]);
    pipeline.run(&mut ctx2).unwrap();
}
