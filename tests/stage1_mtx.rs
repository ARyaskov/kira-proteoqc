use std::fs;
use std::path::Path;

use kira_proteoqc::ctx::Ctx;
use kira_proteoqc::pipeline::Pipeline;
use kira_proteoqc::pipeline::stage1_input::Stage1Input;
use kira_proteoqc::schema::v1::Mode;
use tempfile::TempDir;

fn write_10x(dir: &Path, features: &str, barcodes: &str, mtx: &str) {
    fs::write(dir.join("features.tsv"), features).unwrap();
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
    fs::write(dir.join("matrix.mtx"), mtx).unwrap();
}

#[test]
fn stage1_reads_mtx_summary() {
    let tmp = TempDir::new().unwrap();
    let mtx =
        "%%MatrixMarket matrix coordinate integer general\n%\n3 2 4\n1 1 1\n2 1 2\n3 2 3\n1 2 4\n";
    let features = "g1\tGeneA\n g2\tGeneB\n g3\tGeneC\n".replace(" ", "");
    let barcodes = "cell1\ncell2\n";
    write_10x(tmp.path(), &features, barcodes, mtx);

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

    let pipeline = Pipeline::new(vec![Box::new(Stage1Input::new())]);
    pipeline.run(&mut ctx).unwrap();

    assert_eq!(ctx.genes.len(), 3);
    assert_eq!(ctx.cells.len(), 2);
    assert_eq!(ctx.nnz, 4);
}

#[test]
fn stage1_duplicate_gene_symbol_warns() {
    let tmp = TempDir::new().unwrap();
    let mtx = "%%MatrixMarket matrix coordinate integer general\n2 1 1\n1 1 1\n";
    let features = "g1\tGeneX\ng2\tGeneX\n";
    let barcodes = "cell1\n";
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

    let pipeline = Pipeline::new(vec![Box::new(Stage1Input::new())]);
    pipeline.run(&mut ctx).unwrap();

    assert_eq!(ctx.gene_index.get("GeneX"), Some(&0));
    assert_eq!(ctx.warnings.len(), 1);
}
