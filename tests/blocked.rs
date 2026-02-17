#[cfg(feature = "cache-block")]
mod blocked_tests {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    use kira_proteoqc::expr::layout::{ExprHeaderV1, LAYOUT_CSC, VERSION, write_header};
    use kira_proteoqc::expr::reader::{ExprReader, open_mmap};
    use kira_proteoqc::math::reduce::GeneSetReducer;
    use kira_proteoqc::math::reduce_blocked::per_cell_raw_blocked;
    use tempfile::TempDir;

    fn write_expr(path: &std::path::Path) {
        let header = ExprHeaderV1 {
            version: VERSION,
            n_genes: 1,
            n_cells: 6,
            nnz: 3,
            layout: LAYOUT_CSC,
        };
        let gene_ptr: Vec<u64> = vec![0, 3];
        let cell_idx: Vec<u32> = vec![0, 3, 5];
        let values: Vec<f32> = vec![1.0, 2.0, 3.0];

        let file = File::create(path).unwrap();
        let mut w = BufWriter::new(file);
        write_header(&mut w, &header).unwrap();
        for v in &gene_ptr {
            w.write_all(&v.to_le_bytes()).unwrap();
        }
        for v in &cell_idx {
            w.write_all(&v.to_le_bytes()).unwrap();
        }
        for v in &values {
            w.write_all(&v.to_le_bytes()).unwrap();
        }
        w.flush().unwrap();
    }

    #[test]
    fn blocked_equivalence() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("expr.bin");
        write_expr(&path);

        let (header, mmap) = open_mmap(&path).unwrap();
        let reader = ExprReader::new(&header, &mmap);
        let mut out_block = vec![0.0f32; reader.n_cells()];
        let mut reducer = GeneSetReducer::new(&reader, 1, 0, false);
        let mut out_sc = vec![0.0f32; reader.n_cells()];
        per_cell_raw_blocked(&reader, &[0], &mut out_block, 4, false).unwrap();
        reducer.per_cell_raw(&[0], &mut out_sc).unwrap();
        assert_eq!(out_block, out_sc);
    }
}
