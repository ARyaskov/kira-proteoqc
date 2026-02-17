#[cfg(feature = "fusion")]
mod fusion_tests {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    use kira_proteoqc::ctx::Ctx;
    use kira_proteoqc::expr::layout::{ExprHeaderV1, LAYOUT_CSC, VERSION, write_header};
    use kira_proteoqc::expr::reader::open_mmap;
    use kira_proteoqc::geneset::{GenesetCollection, ResolvedGeneset};
    use kira_proteoqc::schema::v1::Mode;
    use kira_proteoqc::scores::axis_raw::compute_axis_raw_with_mode;
    use tempfile::TempDir;

    fn write_expr(path: &std::path::Path) {
        let header = ExprHeaderV1 {
            version: VERSION,
            n_genes: 2,
            n_cells: 3,
            nnz: 3,
            layout: LAYOUT_CSC,
        };
        let gene_ptr: Vec<u64> = vec![0, 2, 3];
        let cell_idx: Vec<u32> = vec![0, 2, 1];
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

    fn make_genesets() -> GenesetCollection {
        let resolved = vec![
            ResolvedGeneset {
                id: "proteasome_core".into(),
                axis: 'A',
                gene_ids: vec![0],
                missing: vec![],
                total: 1,
            },
            ResolvedGeneset {
                id: "proteasome_regulator".into(),
                axis: 'B',
                gene_ids: vec![1],
                missing: vec![],
                total: 1,
            },
            ResolvedGeneset {
                id: "ubiquitin_axis".into(),
                axis: 'D',
                gene_ids: vec![0],
                missing: vec![],
                total: 1,
            },
            ResolvedGeneset {
                id: "e3_ligases".into(),
                axis: 'D',
                gene_ids: vec![1],
                missing: vec![],
                total: 1,
            },
            ResolvedGeneset {
                id: "dubs".into(),
                axis: 'D',
                gene_ids: vec![0],
                missing: vec![],
                total: 1,
            },
            ResolvedGeneset {
                id: "chaperone_hsp70".into(),
                axis: 'E',
                gene_ids: vec![0],
                missing: vec![],
                total: 1,
            },
            ResolvedGeneset {
                id: "chaperone_hsp90".into(),
                axis: 'E',
                gene_ids: vec![1],
                missing: vec![],
                total: 1,
            },
            ResolvedGeneset {
                id: "chaperone_hsp40".into(),
                axis: 'E',
                gene_ids: vec![0],
                missing: vec![],
                total: 1,
            },
            ResolvedGeneset {
                id: "erad".into(),
                axis: 'F',
                gene_ids: vec![1],
                missing: vec![],
                total: 1,
            },
            ResolvedGeneset {
                id: "ribosome_load".into(),
                axis: 'F',
                gene_ids: vec![0],
                missing: vec![],
                total: 1,
            },
        ];
        GenesetCollection {
            version: "v1".into(),
            defs: vec![],
            resolved,
        }
    }

    #[test]
    fn fusion_equivalence() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("expr.bin");
        write_expr(&path);
        let (header, mmap) = open_mmap(&path).unwrap();

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
        ctx.expr_header = Some(header);
        ctx.expr_mmap = Some(mmap);
        ctx.genes = vec!["G1".into(), "G2".into()];
        ctx.cells = vec!["C1".into(), "C2".into(), "C3".into()];
        ctx.nnz = 3;
        ctx.genesets = Some(make_genesets());

        ctx.fusion = "off".to_string();
        let axis_scalar = compute_axis_raw_with_mode(&mut ctx, Mode::Cell).unwrap();
        ctx.fusion = "proteo".to_string();
        let axis_fused = compute_axis_raw_with_mode(&mut ctx, Mode::Cell).unwrap();
        assert_eq!(axis_scalar.pcs, axis_fused.pcs);
        assert_eq!(axis_scalar.utp, axis_fused.utp);
        assert_eq!(axis_scalar.cls, axis_fused.cls);
        assert_eq!(axis_scalar.erad, axis_fused.erad);
        assert_eq!(axis_scalar.ribo, axis_fused.ribo);
    }
}
