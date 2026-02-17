use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::Context;
use memmap2::Mmap;

use crate::expr::layout::ExprHeaderV1;
use crate::expr::reader;
use crate::expr::reader::ExprReader;
use crate::geneset::GenesetCollection;
use crate::schema::v1::{Mode, ProteoQcV1};
use crate::scores::{AxisRawScores, IntegratedScores, PfsContributions, RiskFlag};
use crate::scores::{TimecourseResult, TimepointSummary};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunMode {
    Standalone,
    Pipeline,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    Mtx10x,
    H5ad,
}

impl InputFormat {
    pub fn detect(path: &Path) -> Self {
        if path.extension().and_then(|s| s.to_str()) == Some("h5ad") {
            Self::H5ad
        } else {
            Self::Mtx10x
        }
    }
}

#[derive(Debug, Clone)]
pub struct InputMeta {
    pub genes: Option<u64>,
    pub cells: Option<u64>,
    pub nnz: Option<u64>,
}

#[derive(Debug, Clone)]
pub struct OutputPaths {
    pub out_dir: PathBuf,
    pub json_path: PathBuf,
    pub tsv_path: PathBuf,
}

#[derive(Debug)]
pub struct Ctx {
    pub input: PathBuf,
    pub input_format: InputFormat,
    pub mode: Mode,
    pub timecourse: bool,
    pub geneset_path: Option<PathBuf>,
    pub log1p: bool,
    pub write_json: bool,
    pub write_tsv: bool,
    pub threads: usize,
    pub cache_block: usize,
    pub prefetch: bool,
    pub fusion: String,
    pub run_mode: RunMode,
    pub cache_override: Option<PathBuf>,
    pub input_prefix: Option<String>,
    pub genes: Vec<String>,
    pub cells: Vec<String>,
    pub nnz: usize,
    pub gene_index: HashMap<String, usize>,
    pub warnings: Vec<String>,
    pub expr_path: PathBuf,
    pub mtx_matrix_path: Option<PathBuf>,
    pub mtx_features_path: Option<PathBuf>,
    pub mtx_barcodes_path: Option<PathBuf>,
    pub shared_cache_path: Option<PathBuf>,
    pub shared_cache_used: bool,
    pub expr_header: Option<ExprHeaderV1>,
    pub expr_mmap: Option<Mmap>,
    pub scratch_cell_buf: Vec<f32>,
    pub scratch_shards: Vec<f32>,
    pub genesets: Option<GenesetCollection>,
    pub axis_raw: Option<AxisRawScores>,
    pub integrated_scores: Option<IntegratedScores>,
    pub pfs_contributions: Option<PfsContributions>,
    pub risk_flags: Vec<RiskFlag>,
    pub timecourse_points: Vec<TimepointSummary>,
    pub timecourse_result: Option<TimecourseResult>,
    pub input_meta: InputMeta,
    pub output: OutputPaths,
    pub report: ProteoQcV1,
}

impl Ctx {
    pub fn new(
        input: PathBuf,
        out_dir: PathBuf,
        mode: Mode,
        timecourse: bool,
        geneset_path: Option<PathBuf>,
        log1p: bool,
        write_json: bool,
        write_tsv: bool,
        tool_version: &str,
    ) -> Self {
        let json_path = out_dir.join("proteoqc.json");
        let tsv_path = out_dir.join("proteoqc.tsv");
        let expr_path = out_dir.join("expr.bin");
        let report = ProteoQcV1::empty(tool_version, mode.clone(), timecourse, log1p);
        Self {
            input_format: InputFormat::detect(&input),
            input,
            mode,
            timecourse,
            geneset_path,
            log1p,
            write_json,
            write_tsv,
            threads: 0,
            cache_block: 4096,
            prefetch: false,
            fusion: "off".to_string(),
            run_mode: RunMode::Standalone,
            cache_override: None,
            input_prefix: None,
            genes: Vec::new(),
            cells: Vec::new(),
            nnz: 0,
            gene_index: HashMap::new(),
            warnings: Vec::new(),
            expr_path,
            mtx_matrix_path: None,
            mtx_features_path: None,
            mtx_barcodes_path: None,
            shared_cache_path: None,
            shared_cache_used: false,
            expr_header: None,
            expr_mmap: None,
            scratch_cell_buf: Vec::new(),
            scratch_shards: Vec::new(),
            genesets: None,
            axis_raw: None,
            integrated_scores: None,
            pfs_contributions: None,
            risk_flags: Vec::new(),
            timecourse_points: Vec::new(),
            timecourse_result: None,
            input_meta: InputMeta {
                genes: None,
                cells: None,
                nnz: None,
            },
            output: OutputPaths {
                out_dir,
                json_path,
                tsv_path,
            },
            report,
        }
    }

    pub fn gene_slice(&self, gene_id: usize) -> anyhow::Result<(&[u32], &[f32])> {
        let header = self.expr_header.as_ref().context("expr header missing")?;
        let mmap = self.expr_mmap.as_ref().context("expr mmap missing")?;
        if gene_id >= header.n_genes as usize {
            anyhow::bail!("gene_id out of range");
        }
        let gene_ptr = reader::gene_ptr_slice(mmap, header);
        let cell_idx = reader::cell_idx_slice(mmap, header);
        let values = reader::values_slice(mmap, header);
        let start = gene_ptr[gene_id] as usize;
        let end = gene_ptr[gene_id + 1] as usize;
        Ok((&cell_idx[start..end], &values[start..end]))
    }

    pub fn expr_reader(&self) -> anyhow::Result<ExprReader<'_>> {
        let header = self.expr_header.as_ref().context("expr header missing")?;
        let mmap = self.expr_mmap.as_ref().context("expr mmap missing")?;
        Ok(ExprReader::new(header, mmap))
    }
}
