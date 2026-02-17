use std::fs::File;
use std::io::{BufWriter, Read};
use std::path::Path;

use anyhow::Result;
use flate2::read::GzDecoder;

use crate::schema::v1::ProteoQcV1;

pub mod barcodes;
pub mod features;
#[cfg(feature = "hdf5")]
pub mod h5ad;
#[cfg(not(feature = "hdf5"))]
pub mod h5ad {
    use anyhow::{Result, bail};
    use std::path::Path;

    #[derive(Debug)]
    pub struct H5adSummary {
        pub genes: Vec<String>,
        pub cells: Vec<String>,
        pub nnz: usize,
        pub nrows: usize,
        pub ncols: usize,
        pub warnings: Vec<String>,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct H5adSparseMeta {
        pub nrows: usize,
        pub ncols: usize,
        pub nnz: usize,
    }

    pub fn read_h5ad_summary(_path: &Path) -> Result<H5adSummary> {
        bail!("H5AD support not enabled. Rebuild with --features hdf5");
    }

    pub fn read_h5ad_sparse(
        _path: &Path,
    ) -> Result<(H5adSparseMeta, String, usize, Vec<u32>, Vec<f32>, Vec<u64>)> {
        bail!("H5AD support not enabled. Rebuild with --features hdf5");
    }
}
pub mod json_writer;
pub mod mtx;
pub mod pipeline_output;
pub mod shared_cache;
pub mod summary;
pub mod tsv_writer;

pub fn write_json(path: &Path, report: &ProteoQcV1) -> Result<()> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, report)?;
    Ok(())
}

pub(crate) fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path)?;
    if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        let decoder = GzDecoder::new(file);
        Ok(Box::new(decoder))
    } else {
        Ok(Box::new(file))
    }
}
