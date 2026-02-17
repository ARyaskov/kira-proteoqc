use std::path::Path;

use anyhow::{Context, Result, bail};
use hdf5::{File, Group};

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

pub fn read_h5ad_summary(path: &Path) -> Result<H5adSummary> {
    let file = File::open(path).context("failed to open H5AD file")?;
    let x = file.group("X").context("missing group: X")?;

    let encoding = read_encoding_type(&x)?;
    if encoding != "csr_matrix" && encoding != "csc_matrix" {
        bail!(
            "unsupported X encoding-type '{}'; supported: csr_matrix, csc_matrix",
            encoding
        );
    }

    let meta = read_sparse_meta(&x)?;

    let (genes, mut warnings) = read_genes(&file, meta.nrows)?;
    let cells = read_cells(&file, meta.ncols, &mut warnings)?;

    Ok(H5adSummary {
        genes,
        cells,
        nnz: meta.nnz,
        nrows: meta.nrows,
        ncols: meta.ncols,
        warnings,
    })
}

pub fn read_h5ad_sparse(
    path: &Path,
) -> Result<(H5adSparseMeta, String, usize, Vec<u32>, Vec<f32>, Vec<u64>)> {
    let file = File::open(path).context("failed to open H5AD file")?;
    let x = file.group("X").context("missing group: X")?;

    let encoding = read_encoding_type(&x)?;
    if encoding != "csr_matrix" && encoding != "csc_matrix" {
        bail!(
            "unsupported X encoding-type '{}'; supported: csr_matrix, csc_matrix",
            encoding
        );
    }

    let meta = read_sparse_meta(&x)?;
    let data_ds = x.dataset("data").context("missing X/data")?;
    let indices_ds = x.dataset("indices").context("missing X/indices")?;
    let indptr_ds = x.dataset("indptr").context("missing X/indptr")?;

    let data = read_f32_vec(&data_ds)?;
    let indices = read_u32_vec(&indices_ds)?;
    let indptr = read_u64_vec(&indptr_ds)?;

    let nnz = data.len();
    if indices.len() != nnz {
        bail!(
            "X/indices length ({}) does not match X/data length ({})",
            indices.len(),
            nnz
        );
    }

    Ok((meta, encoding, indptr.len(), indices, data, indptr))
}

fn read_sparse_meta(x: &Group) -> Result<H5adSparseMeta> {
    let shape = read_shape(x)?;
    let nrows = shape[0];
    let ncols = shape[1];

    let data_ds = x.dataset("data").context("missing X/data")?;
    let indices_ds = x.dataset("indices").context("missing X/indices")?;
    let indptr_ds = x.dataset("indptr").context("missing X/indptr")?;

    let nnz = data_ds.size();
    let indices_len = indices_ds.size();
    let indptr_len = indptr_ds.size();

    if indices_len != nnz {
        bail!(
            "X/indices length ({}) does not match X/data length ({})",
            indices_len,
            nnz
        );
    }

    if indptr_len != nrows + 1 && indptr_len != ncols + 1 {
        bail!(
            "X/indptr length ({}) does not match nrows+1 ({}) or ncols+1 ({})",
            indptr_len,
            nrows + 1,
            ncols + 1
        );
    }

    Ok(H5adSparseMeta { nrows, ncols, nnz })
}

fn read_encoding_type(x: &Group) -> Result<String> {
    let attr = x
        .attr("encoding-type")
        .context("missing X attribute: encoding-type")?;
    let value: String = attr.read_scalar().context("failed to read encoding-type")?;
    Ok(value)
}

fn read_shape(x: &Group) -> Result<[usize; 2]> {
    let shape_ds = x.dataset("shape").context("missing X/shape")?;
    let shape: Vec<u64> = shape_ds.read_1d().context("failed to read X/shape")?;
    if shape.len() != 2 {
        bail!("X/shape must have length 2");
    }
    Ok([shape[0] as usize, shape[1] as usize])
}

fn read_genes(file: &File, nrows: usize) -> Result<(Vec<String>, Vec<String>)> {
    let var = file.group("var").context("missing group: var")?;
    let mut warnings = Vec::new();

    let genes = if let Ok(ds) = var.dataset("gene_symbols") {
        read_string_vec(&ds).context("failed to read var/gene_symbols")?
    } else if let Ok(ds) = var.dataset("gene_names") {
        warnings.push("var/gene_symbols missing; using var/gene_names".to_string());
        read_string_vec(&ds).context("failed to read var/gene_names")?
    } else if let Ok(ds) = var.dataset("_index") {
        warnings.push("var/gene_symbols missing; using var/_index".to_string());
        read_string_vec(&ds).context("failed to read var/_index")?
    } else {
        bail!("missing var/gene_symbols, var/gene_names, and var/_index");
    };

    if genes.len() != nrows {
        bail!(
            "var genes length ({}) does not match X rows ({})",
            genes.len(),
            nrows
        );
    }

    Ok((genes, warnings))
}

fn read_cells(file: &File, ncols: usize, warnings: &mut Vec<String>) -> Result<Vec<String>> {
    let obs = file.group("obs").context("missing group: obs")?;
    let ds = obs.dataset("_index").context("missing obs/_index")?;
    let cells = read_string_vec(&ds).context("failed to read obs/_index")?;

    if cells.len() != ncols {
        bail!(
            "obs cells length ({}) does not match X cols ({})",
            cells.len(),
            ncols
        );
    }

    if cells.is_empty() {
        warnings.push("obs/_index is empty".to_string());
    }

    Ok(cells)
}

fn read_string_vec(ds: &hdf5::Dataset) -> Result<Vec<String>> {
    let values: Vec<String> = ds.read_1d().context("failed to read string vector")?;
    Ok(values)
}

fn read_f32_vec(ds: &hdf5::Dataset) -> Result<Vec<f32>> {
    if let Ok(v) = ds.read_1d::<f32>() {
        return Ok(v);
    }
    if let Ok(v) = ds.read_1d::<f64>() {
        return Ok(v.into_iter().map(|x| x as f32).collect());
    }
    if let Ok(v) = ds.read_1d::<i32>() {
        return Ok(v.into_iter().map(|x| x as f32).collect());
    }
    if let Ok(v) = ds.read_1d::<i64>() {
        return Ok(v.into_iter().map(|x| x as f32).collect());
    }
    bail!("unsupported data type for X/data")
}

fn read_u32_vec(ds: &hdf5::Dataset) -> Result<Vec<u32>> {
    if let Ok(v) = ds.read_1d::<u32>() {
        return Ok(v);
    }
    if let Ok(v) = ds.read_1d::<u64>() {
        let mut out = Vec::with_capacity(v.len());
        for x in v {
            if x > u32::MAX as u64 {
                bail!("index exceeds u32 max");
            }
            out.push(x as u32);
        }
        return Ok(out);
    }
    if let Ok(v) = ds.read_1d::<i64>() {
        let mut out = Vec::with_capacity(v.len());
        for x in v {
            if x < 0 || x > u32::MAX as i64 {
                bail!("index out of range");
            }
            out.push(x as u32);
        }
        return Ok(out);
    }
    bail!("unsupported index type for X/indices")
}

fn read_u64_vec(ds: &hdf5::Dataset) -> Result<Vec<u64>> {
    if let Ok(v) = ds.read_1d::<u64>() {
        return Ok(v);
    }
    if let Ok(v) = ds.read_1d::<u32>() {
        return Ok(v.into_iter().map(|x| x as u64).collect());
    }
    if let Ok(v) = ds.read_1d::<i64>() {
        let mut out = Vec::with_capacity(v.len());
        for x in v {
            if x < 0 {
                bail!("indptr contains negative value");
            }
            out.push(x as u64);
        }
        return Ok(out);
    }
    bail!("unsupported index type for X/indptr")
}
