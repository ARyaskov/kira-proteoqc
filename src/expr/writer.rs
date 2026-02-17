use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use tracing::info;

use crate::ctx::{Ctx, InputFormat};
use crate::expr::layout::{ExprHeaderV1, LAYOUT_CSC, VERSION, write_header};
use crate::expr::reader;
use crate::io::{h5ad, mtx, shared_cache};

pub fn ensure_expr_cache(ctx: &mut Ctx) -> Result<()> {
    let path = &ctx.expr_path;
    if path.exists() {
        if let Ok((header, mmap)) = reader::open_mmap(path) {
            if header.n_genes as usize == ctx.genes.len()
                && header.n_cells as usize == ctx.cells.len()
                && header.nnz as usize == ctx.nnz
            {
                info!(expr = %path.display(), "expr_cache_reuse");
                ctx.expr_header = Some(header);
                ctx.expr_mmap = Some(mmap);
                return Ok(());
            }
        }
        info!(expr = %path.display(), "expr_cache_rebuild");
    }

    build_expr_cache(ctx)?;
    let (header, mmap) = reader::open_mmap(path)?;
    ctx.expr_header = Some(header);
    ctx.expr_mmap = Some(mmap);
    Ok(())
}

fn build_expr_cache(ctx: &Ctx) -> Result<()> {
    fs::create_dir_all(&ctx.output.out_dir)?;
    if let Some(cache_path) = &ctx.shared_cache_path {
        return build_from_shared_cache(ctx, cache_path);
    }
    match ctx.input_format {
        InputFormat::Mtx10x => build_from_mtx(ctx),
        InputFormat::H5ad => build_from_h5ad(ctx),
    }
}

fn build_from_mtx(ctx: &Ctx) -> Result<()> {
    let path = if let Some(path) = &ctx.mtx_matrix_path {
        path.clone()
    } else {
        resolve_mtx_path(&ctx.input)?
    };
    let n_genes = ctx.genes.len();
    let n_cells = ctx.cells.len();

    let mut per_gene: Vec<Vec<(u32, f32)>> = vec![Vec::new(); n_genes];
    let summary = mtx::for_each_entry(&path, |row, col, value| {
        if row >= n_genes || col >= n_cells {
            bail!("MTX index out of bounds after context load");
        }
        per_gene[row].push((col as u32, value));
        Ok(())
    })?;

    if summary.nrows != n_genes || summary.ncols != n_cells {
        bail!(
            "MTX dimensions ({}, {}) do not match context ({}, {})",
            summary.nrows,
            summary.ncols,
            n_genes,
            n_cells
        );
    }
    if ctx.nnz != summary.nnz {
        bail!(
            "MTX nnz ({}) does not match context nnz ({})",
            summary.nnz,
            ctx.nnz
        );
    }

    let nnz = summary.nnz;
    let mut gene_ptr = Vec::with_capacity(n_genes + 1);
    gene_ptr.push(0u64);

    for entries in &mut per_gene {
        entries.sort_by_key(|(cell, _)| *cell);
        let next = gene_ptr.last().copied().unwrap() + entries.len() as u64;
        gene_ptr.push(next);
    }
    if *gene_ptr.last().unwrap_or(&0) as usize != nnz {
        bail!("gene_ptr does not match nnz");
    }

    let header = ExprHeaderV1 {
        version: VERSION,
        n_genes: n_genes as u32,
        n_cells: n_cells as u32,
        nnz: nnz as u64,
        layout: LAYOUT_CSC,
    };

    write_expr_file(&ctx.expr_path, &header, &gene_ptr, &per_gene)?;
    Ok(())
}

fn build_from_shared_cache(ctx: &Ctx, cache_path: &Path) -> Result<()> {
    let cache = shared_cache::SharedCache::open(cache_path).with_context(|| {
        format!(
            "failed to read shared cache while building expr cache: {}",
            cache_path.display()
        )
    })?;
    let n_genes = cache.header.n_genes as usize;
    let n_cells = cache.header.n_cells as usize;
    let nnz = cache.header.nnz as usize;

    if n_genes != ctx.genes.len() || n_cells != ctx.cells.len() || nnz != ctx.nnz {
        bail!(
            "shared cache dimensions ({}, {}, {}) do not match context ({}, {}, {})",
            n_genes,
            n_cells,
            nnz,
            ctx.genes.len(),
            ctx.cells.len(),
            ctx.nnz
        );
    }

    // Shared cache is CSC by cell (col_ptr per cell, row_idx = gene),
    // while expr.bin is CSC by gene. Build gene-major arrays deterministically.
    let mut gene_ptr = vec![0u64; n_genes + 1];
    let col_ptr = cache.col_ptr()?;
    for cell in 0..n_cells {
        let start = col_ptr[cell] as usize;
        let end = col_ptr[cell + 1] as usize;
        for k in start..end {
            let g = cache.row_idx_at(k)? as usize;
            gene_ptr[g + 1] += 1;
        }
    }
    for g in 0..n_genes {
        gene_ptr[g + 1] += gene_ptr[g];
    }

    let mut cell_idx = vec![0u32; nnz];
    let mut values = vec![0f32; nnz];
    let mut offsets = gene_ptr.clone();
    for cell in 0..n_cells {
        let start = col_ptr[cell] as usize;
        let end = col_ptr[cell + 1] as usize;
        for k in start..end {
            let g = cache.row_idx_at(k)? as usize;
            let pos = offsets[g] as usize;
            cell_idx[pos] = cell as u32;
            values[pos] = cache.value_u32_at(k)? as f32;
            offsets[g] += 1;
        }
    }

    let header = ExprHeaderV1 {
        version: VERSION,
        n_genes: n_genes as u32,
        n_cells: n_cells as u32,
        nnz: nnz as u64,
        layout: LAYOUT_CSC,
    };
    write_expr_file_flat(&ctx.expr_path, &header, &gene_ptr, &cell_idx, &values)?;
    Ok(())
}

fn build_from_h5ad(ctx: &Ctx) -> Result<()> {
    let (summary, encoding, indptr_len, indices, data, indptr) =
        h5ad::read_h5ad_sparse(&ctx.input)?;

    if summary.nrows != ctx.genes.len() || summary.ncols != ctx.cells.len() {
        bail!(
            "H5AD shape ({}, {}) does not match context ({}, {})",
            summary.nrows,
            summary.ncols,
            ctx.genes.len(),
            ctx.cells.len()
        );
    }

    let n_genes = ctx.genes.len();
    let n_cells = ctx.cells.len();
    let nnz = summary.nnz;
    if ctx.nnz != nnz {
        bail!(
            "H5AD nnz ({}) does not match context nnz ({})",
            nnz,
            ctx.nnz
        );
    }
    if indptr.is_empty() || *indptr.last().unwrap() as usize != nnz {
        bail!("H5AD indptr last value does not equal nnz");
    }

    let header = ExprHeaderV1 {
        version: VERSION,
        n_genes: n_genes as u32,
        n_cells: n_cells as u32,
        nnz: nnz as u64,
        layout: LAYOUT_CSC,
    };

    let mut gene_ptr = Vec::with_capacity(n_genes + 1);
    let mut cell_idx: Vec<u32> = vec![0; nnz];
    let mut values: Vec<f32> = vec![0.0; nnz];

    if indptr_len == n_genes + 1 {
        gene_ptr.extend(indptr.iter().copied());
        if *gene_ptr.last().unwrap_or(&0) as usize != nnz {
            bail!("H5AD indptr last value does not equal nnz");
        }

        for &cell in &indices {
            if cell as usize >= n_cells {
                bail!("cell index out of bounds in indices");
            }
        }
        cell_idx.copy_from_slice(&indices);
        values.copy_from_slice(&data);
        sort_gene_slices(&gene_ptr, &mut cell_idx, &mut values);
    } else if indptr_len == n_cells + 1 {
        info!(
            encoding = %encoding,
            "transposing sparse layout to gene-major"
        );
        gene_ptr = vec![0u64; n_genes + 1];
        for &gene in &indices {
            let g = gene as usize;
            if g >= n_genes {
                bail!("gene index out of bounds in indices");
            }
            gene_ptr[g + 1] += 1;
        }
        for i in 0..n_genes {
            gene_ptr[i + 1] += gene_ptr[i];
        }
        let mut offsets: Vec<u64> = gene_ptr.clone();
        for cell in 0..n_cells {
            let start = indptr[cell] as usize;
            let end = indptr[cell + 1] as usize;
            for idx in start..end {
                let gene = indices[idx] as usize;
                let pos = offsets[gene] as usize;
                cell_idx[pos] = cell as u32;
                values[pos] = data[idx];
                offsets[gene] += 1;
            }
        }
    } else {
        bail!("H5AD indptr length does not match genes or cells");
    }

    write_expr_file_flat(&ctx.expr_path, &header, &gene_ptr, &cell_idx, &values)?;
    Ok(())
}

fn write_expr_file(
    path: &Path,
    header: &ExprHeaderV1,
    gene_ptr: &[u64],
    per_gene: &[Vec<(u32, f32)>],
) -> Result<()> {
    let file = File::create(path).context("failed to create expr.bin")?;
    let mut writer = BufWriter::new(file);

    write_header(&mut writer, header)?;

    for v in gene_ptr {
        writer.write_all(&v.to_le_bytes())?;
    }

    for entries in per_gene {
        for (cell, _) in entries {
            writer.write_all(&cell.to_le_bytes())?;
        }
    }

    for entries in per_gene {
        for (_, value) in entries {
            writer.write_all(&value.to_le_bytes())?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn write_expr_file_flat(
    path: &Path,
    header: &ExprHeaderV1,
    gene_ptr: &[u64],
    cell_idx: &[u32],
    values: &[f32],
) -> Result<()> {
    let file = File::create(path).context("failed to create expr.bin")?;
    let mut writer = BufWriter::new(file);

    write_header(&mut writer, header)?;
    for v in gene_ptr {
        writer.write_all(&v.to_le_bytes())?;
    }
    for v in cell_idx {
        writer.write_all(&v.to_le_bytes())?;
    }
    for v in values {
        writer.write_all(&v.to_le_bytes())?;
    }
    writer.flush()?;
    Ok(())
}

fn sort_gene_slices(gene_ptr: &[u64], cell_idx: &mut [u32], values: &mut [f32]) {
    let n_genes = gene_ptr.len() - 1;
    for g in 0..n_genes {
        let start = gene_ptr[g] as usize;
        let end = gene_ptr[g + 1] as usize;
        if end - start <= 1 {
            continue;
        }
        let mut pairs: Vec<(u32, f32)> = cell_idx[start..end]
            .iter()
            .copied()
            .zip(values[start..end].iter().copied())
            .collect();
        pairs.sort_by_key(|(cell, _)| *cell);
        for (i, (cell, val)) in pairs.into_iter().enumerate() {
            cell_idx[start + i] = cell;
            values[start + i] = val;
        }
    }
}

fn resolve_mtx_path(input_dir: &Path) -> Result<std::path::PathBuf> {
    let plain = input_dir.join("matrix.mtx");
    if plain.exists() {
        return Ok(plain);
    }
    let gz = input_dir.join("matrix.mtx.gz");
    if gz.exists() {
        return Ok(gz);
    }
    bail!("missing matrix.mtx or matrix.mtx.gz")
}
