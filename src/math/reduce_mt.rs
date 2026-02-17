use anyhow::{Result, bail};

use crate::expr::reader::ExprReader;

#[cfg(feature = "mt")]
use rayon::prelude::*;

const SHARD_LEN: usize = 4096;

pub fn per_cell_raw_mt(
    expr: &ExprReader<'_>,
    genes: &[usize],
    out: &mut [f32],
    threads: usize,
    scratch_shards: &mut Vec<f32>,
) -> Result<()> {
    if out.len() != expr.n_cells() {
        bail!("output buffer length does not match n_cells");
    }
    for v in out.iter_mut() {
        *v = 0.0;
    }
    if genes.is_empty() {
        return Ok(());
    }

    let n_cells = expr.n_cells();
    let shard_len = SHARD_LEN;
    let shard_count = (n_cells + shard_len - 1) / shard_len;
    let total_len = shard_count * shard_len;
    if scratch_shards.len() != total_len {
        *scratch_shards = vec![0.0f32; total_len];
    } else {
        for v in scratch_shards.iter_mut() {
            *v = 0.0;
        }
    }

    #[cfg(feature = "mt")]
    {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.max(1))
            .build()
            .map_err(|e| anyhow::anyhow!("failed to build thread pool: {}", e))?;
        pool.install(|| {
            scratch_shards
                .par_chunks_mut(shard_len)
                .enumerate()
                .for_each(|(shard_id, shard)| {
                    let shard_start = shard_id * shard_len;
                    let shard_end = (shard_start + shard_len).min(n_cells);
                    for &gene_id in genes {
                        if let Ok((cells, values)) = expr.gene_slice(gene_id) {
                            let start = lower_bound(cells, shard_start as u32);
                            let end = lower_bound(cells, shard_end as u32);
                            for idx in start..end {
                                let cell = cells[idx] as usize - shard_start;
                                shard[cell] += values[idx];
                            }
                        }
                    }
                    let denom = genes.len() as f32;
                    for v in &mut shard[..(shard_end - shard_start)] {
                        *v /= denom;
                    }
                });
        });
    }

    #[cfg(not(feature = "mt"))]
    {
        let _ = threads;
        for shard_id in 0..shard_count {
            let shard_start = shard_id * shard_len;
            let shard_end = (shard_start + shard_len).min(n_cells);
            let shard = &mut scratch_shards[shard_id * shard_len..shard_id * shard_len + shard_len];
            for &gene_id in genes {
                let (cells, values) = expr.gene_slice(gene_id)?;
                let start = lower_bound(cells, shard_start as u32);
                let end = lower_bound(cells, shard_end as u32);
                for idx in start..end {
                    let cell = cells[idx] as usize - shard_start;
                    shard[cell] += values[idx];
                }
            }
            let denom = genes.len() as f32;
            for v in &mut shard[..(shard_end - shard_start)] {
                *v /= denom;
            }
        }
    }

    for shard_id in 0..shard_count {
        let shard_start = shard_id * shard_len;
        let shard_end = (shard_start + shard_len).min(n_cells);
        let shard = &scratch_shards[shard_id * shard_len..shard_id * shard_len + shard_len];
        out[shard_start..shard_end].copy_from_slice(&shard[..(shard_end - shard_start)]);
    }

    Ok(())
}

fn lower_bound(slice: &[u32], value: u32) -> usize {
    let mut left = 0usize;
    let mut right = slice.len();
    while left < right {
        let mid = (left + right) / 2;
        if slice[mid] < value {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    left
}
