use anyhow::{Result, bail};
use tracing::warn;

use crate::expr::reader::ExprReader;
#[cfg(feature = "cache-block")]
use crate::math::reduce_blocked;
#[cfg(feature = "mt")]
use crate::math::reduce_mt;
use crate::simd;

pub struct GeneSetReducer<'a> {
    pub expr: &'a ExprReader<'a>,
    pub threads: usize,
    pub cache_block: usize,
    pub prefetch: bool,
    #[cfg(feature = "mt")]
    scratch_shards: Vec<f32>,
}

impl<'a> GeneSetReducer<'a> {
    pub fn new(
        expr: &'a ExprReader<'a>,
        threads: usize,
        cache_block: usize,
        prefetch: bool,
    ) -> Self {
        Self {
            expr,
            threads,
            cache_block,
            prefetch,
            #[cfg(feature = "mt")]
            scratch_shards: Vec::new(),
        }
    }

    pub fn per_cell_raw(&mut self, genes: &[usize], out: &mut [f32]) -> Result<()> {
        if out.len() != self.expr.n_cells() {
            bail!("output buffer length does not match n_cells");
        }
        if genes.is_empty() {
            warn!("empty geneset in per_cell_raw");
            for v in out.iter_mut() {
                *v = 0.0;
            }
            return Ok(());
        }

        #[cfg(feature = "mt")]
        {
            if self.threads > 1 {
                return reduce_mt::per_cell_raw_mt(
                    self.expr,
                    genes,
                    out,
                    self.threads,
                    &mut self.scratch_shards,
                );
            }
        }

        #[cfg(feature = "cache-block")]
        {
            if self.cache_block > 0 {
                return reduce_blocked::per_cell_raw_blocked(
                    self.expr,
                    genes,
                    out,
                    self.cache_block,
                    self.prefetch,
                );
            }
        }

        self.per_cell_raw_scalar(genes, out)
    }

    fn per_cell_raw_scalar(&self, genes: &[usize], out: &mut [f32]) -> Result<()> {
        for v in out.iter_mut() {
            *v = 0.0;
        }
        let n_cells = self.expr.n_cells();
        let mut dense = vec![0.0f32; n_cells];

        for &gene_id in genes {
            let (cells, values) = self.expr.gene_slice(gene_id)?;
            for (&cell, &value) in cells.iter().zip(values.iter()) {
                if value.is_nan() {
                    bail!("NaN encountered in expr.bin values");
                }
                dense[cell as usize] += value;
            }
            simd::add_scaled(out, &dense, 1.0);
            for &cell in cells {
                dense[cell as usize] = 0.0;
            }
        }

        let denom = genes.len() as f32;
        for v in out.iter_mut() {
            *v /= denom;
        }
        Ok(())
    }

    pub fn per_sample_raw(&mut self, genes: &[usize]) -> Result<f32> {
        let mut buf = vec![0.0f32; self.expr.n_cells()];
        self.per_cell_raw(genes, &mut buf)?;
        let sum = simd::sum_f32(&buf);
        if buf.is_empty() {
            Ok(0.0)
        } else {
            Ok(sum / buf.len() as f32)
        }
    }
}
