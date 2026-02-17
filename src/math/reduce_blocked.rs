use anyhow::{Result, bail};

use crate::expr::prefetch::{prefetch_read, prefetch_write};
use crate::expr::reader::ExprReader;

pub fn per_cell_raw_blocked(
    expr: &ExprReader<'_>,
    genes: &[usize],
    out: &mut [f32],
    block_size: usize,
    do_prefetch: bool,
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
    let block = block_size.max(1);

    for &gene_id in genes {
        let (cells, values) = expr.gene_slice(gene_id)?;
        let mut idx = 0usize;
        while idx < cells.len() {
            let cell = cells[idx] as usize;
            let block_start = (cell / block) * block;
            let block_end = (block_start + block).min(n_cells);
            let mut end = idx + 1;
            while end < cells.len() && (cells[end] as usize) < block_end {
                end += 1;
            }
            if do_prefetch {
                let p = out.as_ptr().wrapping_add(block_start);
                prefetch_write(p as *const u8);
                let p2 = values.as_ptr().wrapping_add(end.min(values.len() - 1));
                prefetch_read(p2 as *const u8);
            }
            for i in idx..end {
                let c = cells[i] as usize;
                out[c] += values[i];
            }
            idx = end;
        }
    }

    let denom = genes.len() as f32;
    for v in out.iter_mut() {
        *v /= denom;
    }
    Ok(())
}
