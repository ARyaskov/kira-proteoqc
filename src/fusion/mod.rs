use std::collections::HashMap;

use anyhow::Result;

use crate::expr::reader::ExprReader;
use crate::geneset::ResolvedGeneset;

#[derive(Debug, Clone)]
pub struct TargetRef {
    pub target_id: usize,
    pub scale: f32,
}

#[derive(Debug, Clone)]
pub struct FusionPlan {
    pub targets: Vec<String>,
    pub gene_counts: Vec<usize>,
    pub membership: Vec<Vec<TargetRef>>,
    pub n_cells: usize,
}

pub fn build_plan(
    n_genes: usize,
    n_cells: usize,
    resolved: &[ResolvedGeneset],
    target_ids: &[&str],
) -> FusionPlan {
    let mut target_index = HashMap::new();
    let mut targets = Vec::new();
    let mut gene_counts = Vec::new();

    for &id in target_ids {
        targets.push(id.to_string());
        target_index.insert(id.to_string(), targets.len() - 1);
        gene_counts.push(0);
    }

    let mut membership: Vec<Vec<TargetRef>> = vec![Vec::new(); n_genes];

    for gs in resolved {
        if let Some(&tid) = target_index.get(&gs.id) {
            gene_counts[tid] = gs.gene_ids.len();
            for &gene_id in &gs.gene_ids {
                membership[gene_id].push(TargetRef {
                    target_id: tid,
                    scale: 1.0,
                });
            }
        }
    }

    for targets in &mut membership {
        targets.sort_by_key(|t| t.target_id);
    }

    FusionPlan {
        targets,
        gene_counts,
        membership,
        n_cells,
    }
}

pub fn fused_reduce(expr: &ExprReader<'_>, plan: &FusionPlan, out: &mut [f32]) -> Result<()> {
    let target_count = plan.targets.len();
    let expected = target_count * plan.n_cells;
    if out.len() != expected {
        anyhow::bail!("fused output length mismatch");
    }
    for v in out.iter_mut() {
        *v = 0.0;
    }

    for gene_id in 0..plan.membership.len() {
        let targets = &plan.membership[gene_id];
        if targets.is_empty() {
            continue;
        }
        let (cells, values) = expr.gene_slice(gene_id)?;
        for (cell, value) in cells.iter().zip(values.iter()) {
            let cell_idx = *cell as usize;
            for t in targets {
                let offset = t.target_id * plan.n_cells + cell_idx;
                out[offset] += *value * t.scale;
            }
        }
    }

    Ok(())
}
