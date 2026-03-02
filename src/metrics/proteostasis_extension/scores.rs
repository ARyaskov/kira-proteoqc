use std::collections::{BTreeMap, HashMap};

use anyhow::{Result, bail};

use crate::ctx::Ctx;
use crate::math::stats::{mad, median};

use super::panels::{AGGREGATION_PANEL, CHAPERONE_PANEL, ERAD_PANEL, PROTEASOME_PANEL, UPR_PANEL};

const TRIM_FRACTION: f32 = 0.10;
const MIN_GENES: usize = 2;
const MAD_EPS: f32 = 1e-6;

#[derive(Debug, Clone, Copy)]
pub struct ProteostasisThresholds {
    pub chaperone_high: f32,
    pub proteasome_high: f32,
    pub upr_active: f32,
    pub proteotoxic_high: f32,
    pub imbalance_high: f32,
    pub collapse_risk: f32,
}

impl Default for ProteostasisThresholds {
    fn default() -> Self {
        Self {
            chaperone_high: 2.0,
            proteasome_high: 2.0,
            upr_active: 1.5,
            proteotoxic_high: 2.0,
            imbalance_high: 1.5,
            collapse_risk: 2.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PanelMissingness {
    pub panel_id: String,
    pub total_genes: usize,
    pub mapped_genes: usize,
    pub insufficient_genes: bool,
    pub nan_cells: usize,
}

#[derive(Debug, Clone)]
pub struct ProteostasisMissingness {
    pub panels: Vec<PanelMissingness>,
    pub score_nan_counts: BTreeMap<String, usize>,
}

#[derive(Debug, Clone)]
pub struct ProteostasisScores {
    pub chaperone_core: Vec<f32>,
    pub proteasome_core: Vec<f32>,
    pub upr_core: Vec<f32>,
    pub agg_core: Vec<f32>,
    pub erad_core: Vec<f32>,
    pub cci: Vec<f32>,
    pub pci: Vec<f32>,
    pub upr_a: Vec<f32>,
    pub pls: Vec<f32>,
    pub sci: Vec<f32>,
    pub pcp: Vec<f32>,
    pub chaperone_high: Vec<bool>,
    pub proteasome_high: Vec<bool>,
    pub upr_active: Vec<bool>,
    pub proteotoxic_high: Vec<bool>,
    pub imbalance_high: Vec<bool>,
    pub collapse_risk: Vec<bool>,
    pub thresholds: ProteostasisThresholds,
    pub translation_source: Option<String>,
    pub missingness: ProteostasisMissingness,
}

pub fn compute_scores(ctx: &Ctx) -> Result<ProteostasisScores> {
    let n_cells = ctx.cells.len();
    let thresholds = ProteostasisThresholds::default();
    let upper = uppercase_gene_index(ctx);

    let chaperone_genes = resolve_panel(&ctx.gene_index, &upper, CHAPERONE_PANEL);
    let proteasome_genes = resolve_panel(&ctx.gene_index, &upper, PROTEASOME_PANEL);
    let upr_genes = resolve_panel(&ctx.gene_index, &upper, UPR_PANEL);
    let erad_genes = resolve_panel(&ctx.gene_index, &upper, ERAD_PANEL);
    let agg_genes = resolve_panel(&ctx.gene_index, &upper, AGGREGATION_PANEL);

    let chaperone_dense = build_panel_dense(ctx, &chaperone_genes)?;
    let proteasome_dense = build_panel_dense(ctx, &proteasome_genes)?;
    let upr_dense = build_panel_dense(ctx, &upr_genes)?;
    let erad_dense = build_panel_dense(ctx, &erad_genes)?;
    let agg_dense = build_panel_dense(ctx, &agg_genes)?;

    let mut scratch = Vec::new();
    let chaperone_core = panel_trimmed_mean(&chaperone_dense, n_cells, &mut scratch);
    let proteasome_core = panel_trimmed_mean(&proteasome_dense, n_cells, &mut scratch);
    let upr_core = panel_trimmed_mean(&upr_dense, n_cells, &mut scratch);
    let erad_core = panel_trimmed_mean(&erad_dense, n_cells, &mut scratch);
    let agg_core = panel_trimmed_mean(&agg_dense, n_cells, &mut scratch);

    let cci = robust_z_vec(&chaperone_core);
    let pci = robust_z_vec(&proteasome_core);
    let upr_a = robust_z_vec(&upr_core);
    let agg_z = robust_z_vec(&agg_core);

    let (translation_z, translation_source) = match &ctx.translation_load_z {
        Some(v) if v.len() == n_cells => (Some(v.as_slice()), Some("riboqc".to_string())),
        _ => (None, None),
    };

    let mut pls = vec![f32::NAN; n_cells];
    let mut sci = vec![f32::NAN; n_cells];
    let mut pcp = vec![f32::NAN; n_cells];

    for i in 0..n_cells {
        let upr_v = upr_a[i];
        if upr_v.is_nan() {
            continue;
        }

        let agg_v = agg_z[i];
        let pls_v = if agg_v.is_nan() {
            upr_v
        } else {
            0.6 * upr_v + 0.4 * agg_v.max(0.0)
        };
        pls[i] = pls_v;

        let cap = mean_pair(cci[i], pci[i]);
        if cap.is_nan() {
            continue;
        }
        let sci_v = if let Some(tz) = translation_z {
            let tv = tz[i];
            if tv.is_nan() {
                f32::NAN
            } else {
                (tv - cap).max(0.0)
            }
        } else {
            (upr_v - cap).max(0.0)
        };
        sci[i] = sci_v;
        if sci_v.is_nan() {
            continue;
        }

        let pcp_v = 0.4 * pls_v + 0.3 * sci_v - 0.3 * cap;
        pcp[i] = pcp_v.max(0.0);
    }

    let chaperone_high = threshold_flags(&cci, thresholds.chaperone_high);
    let proteasome_high = threshold_flags(&pci, thresholds.proteasome_high);
    let upr_active = threshold_flags(&upr_a, thresholds.upr_active);
    let proteotoxic_high = threshold_flags(&pls, thresholds.proteotoxic_high);
    let imbalance_high = threshold_flags(&sci, thresholds.imbalance_high);
    let collapse_risk = threshold_flags(&pcp, thresholds.collapse_risk);

    let missingness = ProteostasisMissingness {
        panels: vec![
            panel_missingness(
                "chaperone",
                CHAPERONE_PANEL.len(),
                &chaperone_genes,
                &chaperone_core,
            ),
            panel_missingness(
                "proteasome",
                PROTEASOME_PANEL.len(),
                &proteasome_genes,
                &proteasome_core,
            ),
            panel_missingness("upr", UPR_PANEL.len(), &upr_genes, &upr_core),
            panel_missingness("erad", ERAD_PANEL.len(), &erad_genes, &erad_core),
            panel_missingness(
                "aggregation",
                AGGREGATION_PANEL.len(),
                &agg_genes,
                &agg_core,
            ),
        ],
        score_nan_counts: score_nan_counts(&[
            ("CCI", &cci),
            ("PCI", &pci),
            ("UPR_A", &upr_a),
            ("PLS", &pls),
            ("SCI", &sci),
            ("PCP", &pcp),
        ]),
    };

    Ok(ProteostasisScores {
        chaperone_core,
        proteasome_core,
        upr_core,
        agg_core,
        erad_core,
        cci,
        pci,
        upr_a,
        pls,
        sci,
        pcp,
        chaperone_high,
        proteasome_high,
        upr_active,
        proteotoxic_high,
        imbalance_high,
        collapse_risk,
        thresholds,
        translation_source,
        missingness,
    })
}

fn uppercase_gene_index(ctx: &Ctx) -> HashMap<String, usize> {
    let mut map = HashMap::with_capacity(ctx.gene_index.len());
    for (idx, gene) in ctx.genes.iter().enumerate() {
        let up = gene.to_ascii_uppercase();
        map.entry(up).or_insert(idx);
    }
    map
}

fn resolve_panel(
    direct: &HashMap<String, usize>,
    uppercase: &HashMap<String, usize>,
    panel: &[&str],
) -> Vec<usize> {
    let mut out = Vec::with_capacity(panel.len());
    for gene in panel {
        if let Some(&idx) = direct.get(*gene) {
            out.push(idx);
            continue;
        }
        if let Some(&idx) = uppercase.get(&gene.to_ascii_uppercase()) {
            out.push(idx);
        }
    }
    out
}

fn build_panel_dense(ctx: &Ctx, genes: &[usize]) -> Result<Vec<Vec<f32>>> {
    let n_cells = ctx.cells.len();
    let mut dense = Vec::with_capacity(genes.len());
    for &gene_id in genes {
        let mut row = vec![0.0f32; n_cells];
        let (cells, values) = ctx.gene_slice(gene_id)?;
        for (&cell, &value) in cells.iter().zip(values.iter()) {
            if value.is_nan() {
                bail!("NaN encountered in expression matrix");
            }
            row[cell as usize] = value;
        }
        dense.push(row);
    }
    Ok(dense)
}

fn panel_trimmed_mean(
    panel_dense: &[Vec<f32>],
    n_cells: usize,
    scratch: &mut Vec<f32>,
) -> Vec<f32> {
    if panel_dense.len() < MIN_GENES {
        return vec![f32::NAN; n_cells];
    }
    let mut out = vec![f32::NAN; n_cells];
    scratch.clear();
    scratch.reserve(panel_dense.len());

    for cell in 0..n_cells {
        scratch.clear();
        for gene in panel_dense {
            scratch.push(gene[cell]);
        }
        out[cell] = trimmed_mean_in_place(scratch);
    }
    out
}

fn trimmed_mean_in_place(values: &mut [f32]) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let n = values.len();
    let k = (TRIM_FRACTION * n as f32).floor() as usize;
    let (start, end) = if n < 2 * k + 1 {
        (0usize, n)
    } else {
        (k, n - k)
    };
    let mut sum = 0.0f32;
    let mut count = 0usize;
    for v in values.iter().take(end).skip(start) {
        sum += *v;
        count += 1;
    }
    if count == 0 {
        f32::NAN
    } else {
        sum / count as f32
    }
}

fn robust_z_vec(values: &[f32]) -> Vec<f32> {
    let mut finite = Vec::with_capacity(values.len());
    for &v in values {
        if !v.is_nan() {
            finite.push(v);
        }
    }
    if finite.is_empty() {
        return vec![f32::NAN; values.len()];
    }

    let mut med_buf = finite.clone();
    let med = median(&mut med_buf);
    let mut mad_buf = finite;
    let mad_v = mad(&mut mad_buf, med);

    let mut out = Vec::with_capacity(values.len());
    for &v in values {
        if v.is_nan() {
            out.push(f32::NAN);
            continue;
        }
        if mad_v == 0.0 {
            out.push(0.0);
            continue;
        }
        out.push((v - med) / (1.4826 * mad_v + MAD_EPS));
    }
    out
}

fn threshold_flags(values: &[f32], threshold: f32) -> Vec<bool> {
    let mut out = Vec::with_capacity(values.len());
    for &v in values {
        out.push(!v.is_nan() && v >= threshold);
    }
    out
}

fn panel_missingness(
    panel_id: &str,
    total: usize,
    mapped: &[usize],
    values: &[f32],
) -> PanelMissingness {
    let mut nan_cells = 0usize;
    for v in values {
        if v.is_nan() {
            nan_cells += 1;
        }
    }
    PanelMissingness {
        panel_id: panel_id.to_string(),
        total_genes: total,
        mapped_genes: mapped.len(),
        insufficient_genes: mapped.len() < MIN_GENES,
        nan_cells,
    }
}

fn score_nan_counts(scores: &[(&str, &[f32])]) -> BTreeMap<String, usize> {
    let mut out = BTreeMap::new();
    for (name, values) in scores {
        let mut nan_count = 0usize;
        for v in *values {
            if v.is_nan() {
                nan_count += 1;
            }
        }
        out.insert((*name).to_string(), nan_count);
    }
    out
}

fn mean_pair(a: f32, b: f32) -> f32 {
    if a.is_nan() || b.is_nan() {
        f32::NAN
    } else {
        (a + b) * 0.5
    }
}
