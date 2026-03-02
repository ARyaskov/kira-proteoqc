use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use super::scores::{ProteostasisMissingness, ProteostasisScores, ProteostasisThresholds};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThresholdsOut {
    pub chaperone_high: f64,
    pub proteasome_high: f64,
    pub upr_active: f64,
    pub proteotoxic_high: f64,
    pub imbalance_high: f64,
    pub collapse_risk: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetricStat {
    pub median: f64,
    pub p10: f64,
    pub p90: f64,
    pub mad: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GlobalStatsOut {
    pub metrics: BTreeMap<String, MetricStat>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClusterStatsOut {
    pub cluster: String,
    pub metrics: BTreeMap<String, MetricStat>,
    pub flag_fractions: BTreeMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TopClusterOut {
    pub cluster: String,
    pub collapse_risk_median: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelMissingnessOut {
    pub panel_id: String,
    pub total_genes: u64,
    pub mapped_genes: u64,
    pub insufficient_genes: bool,
    pub nan_cells: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MissingnessOut {
    pub panels: Vec<PanelMissingnessOut>,
    pub score_nan_counts: BTreeMap<String, u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteostasisExtensionSummary {
    pub panel_version: String,
    pub thresholds: ThresholdsOut,
    pub global_stats: GlobalStatsOut,
    pub cluster_stats: Vec<ClusterStatsOut>,
    pub top_clusters_by_collapse_risk: Vec<TopClusterOut>,
    pub missingness: MissingnessOut,
    pub translation_source: Option<String>,
}

pub fn build_summary(
    panel_version: &str,
    scores: &ProteostasisScores,
) -> ProteostasisExtensionSummary {
    let metrics = metric_views(scores);

    let mut global = BTreeMap::new();
    for (name, values) in &metrics {
        global.insert((*name).to_string(), stats(values));
    }

    let mut flag_fractions = BTreeMap::new();
    for (name, flags) in flag_views(scores) {
        flag_fractions.insert(name.to_string(), fraction_true(flags));
    }

    let mut cluster_metrics = BTreeMap::new();
    for (name, values) in &metrics {
        cluster_metrics.insert((*name).to_string(), stats(values));
    }
    let all_cells = ClusterStatsOut {
        cluster: "all_cells".to_string(),
        metrics: cluster_metrics,
        flag_fractions,
    };

    let mut top_clusters = Vec::new();
    let pcp_stat = stats(&scores.pcp);
    if pcp_stat.median.is_finite() {
        top_clusters.push(TopClusterOut {
            cluster: "all_cells".to_string(),
            collapse_risk_median: pcp_stat.median,
        });
    }

    ProteostasisExtensionSummary {
        panel_version: panel_version.to_string(),
        thresholds: thresholds(scores.thresholds),
        global_stats: GlobalStatsOut { metrics: global },
        cluster_stats: vec![all_cells],
        top_clusters_by_collapse_risk: top_clusters,
        missingness: missingness_out(&scores.missingness),
        translation_source: scores.translation_source.clone(),
    }
}

fn thresholds(t: ProteostasisThresholds) -> ThresholdsOut {
    ThresholdsOut {
        chaperone_high: t.chaperone_high as f64,
        proteasome_high: t.proteasome_high as f64,
        upr_active: t.upr_active as f64,
        proteotoxic_high: t.proteotoxic_high as f64,
        imbalance_high: t.imbalance_high as f64,
        collapse_risk: t.collapse_risk as f64,
    }
}

fn metric_views(scores: &ProteostasisScores) -> Vec<(&'static str, &[f32])> {
    vec![
        ("CCI", &scores.cci),
        ("PCI", &scores.pci),
        ("UPR_A", &scores.upr_a),
        ("PLS", &scores.pls),
        ("SCI", &scores.sci),
        ("PCP", &scores.pcp),
    ]
}

fn flag_views(scores: &ProteostasisScores) -> Vec<(&'static str, &[bool])> {
    vec![
        ("chaperone_high", &scores.chaperone_high),
        ("proteasome_high", &scores.proteasome_high),
        ("upr_active", &scores.upr_active),
        ("proteotoxic_high", &scores.proteotoxic_high),
        ("imbalance_high", &scores.imbalance_high),
        ("collapse_risk", &scores.collapse_risk),
    ]
}

fn stats(values: &[f32]) -> MetricStat {
    let mut finite = Vec::with_capacity(values.len());
    for &v in values {
        if !v.is_nan() {
            finite.push(v as f64);
        }
    }
    if finite.is_empty() {
        return MetricStat {
            median: 0.0,
            p10: 0.0,
            p90: 0.0,
            mad: 0.0,
        };
    }
    finite.sort_by(|a, b| a.total_cmp(b));
    let med = percentile_sorted(&finite, 0.50);
    let p10 = percentile_sorted(&finite, 0.10);
    let p90 = percentile_sorted(&finite, 0.90);

    let mut abs_dev = Vec::with_capacity(finite.len());
    for v in &finite {
        abs_dev.push((v - med).abs());
    }
    abs_dev.sort_by(|a, b| a.total_cmp(b));
    let mad = percentile_sorted(&abs_dev, 0.50);

    MetricStat {
        median: round6(med),
        p10: round6(p10),
        p90: round6(p90),
        mad: round6(mad),
    }
}

fn percentile_sorted(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let clamped = q.clamp(0.0, 1.0);
    let pos = clamped * (sorted.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let t = pos - lo as f64;
        sorted[lo] * (1.0 - t) + sorted[hi] * t
    }
}

fn fraction_true(flags: &[bool]) -> f64 {
    if flags.is_empty() {
        return 0.0;
    }
    let mut n_true = 0usize;
    for &f in flags {
        if f {
            n_true += 1;
        }
    }
    round6(n_true as f64 / flags.len() as f64)
}

fn missingness_out(m: &ProteostasisMissingness) -> MissingnessOut {
    let mut panels = Vec::with_capacity(m.panels.len());
    for p in &m.panels {
        panels.push(PanelMissingnessOut {
            panel_id: p.panel_id.clone(),
            total_genes: p.total_genes as u64,
            mapped_genes: p.mapped_genes as u64,
            insufficient_genes: p.insufficient_genes,
            nan_cells: p.nan_cells as u64,
        });
    }

    let mut score_nan_counts = BTreeMap::new();
    for (k, v) in &m.score_nan_counts {
        score_nan_counts.insert(k.clone(), *v as u64);
    }

    MissingnessOut {
        panels,
        score_nan_counts,
    }
}

fn round6(v: f64) -> f64 {
    (v * 1_000_000.0).round() / 1_000_000.0
}
