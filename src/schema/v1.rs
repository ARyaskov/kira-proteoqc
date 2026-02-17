use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Mode {
    Cell,
    Sample,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Normalization {
    pub log1p: bool,
    pub cp10k: bool,
    pub raw: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputMeta {
    pub genes: Option<u64>,
    pub cells: Option<u64>,
    pub nnz: Option<u64>,
    pub normalization: Normalization,
    pub mode: Mode,
    pub timecourse: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerSampleScore {
    pub id: String,
    #[serde(rename = "PCS_raw")]
    pub pcs_raw: Option<f64>,
    #[serde(rename = "UTP_raw")]
    pub utp_raw: Option<f64>,
    #[serde(rename = "CLS_raw")]
    pub cls_raw: Option<f64>,
    #[serde(rename = "ERAD_raw")]
    pub erad_raw: Option<f64>,
    #[serde(rename = "Ribo_raw")]
    pub ribo_raw: Option<f64>,
    #[serde(rename = "Capacity_raw")]
    pub capacity_raw: Option<f64>,
    #[serde(rename = "PII_raw")]
    pub pii_raw: Option<f64>,
    #[serde(rename = "PFS_raw")]
    pub pfs_raw: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Distributions {
    pub placeholder: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Scores {
    pub per_sample: Option<Vec<PerSampleScore>>,
    pub per_cell_tsv_path: Option<String>,
    pub distributions: Option<Distributions>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RiskFlag {
    pub name: String,
    pub fired: bool,
    pub threshold: Option<String>,
    pub details: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComponentContribution {
    pub name: String,
    pub weight: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenesetCoverage {
    pub geneset: String,
    pub found: u64,
    pub total: u64,
    pub fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PfsContributions {
    pub pii: f64,
    pub utp: f64,
    pub ribo: f64,
    pub pcs: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Explainability {
    pub component_contributions: Option<Vec<ComponentContribution>>,
    pub geneset_coverage: Vec<GenesetCoverage>,
    pub pfs_contributions: Option<PfsContributions>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimepointSummary {
    pub label: String,
    pub pfs: f32,
    pub pii: f32,
    pub pcs: f32,
    pub cls: f32,
    pub utp: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeltaSummary {
    pub from: String,
    pub to: String,
    pub delta_pfs: f32,
    pub delta_pii: f32,
    pub delta_pcs: f32,
    pub delta_cls: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimecourseResult {
    pub timepoints: Vec<TimepointSummary>,
    pub deltas: Vec<DeltaSummary>,
    pub trajectory: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteoQcV1 {
    pub tool: String,
    pub version: String,
    pub schema_version: String,
    pub input_meta: InputMeta,
    pub scores: Scores,
    pub risk_flags: Vec<RiskFlag>,
    pub explainability: Explainability,
    pub timecourse: Option<TimecourseResult>,
}

impl ProteoQcV1 {
    pub fn empty(tool_version: &str, mode: Mode, timecourse: bool, log1p: bool) -> Self {
        Self {
            tool: "kira-proteoqc".to_string(),
            version: tool_version.to_string(),
            schema_version: "v1".to_string(),
            input_meta: InputMeta {
                genes: None,
                cells: None,
                nnz: None,
                normalization: Normalization {
                    log1p,
                    cp10k: false,
                    raw: true,
                },
                mode,
                timecourse,
            },
            scores: Scores {
                per_sample: None,
                per_cell_tsv_path: None,
                distributions: None,
            },
            risk_flags: Vec::new(),
            explainability: Explainability {
                component_contributions: None,
                geneset_coverage: Vec::new(),
                pfs_contributions: None,
            },
            timecourse: None,
        }
    }
}
