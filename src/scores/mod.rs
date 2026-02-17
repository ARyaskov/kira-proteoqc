pub mod axis_raw;
pub mod integrated;
pub mod risk;
pub mod timecourse;

#[derive(Debug, Clone)]
pub struct AxisRawScores {
    pub pcs: Vec<f32>,
    pub utp: Vec<f32>,
    pub cls: Vec<f32>,
    pub erad: Vec<f32>,
    pub ribo: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct IntegratedScores {
    pub capacity_raw: Vec<f32>,
    pub pii_raw: Vec<f32>,
    pub pfs_raw: Vec<f32>,
    pub capacity_z: Option<Vec<f32>>,
    pub pii_z: Option<Vec<f32>>,
    pub pfs_z: Option<Vec<f32>>,
}

#[derive(Debug, Clone)]
pub struct PfsContributions {
    pub c_pii: Vec<f32>,
    pub c_utp: Vec<f32>,
    pub c_ribo: Vec<f32>,
    pub c_pcs: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct RiskFlag {
    pub name: String,
    pub fired: bool,
    pub threshold: String,
    pub details: Option<String>,
}

#[derive(Debug, Clone)]
pub struct TimepointSummary {
    pub label: String,
    pub pfs: f32,
    pub pii: f32,
    pub pcs: f32,
    pub cls: f32,
    pub utp: f32,
}

#[derive(Debug, Clone)]
pub struct DeltaSummary {
    pub from: String,
    pub to: String,
    pub delta_pfs: f32,
    pub delta_pii: f32,
    pub delta_pcs: f32,
    pub delta_cls: f32,
}

#[derive(Debug, Clone)]
pub struct TimecourseResult {
    pub timepoints: Vec<TimepointSummary>,
    pub deltas: Vec<DeltaSummary>,
    pub trajectory: String,
}
