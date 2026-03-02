pub mod aggregate;
pub mod panels;
pub mod scores;

use anyhow::Result;

use crate::ctx::Ctx;

use self::aggregate::{ProteostasisExtensionSummary, build_summary};
use self::panels::PROTEO_EXTENSION_PANEL_V1;
use self::scores::{ProteostasisScores, compute_scores};

#[derive(Debug, Clone)]
pub struct ProteostasisExtensionResult {
    pub scores: ProteostasisScores,
    pub summary: ProteostasisExtensionSummary,
}

pub fn compute_extension(ctx: &Ctx) -> Result<ProteostasisExtensionResult> {
    let scores = compute_scores(ctx)?;
    let summary = build_summary(PROTEO_EXTENSION_PANEL_V1, &scores);
    Ok(ProteostasisExtensionResult { scores, summary })
}
