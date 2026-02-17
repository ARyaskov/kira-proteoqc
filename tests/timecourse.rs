use kira_proteoqc::scores::TimepointSummary;
use kira_proteoqc::scores::timecourse::compute_timecourse;

fn tp(label: &str, pfs: f32, pcs: f32, utp: f32, pii: f32, cls: f32) -> TimepointSummary {
    TimepointSummary {
        label: label.to_string(),
        pfs,
        pii,
        pcs,
        cls,
        utp,
    }
}

#[test]
fn deltas_computed() {
    let tps = vec![
        tp("T0", 1.0, 2.0, 0.5, 0.2, 0.1),
        tp("T1", 2.5, 1.5, 0.6, 0.4, 0.2),
    ];
    let tc = compute_timecourse(tps).unwrap();
    assert_eq!(tc.deltas.len(), 1);
    let d = &tc.deltas[0];
    assert!((d.delta_pfs - 1.5).abs() < 1e-6);
    assert!((d.delta_pii - 0.2).abs() < 1e-6);
    assert!((d.delta_pcs + 0.5).abs() < 1e-6);
    assert!((d.delta_cls - 0.1).abs() < 1e-6);
}

#[test]
fn trajectory_adaptive_recovery() {
    let tps = vec![
        tp("T0", 2.0, 1.0, 1.0, 0.0, 0.0),
        tp("T1", 1.0, 1.0, 1.0, 0.0, 0.0),
    ];
    let tc = compute_timecourse(tps).unwrap();
    assert_eq!(tc.trajectory, "adaptive_recovery");
}

#[test]
fn trajectory_collapse_priority() {
    let tps = vec![
        tp("T0", 1.0, 2.0, 1.0, 0.0, 0.0),
        tp("T1", 2.0, 1.0, 2.0, 0.0, 0.0),
    ];
    let tc = compute_timecourse(tps).unwrap();
    assert_eq!(tc.trajectory, "collapse_trajectory");
}

#[test]
fn trajectory_addicted_survival() {
    let tps = vec![
        tp("T0", 1.0, 1.0, 1.0, 0.0, 0.0),
        tp("T1", 1.3, 2.0, 2.0, 0.0, 0.0),
    ];
    let tc = compute_timecourse(tps).unwrap();
    assert_eq!(tc.trajectory, "addicted_survival");
}

#[test]
fn trajectory_priority_order() {
    let tps = vec![
        tp("T0", 1.0, 2.0, 1.0, 0.0, 0.0),
        tp("T1", 2.0, 1.0, 2.0, 0.0, 0.0),
    ];
    let tc = compute_timecourse(tps).unwrap();
    // collapse_trajectory takes priority over addicted_survival
    assert_eq!(tc.trajectory, "collapse_trajectory");
}
