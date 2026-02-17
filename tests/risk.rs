use kira_proteoqc::ctx::Ctx;
use kira_proteoqc::schema::v1::Mode;
use kira_proteoqc::scores::risk::compute_risk_flags;
use kira_proteoqc::scores::{AxisRawScores, IntegratedScores};

fn make_ctx(mode: Mode, axis: AxisRawScores, integrated: IntegratedScores) -> Ctx {
    let mut ctx = Ctx::new(
        std::path::PathBuf::from("input"),
        std::path::PathBuf::from("out"),
        mode,
        false,
        None,
        true,
        false,
        false,
        "0.0.0-test",
    );
    ctx.axis_raw = Some(axis);
    ctx.integrated_scores = Some(integrated);
    ctx
}

fn zeros(n: usize) -> Vec<f32> {
    vec![0.0; n]
}

#[test]
fn fragile_high_cell_only() {
    let axis = AxisRawScores {
        pcs: zeros(3),
        utp: zeros(3),
        cls: zeros(3),
        erad: zeros(3),
        ribo: zeros(3),
    };
    let integrated = IntegratedScores {
        capacity_raw: zeros(3),
        pii_raw: zeros(3),
        pfs_raw: zeros(3),
        capacity_z: Some(zeros(3)),
        pii_z: Some(zeros(3)),
        pfs_z: Some(vec![0.0, 0.0, 2.0]),
    };
    let mut ctx = make_ctx(Mode::Cell, axis, integrated);
    let flags = compute_risk_flags(&mut ctx).unwrap();
    let fired = flags.iter().find(|f| f.name == "fragile_high").unwrap();
    assert!(fired.fired);
    assert!(flags.iter().filter(|f| f.fired).count() == 1);
}

#[test]
fn proteasome_addiction_cell_only() {
    let axis = AxisRawScores {
        pcs: vec![0.0, 1.0, 4.0],
        utp: vec![0.0, 1.0, 4.0],
        cls: zeros(3),
        erad: zeros(3),
        ribo: zeros(3),
    };
    let integrated = IntegratedScores {
        capacity_raw: zeros(3),
        pii_raw: zeros(3),
        pfs_raw: zeros(3),
        capacity_z: Some(zeros(3)),
        pii_z: Some(zeros(3)),
        pfs_z: Some(vec![0.0, 0.0, 0.5]),
    };
    let mut ctx = make_ctx(Mode::Cell, axis, integrated);
    let flags = compute_risk_flags(&mut ctx).unwrap();
    let fired = flags
        .iter()
        .find(|f| f.name == "proteasome_addiction")
        .unwrap();
    assert!(fired.fired);
    assert!(flags.iter().filter(|f| f.fired).count() == 1);
}

#[test]
fn proteotoxic_stress_cell_only() {
    let axis = AxisRawScores {
        pcs: zeros(3),
        utp: zeros(3),
        cls: vec![0.0, 1.0, 4.0],
        erad: zeros(3),
        ribo: zeros(3),
    };
    let integrated = IntegratedScores {
        capacity_raw: zeros(3),
        pii_raw: zeros(3),
        pfs_raw: zeros(3),
        capacity_z: Some(zeros(3)),
        pii_z: Some(vec![0.0, 0.0, 2.0]),
        pfs_z: Some(zeros(3)),
    };
    let mut ctx = make_ctx(Mode::Cell, axis, integrated);
    let flags = compute_risk_flags(&mut ctx).unwrap();
    let fired = flags
        .iter()
        .find(|f| f.name == "proteotoxic_stress")
        .unwrap();
    assert!(fired.fired);
    assert!(flags.iter().filter(|f| f.fired).count() == 1);
}

#[test]
fn er_degradation_overdrive_cell_only() {
    let axis = AxisRawScores {
        pcs: zeros(3),
        utp: zeros(3),
        cls: zeros(3),
        erad: vec![0.0, 1.0, 4.0],
        ribo: zeros(3),
    };
    let integrated = IntegratedScores {
        capacity_raw: zeros(3),
        pii_raw: zeros(3),
        pfs_raw: zeros(3),
        capacity_z: Some(zeros(3)),
        pii_z: Some(vec![0.0, 0.0, 2.0]),
        pfs_z: Some(zeros(3)),
    };
    let mut ctx = make_ctx(Mode::Cell, axis, integrated);
    let flags = compute_risk_flags(&mut ctx).unwrap();
    let fired = flags
        .iter()
        .find(|f| f.name == "er_degradation_overdrive")
        .unwrap();
    assert!(fired.fired);
    assert!(flags.iter().filter(|f| f.fired).count() == 1);
}

#[test]
fn boundary_not_fired() {
    let axis = AxisRawScores {
        pcs: vec![0.0, 1.0, 2.0],
        utp: vec![0.0, 1.0, 2.0],
        cls: zeros(3),
        erad: zeros(3),
        ribo: zeros(3),
    };
    let integrated = IntegratedScores {
        capacity_raw: zeros(3),
        pii_raw: zeros(3),
        pfs_raw: zeros(3),
        capacity_z: Some(zeros(3)),
        pii_z: Some(zeros(3)),
        pfs_z: Some(vec![0.0, 0.0, 1.5]),
    };
    let mut ctx = make_ctx(Mode::Cell, axis, integrated);
    let flags = compute_risk_flags(&mut ctx).unwrap();
    let fragile = flags.iter().find(|f| f.name == "fragile_high").unwrap();
    assert!(!fragile.fired);
}
