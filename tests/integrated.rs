use kira_proteoqc::schema::v1::Mode;
use kira_proteoqc::scores::AxisRawScores;
use kira_proteoqc::scores::integrated::compute_integrated;

fn assert_vec_close(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b.iter()) {
        assert!((x - y).abs() < 1e-6, "{} vs {}", x, y);
    }
}

#[test]
fn integration_math_basic() {
    let axis = AxisRawScores {
        pcs: vec![1.0, 2.0],
        utp: vec![0.5, 1.5],
        cls: vec![2.0, 0.0],
        erad: vec![1.0, 1.0],
        ribo: vec![3.0, 2.0],
    };

    let (integrated, contrib) = compute_integrated(&axis, Mode::Sample).unwrap();

    assert_vec_close(&integrated.capacity_raw, &[1.25, 1.3]);
    assert_vec_close(&integrated.pii_raw, &[1.75, 0.7]);
    assert_vec_close(&integrated.pfs_raw, &[1.275, 0.755]);

    assert_vec_close(&contrib.c_pii, &[0.7, 0.28]);
    assert_vec_close(&contrib.c_utp, &[0.125, 0.375]);
    assert_vec_close(&contrib.c_ribo, &[0.6, 0.4]);
    assert_vec_close(&contrib.c_pcs, &[-0.15, -0.3]);
}

#[test]
fn zscore_known_vector() {
    let axis = AxisRawScores {
        pcs: vec![1.0, 2.0, 3.0],
        utp: vec![0.0, 0.0, 0.0],
        cls: vec![0.0, 0.0, 0.0],
        erad: vec![0.0, 0.0, 0.0],
        ribo: vec![0.0, 0.0, 0.0],
    };

    let (integrated, _contrib) = compute_integrated(&axis, Mode::Cell).unwrap();
    let z = integrated.capacity_z.as_ref().unwrap();
    // capacity_raw == 0.55 * pcs, so [0.55, 1.10, 1.65]
    // median = 1.10, MAD = 0.55
    // z = (x - 1.10) / (1.4826 * 0.55)
    let denom = 1.4826 * 0.55;
    let expected = vec![(-0.55) / denom, 0.0, 0.55 / denom];
    assert_vec_close(z, &expected);
}

#[test]
fn mode_behavior_z() {
    let axis = AxisRawScores {
        pcs: vec![1.0],
        utp: vec![1.0],
        cls: vec![1.0],
        erad: vec![1.0],
        ribo: vec![1.0],
    };

    let (cell_scores, _) = compute_integrated(&axis, Mode::Cell).unwrap();
    assert!(cell_scores.capacity_z.is_some());
    let (sample_scores, _) = compute_integrated(&axis, Mode::Sample).unwrap();
    assert!(sample_scores.capacity_z.is_none());
}
