use kira_proteoqc::math::stats::{mad, median, robust_z, trimmed_mean};

#[test]
fn trimmed_mean_basic() {
    let mut v = vec![1.0, 2.0, 3.0, 4.0, 100.0];
    let t = trimmed_mean(&mut v, 0.2);
    assert!((t - 3.0).abs() < 1e-6);
}

#[test]
fn trimmed_mean_small_fallback() {
    let mut v = vec![1.0, 2.0];
    let t = trimmed_mean(&mut v, 0.4);
    assert!((t - 1.5).abs() < 1e-6);
}

#[test]
fn median_odd_even() {
    let mut v1 = vec![3.0, 1.0, 2.0];
    assert_eq!(median(&mut v1), 2.0);
    let mut v2 = vec![4.0, 1.0, 2.0, 3.0];
    assert_eq!(median(&mut v2), 2.5);
}

#[test]
fn mad_basic() {
    let mut v = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let mut v1 = v.clone();
    let med = median(&mut v1);
    let mut v2 = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let m = mad(&mut v2, med);
    assert!((m - 1.0).abs() < 1e-6);
}

#[test]
fn robust_z_basic() {
    let z = robust_z(2.0, 1.0, 1.0);
    let expected = (2.0 - 1.0) / 1.4826;
    assert!((z - expected).abs() < 1e-6);
}
