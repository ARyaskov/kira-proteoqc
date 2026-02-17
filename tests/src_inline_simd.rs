use kira_proteoqc::simd::{abs_diff, add_scaled, sum_f32};

#[test]
fn sum_is_stable_on_small_input() {
    let v = vec![1.0f32, 2.0, 3.5, 4.0, 5.0];
    let got = sum_f32(&v);
    assert!((got - 15.5).abs() < 1e-6);
}

#[test]
fn add_scaled_matches_expected_tail_behavior() {
    let mut dst = vec![1.0f32, 2.0, 3.0, 4.0, 5.0];
    let src = vec![1.0f32, 1.0, 1.0, 1.0, 1.0];
    add_scaled(&mut dst, &src, 0.5);
    assert_eq!(dst, vec![1.5, 2.5, 3.5, 4.5, 5.5]);
}

#[test]
fn abs_diff_matches_expected_tail_behavior() {
    let mut dst = vec![0.0f32; 5];
    let src = vec![1.0f32, 2.0, 3.0, 4.0, 5.0];
    abs_diff(&mut dst, &src, 3.0);
    assert_eq!(dst, vec![2.0, 1.0, 0.0, 1.0, 2.0]);
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[test]
fn simd_avx2_compiles() {
    let v = vec![1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let sum = sum_f32(&v);
    assert!((sum - 36.0).abs() < 1e-6);
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
#[test]
fn simd_neon_compiles() {
    let v = vec![1.0f32, 2.0, 3.0, 4.0];
    let sum = sum_f32(&v);
    assert!((sum - 10.0).abs() < 1e-6);
}
