#[test]
fn simd_scalar_compiles() {
    // Always compiled; ensures simd module is linked.
    let v = vec![1.0f32, 2.0, 3.0];
    let sum = kira_proteoqc::simd::sum_f32(&v);
    assert!((sum - 6.0).abs() < 1e-6);
}

#[cfg(all(feature = "simd-avx2", target_arch = "x86_64"))]
#[test]
fn simd_avx2_compiles() {
    let v = vec![1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let sum = kira_proteoqc::simd::sum_f32(&v);
    assert!((sum - 36.0).abs() < 1e-6);
}

#[cfg(all(feature = "simd-neon", target_arch = "aarch64"))]
#[test]
fn simd_neon_compiles() {
    let v = vec![1.0f32, 2.0, 3.0, 4.0];
    let sum = kira_proteoqc::simd::sum_f32(&v);
    assert!((sum - 10.0).abs() < 1e-6);
}
