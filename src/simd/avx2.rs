#![cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#![allow(unsafe_op_in_unsafe_fn)]

use core::arch::x86_64::*;

pub unsafe fn sum_f32(slice: &[f32]) -> f32 {
    let mut i = 0usize;
    let mut acc = _mm256_setzero_ps();
    while i + 8 <= slice.len() {
        let v = _mm256_loadu_ps(slice.as_ptr().add(i));
        acc = _mm256_add_ps(acc, v);
        i += 8;
    }
    let mut lanes = [0f32; 8];
    _mm256_storeu_ps(lanes.as_mut_ptr(), acc);
    let mut sum = 0.0f32;
    for j in 0..8 {
        sum += lanes[j];
    }
    while i < slice.len() {
        sum += slice[i];
        i += 1;
    }
    sum
}

pub unsafe fn add_scaled(dst: &mut [f32], src: &[f32], scale: f32) {
    let n = dst.len().min(src.len());
    let mut i = 0usize;
    let scale_v = _mm256_set1_ps(scale);
    while i + 8 <= n {
        let d = _mm256_loadu_ps(dst.as_ptr().add(i));
        let s = _mm256_loadu_ps(src.as_ptr().add(i));
        let mul = _mm256_mul_ps(s, scale_v);
        let out = _mm256_add_ps(d, mul);
        _mm256_storeu_ps(dst.as_mut_ptr().add(i), out);
        i += 8;
    }
    while i < n {
        dst[i] += src[i] * scale;
        i += 1;
    }
}

pub unsafe fn abs_diff(dst: &mut [f32], src: &[f32], center: f32) {
    let n = dst.len().min(src.len());
    let mut i = 0usize;
    let center_v = _mm256_set1_ps(center);
    let sign_mask = _mm256_set1_ps(-0.0);
    while i + 8 <= n {
        let s = _mm256_loadu_ps(src.as_ptr().add(i));
        let diff = _mm256_sub_ps(s, center_v);
        let abs = _mm256_andnot_ps(sign_mask, diff);
        _mm256_storeu_ps(dst.as_mut_ptr().add(i), abs);
        i += 8;
    }
    while i < n {
        dst[i] = (src[i] - center).abs();
        i += 1;
    }
}
