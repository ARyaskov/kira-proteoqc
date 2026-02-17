#![cfg(all(target_arch = "aarch64", target_feature = "neon"))]
#![allow(unsafe_op_in_unsafe_fn)]

use core::arch::aarch64::*;

pub unsafe fn sum_f32(slice: &[f32]) -> f32 {
    let mut i = 0usize;
    let mut acc = vdupq_n_f32(0.0);
    while i + 4 <= slice.len() {
        let v = vld1q_f32(slice.as_ptr().add(i));
        acc = vaddq_f32(acc, v);
        i += 4;
    }
    let mut lanes = [0f32; 4];
    vst1q_f32(lanes.as_mut_ptr(), acc);
    let mut sum = 0.0f32;
    for j in 0..4 {
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
    let scale_v = vdupq_n_f32(scale);
    while i + 4 <= n {
        let d = vld1q_f32(dst.as_ptr().add(i));
        let s = vld1q_f32(src.as_ptr().add(i));
        let mul = vmulq_f32(s, scale_v);
        let out = vaddq_f32(d, mul);
        vst1q_f32(dst.as_mut_ptr().add(i), out);
        i += 4;
    }
    while i < n {
        dst[i] += src[i] * scale;
        i += 1;
    }
}

pub unsafe fn abs_diff(dst: &mut [f32], src: &[f32], center: f32) {
    let n = dst.len().min(src.len());
    let mut i = 0usize;
    let center_v = vdupq_n_f32(center);
    while i + 4 <= n {
        let s = vld1q_f32(src.as_ptr().add(i));
        let diff = vsubq_f32(s, center_v);
        let abs = vabsq_f32(diff);
        vst1q_f32(dst.as_mut_ptr().add(i), abs);
        i += 4;
    }
    while i < n {
        dst[i] = (src[i] - center).abs();
        i += 1;
    }
}
