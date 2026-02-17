#![allow(dead_code)]

pub fn sum_f32(slice: &[f32]) -> f32 {
    let mut sum = 0.0;
    for v in slice {
        sum += *v;
    }
    sum
}

pub fn add_scaled(dst: &mut [f32], src: &[f32], scale: f32) {
    let n = dst.len().min(src.len());
    for i in 0..n {
        dst[i] += src[i] * scale;
    }
}

pub fn abs_diff(dst: &mut [f32], src: &[f32], center: f32) {
    let n = dst.len().min(src.len());
    for i in 0..n {
        dst[i] = (src[i] - center).abs();
    }
}
