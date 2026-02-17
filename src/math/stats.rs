//! Robust statistical primitives.
//!
//! Note: Functions may reorder the input slice.

use crate::simd;

pub fn trimmed_mean(values: &mut [f32], p: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = values.len();
    let k = (p * n as f32).floor() as usize;
    if n < 2 * k + 1 {
        return mean(values);
    }
    mean(&values[k..(n - k)])
}

pub fn median(values: &mut [f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        let a = values[n / 2 - 1];
        let b = values[n / 2];
        (a + b) / 2.0
    }
}

pub fn mad(values: &mut [f32], median_val: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let tmp = values.to_vec();
    simd::abs_diff(values, &tmp, median_val);
    median(values)
}

pub fn robust_z(x: f32, median_val: f32, mad_val: f32) -> f32 {
    if mad_val == 0.0 {
        return 0.0;
    }
    (x - median_val) / (1.4826 * mad_val)
}

fn mean(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    simd::sum_f32(values) / values.len() as f32
}
