#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
mod avx2;
#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
mod neon;
mod scalar;

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
pub fn backend_name() -> &'static str {
    "avx2"
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
pub fn backend_name() -> &'static str {
    "neon"
}

#[cfg(not(any(
    all(target_arch = "x86_64", target_feature = "avx2"),
    all(target_arch = "aarch64", target_feature = "neon")
)))]
pub fn backend_name() -> &'static str {
    "scalar"
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
pub fn sum_f32(slice: &[f32]) -> f32 {
    unsafe { avx2::sum_f32(slice) }
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
pub fn sum_f32(slice: &[f32]) -> f32 {
    unsafe { neon::sum_f32(slice) }
}

#[cfg(not(any(
    all(target_arch = "x86_64", target_feature = "avx2"),
    all(target_arch = "aarch64", target_feature = "neon")
)))]
pub fn sum_f32(slice: &[f32]) -> f32 {
    scalar::sum_f32(slice)
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
pub fn add_scaled(dst: &mut [f32], src: &[f32], scale: f32) {
    unsafe { avx2::add_scaled(dst, src, scale) }
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
pub fn add_scaled(dst: &mut [f32], src: &[f32], scale: f32) {
    unsafe { neon::add_scaled(dst, src, scale) }
}

#[cfg(not(any(
    all(target_arch = "x86_64", target_feature = "avx2"),
    all(target_arch = "aarch64", target_feature = "neon")
)))]
pub fn add_scaled(dst: &mut [f32], src: &[f32], scale: f32) {
    scalar::add_scaled(dst, src, scale)
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
pub fn abs_diff(dst: &mut [f32], src: &[f32], center: f32) {
    unsafe { avx2::abs_diff(dst, src, center) }
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
pub fn abs_diff(dst: &mut [f32], src: &[f32], center: f32) {
    unsafe { neon::abs_diff(dst, src, center) }
}

#[cfg(not(any(
    all(target_arch = "x86_64", target_feature = "avx2"),
    all(target_arch = "aarch64", target_feature = "neon")
)))]
pub fn abs_diff(dst: &mut [f32], src: &[f32], center: f32) {
    scalar::abs_diff(dst, src, center)
}
