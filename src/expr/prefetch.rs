#[inline]
pub fn prefetch_read(_ptr: *const u8) {
    #[cfg(all(feature = "prefetch", target_arch = "x86_64"))]
    unsafe {
        use core::arch::x86_64::_MM_HINT_T0;
        use core::arch::x86_64::_mm_prefetch;
        _mm_prefetch(_ptr as *const i8, _MM_HINT_T0);
    }
}

#[inline]
pub fn prefetch_write(_ptr: *const u8) {
    #[cfg(all(feature = "prefetch", target_arch = "x86_64"))]
    unsafe {
        use core::arch::x86_64::_MM_HINT_T0;
        use core::arch::x86_64::_mm_prefetch;
        _mm_prefetch(_ptr as *const i8, _MM_HINT_T0);
    }
}
