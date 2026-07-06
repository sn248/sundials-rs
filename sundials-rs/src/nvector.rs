use std::ptr::NonNull;
use sundials_rs_sys as sys;
use crate::error::SundialsError;

/// A serial N_Vector wrapping `N_VNew_Serial`.
///
/// Gives safe `&[f64]` / `&mut [f64]` views onto the underlying data.
pub struct NVector {
    ptr: NonNull<sys::_generic_N_Vector>,
    len: usize,
}

impl NVector {
    /// Allocate a new zero-initialised serial N_Vector of length `n`.
    pub(crate) fn new(n: usize, ctx: sys::SUNContext) -> Result<Self, SundialsError> {
        let ptr = unsafe { sys::N_VNew_Serial(n as sys::sunindextype, ctx) };
        let ptr = NonNull::new(ptr).ok_or(SundialsError::NVectorCreate)?;
        Ok(Self { ptr, len: n })
    }

    /// Create from a Rust slice — copies the values in.
    pub(crate) fn from_slice(data: &[f64], ctx: sys::SUNContext) -> Result<Self, SundialsError> {
        let mut v = Self::new(data.len(), ctx)?;
        v.as_mut_slice().copy_from_slice(data);
        Ok(v)
    }

    /// Raw pointer (needed when passing to SUNDIALS functions).
    pub(crate) fn as_ptr(&self) -> sys::N_Vector {
        self.ptr.as_ptr()
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn as_slice(&self) -> &[f64] {
        unsafe {
            let data = sys::N_VGetArrayPointer(self.ptr.as_ptr());
            std::slice::from_raw_parts(data, self.len)
        }
    }

    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        unsafe {
            let data = sys::N_VGetArrayPointer(self.ptr.as_ptr());
            std::slice::from_raw_parts_mut(data, self.len)
        }
    }
}

impl Drop for NVector {
    fn drop(&mut self) {
        unsafe { sys::N_VDestroy(self.ptr.as_ptr()) };
    }
}

unsafe impl Send for NVector {}
