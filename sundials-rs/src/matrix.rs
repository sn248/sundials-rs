use std::ptr::NonNull;
use sundials_rs_sys as sys;
use crate::error::SundialsError;

/// Wraps a dense `SUNMatrix`.
pub struct DenseMatrix {
    pub(crate) ptr: NonNull<sys::_generic_SUNMatrix>,
}

impl DenseMatrix {
    /// Allocate an `m × n` dense matrix.
    pub(crate) fn new(m: usize, n: usize, ctx: sys::SUNContext) -> Result<Self, SundialsError> {
        let ptr = unsafe {
            sys::SUNDenseMatrix(
                m as sys::sunindextype,
                n as sys::sunindextype,
                ctx,
            )
        };
        let ptr = NonNull::new(ptr).ok_or(SundialsError::MatrixCreate)?;
        Ok(Self { ptr })
    }
}

impl Drop for DenseMatrix {
    fn drop(&mut self) {
        unsafe { sys::SUNMatDestroy(self.ptr.as_ptr()) };
    }
}

unsafe impl Send for DenseMatrix {}
