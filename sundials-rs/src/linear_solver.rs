use std::ptr::NonNull;
use sundials_rs_sys as sys;
use crate::{error::SundialsError, matrix::DenseMatrix, nvector::NVector};

/// Wraps a `SUNLinearSolver`.
pub struct LinearSolver {
    pub(crate) ptr: NonNull<sys::_generic_SUNLinearSolver>,
}

impl LinearSolver {
    /// Create a dense direct linear solver paired with `y` and `A`.
    pub(crate) fn dense(y: &NVector, a: &DenseMatrix, ctx: sys::SUNContext) -> Result<Self, SundialsError> {
        let ptr = unsafe {
            sys::SUNLinSol_Dense(y.as_ptr(), a.ptr.as_ptr(), ctx)
        };
        let ptr = NonNull::new(ptr).ok_or(SundialsError::LinearSolverCreate)?;
        Ok(Self { ptr })
    }
}

impl Drop for LinearSolver {
    fn drop(&mut self) {
        unsafe { sys::SUNLinSolFree(self.ptr.as_ptr()) };
    }
}

unsafe impl Send for LinearSolver {}
