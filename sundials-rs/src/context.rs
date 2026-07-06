//! RAII wrapper for `SUNContext` — required by all SUNDIALS 6.x+ creation calls.

use std::ptr::NonNull;
use sundials_rs_sys as sys;
use crate::error::SundialsError;

/// Owns a `SUNContext` and frees it on drop.
///
/// Every solver and vector creation function in SUNDIALS >= 6.0 requires a
/// context.  Create one per solver instance and pass `raw()` wherever needed.
pub struct SunContext {
    ptr: NonNull<sys::SUNContext_>,
}

impl SunContext {
    pub fn new() -> Result<Self, SundialsError> {
        let mut ctx: sys::SUNContext = std::ptr::null_mut();
        // SUN_COMM_NULL == 0 (int) on non-MPI builds, MPI_COMM_NULL on MPI builds.
        // Cast 0 to whatever SUNComm resolves to in the generated bindings.
        let comm: sys::SUNComm = 0;
        let flag = unsafe { sys::SUNContext_Create(comm, &mut ctx) };
        if flag != 0 || ctx.is_null() {
            return Err(SundialsError::Memory("SUNContext_Create"));
        }
        Ok(Self {
            ptr: unsafe { NonNull::new_unchecked(ctx) },
        })
    }

    /// Raw pointer — pass this to SUNDIALS creation functions.
    pub fn raw(&self) -> sys::SUNContext {
        self.ptr.as_ptr()
    }
}

impl Drop for SunContext {
    fn drop(&mut self) {
        let mut ptr = self.ptr.as_ptr();
        unsafe { sys::SUNContext_Free(&mut ptr) };
    }
}

unsafe impl Send for SunContext {}
