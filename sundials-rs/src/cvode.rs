//! Safe wrapper around CVODE — explicit ODE IVP solver.
//!
//! # Method choices
//! - [`Method::Adams`] — non-stiff problems (Adams-Moulton)
//! - [`Method::BDF`]   — stiff problems (Backward Differentiation Formula)
//!
//! # Example
//! ```no_run
//! use sundials_rs::cvode::{CvodeBuilder, Method};
//!
//! // Exponential decay: dy/dt = -0.1 * y,  y(0) = 10.0
//! let y0 = vec![10.0_f64];
//! let mut solver = CvodeBuilder::new(Method::BDF, &y0)
//!     .rtol(1e-8)
//!     .atol(1e-10)
//!     .build(|_t, y, ydot| { ydot[0] = -0.1 * y[0]; Ok(()) })
//!     .unwrap();
//!
//! let mut t = 0.0_f64;
//! for _ in 0..10 {
//!     let (tnew, y) = solver.step(t + 1.0).unwrap();
//!     println!("y({tnew:.1}) = {:.6}", y[0]);
//!     t = tnew;
//! }
//! ```

use std::os::raw::{c_int, c_void};
use sundials_rs_sys as sys;

use crate::{
    context::SunContext,
    error::{check_flag, SundialsError},
    linear_solver::LinearSolver,
    matrix::DenseMatrix,
    nvector::NVector,
};

// ── Integration method ───────────────────────────────────────────────────────

/// ODE integration method.
#[derive(Debug, Clone, Copy)]
pub enum Method {
    /// Adams-Moulton — best for non-stiff problems.
    Adams,
    /// Backward Differentiation Formula — best for stiff problems.
    BDF,
}

impl Method {
    fn as_c_int(self) -> c_int {
        match self {
            Method::Adams => sys::CV_ADAMS as c_int,
            Method::BDF => sys::CV_BDF as c_int,
        }
    }
}

// ── Solver statistics ────────────────────────────────────────────────────────

/// Integration statistics returned by [`CvodeSolver::stats`].
#[derive(Debug, Clone)]
pub struct CvodeStats {
    pub num_steps: i64,
    pub num_rhs_evals: i64,
    pub num_err_test_fails: i64,
    pub last_order: i32,
    pub current_time: f64,
}

// ── User data stored behind the void* ────────────────────────────────────────

struct UserData<F> {
    rhs: F,
    neq: usize,
}

// ── Trampoline (C → Rust) ─────────────────────────────────────────────────────

unsafe extern "C" fn rhs_trampoline<F>(
    t: sys::sunrealtype,
    y: sys::N_Vector,
    ydot: sys::N_Vector,
    user_data: *mut c_void,
) -> c_int
where
    F: Fn(f64, &[f64], &mut [f64]) -> Result<(), i32>,
{
    let ud = &*(user_data as *const UserData<F>);
    let n = ud.neq;
    let y_s    = std::slice::from_raw_parts(sys::N_VGetArrayPointer(y), n);
    let ydot_s = std::slice::from_raw_parts_mut(sys::N_VGetArrayPointer(ydot), n);
    match (ud.rhs)(t, y_s, ydot_s) {
        Ok(()) => 0,
        Err(flag) => flag,
    }
}

// ── Builder ───────────────────────────────────────────────────────────────────

/// Tolerance mode chosen via builder methods.
enum AtolMode {
    Scalar(f64),
    Vector(Vec<f64>),
}

pub struct CvodeBuilder {
    method: Method,
    y0: Vec<f64>,
    t0: f64,
    rtol: f64,
    atol: AtolMode,
    max_steps: Option<i64>,
}

impl CvodeBuilder {
    /// Create a new builder.  Prefer calling [`CvodeBuilder::new`] directly
    /// over `CvodeSolver::new` to avoid type-inference ambiguity.
    pub fn new(method: Method, y0: &[f64]) -> Self {
        Self {
            method,
            y0: y0.to_vec(),
            t0: 0.0,
            rtol: 1e-6,
            atol: AtolMode::Scalar(1e-9),
            max_steps: None,
        }
    }

    pub fn rtol(mut self, rtol: f64) -> Self { self.rtol = rtol; self }

    /// Scalar absolute tolerance (same for all components).
    pub fn atol(mut self, atol: f64) -> Self {
        self.atol = AtolMode::Scalar(atol);
        self
    }

    /// Per-component absolute tolerances (length must equal `y0`).
    pub fn atol_vec(mut self, atol: Vec<f64>) -> Self {
        self.atol = AtolMode::Vector(atol);
        self
    }

    pub fn t0(mut self, t0: f64) -> Self { self.t0 = t0; self }

    /// Maximum number of internal steps before returning an error (default 500).
    pub fn max_steps(mut self, n: i64) -> Self { self.max_steps = Some(n); self }

    /// Finalise the builder with the RHS closure and return a ready solver.
    ///
    /// `rhs(t, y, ydot)` fills `ydot` and returns `Ok(())` on success or
    /// `Err(flag)` (non-zero) for a recoverable (`> 0`) / unrecoverable (`< 0`)
    /// error.
    pub fn build<F>(self, rhs: F) -> Result<CvodeSolver<F>, SundialsError>
    where
        F: Fn(f64, &[f64], &mut [f64]) -> Result<(), i32>,
    {
        let neq = self.y0.len();

        // ── 1. SUNContext ────────────────────────────────────────────────────
        let ctx = SunContext::new()?;

        // ── 2. CVodeMem ──────────────────────────────────────────────────────
        let mem = unsafe { sys::CVodeCreate(self.method.as_c_int(), ctx.raw()) };
        if mem.is_null() {
            return Err(SundialsError::Memory("CVodeCreate"));
        }

        // ── 3. Initial state vector ──────────────────────────────────────────
        let y = NVector::from_slice(&self.y0, ctx.raw())?;

        // ── 4. Pin user_data ─────────────────────────────────────────────────
        let user_data = Box::new(UserData { rhs, neq });
        let ud_ptr = &*user_data as *const UserData<F> as *mut c_void;

        // ── 5. CVodeInit ─────────────────────────────────────────────────────
        check_flag(
            unsafe { sys::CVodeInit(mem, Some(rhs_trampoline::<F>), self.t0, y.as_ptr()) },
            "CVODE", "CVodeInit",
        )?;

        // ── 6. Tolerances ────────────────────────────────────────────────────
        match self.atol {
            AtolMode::Scalar(atol) => {
                check_flag(
                    unsafe { sys::CVodeSStolerances(mem, self.rtol, atol) },
                    "CVODE", "CVodeSStolerances",
                )?;
            }
            AtolMode::Vector(ref v) => {
                let atol_vec = NVector::from_slice(v, ctx.raw())?;
                check_flag(
                    unsafe { sys::CVodeSVtolerances(mem, self.rtol, atol_vec.as_ptr()) },
                    "CVODE", "CVodeSVtolerances",
                )?;
                // atol_vec can be dropped here — CVODE has copied the values.
            }
        }

        // ── 7. user_data ─────────────────────────────────────────────────────
        check_flag(
            unsafe { sys::CVodeSetUserData(mem, ud_ptr) },
            "CVODE", "CVodeSetUserData",
        )?;

        // ── 8. Optional max steps ─────────────────────────────────────────────
        if let Some(n) = self.max_steps {
            check_flag(
                unsafe { sys::CVodeSetMaxNumSteps(mem, n) },
                "CVODE", "CVodeSetMaxNumSteps",
            )?;
        }

        // ── 9. Dense linear solver ────────────────────────────────────────────
        let matrix = DenseMatrix::new(neq, neq, ctx.raw())?;
        let ls = LinearSolver::dense(&y, &matrix, ctx.raw())?;
        check_flag(
            unsafe { sys::CVodeSetLinearSolver(mem, ls.ptr.as_ptr(), matrix.ptr.as_ptr()) },
            "CVODE", "CVodeSetLinearSolver",
        )?;

        Ok(CvodeSolver {
            mem,
            y,
            _matrix: matrix,
            _ls: ls,
            _user_data: user_data,
            _ctx: ctx,
            t: self.t0,
        })
    }
}

// ── Solver ────────────────────────────────────────────────────────────────────

pub struct CvodeSolver<F> {
    mem: *mut c_void,
    y: NVector,
    _matrix: DenseMatrix,
    _ls: LinearSolver,
    _user_data: Box<UserData<F>>,
    // ctx must be dropped LAST — declare it after everything that uses it.
    _ctx: SunContext,
    t: f64,
}

impl<F> CvodeSolver<F>
where
    F: Fn(f64, &[f64], &mut [f64]) -> Result<(), i32>,
{
    /// Start building a solver for the given method and initial condition.
    pub fn new(method: Method, y0: &[f64]) -> CvodeBuilder {
        CvodeBuilder::new(method, y0)
    }

    /// Advance the solution to `tout` (CV_NORMAL mode).
    ///
    /// Returns `(t_reached, y_slice)`.
    pub fn step(&mut self, tout: f64) -> Result<(f64, &[f64]), SundialsError> {
        let mut t_out: sys::sunrealtype = 0.0;
        let flag = unsafe {
            sys::CVode(self.mem, tout, self.y.as_ptr(), &mut t_out, sys::CV_NORMAL as c_int)
        };
        check_flag(flag, "CVODE", "CVode")?;
        self.t = t_out;
        Ok((t_out, self.y.as_slice()))
    }

    /// Reinitialise at `(t0, y0)` without reallocating.
    pub fn reinit(&mut self, t0: f64, y0: &[f64]) -> Result<(), SundialsError> {
        self.y.as_mut_slice().copy_from_slice(y0);
        check_flag(
            unsafe { sys::CVodeReInit(self.mem, t0, self.y.as_ptr()) },
            "CVODE", "CVodeReInit",
        )?;
        self.t = t0;
        Ok(())
    }

    /// Retrieve integration statistics.
    pub fn stats(&self) -> Result<CvodeStats, SundialsError> {
        let mut nsteps: i64     = 0;
        let mut nrhs: i64       = 0;
        let mut netfails: i64   = 0;
        let mut order: c_int    = 0;
        let mut tcur: f64       = 0.0;

        check_flag(
            unsafe { sys::CVodeGetNumSteps(self.mem, &mut nsteps) },
            "CVODE", "CVodeGetNumSteps",
        )?;
        check_flag(
            unsafe { sys::CVodeGetNumRhsEvals(self.mem, &mut nrhs) },
            "CVODE", "CVodeGetNumRhsEvals",
        )?;
        check_flag(
            unsafe { sys::CVodeGetNumErrTestFails(self.mem, &mut netfails) },
            "CVODE", "CVodeGetNumErrTestFails",
        )?;
        check_flag(
            unsafe { sys::CVodeGetLastOrder(self.mem, &mut order) },
            "CVODE", "CVodeGetLastOrder",
        )?;
        check_flag(
            unsafe { sys::CVodeGetCurrentTime(self.mem, &mut tcur) },
            "CVODE", "CVodeGetCurrentTime",
        )?;

        Ok(CvodeStats {
            num_steps: nsteps,
            num_rhs_evals: nrhs,
            num_err_test_fails: netfails,
            last_order: order,
            current_time: tcur,
        })
    }

    pub fn t(&self) -> f64 { self.t }
    pub fn y(&self) -> &[f64] { self.y.as_slice() }
}

impl<F> Drop for CvodeSolver<F> {
    fn drop(&mut self) {
        unsafe { sys::CVodeFree(&mut self.mem) };
    }
}

unsafe impl<F: Send> Send for CvodeSolver<F> {}
