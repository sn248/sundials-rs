//! Safe wrapper around IDAS — IDA with forward sensitivity analysis.
//!
//! IDAS is a superset of IDA.  Use it when you need forward sensitivities
//! `∂y/∂p` in addition to the DAE solution.
//!
//! # DQ requirement
//! When using internal difference-quotient (DQ) sensitivity approximation,
//! the residual closure **must read parameters through the raw pointer**
//! returned by [`IdasBuilder::params_ptr`], not from Rust constants.
//! CVODES/IDAS temporarily perturbs `p[i]` to approximate the sensitivity
//! RHS; if the closure ignores `p`, all sensitivities will be zero.
//!
//! # Example
//! ```no_run
//! use sundials_rs::idas::IdasBuilder;
//!
//! // dy/dt = -k*y  as DAE:  F = y' + k*y = 0
//! let y0  = vec![1.0_f64];
//! let yp0 = vec![-0.5_f64]; // y'(0) = -k*y(0) = -0.5
//! let p   = vec![0.5_f64];  // k
//! let s0  = vec![vec![0.0_f64]]; // ∂y/∂k = 0 at t=0
//!
//! let mut b = IdasBuilder::new(&y0, &yp0)
//!     .rtol(1e-8)
//!     .atol(1e-10)
//!     .with_forward_sensitivity(p, s0);
//! let p_ptr = b.params_ptr();
//!
//! let mut solver = b.build(move |_t, y, yp, res| {
//!     let k = unsafe { *p_ptr };
//!     res[0] = yp[0] + k * y[0];
//!     Ok(())
//! }).unwrap();
//!
//! let (t, y, _yp) = solver.step(1.0).unwrap();
//! let y = y.to_vec(); // copy releases the borrow so sensitivities() can take &mut self
//! let sens = solver.sensitivities().unwrap();
//! println!("y({t:.1}) = {:.6}  ∂y/∂k = {:.6}", y[0], sens[0][0]);
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

struct UserData<F> {
    res: F,
    neq: usize,
}

unsafe extern "C" fn res_trampoline<F>(
    t: sys::sunrealtype,
    y: sys::N_Vector,
    yp: sys::N_Vector,
    rr: sys::N_Vector,
    user_data: *mut c_void,
) -> c_int
where
    F: Fn(f64, &[f64], &[f64], &mut [f64]) -> Result<(), i32>,
{
    let ud = &*(user_data as *const UserData<F>);
    let n = ud.neq;
    let y_s  = std::slice::from_raw_parts(sys::N_VGetArrayPointer(y),  n);
    let yp_s = std::slice::from_raw_parts(sys::N_VGetArrayPointer(yp), n);
    let rr_s = std::slice::from_raw_parts_mut(sys::N_VGetArrayPointer(rr), n);
    match (ud.res)(t, y_s, yp_s, rr_s) {
        Ok(()) => 0,
        Err(f) => f,
    }
}

pub struct IdasBuilder {
    y0: Vec<f64>,
    yp0: Vec<f64>,
    t0: f64,
    rtol: f64,
    atol: f64,
    max_steps: Option<i64>,
    // Box<[f64]> gives a stable heap address so the raw pointer passed to
    // IDASetSensParams stays valid after build() moves it into IdasSolver.
    p: Box<[f64]>,
    s0: Vec<Vec<f64>>,
}

impl IdasBuilder {
    /// Create a new builder.  Prefer calling [`IdasBuilder::new`] directly
    /// over `IdasSolver::new` to avoid type-inference ambiguity.
    pub fn new(y0: &[f64], yp0: &[f64]) -> Self {
        assert_eq!(y0.len(), yp0.len(), "y0 and yp0 must have the same length");
        Self {
            y0: y0.to_vec(),
            yp0: yp0.to_vec(),
            t0: 0.0,
            rtol: 1e-6,
            atol: 1e-9,
            max_steps: None,
            p: Box::new([]),
            s0: Vec::new(),
        }
    }

    pub fn rtol(mut self, rtol: f64) -> Self { self.rtol = rtol; self }
    pub fn atol(mut self, atol: f64) -> Self { self.atol = atol; self }
    pub fn t0(mut self, t0: f64) -> Self { self.t0 = t0; self }

    /// Maximum number of internal steps before returning an error (default 500).
    pub fn max_steps(mut self, n: i64) -> Self { self.max_steps = Some(n); self }

    /// Enable forward sensitivity analysis.
    ///
    /// `p` — initial parameter values.
    /// `s0` — initial sensitivities `∂y/∂p_i` at `t0` (usually all zeros).
    pub fn with_forward_sensitivity(mut self, p: Vec<f64>, s0: Vec<Vec<f64>>) -> Self {
        self.p = p.into_boxed_slice(); // heap-stable allocation
        self.s0 = s0;
        self
    }

    /// Raw pointer to the parameter array.
    ///
    /// Capture this **before** calling `build()` (which consumes the builder).
    /// The heap allocation is moved — not re-allocated — into the resulting
    /// `IdasSolver`, so the pointer remains valid for the solver's lifetime.
    pub fn params_ptr(&self) -> *const f64 {
        self.p.as_ptr()
    }

    pub fn build<F>(self, res: F) -> Result<IdasSolver<F>, SundialsError>
    where
        F: Fn(f64, &[f64], &[f64], &mut [f64]) -> Result<(), i32>,
    {
        let neq = self.y0.len();

        let ctx = SunContext::new()?;

        let mem = unsafe { sys::IDACreate(ctx.raw()) };
        if mem.is_null() {
            return Err(SundialsError::Memory("IDACreate (IDAS)"));
        }

        let y  = NVector::from_slice(&self.y0,  ctx.raw())?;
        let yp = NVector::from_slice(&self.yp0, ctx.raw())?;

        let user_data = Box::new(UserData { res, neq });
        let ud_ptr = &*user_data as *const UserData<F> as *mut c_void;

        check_flag(
            unsafe { sys::IDAInit(mem, Some(res_trampoline::<F>), self.t0, y.as_ptr(), yp.as_ptr()) },
            "IDAS", "IDAInit",
        )?;
        check_flag(
            unsafe { sys::IDASStolerances(mem, self.rtol, self.atol) },
            "IDAS", "IDASStolerances",
        )?;
        check_flag(
            unsafe { sys::IDASetUserData(mem, ud_ptr) },
            "IDAS", "IDASetUserData",
        )?;

        let matrix = DenseMatrix::new(neq, neq, ctx.raw())?;
        let ls = LinearSolver::dense(&y, &matrix, ctx.raw())?;
        check_flag(
            unsafe { sys::IDASetLinearSolver(mem, ls.ptr.as_ptr(), matrix.ptr.as_ptr()) },
            "IDAS", "IDASetLinearSolver",
        )?;

        if let Some(n) = self.max_steps {
            check_flag(
                unsafe { sys::IDASetMaxNumSteps(mem, n) },
                "IDAS", "IDASetMaxNumSteps",
            )?;
        }

        let mut sens_vecs: Vec<NVector> = Vec::new();
        let mut sens_vecs_yp: Vec<NVector> = Vec::new();
        let ns = self.p.len();

        if ns > 0 {
            for s in &self.s0 {
                sens_vecs.push(NVector::from_slice(s, ctx.raw())?);
                sens_vecs_yp.push(NVector::new(neq, ctx.raw())?);
            }
            let mut sv_ptrs: Vec<sys::N_Vector> =
                sens_vecs.iter().map(|v| v.as_ptr()).collect();
            let mut svp_ptrs: Vec<sys::N_Vector> =
                sens_vecs_yp.iter().map(|v| v.as_ptr()).collect();

            check_flag(
                unsafe {
                    sys::IDASensInit(
                        mem, ns as c_int, sys::IDA_SIMULTANEOUS as c_int,
                        None, sv_ptrs.as_mut_ptr(), svp_ptrs.as_mut_ptr(),
                    )
                },
                "IDAS", "IDASensInit",
            )?;
            check_flag(
                unsafe { sys::IDASensEEtolerances(mem) },
                "IDAS", "IDASensEEtolerances",
            )?;
        }

        // Move the heap-stable Box into the solver first — its address is
        // unchanged.  Then pass that address to IDASetSensParams so the pointer
        // outlives this call.
        let p = self.p;

        let solver = IdasSolver {
            mem,
            y,
            yp,
            _matrix: matrix,
            _ls: ls,
            _user_data: user_data,
            sens_vecs,
            sens_vecs_yp,
            p,
            _ctx: ctx,
            t: self.t0,
        };

        if solver.ns() > 0 {
            check_flag(
                unsafe {
                    sys::IDASetSensParams(
                        solver.mem,
                        solver.p.as_ptr() as *mut f64,
                        std::ptr::null_mut(),
                        std::ptr::null_mut(),
                    )
                },
                "IDAS", "IDASetSensParams",
            )?;
        }

        Ok(solver)
    }
}

pub struct IdasSolver<F> {
    mem: *mut c_void,
    y: NVector,
    yp: NVector,
    _matrix: DenseMatrix,
    _ls: LinearSolver,
    _user_data: Box<UserData<F>>,
    sens_vecs: Vec<NVector>,
    #[allow(dead_code)]
    sens_vecs_yp: Vec<NVector>,
    /// Parameter array — heap-stable Box so IDASetSensParams' raw ptr stays valid.
    p: Box<[f64]>,
    // _ctx must be dropped LAST — declare after everything that uses it.
    _ctx: SunContext,
    t: f64,
}

impl<F> IdasSolver<F> {
    fn ns(&self) -> usize { self.sens_vecs.len() }
}

impl<F> IdasSolver<F>
where
    F: Fn(f64, &[f64], &[f64], &mut [f64]) -> Result<(), i32>,
{
    /// Start building an IDAS solver.
    pub fn new(y0: &[f64], yp0: &[f64]) -> IdasBuilder {
        IdasBuilder::new(y0, yp0)
    }

    /// Ask IDAS to correct initial conditions to satisfy `F(t0, y0, yp0) = 0`.
    ///
    /// `tout1` is the first output time — used only to determine the direction
    /// of integration, not as a stopping point.
    pub fn calc_ic(&mut self, tout1: f64) -> Result<(), SundialsError> {
        check_flag(
            unsafe { sys::IDACalcIC(self.mem, sys::IDA_YA_YDP_INIT as c_int, tout1) },
            "IDAS", "IDACalcIC",
        )
    }

    /// Advance to `tout`.  Returns `(t_reached, y, yp)`.
    pub fn step(&mut self, tout: f64) -> Result<(f64, &[f64], &[f64]), SundialsError> {
        let mut t_out: sys::sunrealtype = 0.0;
        let flag = unsafe {
            sys::IDASolve(
                self.mem, tout, &mut t_out,
                self.y.as_ptr(), self.yp.as_ptr(),
                sys::IDA_NORMAL as c_int,
            )
        };
        check_flag(flag, "IDAS", "IDASolve")?;
        self.t = t_out;
        Ok((t_out, self.y.as_slice(), self.yp.as_slice()))
    }

    /// Retrieve forward sensitivities `∂y/∂p_i` after the most recent step.
    ///
    /// Call immediately after `step()`.  Copy `y` to a `Vec` first if you also
    /// need the state (both methods require `&mut self`).
    pub fn sensitivities(&mut self) -> Result<Vec<&[f64]>, SundialsError> {
        if self.sens_vecs.is_empty() {
            return Ok(Vec::new());
        }
        let mut sv_ptrs: Vec<sys::N_Vector> =
            self.sens_vecs.iter().map(|v| v.as_ptr()).collect();
        check_flag(
            unsafe { sys::IDAGetSens(self.mem, &mut self.t, sv_ptrs.as_mut_ptr()) },
            "IDAS", "IDAGetSens",
        )?;
        Ok(self.sens_vecs.iter().map(|v| v.as_slice()).collect())
    }

    pub fn t(&self) -> f64 { self.t }
    pub fn y(&self) -> &[f64] { self.y.as_slice() }
    pub fn yp(&self) -> &[f64] { self.yp.as_slice() }
}

impl<F> Drop for IdasSolver<F> {
    fn drop(&mut self) {
        unsafe { sys::IDAFree(&mut self.mem) };
    }
}

unsafe impl<F: Send> Send for IdasSolver<F> {}
