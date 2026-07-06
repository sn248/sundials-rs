//! Safe wrapper around IDA — DAE (Differential-Algebraic Equation) solver.
//!
//! IDA solves **implicit** systems of the form:
//!
//! ```text
//! F(t, y, y') = 0,    y(t₀) = y₀,   y'(t₀) = y'₀
//! ```
//!
//! This covers both:
//! - **Implicit ODEs** — every component has a `y'` term (no algebraic variables).
//! - **True DAEs** — some components are purely algebraic (no `y'` term),
//!   enforcing constraints between state variables.
//!
//! If you only need an explicit ODE `y' = f(t, y)`, prefer [`crate::cvode`].
//! Use IDA when the equations are naturally in implicit form, or when algebraic
//! constraints must be maintained exactly.
//!
//! # Initial conditions
//!
//! IDA requires both `y₀` and `y'₀` to be **consistent**, meaning
//! `F(t₀, y₀, y'₀) = 0`.  For implicit ODEs you can compute `y'₀` directly
//! from the equations.  For DAEs with algebraic variables, call
//! [`IdaSolver::calc_ic`] to let IDA correct the initial conditions.
//!
//! # Example — harmonic oscillator as implicit ODE
//! ```no_run
//! use sundials_rs::ida::IdaBuilder;
//!
//! // d²x/dt² = -x  written as first-order system:
//! //   y[0] = x,  y[1] = x'
//! // Residual form:  F = [y'[0] - y[1],  y'[1] + y[0]]
//! let y0  = vec![1.0_f64, 0.0];   // x(0)=1, x'(0)=0
//! let yp0 = vec![0.0_f64, -1.0];  // consistent: y'[0]=y[1]=0, y'[1]=-y[0]=-1
//!
//! let mut solver = IdaBuilder::new(&y0, &yp0)
//!     .rtol(1e-8)
//!     .atol(1e-10)
//!     .build(|_t, y, yp, res| {
//!         res[0] = yp[0] - y[1];
//!         res[1] = yp[1] + y[0];
//!         Ok(())
//!     })
//!     .unwrap();
//!
//! let (t, y, _yp) = solver.step(1.0).unwrap();
//! println!("y({t:.3}) = [{:.6}, {:.6}]  (exact: [{:.6}, {:.6}])",
//!          y[0], y[1], t.cos(), -t.sin());
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

pub struct IdaBuilder {
    y0: Vec<f64>,
    yp0: Vec<f64>,
    t0: f64,
    rtol: f64,
    atol: f64,
    max_steps: Option<i64>,
}

impl IdaBuilder {
    /// Create a new builder.  Prefer calling [`IdaBuilder::new`] directly
    /// over `IdaSolver::new` to avoid type-inference ambiguity.
    pub fn new(y0: &[f64], yp0: &[f64]) -> Self {
        assert_eq!(y0.len(), yp0.len(), "y0 and yp0 must have the same length");
        Self {
            y0: y0.to_vec(),
            yp0: yp0.to_vec(),
            t0: 0.0,
            rtol: 1e-6,
            atol: 1e-9,
            max_steps: None,
        }
    }

    pub fn rtol(mut self, rtol: f64) -> Self { self.rtol = rtol; self }
    pub fn atol(mut self, atol: f64) -> Self { self.atol = atol; self }
    pub fn t0(mut self, t0: f64) -> Self { self.t0 = t0; self }

    /// Maximum number of internal steps before returning an error (default 500).
    pub fn max_steps(mut self, n: i64) -> Self { self.max_steps = Some(n); self }

    pub fn build<F>(self, res: F) -> Result<IdaSolver<F>, SundialsError>
    where
        F: Fn(f64, &[f64], &[f64], &mut [f64]) -> Result<(), i32>,
    {
        let neq = self.y0.len();

        let ctx = SunContext::new()?;

        let mem = unsafe { sys::IDACreate(ctx.raw()) };
        if mem.is_null() {
            return Err(SundialsError::Memory("IDACreate"));
        }

        let y  = NVector::from_slice(&self.y0,  ctx.raw())?;
        let yp = NVector::from_slice(&self.yp0, ctx.raw())?;

        let user_data = Box::new(UserData { res, neq });
        let ud_ptr = &*user_data as *const UserData<F> as *mut c_void;

        check_flag(
            unsafe { sys::IDAInit(mem, Some(res_trampoline::<F>), self.t0, y.as_ptr(), yp.as_ptr()) },
            "IDA", "IDAInit",
        )?;
        check_flag(
            unsafe { sys::IDASStolerances(mem, self.rtol, self.atol) },
            "IDA", "IDASStolerances",
        )?;
        check_flag(
            unsafe { sys::IDASetUserData(mem, ud_ptr) },
            "IDA", "IDASetUserData",
        )?;

        let matrix = DenseMatrix::new(neq, neq, ctx.raw())?;
        let ls = LinearSolver::dense(&y, &matrix, ctx.raw())?;
        check_flag(
            unsafe { sys::IDASetLinearSolver(mem, ls.ptr.as_ptr(), matrix.ptr.as_ptr()) },
            "IDA", "IDASetLinearSolver",
        )?;

        if let Some(n) = self.max_steps {
            check_flag(
                unsafe { sys::IDASetMaxNumSteps(mem, n) },
                "IDA", "IDASetMaxNumSteps",
            )?;
        }

        Ok(IdaSolver {
            mem,
            y,
            yp,
            _matrix: matrix,
            _ls: ls,
            _user_data: user_data,
            _ctx: ctx,
            t: self.t0,
        })
    }
}

pub struct IdaSolver<F> {
    mem: *mut c_void,
    y: NVector,
    yp: NVector,
    _matrix: DenseMatrix,
    _ls: LinearSolver,
    _user_data: Box<UserData<F>>,
    _ctx: SunContext,
    t: f64,
}

impl<F> IdaSolver<F>
where
    F: Fn(f64, &[f64], &[f64], &mut [f64]) -> Result<(), i32>,
{
    pub fn new(y0: &[f64], yp0: &[f64]) -> IdaBuilder {
        IdaBuilder::new(y0, yp0)
    }

    /// Ask IDA to correct the initial conditions to satisfy the residual.
    pub fn calc_ic(&mut self, tout1: f64) -> Result<(), SundialsError> {
        check_flag(
            unsafe { sys::IDACalcIC(self.mem, sys::IDA_YA_YDP_INIT as c_int, tout1) },
            "IDA", "IDACalcIC",
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
        check_flag(flag, "IDA", "IDASolve")?;
        self.t = t_out;
        Ok((t_out, self.y.as_slice(), self.yp.as_slice()))
    }

    pub fn reinit(&mut self, t0: f64, y0: &[f64], yp0: &[f64]) -> Result<(), SundialsError> {
        self.y.as_mut_slice().copy_from_slice(y0);
        self.yp.as_mut_slice().copy_from_slice(yp0);
        check_flag(
            unsafe { sys::IDAReInit(self.mem, t0, self.y.as_ptr(), self.yp.as_ptr()) },
            "IDA", "IDAReInit",
        )?;
        self.t = t0;
        Ok(())
    }

    pub fn t(&self) -> f64 { self.t }
    pub fn y(&self) -> &[f64] { self.y.as_slice() }
    pub fn yp(&self) -> &[f64] { self.yp.as_slice() }
}

impl<F> Drop for IdaSolver<F> {
    fn drop(&mut self) {
        unsafe { sys::IDAFree(&mut self.mem) };
    }
}

unsafe impl<F: Send> Send for IdaSolver<F> {}
