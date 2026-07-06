//! Safe wrapper around CVODES — CVODE with forward/adjoint sensitivity analysis.

use std::os::raw::{c_int, c_void};
use sundials_rs_sys as sys;

use crate::{
    context::SunContext,
    error::{check_flag, SundialsError},
    linear_solver::LinearSolver,
    matrix::DenseMatrix,
    nvector::NVector,
};

#[derive(Debug, Clone, Copy)]
pub enum SensMethod {
    Simultaneous,
    Staggered,
}

impl SensMethod {
    fn as_c_int(self) -> c_int {
        match self {
            SensMethod::Simultaneous => sys::CV_SIMULTANEOUS as c_int,
            SensMethod::Staggered    => sys::CV_STAGGERED    as c_int,
        }
    }
}

struct UserData<F> {
    rhs: F,
    neq: usize,
}

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
        Err(f) => f,
    }
}

enum AtolMode {
    Scalar(f64),
    Vector(Vec<f64>),
}

pub struct CvodesBuilder {
    y0: Vec<f64>,
    t0: f64,
    rtol: f64,
    atol: AtolMode,
    ns: Option<usize>,
    sens_method: SensMethod,
    // Box<[f64]> gives a stable heap address that won't change when the
    // builder is consumed by build() — safe to hand out a raw pointer to it.
    p: Box<[f64]>,
    s0: Vec<Vec<f64>>,
}

impl CvodesBuilder {
    pub fn new(y0: &[f64]) -> Self {
        Self {
            y0: y0.to_vec(),
            t0: 0.0,
            rtol: 1e-6,
            atol: AtolMode::Scalar(1e-9),
            ns: None,
            sens_method: SensMethod::Simultaneous,
            p: Box::new([]),
            s0: Vec::new(),
        }
    }

    pub fn rtol(mut self, rtol: f64) -> Self { self.rtol = rtol; self }
    pub fn atol(mut self, atol: f64) -> Self { self.atol = AtolMode::Scalar(atol); self }
    pub fn atol_vec(mut self, atol: Vec<f64>) -> Self { self.atol = AtolMode::Vector(atol); self }
    pub fn t0(mut self, t0: f64) -> Self { self.t0 = t0; self }

    pub fn with_forward_sensitivity(
        mut self, p: Vec<f64>, s0: Vec<Vec<f64>>, method: SensMethod,
    ) -> Self {
        self.ns = Some(p.len());
        self.p = p.into_boxed_slice(); // heap-stable
        self.s0 = s0;
        self.sens_method = method;
        self
    }

    /// Raw pointer to the parameter array.
    ///
    /// Call this **after** [`CvodesBuilder::with_forward_sensitivity`] and capture the pointer
    /// in the RHS closure so that CVODES's internal difference-quotient
    /// perturbations of `p` are visible to the integrator:
    ///
    /// ```no_run
    /// # use sundials_rs::cvodes::{CvodesBuilder, SensMethod};
    /// let mut b = CvodesBuilder::new(&[1.0, 0.0, 0.0])
    ///     .with_forward_sensitivity(vec![0.04], vec![vec![0.0; 3]], SensMethod::Simultaneous);
    /// let p_ptr = b.params_ptr();
    /// let solver = b.build(move |_t, y, ydot| {
    ///     let k1 = unsafe { *p_ptr };
    ///     ydot[0] = -k1 * y[0];
    ///     Ok(())
    /// });
    /// ```
    ///
    /// # Safety
    /// The pointer is valid for the entire lifetime of the resulting
    /// `CvodesSolver` — `build()` moves the allocation into the solver without
    /// changing its heap address.
    pub fn params_ptr(&self) -> *const f64 {
        self.p.as_ptr()
    }

    pub fn build<F>(self, rhs: F) -> Result<CvodesSolver<F>, SundialsError>
    where
        F: Fn(f64, &[f64], &mut [f64]) -> Result<(), i32>,
    {
        let neq = self.y0.len();

        let ctx = SunContext::new()?;

        let mem = unsafe { sys::CVodeCreate(sys::CV_BDF as c_int, ctx.raw()) };
        if mem.is_null() {
            return Err(SundialsError::Memory("CVodeCreate (CVODES)"));
        }

        let y = NVector::from_slice(&self.y0, ctx.raw())?;
        let user_data = Box::new(UserData { rhs, neq });
        let ud_ptr = &*user_data as *const UserData<F> as *mut c_void;

        check_flag(
            unsafe { sys::CVodeInit(mem, Some(rhs_trampoline::<F>), self.t0, y.as_ptr()) },
            "CVODES", "CVodeInit",
        )?;
        match self.atol {
            AtolMode::Scalar(atol) => {
                check_flag(
                    unsafe { sys::CVodeSStolerances(mem, self.rtol, atol) },
                    "CVODES", "CVodeSStolerances",
                )?;
            }
            AtolMode::Vector(ref v) => {
                let atol_vec = NVector::from_slice(v, ctx.raw())?;
                check_flag(
                    unsafe { sys::CVodeSVtolerances(mem, self.rtol, atol_vec.as_ptr()) },
                    "CVODES", "CVodeSVtolerances",
                )?;
            }
        }
        check_flag(
            unsafe { sys::CVodeSetUserData(mem, ud_ptr) },
            "CVODES", "CVodeSetUserData",
        )?;

        let matrix = DenseMatrix::new(neq, neq, ctx.raw())?;
        let ls = LinearSolver::dense(&y, &matrix, ctx.raw())?;
        check_flag(
            unsafe { sys::CVodeSetLinearSolver(mem, ls.ptr.as_ptr(), matrix.ptr.as_ptr()) },
            "CVODES", "CVodeSetLinearSolver",
        )?;

        let mut sens_vecs: Vec<NVector> = Vec::new();

        if let Some(ns) = self.ns {
            for s in &self.s0 {
                sens_vecs.push(NVector::from_slice(s, ctx.raw())?);
            }
            let mut sv_ptrs: Vec<sys::N_Vector> =
                sens_vecs.iter().map(|v| v.as_ptr()).collect();

            check_flag(
                unsafe {
                    sys::CVodeSensInit(
                        mem, ns as c_int, self.sens_method.as_c_int(),
                        None, sv_ptrs.as_mut_ptr(),
                    )
                },
                "CVODES", "CVodeSensInit",
            )?;
            // NOTE: CVodeSetSensParams stores the raw pointer — p must live as
            // long as the solver.  We move self.p into the solver struct below.
            check_flag(
                unsafe { sys::CVodeSensEEtolerances(mem) },
                "CVODES", "CVodeSensEEtolerances",
            )?;
        }

        // Move the heap-stable Box into the solver — address doesn't change.
        let p = self.p;

        let solver = CvodesSolver {
            mem,
            y,
            _matrix: matrix,
            _ls: ls,
            _user_data: user_data,
            sens_vecs,
            p,
            _ctx: ctx,
            t: self.t0,
        };

        // Set sens params now that solver.p has its final address in the heap.
        if solver.ns() > 0 {
            check_flag(
                unsafe {
                    sys::CVodeSetSensParams(
                        solver.mem,
                        solver.p.as_ptr() as *mut f64,
                        std::ptr::null_mut(),
                        std::ptr::null_mut(),
                    )
                },
                "CVODES", "CVodeSetSensParams",
            )?;
        }

        Ok(solver)
    }
}

pub struct CvodesSolver<F> {
    mem: *mut c_void,
    y: NVector,
    _matrix: DenseMatrix,
    _ls: LinearSolver,
    _user_data: Box<UserData<F>>,
    sens_vecs: Vec<NVector>,
    /// Parameter array — heap-stable Box so CVodeSetSensParams' raw ptr stays valid.
    p: Box<[f64]>,
    _ctx: SunContext,
    t: f64,
}

impl<F> CvodesSolver<F> {
    fn ns(&self) -> usize { self.sens_vecs.len() }
}

impl<F> CvodesSolver<F>
where
    F: Fn(f64, &[f64], &mut [f64]) -> Result<(), i32>,
{
    pub fn new(y0: &[f64]) -> CvodesBuilder {
        CvodesBuilder::new(y0)
    }

    pub fn step(&mut self, tout: f64) -> Result<(f64, &[f64]), SundialsError> {
        let mut t_out: sys::sunrealtype = 0.0;
        let flag = unsafe {
            sys::CVode(self.mem, tout, self.y.as_ptr(), &mut t_out, sys::CV_NORMAL as c_int)
        };
        check_flag(flag, "CVODES", "CVode")?;
        self.t = t_out;
        Ok((t_out, self.y.as_slice()))
    }

    pub fn sensitivities(&mut self) -> Result<Vec<&[f64]>, SundialsError> {
        if self.sens_vecs.is_empty() {
            return Ok(Vec::new());
        }
        let mut sv_ptrs: Vec<sys::N_Vector> =
            self.sens_vecs.iter().map(|v| v.as_ptr()).collect();
        check_flag(
            unsafe { sys::CVodeGetSens(self.mem, &mut self.t, sv_ptrs.as_mut_ptr()) },
            "CVODES", "CVodeGetSens",
        )?;
        Ok(self.sens_vecs.iter().map(|v| v.as_slice()).collect())
    }

    pub fn t(&self) -> f64 { self.t }
    pub fn y(&self) -> &[f64] { self.y.as_slice() }
}

impl<F> Drop for CvodesSolver<F> {
    fn drop(&mut self) {
        unsafe { sys::CVodeFree(&mut self.mem) };
    }
}

unsafe impl<F: Send> Send for CvodesSolver<F> {}
