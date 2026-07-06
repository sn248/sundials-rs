#![doc = include_str!("../GUIDE.md")]
//! Safe, idiomatic Rust bindings to the [SUNDIALS](https://computing.llnl.gov/projects/sundials)
//! solver library.
//!
//! # Solvers
//!
//! | Module | Underlying library | What it solves |
//! |---|---|---|
//! | [`cvode`] | CVODE | Explicit ODE IVP: `y' = f(t, y)` |
//! | [`cvodes`] | CVODES | CVODE + forward/adjoint sensitivity |
//! | [`ida`] | IDA | Implicit DAE IVP: `F(t, y, y') = 0` |
//! | [`idas`] | IDAS | IDA + forward sensitivity |
//!
//! All four solvers use SUNDIALS's variable-order BDF (or Adams for CVODE)
//! method with an adaptive step-size controller — well-suited for stiff and
//! mildly stiff problems.
//!
//! # Choosing a solver
//!
//! - **Explicit ODE** `y' = f(t, y)` → [`cvode`] (no sensitivities) or [`cvodes`]
//! - **DAE / implicit ODE** `F(t, y, y') = 0` → [`ida`] (no sensitivities) or [`idas`]
//! - **Need `∂y/∂p`** (parameter gradients, UQ) → [`cvodes`] or [`idas`]
//!
//! # Quick start — CVODE BDF
//! ```no_run
//! use sundials_rs::cvode::{CvodeBuilder, Method};
//!
//! // dy/dt = -y,  y(0) = 1.0  (exact: exp(-t))
//! let y0 = vec![1.0_f64];
//! let mut solver = CvodeBuilder::new(Method::BDF, &y0)
//!     .rtol(1e-8)
//!     .atol(1e-10)
//!     .build(|_t, y, ydot| { ydot[0] = -y[0]; Ok(()) })
//!     .unwrap();
//!
//! let (t, y) = solver.step(1.0).unwrap();
//! println!("y({t:.3}) = {:.8}  (exact {:.8})", y[0], (-t).exp());
//! ```
//!
//! # Quick start — IDA (implicit ODE / DAE)
//! ```no_run
//! use sundials_rs::ida::IdaBuilder;
//!
//! // Same ODE in residual form:  F = y' + y = 0
//! let y0  = vec![1.0_f64];
//! let yp0 = vec![-1.0_f64]; // consistent: y'(0) = -y(0)
//! let mut solver = IdaBuilder::new(&y0, &yp0)
//!     .rtol(1e-8)
//!     .atol(1e-10)
//!     .build(|_t, y, yp, res| { res[0] = yp[0] + y[0]; Ok(()) })
//!     .unwrap();
//!
//! let (t, y, _yp) = solver.step(1.0).unwrap();
//! println!("y({t:.3}) = {:.8}  (exact {:.8})", y[0], (-t).exp());
//! ```

pub mod context;
pub mod error;
pub mod nvector;
pub mod matrix;
pub mod linear_solver;
pub mod cvode;
pub mod cvodes;
pub mod ida;
pub mod idas;

pub use error::SundialsError;
