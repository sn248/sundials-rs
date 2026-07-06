//! Robertson chemical kinetics — CVODES forward sensitivity analysis (FSA).
//!
//! Mirrors the SUNDIALS example `cvsRoberts_FSA_dns.c`.
//!
//! Problem
//! -------
//! Same 3-species Robertson ODE as `robertson_cvode.rs`.  In addition, we
//! compute the sensitivities of each state variable with respect to the three
//! rate constants:
//!
//!   s_i(t) = ∂y(t)/∂p_i   for  p = [k1, k2, k3]
//!
//! The sensitivity system is integrated simultaneously using CVODES internal
//! difference-quotient (DQ) approximation.
//!
//! KEY REQUIREMENT for DQ:
//! The RHS closure MUST read rate constants from the `p` array (via a raw
//! pointer captured from `builder.params_ptr()`), NOT from Rust constants.
//! CVODES temporarily perturbs `p[i]` when computing the DQ; if the RHS
//! ignores `p`, all sensitivities will be zero.
//!
//! Run with:
//!   cargo run --example robertson_cvodes_fsa

use sundials_rs::cvodes::{CvodesBuilder, SensMethod};

const NEQ: usize = 3;
const NS:  usize = 3; // number of sensitivity parameters

fn main() {
    let y0 = [1.0_f64, 0.0, 0.0];

    // Initial parameter values: k1, k2, k3
    let p  = vec![0.04_f64, 1.0e4, 3.0e7];

    // Initial sensitivities: ∂y/∂p_i = 0 at t = 0
    let s0: Vec<Vec<f64>> = vec![vec![0.0; NEQ]; NS];

    // ── Build the builder (consuming chain stops here so we can call params_ptr)
    let builder = CvodesBuilder::new(&y0)
        .rtol(1.0e-4)
        .atol_vec(vec![1.0e-8, 1.0e-14, 1.0e-6])
        .with_forward_sensitivity(p, s0, SensMethod::Simultaneous);

    // Capture a raw pointer to the parameter array BEFORE build() consumes
    // the builder.  The heap allocation is moved (not copied) into the solver,
    // so the pointer remains valid for the solver's entire lifetime.
    let p_ptr = builder.params_ptr();

    let mut solver = builder
        .build(move |_t, y, ydot| {
            // SAFETY: p_ptr points into solver.p which lives as long as the solver.
            // CVODES may temporarily change p[i] here for DQ — reading via the
            // pointer ensures those perturbations are visible.
            let (k1, k2, k3) = unsafe { (*p_ptr, *p_ptr.add(1), *p_ptr.add(2)) };

            ydot[0] = -k1 * y[0] + k2 * y[1] * y[2];
            ydot[1] =  k1 * y[0] - k2 * y[1] * y[2] - k3 * y[1] * y[1];
            ydot[2] =                                   k3 * y[1] * y[1];
            Ok(())
        })
        .expect("failed to create CVODES solver");

    // ── Print header ──────────────────────────────────────────────────────────
    println!("\nRobertson Chemical Kinetics — CVODES BDF + Forward Sensitivity");
    println!("{:-<72}", "");
    println!("{:>12}  {:>12}  {:>12}  {:>12}", "t", "y1", "y2×1e4", "y3");
    println!("{:-<72}", "");

    // ── Integrate and print at each decade ────────────────────────────────────
    let mut tout = 0.4_f64;
    for _ in 0..12 {
        let (t, y_ref) = solver.step(tout).expect("step failed");
        let y: Vec<f64> = y_ref.to_vec();
        let sens = solver.sensitivities().expect("sens failed");

        println!(
            "{:12.4e}  {:12.4e}  {:12.4e}  {:12.4e}",
            t, y[0], y[1] * 1.0e4, y[2]
        );
        for (i, name) in ["k1", "k2", "k3"].iter().enumerate() {
            println!(
                "  ∂y/∂{} = [{:+.3e}  {:+.3e}  {:+.3e}]",
                name,
                sens[i][0],
                sens[i][1] * 1.0e4,
                sens[i][2],
            );
        }
        println!();

        tout *= 10.0;
    }
}
