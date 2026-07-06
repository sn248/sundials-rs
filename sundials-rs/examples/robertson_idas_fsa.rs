//! Robertson chemical kinetics — IDAS forward sensitivity analysis (FSA).
//!
//! Mirrors `robertson_cvodes_fsa.rs` but uses the IDA-family solver, which
//! works with residual equations F(t, y, y') = 0 rather than explicit
//! right-hand sides y' = f(t, y).
//!
//! Problem
//! -------
//! Same Robertson residual system as `robertson_ida.rs`.  In addition, we
//! compute the sensitivities of each state variable with respect to the three
//! rate constants:
//!
//!   s_i(t) = ∂y(t)/∂p_i   for  p = [k1, k2, k3]
//!
//! KEY REQUIREMENT for DQ:
//! The residual closure MUST read rate constants from the `p` array via the
//! raw pointer from `builder.params_ptr()`.  IDAS temporarily perturbs `p[i]`
//! when computing the difference-quotient sensitivity; if the closure uses
//! hardcoded constants, all sensitivities will be zero.
//!
//! Run with:
//!   cargo run --example robertson_idas_fsa

use sundials_rs::idas::IdasBuilder;

const NEQ: usize = 3;
const NS:  usize = 3; // number of sensitivity parameters (k1, k2, k3)

fn main() {
    let y0  = [1.0_f64, 0.0, 0.0];
    let p   = vec![0.04_f64, 1.0e4, 3.0e7];   // [k1, k2, k3]

    // Consistent initial derivatives computed from p
    let yp0 = [-p[0], p[0], 0.0_f64];

    // Initial sensitivities: ∂y/∂p_i = 0 at t = 0
    let s0: Vec<Vec<f64>> = vec![vec![0.0; NEQ]; NS];

    // ── Build the builder (chain stops here so we can call params_ptr) ────────
    let builder = IdasBuilder::new(&y0, &yp0)
        .rtol(1.0e-4)
        .atol(1.0e-8)
        .max_steps(100_000)
        .with_forward_sensitivity(p, s0);

    // Capture a raw pointer to the parameter array BEFORE build() consumes the
    // builder.  The heap allocation is moved (not re-allocated) into the solver,
    // so the pointer remains valid for the solver's entire lifetime.
    let p_ptr = builder.params_ptr();

    let mut solver = builder
        .build(move |_t, y, yp, res| {
            // SAFETY: p_ptr points into solver.p which lives as long as the solver.
            // IDAS may temporarily change p[i] here for DQ — reading via the
            // pointer ensures those perturbations are visible.
            let (k1, k2, k3) = unsafe { (*p_ptr, *p_ptr.add(1), *p_ptr.add(2)) };

            res[0] = yp[0] + k1 * y[0] - k2 * y[1] * y[2];
            res[1] = yp[1] - k1 * y[0] + k2 * y[1] * y[2] + k3 * y[1] * y[1];
            res[2] = yp[2]                                   - k3 * y[1] * y[1];
            Ok(())
        })
        .expect("failed to create IDAS solver");

    // ── Print header ──────────────────────────────────────────────────────────
    println!("\nRobertson Chemical Kinetics — IDAS BDF + Forward Sensitivity");
    println!("{:-<72}", "");
    println!("{:>12}  {:>12}  {:>12}  {:>12}", "t", "y1", "y2×1e4", "y3");
    println!("{:-<72}", "");

    // ── Integrate and print at each decade ────────────────────────────────────
    let mut tout = 0.4_f64;
    for _ in 0..12 {
        let (t, y_ref, _yp) = solver.step(tout).expect("step failed");
        let y: Vec<f64> = y_ref.to_vec(); // copy before calling sensitivities()

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
