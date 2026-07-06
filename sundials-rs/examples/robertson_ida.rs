//! Robertson chemical kinetics — IDA implicit solver.
//!
//! Background
//! ----------
//! The Robertson problem is normally written as an explicit ODE (see
//! `robertson_cvode.rs`).  Here we solve the identical system using IDA
//! by expressing every equation in *residual* form:
//!
//!   F₁ = y₁' + k₁·y₁ − k₂·y₂·y₃ = 0
//!   F₂ = y₂' − k₁·y₁ + k₂·y₂·y₃ + k₃·y₂² = 0
//!   F₃ = y₃' − k₃·y₂² = 0
//!
//! IDA treats every equation as F(t, y, y') = 0, which makes it the natural
//! choice for true differential-algebraic systems (e.g. Robertson with the
//! conservation constraint y₁ + y₂ + y₃ = 1 made algebraic).  Here all three
//! equations are differential, so IDA and CVODE give identical results.
//!
//! Key differences vs CVODE
//! ------------------------
//! - The residual closure takes FOUR arguments: (t, y, y', res).
//! - Both y₀ and y'₀ must be provided and must be *consistent*:
//!   F(t₀, y₀, y'₀) = 0.  For this problem y'₀ is computed directly.
//! - Use IDASolve / IDA_NORMAL instead of CVode.
//!
//! Run with:
//!   cargo run --example robertson_ida

use sundials_rs::ida::IdaBuilder;

const K1: f64 = 0.04;
const K2: f64 = 1.0e4;
const K3: f64 = 3.0e7;

fn main() {
    let y0 = [1.0_f64, 0.0, 0.0];

    // Consistent initial derivatives:  yp[i] = dy[i]/dt at t=0
    //   y1'(0) = -k1*y1(0) + k2*y2(0)*y3(0) = -0.04
    //   y2'(0) =  k1*y1(0) - k2*y2(0)*y3(0) - k3*y2(0)² = 0.04
    //   y3'(0) =  k3*y2(0)² = 0
    let yp0 = [-K1, K1, 0.0_f64];

    let mut solver = IdaBuilder::new(&y0, &yp0)
        .rtol(1.0e-4)
        .atol(1.0e-8)   // scalar; Robertson ideally uses per-component atol,
                        // but scalar works here for a tutorial example
        .max_steps(20_000)
        .build(|_t, y, yp, res| {
            res[0] = yp[0] + K1 * y[0] - K2 * y[1] * y[2];
            res[1] = yp[1] - K1 * y[0] + K2 * y[1] * y[2] + K3 * y[1] * y[1];
            res[2] = yp[2]                                   - K3 * y[1] * y[1];
            Ok(())
        })
        .expect("solver creation failed");

    println!("\nRobertson Chemical Kinetics — IDA BDF (residual form)");
    println!("{:-<65}", "");
    println!("{:>12}  {:>12}  {:>12}  {:>12}", "t", "y1", "y2×1e4", "y3");
    println!("{:-<65}", "");

    let mut tout = 0.4_f64;
    for _ in 0..12 {
        let (t, y, _yp) = solver.step(tout).expect("step failed");
        println!(
            "{:12.4e}  {:12.4e}  {:12.4e}  {:12.4e}",
            t, y[0], y[1] * 1.0e4, y[2]
        );
        tout *= 10.0;
    }

    println!("{:-<65}", "");
    println!("(same result as robertson_cvode — IDA solving in residual form)");
}
