//! Robertson chemical kinetics — serial CVODE with dense Jacobian.
//!
//! Mirrors the SUNDIALS example `cvRoberts_dns.c`.
//!
//! Problem
//! -------
//! A stiff 3-species autocatalytic reaction:
//!
//!   dy1/dt = -k1·y1 + k2·y2·y3
//!   dy2/dt =  k1·y1 - k2·y2·y3 - k3·y2²
//!   dy3/dt =  k3·y2²
//!
//! with rate constants  k1 = 0.04,  k2 = 1e4,  k3 = 3e7
//! and initial conditions  y(0) = [1, 0, 0].
//!
//! The problem is solved on [0, 4×10¹¹] and the solution is printed at
//! each decade of time just as in the original C example.
//!
//! Run with:
//!   cargo run --example robertson_cvode

use sundials_rs::cvode::{CvodeBuilder, Method};

// Rate constants
const K1: f64 = 0.04;
const K2: f64 = 1.0e4;
const K3: f64 = 3.0e7;

fn main() {
    // ── Initial conditions ────────────────────────────────────────────────────
    let y0 = [1.0_f64, 0.0, 0.0];
    let t0 = 0.0_f64;

    // ── Build solver ──────────────────────────────────────────────────────────
    // BDF method — the standard choice for stiff chemistry.
    // Per-component absolute tolerances: y2 can reach ~3e-13 so we set a very
    // tight atol on that component (matches the original C example).
    let mut solver = CvodeBuilder::new(Method::BDF, &y0)
        .t0(t0)
        .rtol(1.0e-4)
        .atol_vec(vec![1.0e-8, 1.0e-14, 1.0e-6])
        .build(|_t, y, ydot| {
            ydot[0] = -K1 * y[0] + K2 * y[1] * y[2];
            ydot[1] =  K1 * y[0] - K2 * y[1] * y[2] - K3 * y[1] * y[1];
            ydot[2] =                                   K3 * y[1] * y[1];
            Ok(())
        })
        .expect("failed to create CVODE solver");

    // ── Print header ──────────────────────────────────────────────────────────
    println!("\nRobertson Chemical Kinetics — CVODE BDF / Dense");
    println!("{:-<60}", "");
    println!("{:>12}  {:>14}  {:>14}  {:>14}", "t", "y1", "y2", "y3");
    println!("{:-<60}", "");

    // ── Integrate and print at each decade ────────────────────────────────────
    let mut tout = 0.4_f64;
    for _ in 0..12 {
        let (t, y) = solver.step(tout).expect("step failed");

        // y2 is tiny — scale it by 1e4 for display (same as SUNDIALS C example)
        println!(
            "{:12.4e}  {:14.6e}  {:14.6e}  {:14.6e}",
            t, y[0], y[1] * 1.0e4, y[2]
        );

        tout *= 10.0;
    }

    // ── Final statistics ──────────────────────────────────────────────────────
    let stats = solver.stats().expect("failed to get stats");
    println!("{:-<60}", "");
    println!("Steps taken         : {}", stats.num_steps);
    println!("RHS evaluations     : {}", stats.num_rhs_evals);
    println!("Error test failures : {}", stats.num_err_test_fails);
    println!("Last order used     : {}", stats.last_order);
}
