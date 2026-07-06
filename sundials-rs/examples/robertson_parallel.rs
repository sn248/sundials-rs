//! Robertson chemical kinetics — parallel parameter sweep using Rayon.
//!
//! Background
//! ----------
//! The SUNDIALS C library ships MPI-based parallel examples that distribute a
//! *single* large ODE across many nodes using the parallel N_Vector.  That
//! approach requires linking against MPI and is not yet covered by this crate.
//!
//! For the Robertson problem — which has only 3 equations — the practical
//! "parallel" use case is different: running many independent solves
//! simultaneously while varying a parameter.  This is extremely common in
//! pharmacokinetics, parameter estimation, and uncertainty quantification.
//!
//! This example sweeps k1 ∈ [0.02, 0.04, 0.06, 0.08, 0.10] and solves each
//! Robertson problem to t = 4×10¹¹ in parallel across CPU threads using Rayon.
//!
//! Key insight
//! -----------
//! Each `CvodeSolver` owns its SUNDIALS memory independently, so it is safe to
//! send across threads (`Send`).  Multiple threads never share the same solver.
//!
//! Run with:
//!   cargo run --example robertson_parallel

use rayon::prelude::*;
use sundials_rs::cvode::{CvodeBuilder, Method};

const K2: f64 = 1.0e4;
const K3: f64 = 3.0e7;

/// Solve Robertson to t_end and return the final state.
fn solve_robertson(k1: f64, t_end: f64) -> (f64, [f64; 3]) {
    let y0 = [1.0_f64, 0.0, 0.0];

    let mut solver = CvodeBuilder::new(Method::BDF, &y0)
        .rtol(1.0e-4)
        .atol_vec(vec![1.0e-8, 1.0e-14, 1.0e-6])
        .max_steps(5000)
        .build(move |_t, y, ydot| {
            ydot[0] = -k1 * y[0] + K2 * y[1] * y[2];
            ydot[1] =  k1 * y[0] - K2 * y[1] * y[2] - K3 * y[1] * y[1];
            ydot[2] =                                   K3 * y[1] * y[1];
            Ok(())
        })
        .expect("solver creation failed");

    let (t, y_slice) = solver.step(t_end).expect("step failed");
    let y_arr = [y_slice[0], y_slice[1], y_slice[2]];
    (t, y_arr)
}

fn main() {
    // Parameter sweep: k1 varies, k2/k3 fixed.
    let k1_values = [0.02_f64, 0.04, 0.06, 0.08, 0.10];
    let t_end = 4.0e11_f64;

    println!("\nRobertson — Parallel parameter sweep over k1");
    println!("{:-<65}", "");
    println!("{:>8}  {:>14}  {:>14}  {:>14}", "k1", "y1(T)", "y2(T)×1e4", "y3(T)");
    println!("{:-<65}", "");

    // Collect results in parallel, then sort by k1 for deterministic output.
    let mut results: Vec<(f64, f64, [f64; 3])> = k1_values
        .par_iter()                      // Rayon parallel iterator
        .map(|&k1| {
            let (t, y) = solve_robertson(k1, t_end);
            (k1, t, y)
        })
        .collect();

    results.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    for (k1, _t, y) in &results {
        println!(
            "{:8.4}  {:14.6e}  {:14.6e}  {:14.6e}",
            k1, y[0], y[1] * 1.0e4, y[2]
        );
    }

    println!("{:-<65}", "");
    println!("(Each row solved independently on a separate thread)");
}
