# Using CVODE, CVODES, IDA, and IDAS in Rust

A practical guide to integrating ODE and DAE systems and computing parameter
sensitivities using the `sundials-rs` crate — safe Rust bindings to the
[SUNDIALS](https://computing.llnl.gov/projects/sundials) library.

---

## Table of contents

1. [Prerequisites](#prerequisites)
2. [Core concepts](#core-concepts)
3. [Solving an ODE with CVODE](#solving-an-ode-with-cvode)
   - [Minimal example](#minimal-example)
   - [Choosing a method](#choosing-a-method)
   - [Setting tolerances](#setting-tolerances)
   - [Stepping through time](#stepping-through-time)
   - [Reinitialising the solver](#reinitialising-the-solver)
   - [Reading integration statistics](#reading-integration-statistics)
4. [Forward sensitivity analysis with CVODES](#forward-sensitivity-analysis-with-cvodes)
   - [What sensitivities mean](#what-sensitivities-mean)
   - [The DQ requirement](#the-dq-requirement)
   - [Full CVODES example](#full-cvodes-example)
   - [Reading sensitivities](#reading-sensitivities)
5. [Solving a DAE with IDA](#solving-a-dae-with-ida)
   - [What IDA solves](#what-ida-solves)
   - [Consistent initial conditions](#consistent-initial-conditions)
   - [Full IDA example](#full-ida-example)
6. [Forward sensitivity analysis with IDAS](#forward-sensitivity-analysis-with-idas)
   - [Full IDAS example](#full-idas-example)
7. [The Robertson problem — worked example](#the-robertson-problem--worked-example)
8. [Common pitfalls](#common-pitfalls)
9. [API reference summary](#api-reference-summary)

---

## Prerequisites

Install SUNDIALS (>= 6.0) on your system:

```bash
# Ubuntu / Debian
sudo apt install libsundials-dev

# Fedora / RHEL
sudo dnf install sundials-devel

# macOS (Homebrew)
brew install sundials
```

If SUNDIALS is installed in a non-standard location, set:

```bash
export SUNDIALS_DIR=/path/to/sundials/install
```

Add the crate to your `Cargo.toml`:

```toml
[dependencies]
sundials-rs = "0.1"

# Or build SUNDIALS from source (requires cmake on PATH — no system library needed):
# sundials-rs = { version = "0.1", features = ["vendored"] }
```

---

## Core concepts

### What CVODE solves

CVODE solves an **explicit ODE initial value problem** (IVP):

```text
dy/dt = f(t, y),    y(t₀) = y₀
```

where `y` is a vector of state variables and `f` is your right-hand side (RHS) function.

### What CVODES adds

CVODES is a superset of CVODE.  In addition to solving the ODE, it
simultaneously computes **forward sensitivities**:

```text
s_i(t) = ∂y(t) / ∂p_i
```

i.e. how each state variable changes with respect to each parameter `p_i`.
This is essential for parameter estimation, uncertainty quantification, and
gradient-based optimisation.

### Builder pattern

Both solvers use a builder to configure options before the first step:

```rust,ignore
let mut solver = CvodeBuilder::new(Method::BDF, &y0)
    .rtol(1e-6)
    .atol(1e-9)
    .build(|t, y, ydot| { /* fill ydot */ Ok(()) })
    .unwrap();
```

---

## Solving an ODE with CVODE

### Minimal example

```rust
use sundials_rs::cvode::{CvodeBuilder, Method};

fn main() {
    // ODE: dy/dt = -y,  y(0) = 1.0  (exact solution: y(t) = exp(-t))
    let y0 = vec![1.0_f64];

    let mut solver = CvodeBuilder::new(Method::BDF, &y0)
        .rtol(1e-8)
        .atol(1e-10)
        .build(|_t, y, ydot| {
            ydot[0] = -y[0];
            Ok(())
        })
        .unwrap();

    let (t, y) = solver.step(1.0).unwrap();
    println!("y({t}) = {:.8}  (exact: {:.8})", y[0], (-t).exp());
}
```

### Choosing a method

| Method | When to use |
|---|---|
| `Method::BDF` | Stiff ODEs (chemistry, pharmacokinetics, electrical circuits) |
| `Method::Adams` | Non-stiff ODEs (simple mechanics, low-frequency oscillators) |

**If in doubt, use BDF.** It is more robust for systems with widely separated time scales.

```rust,ignore
CvodeBuilder::new(Method::BDF, &y0)   // stiff — recommended default
CvodeBuilder::new(Method::Adams, &y0) // non-stiff
```

### Setting tolerances

Tolerances control accuracy vs speed.  Two types:

#### Scalar absolute tolerance (same for every component)

```rust,ignore
CvodeBuilder::new(Method::BDF, &y0)
    .rtol(1e-6)   // relative — scales with the solution magnitude
    .atol(1e-9)   // absolute — floor when the solution is near zero
```

#### Per-component absolute tolerance (recommended for mixed-scale systems)

Use this when different state variables have very different magnitudes.
For example, if `y[1]` drops to `1e-13` during the solve, a scalar `atol`
of `1e-9` would waste steps trying to resolve it unnecessarily:

```rust,ignore
CvodeBuilder::new(Method::BDF, &y0)
    .rtol(1e-4)
    .atol_vec(vec![1e-8, 1e-14, 1e-6]) // one entry per state variable
```

**Rule of thumb:** set `atol[i]` to the smallest value of `y[i]` that is
physically significant.  The Robertson chemistry example uses `[1e-8, 1e-14, 1e-6]`
because species 2 (`y[1]`) reaches concentrations around `3e-13`.

### Stepping through time

`solver.step(tout)` advances the solution to *at least* `tout` and returns
the actual time reached along with a slice of the state:

```rust,ignore
let (t, y) = solver.step(tout).unwrap();
//    ^         ^
//    actual    &[f64] view into solver's internal state
//    time
```

To march through a sequence of output times:

```rust,ignore
let output_times = [0.1, 1.0, 10.0, 100.0];

for &tout in &output_times {
    let (t, y) = solver.step(tout).unwrap();
    println!("t = {t:.2e}  y = {:?}", y);
}
```

> **Note:** `y` is a borrowed slice into the solver.  If you need to keep
> the values after the next call to `step`, copy them first:
> ```rust,ignore
> let y_saved: Vec<f64> = y.to_vec();
> ```

### Reinitialising the solver

To restart the integration from a new initial condition without allocating
a new solver (useful in parameter sweeps):

```rust,ignore
solver.reinit(t_new, &new_y0).unwrap();
```

### Reading integration statistics

After one or more steps, retrieve diagnostic information:

```rust,ignore
let stats = solver.stats().unwrap();
println!("Steps:              {}", stats.num_steps);
println!("RHS evaluations:    {}", stats.num_rhs_evals);
println!("Error test failures:{}", stats.num_err_test_fails);
println!("Last BDF order:     {}", stats.last_order);
```

A high `num_err_test_fails` (more than ~5–10% of `num_steps`) suggests
the tolerances are too tight for the step-size control, or the problem is
exceptionally stiff.

---

## Forward sensitivity analysis with CVODES

### What sensitivities mean

Given parameters `p = [p₀, p₁, ..., p_{ns-1}]`, forward sensitivity analysis
computes the `ns × neq` matrix:

```text
S[i][j] = ∂y_j(t) / ∂p_i
```

at every time point alongside the state `y(t)`.  This tells you: *if I
change parameter `p_i` by a small amount, how much does state `y_j` change?*

Common applications:
- **Parameter estimation** — gradient of the objective w.r.t. parameters
- **Uncertainty quantification** — propagate parameter uncertainty to states
- **Identifiability** — detect which parameters can be determined from data

### The DQ requirement

CVODES can approximate the sensitivity RHS automatically using **internal
difference quotients (DQ)** — no analytic derivatives needed.  However, DQ
works by **temporarily perturbing `p[i]`** and re-evaluating the RHS.

**This only works if your RHS closure reads the rate constants from the `p`
array rather than using hardcoded Rust constants.**

The `CvodesBuilder` provides a `params_ptr()` method for this purpose:

```rust,ignore
//                              ┌─ build up to here first
let builder = CvodesBuilder::new(&y0)
    .rtol(1e-6)
    .atol(1e-9)
    .with_forward_sensitivity(p_values, s0, SensMethod::Simultaneous);

// Get a raw pointer to the parameter array BEFORE build() consumes the builder.
// The heap allocation is moved (not re-allocated) into the solver, so this
// pointer remains valid for the entire lifetime of the solver.
let p_ptr = builder.params_ptr();

let mut solver = builder
    .build(move |_t, y, ydot| {
        // Read parameters via the raw pointer — CVODES may temporarily
        // change these during DQ perturbation.
        let p = unsafe { std::slice::from_raw_parts(p_ptr, num_params) };

        ydot[0] = -p[0] * y[0] + p[1] * y[1] * y[2];
        // ...
        Ok(())
    })
    .unwrap();
```

> **Why a raw pointer?**  
> `CvodeBuilder` is consumed by `build()`.  The parameter array lives inside
> the resulting `CvodesSolver` struct on the heap.  Moving a `Box<[f64]>` does
> not change the address of the underlying data, so a raw pointer captured
> before `build()` stays valid after.

### Full CVODES example

```rust
use sundials_rs::cvodes::{CvodesBuilder, SensMethod};

fn main() {
    // ODE: dy/dt = -k * y,  y(0) = 1
    // Parameters: p = [k]
    // Sensitivity: ds/dt = ∂(dy/dt)/∂k = -y   (exact: s(t) = -t * exp(-k*t))

    let y0 = vec![1.0_f64];
    let p  = vec![0.5_f64];           // k = 0.5
    let s0 = vec![vec![0.0_f64]];    // ds/dk = 0 at t = 0

    let builder = CvodesBuilder::new(&y0)
        .rtol(1e-8)
        .atol(1e-10)
        .with_forward_sensitivity(p, s0, SensMethod::Simultaneous);

    let p_ptr = builder.params_ptr();  // stable for the solver's lifetime

    let mut solver = builder
        .build(move |_t, y, ydot| {
            let k = unsafe { *p_ptr };  // read k from the parameter array
            ydot[0] = -k * y[0];
            Ok(())
        })
        .unwrap();

    let t_end = 2.0_f64;
    let k     = 0.5_f64;

    let (t, y_ref) = solver.step(t_end).unwrap();
    let y = y_ref.to_vec();   // copy before calling sensitivities()

    let sens = solver.sensitivities().unwrap();

    let y_exact = (-k * t).exp();
    let s_exact = -t * (-k * t).exp();

    println!("t = {t}");
    println!("y   = {:.8}  (exact {:.8})", y[0], y_exact);
    println!("s₀  = {:.8}  (exact {:.8})", sens[0][0], s_exact);
}
```

### Reading sensitivities

Call `solver.sensitivities()` immediately after `solver.step()`:

```rust,ignore
// step() borrows solver mutably and returns a &[f64] into it.
// Copy y first so the borrow is released before calling sensitivities().
let (t, y_ref) = solver.step(tout).unwrap();
let y: Vec<f64> = y_ref.to_vec();

let sens = solver.sensitivities().unwrap();
// sens[i] = ∂y/∂p_i  as &[f64] of length neq
```

`sens[i][j]` is `∂y_j / ∂p_i`.

> **Why copy `y` first?**  
> Both `step()` and `sensitivities()` take `&mut self`.  Rust's borrow
> checker prevents holding the `&[f64]` from `step()` while calling
> `sensitivities()`.  Copying `y` to a `Vec` releases the borrow.

---

## Solving a DAE with IDA

### What IDA solves

IDA solves **implicit** initial value problems:

```text
F(t, y, y') = 0,    y(t₀) = y₀,   y'(t₀) = y'₀
```

This is more general than CVODE's explicit form.  Use IDA when:

- Your equations are naturally implicit (e.g. circuit equations, chemical
  equilibria).
- Some components are **algebraic** — they have no `y'` term and instead
  enforce a constraint such as `y₁ + y₂ + y₃ = 1`.
- You are computing a constrained mechanical system (pendulum, contact
  dynamics).

For a pure explicit ODE `y' = f(t, y)` you *can* use IDA by writing
`F = y' − f(t, y)`, but CVODE is simpler and slightly more efficient.

### Consistent initial conditions

IDA requires `F(t₀, y₀, y'₀) = 0`.  Two strategies:

1. **Compute `y'₀` analytically** — for an implicit ODE, evaluate `y'` from
   the equations at `t₀`.  This is the simplest approach when there are no
   algebraic variables.
2. **Call `calc_ic`** — for DAEs with algebraic variables, provide best-guess
   initial values and let IDA correct them:

```rust,ignore
solver.calc_ic(tout1)?;  // tout1 is the first output time — sets direction only
```

### Full IDA example

Harmonic oscillator written as an implicit first-order system:

```rust
use sundials_rs::ida::IdaBuilder;

// d²x/dt² = -x  →  y[0]=x, y[1]=x'
// Residuals: F[0] = y'[0] - y[1] = 0
//            F[1] = y'[1] + y[0] = 0
let y0  = vec![1.0_f64, 0.0];   // x(0)=1, x'(0)=0
let yp0 = vec![0.0_f64, -1.0];  // y'[0]=y[1]=0, y'[1]=-y[0]=-1

let mut solver = IdaBuilder::new(&y0, &yp0)
    .rtol(1e-8)
    .atol(1e-10)
    .build(|_t, y, yp, res| {
        res[0] = yp[0] - y[1];
        res[1] = yp[1] + y[0];
        Ok(())
    })
    .unwrap();

let output_times = [std::f64::consts::PI / 2.0,
                    std::f64::consts::PI,
                    3.0 * std::f64::consts::PI / 2.0,
                    2.0 * std::f64::consts::PI];

for &tout in &output_times {
    let (t, y, _yp) = solver.step(tout).unwrap();
    println!("t={t:.4}  x={:.6}  (exact {:.6})", y[0], t.cos());
}
```

The residual closure signature is `|t, y, yp, res|` — note the extra `yp`
argument compared to CVODE's `|t, y, ydot|`.

---

## Forward sensitivity analysis with IDAS

IDAS is to IDA what CVODES is to CVODE — it augments the DAE solve with
forward sensitivities `∂y/∂p_i` using the same internal difference-quotient
(DQ) approximation.

The **DQ requirement** is identical: the residual closure must read parameters
through the raw pointer from `IdasSolver::params_ptr()`, not from hardcoded
Rust constants.

### Full IDAS example

```rust
use sundials_rs::idas::IdasBuilder;

// DAE residual:  F = y' + k*y = 0    (implicit decay)
// Parameter:     p = [k]
let y0  = vec![1.0_f64];
let p   = vec![0.5_f64];         // k = 0.5
let yp0 = vec![-p[0] * y0[0]];  // y'(0) = -k*y(0) — consistent
let s0  = vec![vec![0.0_f64]];   // ∂y/∂k = 0 at t=0

let mut builder = IdasBuilder::new(&y0, &yp0)
    .rtol(1e-8)
    .atol(1e-10)
    .with_forward_sensitivity(p, s0);

let p_ptr = builder.params_ptr(); // capture BEFORE build()

let mut solver = builder
    .build(move |_t, y, yp, res| {
        let k = unsafe { *p_ptr }; // read (possibly perturbed) parameter
        res[0] = yp[0] + k * y[0];
        Ok(())
    })
    .unwrap();

let t_end = 2.0_f64;
let (t, y_ref, _yp) = solver.step(t_end).unwrap();
let y = y_ref.to_vec(); // copy y before calling sensitivities()

let sens = solver.sensitivities().unwrap();

let k = 0.5_f64;
println!("t = {t}");
println!("y        = {:.8}  (exact {:.8})", y[0], (-k * t).exp());
println!("∂y/∂k    = {:.8}  (exact {:.8})", sens[0][0], -t * (-k * t).exp());
```

---

## The Robertson problem — worked example

The Robertson chemical kinetics system is the canonical benchmark for stiff
ODE solvers:

```text
dy₁/dt = −k₁·y₁ + k₂·y₂·y₃
dy₂/dt =  k₁·y₁ − k₂·y₂·y₃ − k₃·y₂²
dy₃/dt =                        k₃·y₂²

k₁ = 0.04,  k₂ = 1×10⁴,  k₃ = 3×10⁷
y(0) = [1, 0, 0]
```

The system is **extremely stiff**: k₃/k₁ = 7.5×10⁸.  Only BDF is practical.

### CVODE (integration only)

```rust
use sundials_rs::cvode::{CvodeBuilder, Method};

const K1: f64 = 0.04;
const K2: f64 = 1.0e4;
const K3: f64 = 3.0e7;

let y0 = [1.0_f64, 0.0, 0.0];

let mut solver = CvodeBuilder::new(Method::BDF, &y0)
    .rtol(1.0e-4)
    // Per-component atol: y[1] gets as small as ~3e-13
    .atol_vec(vec![1.0e-8, 1.0e-14, 1.0e-6])
    .max_steps(5000)
    .build(|_t, y, ydot| {
        ydot[0] = -K1 * y[0] + K2 * y[1] * y[2];
        ydot[1] =  K1 * y[0] - K2 * y[1] * y[2] - K3 * y[1] * y[1];
        ydot[2] =                                   K3 * y[1] * y[1];
        Ok(())
    })
    .unwrap();

let mut tout = 0.4_f64;
for _ in 0..12 {
    let (t, y) = solver.step(tout).unwrap();
    println!("t={:.1e}  y=[{:.4e}  {:.4e}  {:.4e}]",
             t, y[0], y[1] * 1e4, y[2]);
    tout *= 10.0;
}
```

### CVODES (integration + sensitivity w.r.t. k1, k2, k3)

```rust
use sundials_rs::cvodes::{CvodesBuilder, SensMethod};

let y0 = [1.0_f64, 0.0, 0.0];
let p  = vec![0.04_f64, 1.0e4, 3.0e7];   // [k1, k2, k3]
let s0 = vec![vec![0.0; 3]; 3];           // ∂y/∂pᵢ = 0 at t=0

let builder = CvodesBuilder::new(&y0)
    .rtol(1.0e-4)
    .atol_vec(vec![1.0e-8, 1.0e-14, 1.0e-6])
    .with_forward_sensitivity(p, s0, SensMethod::Simultaneous);

let p_ptr = builder.params_ptr();  // capture before build() consumes builder

let mut solver = builder
    .build(move |_t, y, ydot| {
        // Must use p_ptr — NOT hardcoded constants — so DQ perturbations work
        let p = unsafe { std::slice::from_raw_parts(p_ptr, 3) };
        ydot[0] = -p[0] * y[0] + p[1] * y[1] * y[2];
        ydot[1] =  p[0] * y[0] - p[1] * y[1] * y[2] - p[2] * y[1] * y[1];
        ydot[2] =                                        p[2] * y[1] * y[1];
        Ok(())
    })
    .unwrap();

let mut tout = 0.4_f64;
for _ in 0..12 {
    let (t, y_ref) = solver.step(tout).unwrap();
    let y = y_ref.to_vec();
    let sens = solver.sensitivities().unwrap();

    println!("t = {:.1e}  y = [{:.4e}  {:.4e}  {:.4e}]",
             t, y[0], y[1] * 1e4, y[2]);

    for (i, name) in ["k1", "k2", "k3"].iter().enumerate() {
        println!("  ∂y/∂{} = [{:+.3e}  {:+.3e}  {:+.3e}]",
                 name, sens[i][0], sens[i][1] * 1e4, sens[i][2]);
    }

    tout *= 10.0;
}
```

### IDA (residual form)

The same Robertson system written for IDA — just rearrange each equation
so the right-hand side is zero:

```rust
use sundials_rs::ida::IdaBuilder;

const K1: f64 = 0.04;
const K2: f64 = 1.0e4;
const K3: f64 = 3.0e7;

let y0  = [1.0_f64, 0.0, 0.0];
let yp0 = [-K1, K1, 0.0_f64]; // consistent derivatives at t=0

let mut solver = IdaBuilder::new(&y0, &yp0)
    .rtol(1.0e-4)
    .atol(1.0e-8)
    .build(|_t, y, yp, res| {
        res[0] = yp[0] + K1 * y[0] - K2 * y[1] * y[2];
        res[1] = yp[1] - K1 * y[0] + K2 * y[1] * y[2] + K3 * y[1] * y[1];
        res[2] = yp[2]                                   - K3 * y[1] * y[1];
        Ok(())
    })
    .unwrap();
```

The output is numerically identical to the CVODE version — IDA is solving
the same physics in a different representation.

### IDAS (residual form + sensitivity w.r.t. k1, k2, k3)

```rust
use sundials_rs::idas::IdasBuilder;

let y0  = [1.0_f64, 0.0, 0.0];
let p   = vec![0.04_f64, 1.0e4, 3.0e7]; // [k1, k2, k3]
let yp0 = [-p[0], p[0], 0.0_f64];       // consistent
let s0  = vec![vec![0.0; 3]; 3];

let mut builder = IdasBuilder::new(&y0, &yp0)
    .rtol(1.0e-4)
    .atol(1.0e-8)
    .with_forward_sensitivity(p, s0);

let p_ptr = builder.params_ptr(); // capture BEFORE build()

let mut solver = builder
    .build(move |_t, y, yp, res| {
        let (k1, k2, k3) = unsafe { (*p_ptr, *p_ptr.add(1), *p_ptr.add(2)) };
        res[0] = yp[0] + k1 * y[0] - k2 * y[1] * y[2];
        res[1] = yp[1] - k1 * y[0] + k2 * y[1] * y[2] + k3 * y[1] * y[1];
        res[2] = yp[2]                                   - k3 * y[1] * y[1];
        Ok(())
    })
    .unwrap();
```

> Full runnable versions are in `examples/robertson_ida.rs` and
> `examples/robertson_idas_fsa.rs`.

---

## Common pitfalls

### 1. Using hardcoded constants in the CVODES / IDAS residual

**Wrong** — sensitivities will all be zero:
```rust,ignore
let builder = CvodesBuilder::new(&y0)
    .with_forward_sensitivity(vec![0.04, 1e4], s0, SensMethod::Simultaneous);
let p_ptr = builder.params_ptr();

builder.build(|_t, y, ydot| {
    ydot[0] = -0.04 * y[0];  // ← hardcoded: DQ perturbs p[0] but this ignores it
    Ok(())
})
```

**Correct** — read rate constants through `p_ptr`:
```rust,ignore
builder.build(move |_t, y, ydot| {
    let k = unsafe { *p_ptr };  // ← reads the (possibly perturbed) value
    ydot[0] = -k * y[0];
    Ok(())
})
```

### 2. Borrow conflict between `step()` and `sensitivities()`

**Wrong** — borrow checker error:
```rust,ignore
let (t, y) = solver.step(tout).unwrap();  // borrows solver
let sens    = solver.sensitivities()?;    // second borrow — compile error
println!("{}", y[0]);                     // first borrow used here
```

**Correct** — copy `y` to release the borrow first:
```rust,ignore
let (t, y_ref) = solver.step(tout).unwrap();
let y = y_ref.to_vec();          // copy releases the borrow on solver
let sens = solver.sensitivities().unwrap();
println!("{}", y[0]);            // fine — using the owned Vec
```

### 3. Using scalar `atol` for mixed-magnitude systems

If any state variable drops many orders of magnitude below 1, a scalar
absolute tolerance will either waste steps (tight `atol`) or miss the
small-scale dynamics (loose `atol`).  Use `atol_vec` with per-component values.

### 4. Hitting the default step limit

CVODE defaults to 500 internal steps between output times.  Very stiff
problems or large output intervals will exceed this.  Increase it with
`.max_steps(n)`:

```rust,ignore
CvodeBuilder::new(Method::BDF, &y0)
    .max_steps(10_000)
    // ...
```

### 5. Forgetting `params_ptr()` must be called before `build()`

`build()` consumes the builder.  Get the parameter pointer first:

```rust,ignore
// ✓ correct order
let p_ptr = builder.params_ptr();
let solver = builder.build(move |_t, y, ydot| { /* uses p_ptr */ }).unwrap();

// ✗ wrong — builder is moved, p_ptr call is unreachable
let solver = builder.build(...).unwrap();
let p_ptr  = builder.params_ptr();  // compile error: builder moved
```

---

## API reference summary

### `CvodeBuilder` (from `sundials_rs::cvode`)

| Method | Description |
|---|---|
| `CvodeBuilder::new(method, y0)` | Create builder with integration method and initial conditions |
| `.t0(f64)` | Initial time (default `0.0`) |
| `.rtol(f64)` | Relative tolerance (default `1e-6`) |
| `.atol(f64)` | Scalar absolute tolerance (default `1e-9`) |
| `.atol_vec(Vec<f64>)` | Per-component absolute tolerance |
| `.max_steps(i64)` | Max internal steps per output interval (default `500`) |
| `.build(rhs)` | Consume builder, return `CvodeSolver<F>` |

### `CvodeSolver<F>` (from `sundials_rs::cvode`)

| Method | Description |
|---|---|
| `.step(tout) -> (f64, &[f64])` | Advance to `tout`, return `(t, y)` |
| `.reinit(t0, y0)` | Restart from new initial condition |
| `.stats()` | Return `CvodeStats` (step count, RHS evals, …) |
| `.t()` | Current time |
| `.y()` | Current state slice |

### `CvodesBuilder` (from `sundials_rs::cvodes`)

| Method | Description |
|---|---|
| `CvodesBuilder::new(y0)` | Create builder (BDF method; CVODES is always BDF) |
| `.t0(f64)` | Initial time |
| `.rtol(f64)` | Relative tolerance |
| `.atol(f64)` / `.atol_vec(Vec<f64>)` | Absolute tolerance |
| `.with_forward_sensitivity(p, s0, method)` | Enable FSA with parameters `p`, initial sensitivities `s0` |
| `.params_ptr()` | Raw pointer to the parameter array — capture before `build()` |
| `.build(rhs)` | Consume builder, return `CvodesSolver<F>` |

### `CvodesSolver<F>` (from `sundials_rs::cvodes`)

| Method | Description |
|---|---|
| `.step(tout) -> (f64, &[f64])` | Advance to `tout`, return `(t, y)` |
| `.sensitivities() -> Vec<&[f64]>` | Return `sens[i] = ∂y/∂pᵢ` after a step |
| `.t()` | Current time |
| `.y()` | Current state slice |

### Sensitivity methods (`SensMethod`)

| Variant | When to use |
|---|---|
| `SensMethod::Simultaneous` | Default — corrects state and sensitivities together |
| `SensMethod::Staggered` | Useful when the sensitivity RHS dominates cost |

---

### `IdaBuilder` (from `sundials_rs::ida`)

| Method | Description |
|---|---|
| `IdaBuilder::new(y0, yp0)` | Create builder with initial state and its derivative |
| `.t0(f64)` | Initial time (default `0.0`) |
| `.rtol(f64)` | Relative tolerance (default `1e-6`) |
| `.atol(f64)` | Scalar absolute tolerance (default `1e-9`) |
| `.build(res)` | Consume builder, return `IdaSolver<F>` |

The residual closure signature is `|t, y, yp, res| -> Result<(), i32>`.

### `IdaSolver<F>` (from `sundials_rs::ida`)

| Method | Description |
|---|---|
| `.calc_ic(tout1)` | Correct initial conditions so `F(t₀, y₀, yp₀) = 0` |
| `.step(tout) -> (f64, &[f64], &[f64])` | Advance to `tout`, return `(t, y, yp)` |
| `.reinit(t0, y0, yp0)` | Restart from new initial condition |
| `.t()` | Current time |
| `.y()` | Current state slice |
| `.yp()` | Current derivative slice |

---

### `IdasBuilder` (from `sundials_rs::idas`)

| Method | Description |
|---|---|
| `IdasBuilder::new(y0, yp0)` | Create builder with initial state and derivative |
| `.t0(f64)` | Initial time |
| `.rtol(f64)` | Relative tolerance |
| `.atol(f64)` | Scalar absolute tolerance |
| `.with_forward_sensitivity(p, s0)` | Enable FSA with parameters `p`, initial sensitivities `s0` |
| `.params_ptr()` | Raw pointer to the parameter array — capture before `build()` |
| `.build(res)` | Consume builder, return `IdasSolver<F>` |

### `IdasSolver<F>` (from `sundials_rs::idas`)

| Method | Description |
|---|---|
| `.calc_ic(tout1)` | Correct initial conditions to satisfy the DAE residual |
| `.step(tout) -> (f64, &[f64], &[f64])` | Advance to `tout`, return `(t, y, yp)` |
| `.sensitivities() -> Vec<&[f64]>` | Return `sens[i] = ∂y/∂pᵢ` after a step |
| `.t()` | Current time |
| `.y()` | Current state slice |
| `.yp()` | Current derivative slice |
