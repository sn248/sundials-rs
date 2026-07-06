# sundials-rs

Safe, idiomatic Rust bindings to the [SUNDIALS](https://computing.llnl.gov/projects/sundials)
ODE/DAE solver library.

## Solvers

| Crate module | C library | Solves |
|---|---|---|
| `cvode` | CVODE | Explicit ODE IVP: `y' = f(t, y)` |
| `cvodes` | CVODES | CVODE + forward sensitivity `∂y/∂p` |
| `ida` | IDA | Implicit DAE IVP: `F(t, y, y') = 0` |
| `idas` | IDAS | IDA + forward sensitivity |

All four use SUNDIALS's variable-order BDF (stiff) or Adams (non-stiff) method
with adaptive step-size control.

## Quick start

```toml
[dependencies]
sundials-rs = "0.1"

# Build SUNDIALS from source (requires cmake on PATH — no system library needed):
# sundials-rs = { version = "0.1", features = ["vendored"] }
```

```rust
use sundials_rs::cvode::{CvodeBuilder, Method};

let y0 = vec![1.0_f64];
let mut solver = CvodeBuilder::new(Method::BDF, &y0)
    .rtol(1e-8).atol(1e-10)
    .build(|_t, y, ydot| { ydot[0] = -y[0]; Ok(()) })
    .unwrap();

let (t, y) = solver.step(1.0).unwrap();
println!("y({t:.3}) = {:.8}  (exact {:.8})", y[0], (-t).exp());
```

## System requirements

**Without `vendored` feature** — SUNDIALS >= 6.0 must be installed:

```bash
# Ubuntu / Debian
sudo apt install libsundials-dev

# macOS
brew install sundials
```

**With `vendored` feature** — only `cmake` is required; SUNDIALS is downloaded
and built automatically.

## Documentation

Full API documentation and a worked guide (CVODE, CVODES, IDA, IDAS, Robertson
examples, common pitfalls) is available at
[docs.rs/sundials-rs](https://docs.rs/sundials-rs).

## License

Licensed under either of [MIT](LICENSE-MIT) or
[Apache-2.0](LICENSE-APACHE) at your option.

SUNDIALS itself is distributed under the
[BSD-3-Clause](https://github.com/LLNL/sundials/blob/main/LICENSE) license.
