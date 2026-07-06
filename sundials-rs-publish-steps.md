# sundials-rs — Manual steps before publishing to crates.io

## 1. Verify the repository URL

The `repository` field in both `Cargo.toml` files points to
`https://github.com/sn248/sundials-rs`. Make sure the repository exists
(push this workspace to GitHub first).

## 2. Create a crates.io account and API token

```bash
# Log in (opens browser for crates.io OAuth)
cargo login
```

## 3. Publish `sundials-rs-sys` first

```bash
cd sundials-rs-sys
cargo publish
```

## 4. Publish `sundials-rs`

`sundials-rs/Cargo.toml` declares the dependency with both `version` and
`path`, so no edit is needed: the workspace uses the path, and the published
package uses the version.

```bash
cd ../sundials-rs
cargo publish
```

Note: `sundials-rs-sys` must be indexed by crates.io before `sundials-rs`
can be published (usually takes a few minutes).

## 5. Verify on docs.rs

After publishing, docs.rs automatically builds and hosts the documentation.
Check that the GUIDE.md renders correctly at:

```
https://docs.rs/sundials-rs/0.1.0
```

---

## Notes

- The crate names are `sundials-rs` / `sundials-rs-sys` because `sundials`
  and `sundials-sys` were already taken on crates.io by other maintainers.
- The `vendored` feature requires `cmake` on the user's `PATH` to build
  SUNDIALS from source.
- The `SUNDIALS_DIR` environment variable always takes priority over both
  the system search and the vendored build.
