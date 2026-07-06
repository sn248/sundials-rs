use std::env;
use std::path::{Path, PathBuf};

// ── Version pinned here; bump to upgrade ──────────────────────────────────────
#[cfg(feature = "vendored")]
const SUNDIALS_VERSION: &str = "7.4.0";
#[cfg(feature = "vendored")]
const SUNDIALS_URL: &str = concat!(
    "https://github.com/LLNL/sundials/releases/download/",
    "v7.4.0/sundials-7.4.0.tar.gz"
);

fn main() {
    // Re-run only when this file or the header wrapper changes, or when
    // SUNDIALS_DIR is set/cleared.  Avoids re-running (and re-downloading) on
    // every source change.
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=wrapper.h");
    println!("cargo:rerun-if-env-changed=SUNDIALS_DIR");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    // ── 1. Locate (or build) SUNDIALS ─────────────────────────────────────────
    let (lib_dir, include_dir) = if let Ok(dir) = env::var("SUNDIALS_DIR") {
        // Explicit user-provided install — always takes priority.
        let root = PathBuf::from(&dir);
        (root.join("lib"), root.join("include"))
    } else if cfg!(feature = "vendored") {
        // Download SUNDIALS source and build a private static copy.
        build_vendored(&out_dir)
    } else {
        // Use whatever is installed on the system.
        locate_system()
    };

    // ── 2. Link ───────────────────────────────────────────────────────────────
    println!("cargo:rustc-link-search=native={}", lib_dir.display());

    if cfg!(feature = "vendored") {
        // Vendored build: we control exactly what was built, so link statically
        // against the precise set of libraries we asked CMake to produce.
        for lib in &[
            "sundials_cvodes",
            "sundials_idas",
            "sundials_nvecserial",
            "sundials_sunmatrixdense",
            "sundials_sunlinsoldense",
            "sundials_core",
        ] {
            println!("cargo:rustc-link-lib=static={lib}");
        }
        // SUNDIALS uses libm on Unix for math functions.
        #[cfg(unix)]
        println!("cargo:rustc-link-lib=m");
    } else {
        // System install: prefer cvodes over cvode, idas over ida, to avoid
        // duplicate CVodeCreate symbols (see comment in locate_system).
        let has = |name: &str| -> bool {
            lib_dir.join(format!("lib{name}.so")).exists()
                || lib_dir.join(format!("lib{name}.a")).exists()
        };

        if has("sundials_core") {
            println!("cargo:rustc-link-lib=sundials_core");
        }
        if has("sundials_cvodes") {
            println!("cargo:rustc-link-lib=sundials_cvodes");
        } else if has("sundials_cvode") {
            println!("cargo:rustc-link-lib=sundials_cvode");
        }
        if has("sundials_idas") {
            println!("cargo:rustc-link-lib=sundials_idas");
        } else if has("sundials_ida") {
            println!("cargo:rustc-link-lib=sundials_ida");
        }
        for lib in &[
            "sundials_nvecserial",
            "sundials_sunmatrixdense",
            "sundials_sunlinsoldense",
        ] {
            if has(lib) {
                println!("cargo:rustc-link-lib={lib}");
            }
        }
    }

    // ── 3. bindgen ────────────────────────────────────────────────────────────
    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .clang_arg(format!("-I{}", include_dir.display()))
        .allowlist_function("CV.*|IDA.*|N_V.*|SUN.*|Sundials.*")
        .allowlist_type(
            "CV.*|IDA.*|N_V.*|SUN.*|_generic.*|_sundials.*\
             |realtype|sunindextype|booleantype|SUNComm|SUNErrCode",
        )
        .allowlist_var("CV_.*|IDA_.*|SUN_.*|SUNTRUE|SUNFALSE")
        .derive_debug(true)
        .derive_default(true)
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("bindgen failed to generate bindings");

    bindings
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("failed to write bindings.rs");
}

// ── Vendored path ─────────────────────────────────────────────────────────────

/// Download, extract, and build SUNDIALS.  Returns (lib_dir, include_dir).
///
/// The cmake crate installs into `OUT_DIR/` so Cargo's caching is automatic:
/// if `OUT_DIR/lib/libsundials_cvodes.a` already exists (e.g. from a previous
/// build invocation that did not need to re-run this script), the cmake step
/// is skipped entirely.
#[cfg(feature = "vendored")]
fn build_vendored(out_dir: &Path) -> (PathBuf, PathBuf) {
    use flate2::read::GzDecoder;
    use std::fs;
    use tar::Archive;

    let tarball = out_dir.join(format!("sundials-{SUNDIALS_VERSION}.tar.gz"));
    let src_dir = out_dir.join(format!("sundials-{SUNDIALS_VERSION}"));
    // cmake crate installs to OUT_DIR by default.
    let lib_dir     = out_dir.join("lib");
    let include_dir = out_dir.join("include");

    // ── a) Download ───────────────────────────────────────────────────────────
    if !tarball.exists() {
        eprintln!("[sundials-rs-sys] Downloading SUNDIALS {SUNDIALS_VERSION}...");
        let resp = ureq::get(SUNDIALS_URL)
            .call()
            .expect("failed to download SUNDIALS; check network / proxy settings");

        let mut buf = Vec::new();
        use std::io::Read;
        resp.into_reader().read_to_end(&mut buf)
            .expect("failed to read SUNDIALS tarball");

        fs::write(&tarball, buf).expect("failed to write SUNDIALS tarball");
        eprintln!("[sundials-rs-sys] Download complete ({} bytes)", tarball.metadata().unwrap().len());
    }

    // ── b) Extract ────────────────────────────────────────────────────────────
    if !src_dir.exists() {
        eprintln!("[sundials-rs-sys] Extracting...");
        let f   = fs::File::open(&tarball).expect("failed to open tarball");
        let gz  = GzDecoder::new(f);
        let mut archive = Archive::new(gz);
        archive.unpack(out_dir).expect("failed to extract tarball");
    }

    // ── c) CMake configure + build + install ──────────────────────────────────
    if !lib_dir.join("libsundials_cvodes.a").exists() {
        eprintln!("[sundials-rs-sys] Building SUNDIALS with CMake (this takes ~1 minute the first time)...");

        cmake::Config::new(&src_dir)
            // Build type
            .define("CMAKE_BUILD_TYPE", "Release")
            // Static libraries only — no .so/.dylib produced.
            .define("BUILD_SHARED_LIBS", "OFF")
            .define("BUILD_STATIC_LIBS", "ON")
            // No MPI.
            .define("ENABLE_MPI", "OFF")
            // Skip tests and examples — they are slow to compile and not needed.
            .define("BUILD_TESTING",     "OFF")
            .define("EXAMPLES_ENABLE_C", "OFF")
            .define("EXAMPLES_ENABLE_CXX", "OFF")
            // Build only the solvers this crate exposes.
            // CVODES is a strict superset of CVODE; IDAS is a strict superset of
            // IDA.  Building both the plain and -S variants would create duplicate
            // symbols, so we build only the -S variants.
            .define("BUILD_CVODE",   "OFF")
            .define("BUILD_CVODES",  "ON")
            .define("BUILD_IDA",     "OFF")
            .define("BUILD_IDAS",    "ON")
            .define("BUILD_ARKODE",  "OFF")
            .define("BUILD_KINSOL",  "OFF")
            .build(); // returns OUT_DIR — the install prefix

        eprintln!("[sundials-rs-sys] SUNDIALS build complete.");
    }

    (lib_dir, include_dir)
}

// Stub so the file compiles without the feature — the real function is never
// called in this path, but Rust still needs to see a consistent return type.
#[cfg(not(feature = "vendored"))]
fn build_vendored(_out_dir: &Path) -> (PathBuf, PathBuf) {
    unreachable!("vendored feature is not enabled")
}

// ── System-library path ───────────────────────────────────────────────────────

fn locate_system() -> (PathBuf, PathBuf) {
    // a) SUNDIALS_DIR is handled in main() before this is called.

    // b) pkg-config (works when .pc files are present).
    if let Ok(lib_info) = pkg_config::Config::new()
        .atleast_version("6.0")
        .cargo_metadata(false)
        .probe("sundials_cvode")
    {
        if let (Some(lib_dir), Some(inc_dir)) = (
            lib_info.link_paths.into_iter().next(),
            lib_info.include_paths.into_iter().next(),
        ) {
            return (lib_dir, inc_dir);
        }
    }

    // c) Well-known filesystem locations.
    let candidates: &[(&str, &str)] = &[
        ("/usr/local/lib",             "/usr/local/include"),
        ("/usr/lib/x86_64-linux-gnu",  "/usr/include"),
        ("/usr/lib",                    "/usr/include"),
        ("/opt/sundials/lib",           "/opt/sundials/include"),
        ("/usr/local/lib64",           "/usr/local/include"),
    ];

    for (lib_str, inc_str) in candidates {
        let lib = PathBuf::from(lib_str);
        let inc = PathBuf::from(inc_str);
        let header = inc.join("cvode/cvode.h");
        if header.exists()
            && (lib.join("libsundials_cvode.so").exists()
                || lib.join("libsundials_cvode.a").exists())
        {
            return (lib, inc);
        }
    }

    panic!(
        "Could not locate SUNDIALS.\n\
         Options:\n\
         1. Install it:  sudo apt install libsundials-dev\n\
         2. Set env var: SUNDIALS_DIR=/path/to/install\n\
         3. Build from source: add features = [\"vendored\"] to sundials-rs-sys\n"
    );
}
