//! Build script — embeds git/build metadata as compile-time environment
//! variables so the binary can report exactly which commit and target it
//! was built from.
//!
//! Embedded variables:
//! - `GIT_SHORT_HASH` — short commit hash, or `unknown`
//! - `BUILD_TIMESTAMP` — UTC timestamp of the build (ISO-8601)
//! - `VERSION_BODY` — pre-formatted `--version` detail line, e.g.
//!   `c0892fb - linux/aarch64 (NEON) - built 2026-04-16T21:36:44Z`

use std::process::Command;

fn main() {
    // Prefer the GIT_SHORT_HASH env var (set via Docker build arg or CI),
    // fall back to running `git rev-parse --short HEAD`.
    let git_hash = std::env::var("GIT_SHORT_HASH")
        .ok()
        .filter(|s| !s.is_empty() && s != "unknown")
        .unwrap_or_else(|| {
            Command::new("git")
                .args(["rev-parse", "--short", "HEAD"])
                .output()
                .ok()
                .filter(|o| o.status.success())
                .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
                .unwrap_or_else(|| "unknown".to_string())
        });
    println!("cargo:rustc-env=GIT_SHORT_HASH={git_hash}");

    let build_ts = chrono::Utc::now().format("%Y-%m-%dT%H:%M:%SZ").to_string();
    println!("cargo:rustc-env=BUILD_TIMESTAMP={build_ts}");

    // Derive the SIMD level and label from the TARGET being compiled for,
    // using the `CARGO_CFG_*` env vars Cargo exposes to build scripts.
    let target_arch = std::env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let target_features = std::env::var("CARGO_CFG_TARGET_FEATURE").unwrap_or_default();
    let has_feature = |f: &str| target_features.split(',').any(|x| x == f);

    let (target, label) = match target_arch.as_str() {
        "x86_64" => {
            if has_feature("avx512f") {
                ("x86-64-v4", "AVX-512")
            } else if has_feature("avx2") {
                ("x86-64-v3", "AVX2")
            } else {
                ("x86-64", "baseline")
            }
        }
        "aarch64" => {
            if has_feature("sve") {
                ("aarch64-neoverse-v1", "SVE")
            } else {
                ("aarch64", "NEON")
            }
        }
        _ => ("unknown", "baseline"),
    };

    let os = match target_os.as_str() {
        "linux" | "macos" | "windows" => target_os.as_str(),
        _ => "unknown",
    };

    println!(
        "cargo:rustc-env=VERSION_BODY={git_hash} - {os}/{target} ({label}) - built {build_ts}"
    );

    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-env-changed=GIT_SHORT_HASH");
}
