//! Build script — embeds git/build metadata as compile-time environment
//! variables so the binary can report exactly which commit it was built from.
//!
//! Embedded variables:
//! - `GIT_SHORT_HASH`   — short commit hash, or `unknown`
//! - `BUILD_TIMESTAMP`  — UTC timestamp of the build (ISO-8601)

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

    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-env-changed=GIT_SHORT_HASH");
}
