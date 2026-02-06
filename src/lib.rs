#![allow(non_snake_case)]

pub mod error;
pub mod params;

// Future module stubs — uncomment as implemented:
// pub mod genome;
// pub mod index;
// pub mod align;
// pub mod junction;
// pub mod io;
// pub mod threading;
// pub mod chimeric;
// pub mod mapq;
// pub mod stats;

use log::info;

use crate::params::{Parameters, RunMode};

/// Top-level dispatcher. Called from `main()` after CLI parsing.
pub fn run(params: &Parameters) -> anyhow::Result<()> {
    params.validate()?;

    info!("ruSTAR v{}", env!("CARGO_PKG_VERSION"));
    info!("runMode: {}", params.run_mode);
    info!("runThreadN: {}", params.run_thread_n);

    match params.run_mode {
        RunMode::GenomeGenerate => genome_generate(params),
        RunMode::AlignReads => align_reads(params),
    }
}

fn genome_generate(params: &Parameters) -> anyhow::Result<()> {
    info!("genomeDir: {}", params.genome_dir.display());
    info!(
        "genomeFastaFiles: {:?}",
        params
            .genome_fasta_files
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
    );

    anyhow::bail!("genomeGenerate is not yet implemented (Phase 2+3)")
}

fn align_reads(params: &Parameters) -> anyhow::Result<()> {
    info!("genomeDir: {}", params.genome_dir.display());
    info!(
        "readFilesIn: {:?}",
        params
            .read_files_in
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
    );

    anyhow::bail!("alignReads is not yet implemented (Phase 4+)")
}
