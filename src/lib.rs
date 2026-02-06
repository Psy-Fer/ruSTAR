#![allow(non_snake_case)]

pub mod error;
pub mod params;

pub mod align;
pub mod genome;
pub mod index;

// Future module stubs — uncomment as implemented:
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
    use index::GenomeIndex;

    info!("genomeDir: {}", params.genome_dir.display());
    info!(
        "genomeFastaFiles: {:?}",
        params
            .genome_fasta_files
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
    );

    info!("Building genome index...");
    let index = GenomeIndex::build(params)?;

    info!(
        "Writing index files to {}...",
        params.genome_dir.display()
    );
    index.write(&params.genome_dir, params)?;

    info!("Genome generation complete!");
    Ok(())
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
