#![allow(non_snake_case)]

pub mod error;
pub mod params;

pub mod align;
pub mod genome;
pub mod index;
pub mod io;
pub mod mapq;
pub mod stats;

// Future module stubs — uncomment as implemented:
// pub mod junction;
// pub mod threading;
// pub mod chimeric;

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

    info!("Writing index files to {}...", params.genome_dir.display());
    index.write(&params.genome_dir, params)?;

    info!("Genome generation complete!");
    Ok(())
}

fn align_reads(params: &Parameters) -> anyhow::Result<()> {
    use crate::align::read_align::align_read;
    use crate::index::GenomeIndex;
    use crate::io::fastq::{FastqReader, clip_read};
    use crate::io::sam::SamWriter;
    use crate::stats::AlignmentStats;

    info!("Starting read alignment...");

    // 1. Load genome index
    info!("Loading genome index from {}", params.genome_dir.display());
    let index = GenomeIndex::load(&params.genome_dir, params)?;
    info!(
        "Loaded {} chromosomes, {} bases",
        index.genome.n_chr_real, index.genome.n_genome
    );

    // 2. Open FASTQ reader
    if params.read_files_in.is_empty() {
        anyhow::bail!("No read files specified (--readFilesIn)");
    }
    let read_file = &params.read_files_in[0];
    info!("Reading from {}", read_file.display());

    let mut reader = FastqReader::open(read_file, params.read_files_command.as_deref())?;

    // 3. Open SAM writer
    let output_path = params.out_file_name_prefix.join("Aligned.out.sam");
    info!("Writing to {}", output_path.display());

    // Create output directory if it doesn't exist
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut writer = SamWriter::create(&output_path, &index.genome, params)?;

    // 4. Process reads
    let mut stats = AlignmentStats::default();
    let mut read_count = 0u64;
    let max_reads = if params.read_map_number < 0 {
        u64::MAX
    } else {
        params.read_map_number as u64
    };

    info!("Aligning reads...");
    while let Some(read) = reader.next_encoded()? {
        if read_count >= max_reads {
            break;
        }
        read_count += 1;

        if read_count % 100000 == 0 {
            info!("Processed {} reads...", read_count);
        }

        // Apply clipping
        let (clipped_seq, clipped_qual) = clip_read(
            &read.sequence,
            &read.quality,
            params.clip5p_nbases as usize,
            params.clip3p_nbases as usize,
        );

        // Skip if read is too short after clipping
        if clipped_seq.is_empty() {
            stats.record_alignment(0, params.out_filter_multimap_nmax as usize);
            if params.out_sam_unmapped != params::OutSamUnmapped::None {
                writer.write_unmapped(&read.name, &clipped_seq, &clipped_qual)?;
            }
            continue;
        }

        // Align read
        let transcripts = align_read(&clipped_seq, &index, params)?;

        // Record stats
        stats.record_alignment(transcripts.len(), params.out_filter_multimap_nmax as usize);

        // Write to SAM
        if transcripts.is_empty() {
            // Unmapped
            if params.out_sam_unmapped != params::OutSamUnmapped::None {
                writer.write_unmapped(&read.name, &clipped_seq, &clipped_qual)?;
            }
        } else if transcripts.len() <= params.out_filter_multimap_nmax as usize {
            // Mapped (within multimap limit)
            writer.write_alignment(
                &read.name,
                &clipped_seq,
                &clipped_qual,
                &transcripts,
                &index.genome,
                params,
            )?;
        }
        // else: too many loci, skip output
    }

    // 5. Print summary
    info!("Alignment complete!");
    stats.print_summary();

    Ok(())
}
