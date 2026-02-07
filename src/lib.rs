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
    use crate::io::sam::{BufferedSamRecords, SamWriter};
    use crate::stats::AlignmentStats;
    use rayon::prelude::*;
    use std::sync::Arc;

    info!("Starting read alignment...");

    // Configure Rayon thread pool based on --runThreadN
    if params.run_thread_n > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(params.run_thread_n)
            .build_global()
            .map_err(|e| {
                error::Error::Parameter(format!("Failed to configure thread pool: {}", e))
            })?;
        info!("Using {} threads for alignment", params.run_thread_n);
    } else {
        info!("Using single-threaded mode");
    }

    // 1. Load genome index
    info!("Loading genome index from {}", params.genome_dir.display());
    let index = Arc::new(GenomeIndex::load(&params.genome_dir, params)?);
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

    // 4. Process reads in batches
    let stats = Arc::new(AlignmentStats::new());
    let mut read_count = 0u64;
    let max_reads = if params.read_map_number < 0 {
        u64::MAX
    } else {
        params.read_map_number as u64
    };

    let batch_size = 10000;
    let clip5p = params.clip5p_nbases as usize;
    let clip3p = params.clip3p_nbases as usize;
    let max_multimaps = params.out_filter_multimap_nmax as usize;
    let output_unmapped = params.out_sam_unmapped != params::OutSamUnmapped::None;

    info!("Aligning reads...");
    loop {
        // Sequential FASTQ reading (unavoidable bottleneck)
        let batch = reader.read_batch(batch_size)?;
        if batch.is_empty() {
            break;
        }

        // Check max reads limit
        let reads_to_process = if read_count + batch.len() as u64 > max_reads {
            (max_reads - read_count) as usize
        } else {
            batch.len()
        };

        let batch_to_process = &batch[..reads_to_process];

        // Parallel alignment processing
        let sam_buffers: Vec<Result<BufferedSamRecords, error::Error>> = batch_to_process
            .par_iter()
            .map(|read| {
                let index = Arc::clone(&index);
                let stats = Arc::clone(&stats);

                // Apply clipping
                let (clipped_seq, clipped_qual) =
                    clip_read(&read.sequence, &read.quality, clip5p, clip3p);

                let mut buffer = BufferedSamRecords::new();

                // Skip if read is too short after clipping
                if clipped_seq.is_empty() {
                    stats.record_alignment(0, max_multimaps);
                    if output_unmapped {
                        let record = SamWriter::build_unmapped_record(
                            &read.name,
                            &clipped_seq,
                            &clipped_qual,
                        )?;
                        buffer.push(record);
                    }
                    return Ok(buffer);
                }

                // Align read (CPU-intensive, pure function)
                let transcripts = align_read(&clipped_seq, &index, params)?;

                // Record stats (atomic, lock-free)
                stats.record_alignment(transcripts.len(), max_multimaps);

                // Build SAM records (no I/O, just construction)
                if transcripts.is_empty() {
                    // Unmapped
                    if output_unmapped {
                        let record = SamWriter::build_unmapped_record(
                            &read.name,
                            &clipped_seq,
                            &clipped_qual,
                        )?;
                        buffer.push(record);
                    }
                } else if transcripts.len() <= max_multimaps {
                    // Mapped (within multimap limit)
                    let records = SamWriter::build_alignment_records(
                        &read.name,
                        &clipped_seq,
                        &clipped_qual,
                        &transcripts,
                        &index.genome,
                        params,
                    )?;
                    for record in records {
                        buffer.push(record);
                    }
                }
                // else: too many loci, skip output

                Ok(buffer)
            })
            .collect();

        // Sequential SAM writing (merge buffers in chunk order)
        for buffer_result in sam_buffers {
            let buffer = buffer_result?;
            writer.write_batch(&buffer.records)?;
        }

        read_count += reads_to_process as u64;

        // Progress logging
        if read_count % 100000 < batch_size as u64 {
            info!("Processed {} reads...", read_count);
        }

        if read_count >= max_reads {
            break;
        }
    }

    // 5. Print summary
    info!("Alignment complete!");
    stats.print_summary();

    Ok(())
}
