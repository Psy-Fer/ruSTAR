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
    use crate::index::GenomeIndex;
    use crate::io::sam::SamWriter;
    use crate::stats::AlignmentStats;
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

    // 2. Detect single-end vs paired-end mode
    if params.read_files_in.is_empty() {
        anyhow::bail!("No read files specified (--readFilesIn)");
    }

    // 3. Open SAM writer
    let output_path = params.out_file_name_prefix.join("Aligned.out.sam");
    info!("Writing to {}", output_path.display());

    // Create output directory if it doesn't exist
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut writer = SamWriter::create(&output_path, &index.genome, params)?;

    // 4. Route to single-end or paired-end mode
    let stats = Arc::new(AlignmentStats::new());

    match params.read_files_in.len() {
        1 => align_reads_single_end(params, &index, &mut writer, &stats),
        2 => align_reads_paired_end(params, &index, &mut writer, &stats),
        n => anyhow::bail!("Invalid number of read files: {} (expected 1 or 2)", n),
    }?;

    // 5. Print summary
    info!("Alignment complete!");
    stats.print_summary();

    Ok(())
}

/// Align single-end reads
fn align_reads_single_end(
    params: &Parameters,
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    writer: &mut crate::io::sam::SamWriter,
    stats: &std::sync::Arc<crate::stats::AlignmentStats>,
) -> anyhow::Result<()> {
    use crate::align::read_align::align_read;
    use crate::io::fastq::{FastqReader, clip_read};
    use crate::io::sam::{BufferedSamRecords, SamWriter};
    use rayon::prelude::*;
    use std::sync::Arc;

    let read_file = &params.read_files_in[0];
    info!("Reading single-end from {}", read_file.display());

    let mut reader = FastqReader::open(read_file, params.read_files_command.as_deref())?;

    let stats = Arc::clone(stats);
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
                #[allow(clippy::needless_borrow)]
                let index = Arc::clone(&index);
                #[allow(clippy::needless_borrow)]
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

    Ok(())
}

/// Align paired-end reads
fn align_reads_paired_end(
    params: &Parameters,
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    writer: &mut crate::io::sam::SamWriter,
    stats: &std::sync::Arc<crate::stats::AlignmentStats>,
) -> anyhow::Result<()> {
    use crate::align::read_align::align_paired_read;
    use crate::io::fastq::{PairedFastqReader, clip_read};
    use crate::io::sam::{BufferedSamRecords, SamWriter};
    use rayon::prelude::*;
    use std::sync::Arc;

    info!(
        "Reading paired-end from {} and {}",
        params.read_files_in[0].display(),
        params.read_files_in[1].display()
    );

    let mut reader = PairedFastqReader::open(
        &params.read_files_in[0],
        &params.read_files_in[1],
        params.read_files_command.as_deref(),
    )?;

    let stats = Arc::clone(stats);
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

    info!("Aligning paired-end reads...");
    loop {
        // Sequential FASTQ reading
        let batch = reader.read_paired_batch(batch_size)?;
        if batch.is_empty() {
            break;
        }

        // Check max reads limit (pairs, not individual reads)
        let pairs_to_process = if read_count + batch.len() as u64 > max_reads {
            (max_reads - read_count) as usize
        } else {
            batch.len()
        };

        let batch_to_process = &batch[..pairs_to_process];

        // Parallel alignment processing
        let sam_buffers: Vec<Result<BufferedSamRecords, error::Error>> = batch_to_process
            .par_iter()
            .map(|paired_read| {
                #[allow(clippy::needless_borrow)]
                let index = Arc::clone(&index);
                #[allow(clippy::needless_borrow)]
                let stats = Arc::clone(&stats);

                // Apply clipping to both mates
                let (m1_seq, m1_qual) = clip_read(
                    &paired_read.mate1.sequence,
                    &paired_read.mate1.quality,
                    clip5p,
                    clip3p,
                );
                let (m2_seq, m2_qual) = clip_read(
                    &paired_read.mate2.sequence,
                    &paired_read.mate2.quality,
                    clip5p,
                    clip3p,
                );

                let mut buffer = BufferedSamRecords::new();

                // Skip if either mate is too short after clipping
                if m1_seq.is_empty() || m2_seq.is_empty() {
                    stats.record_alignment(0, max_multimaps);
                    if output_unmapped {
                        let records = SamWriter::build_paired_unmapped_records(
                            &paired_read.name,
                            &m1_seq,
                            &m1_qual,
                            &m2_seq,
                            &m2_qual,
                        )?;
                        for record in records {
                            buffer.push(record);
                        }
                    }
                    return Ok(buffer);
                }

                // Align paired read (CPU-intensive)
                let paired_alns = align_paired_read(&m1_seq, &m2_seq, &index, params)?;

                // Record stats (count pairs, not individual reads)
                stats.record_alignment(paired_alns.len(), max_multimaps);

                // Build SAM records (2 per pair)
                if paired_alns.is_empty() {
                    // Unmapped pair
                    if output_unmapped {
                        let records = SamWriter::build_paired_unmapped_records(
                            &paired_read.name,
                            &m1_seq,
                            &m1_qual,
                            &m2_seq,
                            &m2_qual,
                        )?;
                        for record in records {
                            buffer.push(record);
                        }
                    }
                } else if paired_alns.len() <= max_multimaps {
                    // Mapped pair (within multimap limit)
                    let records = SamWriter::build_paired_records(
                        &paired_read.name,
                        &m1_seq,
                        &m1_qual,
                        &m2_seq,
                        &m2_qual,
                        &paired_alns,
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

        // Sequential SAM writing
        for buffer_result in sam_buffers {
            let buffer = buffer_result?;
            writer.write_batch(&buffer.records)?;
        }

        read_count += pairs_to_process as u64;

        // Progress logging
        if read_count % 100000 < batch_size as u64 {
            info!("Processed {} pairs...", read_count);
        }

        if read_count >= max_reads {
            break;
        }
    }

    Ok(())
}
