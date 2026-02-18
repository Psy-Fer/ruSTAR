#![allow(non_snake_case)]

pub mod error;
pub mod params;

pub mod align;
pub mod chimeric;
pub mod genome;
pub mod index;
pub mod io;
pub mod junction;
pub mod mapq;
pub mod stats;

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

/// Trait for alignment output writers (SAM or BAM)
trait AlignmentWriter {
    fn write_batch(
        &mut self,
        batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error>;
}

/// Null writer that discards all output (for two-pass mode pass 1)
struct NullWriter;

impl AlignmentWriter for NullWriter {
    fn write_batch(
        &mut self,
        _batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error> {
        Ok(()) // Discard all records
    }
}

impl AlignmentWriter for crate::io::sam::SamWriter {
    fn write_batch(
        &mut self,
        batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error> {
        self.write_batch(batch)
    }
}

impl AlignmentWriter for crate::io::bam::BamWriter {
    fn write_batch(
        &mut self,
        batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error> {
        self.write_batch(batch)
    }
}

fn align_reads(params: &Parameters) -> anyhow::Result<()> {
    use crate::index::GenomeIndex;

    use crate::params::TwopassMode;

    use std::sync::Arc;

    let time_start = chrono::Local::now();

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

    // Validate read files
    if params.read_files_in.is_empty() {
        anyhow::bail!("No read files specified (--readFilesIn)");
    }

    // 1. Load genome index
    info!("Loading genome index from {}", params.genome_dir.display());
    let index = Arc::new(GenomeIndex::load(&params.genome_dir, params)?);
    info!(
        "Loaded {} chromosomes, {} bases",
        index.genome.n_chr_real, index.genome.n_genome
    );

    let time_map_start = chrono::Local::now();

    // 2. Dispatch based on two-pass mode
    let stats = match params.twopass_mode {
        TwopassMode::None => {
            info!("Running single-pass alignment");
            run_single_pass(&index, params)?
        }
        TwopassMode::Basic => {
            info!("Running two-pass alignment mode");
            run_two_pass(&index, params)?
        }
    };

    let time_finish = chrono::Local::now();

    // Write Log.final.out
    let log_path = params.out_file_name_prefix.join("Log.final.out");
    if let Some(parent) = log_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    stats.write_log_final(&log_path, time_start, time_map_start, time_finish)?;
    info!("Wrote {}", log_path.display());

    info!("Alignment complete!");
    Ok(())
}

/// Run single-pass alignment (original logic)
fn run_single_pass(
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    params: &Parameters,
) -> anyhow::Result<std::sync::Arc<crate::stats::AlignmentStats>> {
    use crate::io::bam::BamWriter;
    use crate::io::sam::SamWriter;
    use crate::params::OutSamFormat;
    use std::sync::Arc;

    // Initialize statistics collectors
    let stats = Arc::new(crate::stats::AlignmentStats::new());
    let sj_stats = Arc::new(crate::junction::SpliceJunctionStats::new());

    // 4. Route to SAM or BAM output based on --outSAMtype
    let out_type = params
        .out_sam_type()
        .map_err(|e| anyhow::anyhow!("Invalid --outSAMtype: {}", e))?;

    match out_type.format {
        OutSamFormat::Sam => {
            let output_path = params.out_file_name_prefix.join("Aligned.out.sam");
            info!("Writing SAM to {}", output_path.display());

            // Create output directory if it doesn't exist
            if let Some(parent) = output_path.parent() {
                std::fs::create_dir_all(parent)?;
            }

            let mut writer = SamWriter::create(&output_path, &index.genome, params)?;

            // Route to single-end or paired-end mode
            match params.read_files_in.len() {
                1 => align_reads_single_end(params, index, &mut writer, &stats, &sj_stats),
                2 => align_reads_paired_end(params, index, &mut writer, &stats, &sj_stats),
                n => anyhow::bail!("Invalid number of read files: {} (expected 1 or 2)", n),
            }?;
        }
        OutSamFormat::Bam => {
            let output_path = params.out_file_name_prefix.join("Aligned.out.bam");
            info!("Writing BAM to {}", output_path.display());

            // Create output directory if it doesn't exist
            if let Some(parent) = output_path.parent() {
                std::fs::create_dir_all(parent)?;
            }

            let mut writer = BamWriter::create(&output_path, &index.genome, params)?;

            // Route to single-end or paired-end mode (same functions as SAM, generic!)
            match params.read_files_in.len() {
                1 => align_reads_single_end(params, index, &mut writer, &stats, &sj_stats),
                2 => align_reads_paired_end(params, index, &mut writer, &stats, &sj_stats),
                n => anyhow::bail!("Invalid number of read files: {} (expected 1 or 2)", n),
            }?;

            // Finish BAM file (flush BGZF buffers)
            writer.finish()?;
        }
        OutSamFormat::None => {
            info!("Output format set to None, skipping alignment output");
            anyhow::bail!("Output format 'None' not yet implemented");
        }
    }

    // 5. Write SJ.out.tab file
    let sj_output_path = params.out_file_name_prefix.join("SJ.out.tab");
    if !sj_stats.is_empty() {
        info!(
            "Writing splice junction statistics to {}",
            sj_output_path.display()
        );
        sj_stats.write_output(&sj_output_path, &index.genome, params)?;
    }

    // 6. Print summary
    stats.print_summary();

    Ok(stats)
}

/// Run two-pass alignment mode
fn run_two_pass(
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    params: &Parameters,
) -> anyhow::Result<std::sync::Arc<crate::stats::AlignmentStats>> {
    use std::sync::Arc;

    // PASS 1: Junction discovery
    info!("Two-pass mode: Pass 1 - Junction discovery");
    let (sj_stats_pass1, novel_junctions) = run_pass1(index, params)?;

    // Write SJ.pass1.out.tab
    let pass1_path = params.out_file_name_prefix.join("SJ.pass1.out.tab");

    // Create output directory if it doesn't exist
    if let Some(parent) = pass1_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    info!("Writing pass 1 junctions to {}", pass1_path.display());
    sj_stats_pass1.write_output(&pass1_path, &index.genome, params)?;
    info!(
        "Pass 1 discovered {} novel junctions",
        novel_junctions.len()
    );

    // Insert novel junctions into DB
    let mut merged_index = (**index).clone();
    merged_index
        .junction_db
        .insert_novel(novel_junctions.clone());
    info!(
        "Merged junction DB: {} total junctions",
        merged_index.junction_db.len()
    );

    // PASS 2: Re-alignment with merged DB
    info!("Two-pass mode: Pass 2 - Re-alignment");
    let stats = run_single_pass(&Arc::new(merged_index), params)?;

    Ok(stats)
}

/// Run pass 1 of two-pass mode (junction discovery)
fn run_pass1(
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    params: &Parameters,
) -> anyhow::Result<(
    crate::junction::SpliceJunctionStats,
    Vec<(
        crate::junction::NovelJunctionKey,
        crate::junction::JunctionInfo,
    )>,
)> {
    use std::sync::Arc;

    let stats = Arc::new(crate::stats::AlignmentStats::new());
    let sj_stats = Arc::new(crate::junction::SpliceJunctionStats::new());

    // Modify params to limit reads for pass 1
    let mut params_pass1 = params.clone();
    if params.twopass1_reads_n >= 0 {
        params_pass1.read_map_number = params.twopass1_reads_n;
        info!("Pass 1 will align {} reads", params.twopass1_reads_n);
    } else {
        info!("Pass 1 will align all reads");
    }

    // Create NullWriter (discard SAM/BAM output in pass 1)
    let mut null_writer = NullWriter;

    // Align reads (single-end or paired-end)
    match params.read_files_in.len() {
        1 => align_reads_single_end(&params_pass1, index, &mut null_writer, &stats, &sj_stats)?,
        2 => align_reads_paired_end(&params_pass1, index, &mut null_writer, &stats, &sj_stats)?,
        n => anyhow::bail!("Invalid number of read files: {} (expected 1 or 2)", n),
    }

    info!("Pass 1 aligned {} reads", stats.total_reads());

    // Filter novel junctions
    let novel_junctions = crate::junction::filter_novel_junctions(&sj_stats, params);

    // Return ownership of sj_stats
    let sj_stats = Arc::try_unwrap(sj_stats).unwrap_or_else(|arc| (*arc).clone());

    Ok((sj_stats, novel_junctions))
}

/// Helper struct to hold alignment results from parallel processing
struct AlignmentBatchResults {
    sam_records: crate::io::sam::BufferedSamRecords,
    chimeric_alns: Vec<crate::chimeric::ChimericAlignment>,
    /// Junction keys from the primary (best) alignment for BySJout filtering.
    /// Empty if unmapped or no junctions.
    primary_junction_keys: Vec<crate::junction::SjKey>,
}

/// Extract SjKey junction identifiers from a transcript's CIGAR.
/// Used to check if a read's junctions survive outSJfilter* for BySJout mode.
fn extract_junction_keys(
    transcript: &crate::align::transcript::Transcript,
    index: &crate::index::GenomeIndex,
) -> Vec<crate::junction::SjKey> {
    use crate::align::score::AlignmentScorer;
    use crate::align::transcript::CigarOp;

    let scorer = AlignmentScorer::from_params_minimal();
    let mut keys = Vec::new();
    let mut genome_pos = transcript.genome_start;

    for op in &transcript.cigar {
        match op {
            CigarOp::RefSkip(len) => {
                let intron_len = *len;
                let intron_start = genome_pos + 1;
                let intron_end = genome_pos + intron_len as u64;

                let motif = scorer.detect_splice_motif(genome_pos, intron_len, &index.genome);
                let strand = match motif.implied_strand() {
                    Some('+') => 1u8,
                    Some('-') => 2u8,
                    _ => 0u8,
                };
                let encoded_motif = crate::junction::encode_motif(motif);

                keys.push(crate::junction::SjKey {
                    chr_idx: transcript.chr_idx,
                    intron_start,
                    intron_end,
                    strand,
                    motif: encoded_motif,
                });

                genome_pos += intron_len as u64;
            }
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::Del(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) => {}
        }
    }

    keys
}

/// Align single-end reads
fn align_reads_single_end<W: AlignmentWriter>(
    params: &Parameters,
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    writer: &mut W,
    stats: &std::sync::Arc<crate::stats::AlignmentStats>,
    sj_stats: &std::sync::Arc<crate::junction::SpliceJunctionStats>,
) -> anyhow::Result<()> {
    use crate::align::read_align::align_read;
    use crate::io::fastq::{FastqReader, clip_read};
    use crate::io::sam::{BufferedSamRecords, SamWriter};
    use crate::params::OutFilterType;
    use rayon::prelude::*;
    use std::sync::Arc;

    let read_file = &params.read_files_in[0];
    info!("Reading single-end from {}", read_file.display());

    let mut reader = FastqReader::open(read_file, params.read_files_command.as_deref())?;

    // Create chimeric output writer if enabled
    let mut chimeric_writer = if params.chim_segment_min > 0 {
        use crate::chimeric::ChimericJunctionWriter;
        let prefix = params.out_file_name_prefix.to_str().unwrap_or(".");
        info!(
            "Chimeric detection enabled (chimSegmentMin={})",
            params.chim_segment_min
        );
        Some(ChimericJunctionWriter::new(prefix)?)
    } else {
        None
    };

    let stats = Arc::clone(stats);
    let sj_stats = Arc::clone(sj_stats);
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
    let by_sjout = params.out_filter_type == OutFilterType::BySJout;

    // Buffer for BySJout mode: accumulate all results before filtering
    let mut bysj_buffer: Vec<AlignmentBatchResults> = Vec::new();

    if by_sjout {
        info!("outFilterType=BySJout: buffering reads for post-alignment junction filtering");
    }

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
        let batch_results: Vec<Result<AlignmentBatchResults, error::Error>> = batch_to_process
            .par_iter()
            .map(|read| {
                #[allow(clippy::needless_borrow)]
                let index = Arc::clone(&index);
                #[allow(clippy::needless_borrow)]
                let stats = Arc::clone(&stats);
                #[allow(clippy::needless_borrow)]
                let sj_stats = Arc::clone(&sj_stats);

                // Apply clipping
                let (clipped_seq, clipped_qual) =
                    clip_read(&read.sequence, &read.quality, clip5p, clip3p);

                let mut buffer = BufferedSamRecords::new();
                let mut chimeric_alns = Vec::new();

                // Record read bases for Log.final.out
                stats.record_read_bases(clipped_seq.len() as u64);

                // Skip if read is too short after clipping
                if clipped_seq.is_empty() {
                    stats.record_alignment(0, max_multimaps);
                    stats.record_unmapped_reason(crate::stats::UnmappedReason::Other);
                    if output_unmapped {
                        let record = SamWriter::build_unmapped_record(
                            &read.name,
                            &clipped_seq,
                            &clipped_qual,
                        )?;
                        buffer.push(record);
                    }
                    return Ok(AlignmentBatchResults {
                        sam_records: buffer,
                        chimeric_alns,
                        primary_junction_keys: Vec::new(),
                    });
                }

                // Align read (CPU-intensive, pure function)
                let (transcripts, chimeric_results, n_for_mapq, unmapped_reason) =
                    align_read(&clipped_seq, &read.name, &index, params)?;

                // Collect chimeric alignments if enabled
                if params.chim_segment_min > 0 {
                    chimeric_alns.extend(chimeric_results);
                    if !chimeric_alns.is_empty() {
                        stats.record_chimeric();
                    }
                }

                // Record stats (atomic, lock-free)
                let n = transcripts.len();
                stats.record_alignment(n, max_multimaps);
                if n == 0 {
                    stats.record_unmapped_reason(
                        unmapped_reason.unwrap_or(crate::stats::UnmappedReason::Other),
                    );
                } else if n == 1 {
                    stats.record_transcript_stats(&transcripts[0]);
                }

                // Record junction statistics
                let is_unique = transcripts.len() == 1;
                for transcript in &transcripts {
                    record_transcript_junctions(transcript, &index, &sj_stats, is_unique);
                }

                // Extract junction keys from primary alignment for BySJout filtering
                let primary_junction_keys =
                    if by_sjout && !transcripts.is_empty() && transcripts[0].n_junction > 0 {
                        extract_junction_keys(&transcripts[0], &index)
                    } else {
                        Vec::new()
                    };

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
                        n_for_mapq,
                    )?;
                    for record in records {
                        buffer.push(record);
                    }
                }
                // else: too many loci, skip output

                Ok(AlignmentBatchResults {
                    sam_records: buffer,
                    chimeric_alns,
                    primary_junction_keys,
                })
            })
            .collect();

        if by_sjout {
            // Buffer all results for post-alignment filtering
            for result in batch_results {
                bysj_buffer.push(result?);
            }
        } else {
            // Normal mode: sequential writing (merge buffers in chunk order)
            for result in batch_results {
                let batch = result?;

                // Write SAM/BAM records
                writer.write_batch(&batch.sam_records.records)?;

                // Write chimeric alignments
                if let Some(ref mut chim_writer) = chimeric_writer {
                    for chim_aln in &batch.chimeric_alns {
                        chim_writer.write_alignment(
                            chim_aln,
                            &index.genome.chr_name,
                            &chim_aln.read_name,
                        )?;
                    }
                }
            }
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

    // BySJout post-alignment filtering
    if by_sjout {
        let surviving_junctions = sj_stats.compute_surviving_junctions(params);
        info!(
            "BySJout filtering: {} surviving junctions from {} total",
            surviving_junctions.len(),
            sj_stats.len()
        );

        let mut filtered_count = 0u64;
        for batch in &bysj_buffer {
            if !batch.primary_junction_keys.is_empty() {
                // Read has junctions — check if ALL survive
                let all_survive = batch
                    .primary_junction_keys
                    .iter()
                    .all(|key| surviving_junctions.contains(key));

                if !all_survive {
                    // Filter this read: skip writing, adjust stats
                    filtered_count += 1;
                    stats.undo_mapped_record_bysj();
                    continue;
                }
            }

            // No junctions or all junctions survive — write normally
            writer.write_batch(&batch.sam_records.records)?;

            // Write chimeric alignments
            if let Some(ref mut chim_writer) = chimeric_writer {
                for chim_aln in &batch.chimeric_alns {
                    chim_writer.write_alignment(
                        chim_aln,
                        &index.genome.chr_name,
                        &chim_aln.read_name,
                    )?;
                }
            }
        }

        info!(
            "BySJout: filtered {} reads with non-surviving junctions",
            filtered_count
        );
    }

    // Flush chimeric output if enabled
    if let Some(ref mut chim_writer) = chimeric_writer {
        chim_writer.flush()?;
        info!("Chimeric junction output complete");
    }

    Ok(())
}

/// Align paired-end reads
fn align_reads_paired_end<W: AlignmentWriter>(
    params: &Parameters,
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    writer: &mut W,
    stats: &std::sync::Arc<crate::stats::AlignmentStats>,
    sj_stats: &std::sync::Arc<crate::junction::SpliceJunctionStats>,
) -> anyhow::Result<()> {
    use crate::align::read_align::{PairedAlignment, PairedAlignmentResult, align_paired_read};
    use crate::io::fastq::{PairedFastqReader, clip_read};
    use crate::io::sam::{BufferedSamRecords, SamWriter};
    use crate::params::OutFilterType;
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

    // Create chimeric output writer if enabled (paired-end chimeric detection not yet implemented)
    let mut chimeric_writer = if params.chim_segment_min > 0 {
        use crate::chimeric::ChimericJunctionWriter;
        let prefix = params.out_file_name_prefix.to_str().unwrap_or(".");
        info!(
            "Chimeric detection enabled (chimSegmentMin={}) - paired-end chimeric detection not yet implemented",
            params.chim_segment_min
        );
        Some(ChimericJunctionWriter::new(prefix)?)
    } else {
        None
    };

    let stats = Arc::clone(stats);
    let sj_stats = Arc::clone(sj_stats);
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
    let by_sjout = params.out_filter_type == OutFilterType::BySJout;

    // Buffer for BySJout mode
    let mut bysj_buffer: Vec<AlignmentBatchResults> = Vec::new();

    if by_sjout {
        info!("outFilterType=BySJout: buffering pairs for post-alignment junction filtering");
    }

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
        let batch_results: Vec<Result<AlignmentBatchResults, error::Error>> = batch_to_process
            .par_iter()
            .map(|paired_read| {
                #[allow(clippy::needless_borrow)]
                let index = Arc::clone(&index);
                #[allow(clippy::needless_borrow)]
                let stats = Arc::clone(&stats);
                #[allow(clippy::needless_borrow)]
                let sj_stats = Arc::clone(&sj_stats);

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

                // Record read bases for Log.final.out (both mates)
                stats.record_read_bases(m1_seq.len() as u64 + m2_seq.len() as u64);

                // Skip if either mate is too short after clipping
                if m1_seq.is_empty() || m2_seq.is_empty() {
                    stats.record_alignment(0, max_multimaps);
                    stats.record_unmapped_reason(crate::stats::UnmappedReason::Other);
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
                    return Ok(AlignmentBatchResults {
                        sam_records: buffer,
                        chimeric_alns: Vec::new(),
                        primary_junction_keys: Vec::new(),
                    });
                }

                // Align paired read (CPU-intensive)
                let (results, n_for_mapq, unmapped_reason) =
                    align_paired_read(&m1_seq, &m2_seq, &index, params)?;

                // Classify the result for stats and SAM output
                let has_half_mapped = results
                    .iter()
                    .any(|r| matches!(r, PairedAlignmentResult::HalfMapped { .. }));
                let both_mapped: Vec<_> = results
                    .iter()
                    .filter_map(|r| {
                        if let PairedAlignmentResult::BothMapped(pa) = r {
                            Some(pa)
                        } else {
                            None
                        }
                    })
                    .collect();

                if results.is_empty() {
                    // Both mates unmapped
                    stats.record_alignment(0, max_multimaps);
                    stats.record_unmapped_reason(
                        unmapped_reason.unwrap_or(crate::stats::UnmappedReason::Other),
                    );
                } else if has_half_mapped {
                    // Half-mapped: count as mapped for the mapped mate
                    stats.record_alignment(1, max_multimaps);
                    stats.record_half_mapped();
                    // Record transcript stats from the mapped mate only
                    if let Some(PairedAlignmentResult::HalfMapped {
                        mapped_transcript, ..
                    }) = results.first()
                    {
                        stats.record_transcript_stats(mapped_transcript);
                    }
                } else {
                    // Both-mapped pairs
                    let n = both_mapped.len();
                    stats.record_alignment(n, max_multimaps);
                    if n == 1 {
                        stats.record_transcript_stats(&both_mapped[0].mate1_transcript);
                        stats.record_transcript_stats(&both_mapped[0].mate2_transcript);
                    }
                }

                // Record junction statistics
                let is_unique = both_mapped.len() == 1 || (has_half_mapped && results.len() == 1);
                for result in &results {
                    match result {
                        PairedAlignmentResult::BothMapped(pair) => {
                            record_transcript_junctions(
                                &pair.mate1_transcript,
                                &index,
                                &sj_stats,
                                is_unique,
                            );
                            record_transcript_junctions(
                                &pair.mate2_transcript,
                                &index,
                                &sj_stats,
                                is_unique,
                            );
                        }
                        PairedAlignmentResult::HalfMapped {
                            mapped_transcript, ..
                        } => {
                            record_transcript_junctions(
                                mapped_transcript,
                                &index,
                                &sj_stats,
                                is_unique,
                            );
                        }
                    }
                }

                // Extract junction keys from primary alignment for BySJout
                let primary_junction_keys = if by_sjout && !results.is_empty() {
                    let mut keys = Vec::new();
                    match &results[0] {
                        PairedAlignmentResult::BothMapped(pair) => {
                            if pair.mate1_transcript.n_junction > 0 {
                                keys.extend(extract_junction_keys(&pair.mate1_transcript, &index));
                            }
                            if pair.mate2_transcript.n_junction > 0 {
                                keys.extend(extract_junction_keys(&pair.mate2_transcript, &index));
                            }
                        }
                        PairedAlignmentResult::HalfMapped {
                            mapped_transcript, ..
                        } => {
                            if mapped_transcript.n_junction > 0 {
                                keys.extend(extract_junction_keys(mapped_transcript, &index));
                            }
                        }
                    }
                    keys
                } else {
                    Vec::new()
                };

                // Build SAM records
                if results.is_empty() {
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
                } else if has_half_mapped {
                    // Half-mapped pair
                    if let Some(PairedAlignmentResult::HalfMapped {
                        mapped_transcript,
                        mate1_is_mapped,
                    }) = results.first()
                    {
                        let records = SamWriter::build_half_mapped_records(
                            &paired_read.name,
                            &m1_seq,
                            &m1_qual,
                            &m2_seq,
                            &m2_qual,
                            mapped_transcript,
                            *mate1_is_mapped,
                            &index.genome,
                            params,
                            n_for_mapq,
                        )?;
                        for record in records {
                            buffer.push(record);
                        }
                    }
                } else if both_mapped.len() <= max_multimaps {
                    // Both-mapped pairs (within multimap limit)
                    // Extract PairedAlignments for the existing build_paired_records
                    let paired_alns: Vec<PairedAlignment> = both_mapped
                        .iter()
                        .map(|pa| PairedAlignment::clone(pa))
                        .collect();
                    let records = SamWriter::build_paired_records(
                        &paired_read.name,
                        &m1_seq,
                        &m1_qual,
                        &m2_seq,
                        &m2_qual,
                        &paired_alns,
                        &index.genome,
                        params,
                        n_for_mapq,
                    )?;
                    for record in records {
                        buffer.push(record);
                    }
                }
                // else: too many loci, skip output

                Ok(AlignmentBatchResults {
                    sam_records: buffer,
                    chimeric_alns: Vec::new(),
                    primary_junction_keys,
                })
            })
            .collect();

        if by_sjout {
            // Buffer all results for post-alignment filtering
            for result in batch_results {
                bysj_buffer.push(result?);
            }
        } else {
            // Normal mode: sequential SAM writing
            for result in batch_results {
                let batch = result?;
                writer.write_batch(&batch.sam_records.records)?;
            }
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

    // BySJout post-alignment filtering
    if by_sjout {
        let surviving_junctions = sj_stats.compute_surviving_junctions(params);
        info!(
            "BySJout filtering: {} surviving junctions from {} total",
            surviving_junctions.len(),
            sj_stats.len()
        );

        let mut filtered_count = 0u64;
        for batch in &bysj_buffer {
            if !batch.primary_junction_keys.is_empty() {
                let all_survive = batch
                    .primary_junction_keys
                    .iter()
                    .all(|key| surviving_junctions.contains(key));

                if !all_survive {
                    filtered_count += 1;
                    stats.undo_mapped_record_bysj();
                    continue;
                }
            }

            writer.write_batch(&batch.sam_records.records)?;
        }

        info!(
            "BySJout: filtered {} pairs with non-surviving junctions",
            filtered_count
        );
    }

    // Flush chimeric output if enabled (currently no chimeric alignments from paired-end)
    if let Some(ref mut chim_writer) = chimeric_writer {
        chim_writer.flush()?;
    }

    Ok(())
}

/// Record junctions from a transcript into SJ statistics
fn record_transcript_junctions(
    transcript: &crate::align::transcript::Transcript,
    index: &crate::index::GenomeIndex,
    sj_stats: &crate::junction::SpliceJunctionStats,
    is_unique: bool,
) {
    use crate::align::score::AlignmentScorer;
    use crate::align::transcript::CigarOp;

    // First pass: compute exon segment lengths (query-consuming bases between N operations)
    // An "exon segment" is the query bases on each side of a splice junction.
    let mut exon_lengths: Vec<u32> = Vec::new();
    let mut current_exon_len = 0u32;

    for op in &transcript.cigar {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                current_exon_len += *len;
            }
            CigarOp::Ins(len) => {
                current_exon_len += *len;
            }
            CigarOp::SoftClip(len) => {
                current_exon_len += *len;
            }
            CigarOp::RefSkip(_) => {
                exon_lengths.push(current_exon_len);
                current_exon_len = 0;
            }
            CigarOp::Del(_) | CigarOp::HardClip(_) => {
                // Do not consume query bases
            }
        }
    }
    exon_lengths.push(current_exon_len); // Final exon segment

    // Second pass: record junctions with computed overhangs
    let mut genome_pos = transcript.genome_start;
    let mut junction_idx = 0usize;

    let scorer = AlignmentScorer::from_params_minimal();

    for op in &transcript.cigar {
        match op {
            CigarOp::RefSkip(len) => {
                // This is a splice junction
                let intron_len = *len;
                let intron_start = genome_pos + 1; // 1-based, first intronic base
                let intron_end = genome_pos + intron_len as u64; // 1-based, last intronic base

                // Detect splice motif
                let motif = scorer.detect_splice_motif(genome_pos, intron_len, &index.genome);

                // Compute overhang: min(left_exon_length, right_exon_length)
                let left_exon = exon_lengths[junction_idx];
                let right_exon = exon_lengths[junction_idx + 1];
                let overhang = left_exon.min(right_exon);

                // Derive strand from splice motif (STAR convention)
                let strand = match motif.implied_strand() {
                    Some('+') => 1u8,
                    Some('-') => 2u8,
                    _ => 0u8, // non-canonical: unknown strand
                };
                let annotated = index.junction_db.is_annotated(
                    transcript.chr_idx,
                    intron_start,
                    intron_end,
                    strand,
                );

                // Record junction
                sj_stats.record_junction(
                    transcript.chr_idx,
                    intron_start,
                    intron_end,
                    strand,
                    motif,
                    is_unique,
                    overhang,
                    annotated,
                );

                // Advance genome position past the intron
                genome_pos += intron_len as u64;
                junction_idx += 1;
            }
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::Ins(_) => {}
            CigarOp::Del(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::SoftClip(_) | CigarOp::HardClip(_) => {}
        }
    }
}
