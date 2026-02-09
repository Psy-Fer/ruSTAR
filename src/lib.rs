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

    // 2. Dispatch based on two-pass mode
    match params.twopass_mode {
        TwopassMode::None => {
            info!("Running single-pass alignment");
            run_single_pass(&index, params)?;
        }
        TwopassMode::Basic => {
            info!("Running two-pass alignment mode");
            run_two_pass(&index, params)?;
        }
    }

    info!("Alignment complete!");
    Ok(())
}

/// Run single-pass alignment (original logic)
fn run_single_pass(
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    params: &Parameters,
) -> anyhow::Result<()> {
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
        sj_stats.write_output(&sj_output_path, &index.genome)?;
    }

    // 6. Print summary
    stats.print_summary();

    Ok(())
}

/// Run two-pass alignment mode
fn run_two_pass(
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    params: &Parameters,
) -> anyhow::Result<()> {
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
    sj_stats_pass1.write_output(&pass1_path, &index.genome)?;
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
    run_single_pass(&Arc::new(merged_index), params)?;

    Ok(())
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
                    return Ok(AlignmentBatchResults {
                        sam_records: buffer,
                        chimeric_alns,
                    });
                }

                // Align read (CPU-intensive, pure function)
                let (transcripts, chimeric_results) =
                    align_read(&clipped_seq, &read.name, &index, params)?;

                // Collect chimeric alignments if enabled
                if params.chim_segment_min > 0 {
                    chimeric_alns.extend(chimeric_results);
                }

                // Record stats (atomic, lock-free)
                stats.record_alignment(transcripts.len(), max_multimaps);

                // Record junction statistics
                let is_unique = transcripts.len() == 1;
                for transcript in &transcripts {
                    record_transcript_junctions(transcript, &index, &sj_stats, is_unique);
                }

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

                Ok(AlignmentBatchResults {
                    sam_records: buffer,
                    chimeric_alns,
                })
            })
            .collect();

        // Sequential writing (merge buffers in chunk order)
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

        read_count += reads_to_process as u64;

        // Progress logging
        if read_count % 100000 < batch_size as u64 {
            info!("Processed {} reads...", read_count);
        }

        if read_count >= max_reads {
            break;
        }
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

                // Record junction statistics (unified transcript covers both mates)
                let is_unique = paired_alns.len() == 1;
                for pair in &paired_alns {
                    record_transcript_junctions(&pair.transcript, &index, &sj_stats, is_unique);
                }

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

    // Track genome position as we traverse CIGAR
    let mut genome_pos = transcript.genome_start;
    let mut _read_pos = 0usize;

    for op in &transcript.cigar {
        match op {
            CigarOp::RefSkip(len) => {
                // This is a splice junction
                let intron_len = *len;
                let intron_start = genome_pos + 1; // 1-based, first intronic base
                let intron_end = genome_pos + intron_len as u64; // 1-based, last intronic base

                // Detect splice motif
                let scorer = AlignmentScorer {
                    score_gap: 0,
                    score_gap_noncan: -8,
                    score_gap_gcag: -4,
                    score_gap_atac: -8,
                    score_del_open: -2,
                    score_del_base: -2,
                    score_ins_open: -2,
                    score_ins_base: -2,
                    align_intron_min: 21,
                    sjdb_score: 2,
                    align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
                    n_mm_max: 10,
                    p_mm_max: 0.3,
                    align_sj_overhang_min: 5,
                    align_sjdb_overhang_min: 3,
                };
                let motif = scorer.detect_splice_motif(genome_pos, intron_len, &index.genome);

                // Calculate overhang (simplified: use exon lengths)
                // TODO: More accurate overhang calculation from read/exon boundaries
                let overhang = 5u32; // Placeholder - actual overhang needs exon boundary calculation

                // Check if annotated
                let strand = if transcript.is_reverse { 2 } else { 1 };
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
            }
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                genome_pos += *len as u64;
                _read_pos += *len as usize;
            }
            CigarOp::Ins(len) => {
                _read_pos += *len as usize;
            }
            CigarOp::Del(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::SoftClip(len) | CigarOp::HardClip(len) => {
                _read_pos += *len as usize;
            }
        }
    }
}
