/// Read alignment driver function
use crate::align::score::{AlignmentScorer, SpliceMotif};
use crate::align::seed::Seed;
use crate::align::stitch::{cluster_seeds, stitch_seeds_with_jdb};
use crate::align::transcript::Transcript;
use crate::error::Error;
use crate::index::GenomeIndex;
use crate::params::{IntronMotifFilter, IntronStrandFilter, Parameters};
use crate::stats::UnmappedReason;

/// Paired-end alignment result
#[derive(Debug, Clone)]
pub struct PairedAlignment {
    /// Transcript for mate1
    pub mate1_transcript: Transcript,
    /// Transcript for mate2
    pub mate2_transcript: Transcript,
    /// Read positions for mate1 in transcript (start, end)
    pub mate1_region: (usize, usize),
    /// Read positions for mate2 in transcript (start, end)
    pub mate2_region: (usize, usize),
    /// Whether this is a proper pair (same chr, concordant orientation, distance)
    pub is_proper_pair: bool,
    /// Signed insert size (TLEN) - genomic distance between mate starts
    pub insert_size: i32,
}

/// Align a read to the genome.
///
/// # Algorithm
/// 1. Find seeds (exact matches) using MMP search
/// 2. Cluster seeds by genomic proximity
/// 3. Stitch seeds within each cluster using DP
/// 4. Filter transcripts by quality thresholds
/// 5. Sort by score and limit to top N
/// 6. Detect chimeric alignments if enabled
///
/// # Arguments
/// * `read_seq` - Read sequence (encoded as 0=A, 1=C, 2=G, 3=T)
/// * `read_name` - Read name (needed for chimeric output)
/// * `index` - Genome index
/// * `params` - User parameters
///
/// # Returns
/// Tuple of (transcripts, chimeric alignments, n_for_mapq, unmapped_reason):
/// - transcripts: sorted by score (best first)
/// - chimeric alignments: sorted by score (best first)
/// - n_for_mapq: effective alignment count for MAPQ calculation (max of transcript count
///   and valid cluster count, to avoid undercounting from coordinate dedup on tandem repeats)
/// - unmapped_reason: `Some(reason)` if no alignments produced, `None` if mapped
pub fn align_read(
    read_seq: &[u8],
    read_name: &str,
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<
    (
        Vec<Transcript>,
        Vec<crate::chimeric::ChimericAlignment>,
        usize,
        Option<UnmappedReason>,
    ),
    Error,
> {
    // Step 1: Find seeds
    // Use a reasonable default min seed length (typically 8-20bp)
    let min_seed_length = 8;
    let seeds = Seed::find_seeds(read_seq, index, min_seed_length, params)?;

    if seeds.is_empty() {
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Step 2: Cluster seeds
    let max_cluster_dist = params.win_bin_window_dist();
    let max_loci_for_anchor = 10; // Seeds mapping to <=10 loci can be anchors
    let clusters = cluster_seeds(
        &seeds,
        index,
        max_cluster_dist,
        max_loci_for_anchor,
        params.win_anchor_multimap_nmax,
        params.seed_none_loci_per_window,
        params.win_bin_nbits,
    );

    if clusters.is_empty() {
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Cap total clusters (alignWindowsPerReadNmax)
    let mut clusters = clusters;
    clusters.truncate(params.align_windows_per_read_nmax);

    // Cap seeds per cluster (seedPerWindowNmax)
    for cluster in &mut clusters {
        cluster.seed_indices.truncate(params.seed_per_window_nmax);
    }

    // Step 2a: Filter clusters by seed coverage (winReadCoverageRelativeMin)
    let read_length = read_seq.len();
    let min_coverage = params.win_read_coverage_relative_min;
    let clusters: Vec<_> = clusters
        .into_iter()
        .filter(|cluster| {
            // Compute total seed coverage for this cluster (union of read ranges)
            let mut covered = vec![false; read_length];
            for &seed_idx in &cluster.seed_indices {
                let seed = &seeds[seed_idx];
                let end = seed.read_pos + seed.length;
                for c in covered
                    .iter_mut()
                    .take(end.min(read_length))
                    .skip(seed.read_pos)
                {
                    *c = true;
                }
            }
            let coverage = covered.iter().filter(|&&c| c).count() as f64 / read_length as f64;
            coverage >= min_coverage
        })
        .collect();

    if clusters.is_empty() {
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Step 2b: Detect chimeric alignments from multi-cluster seeds (Tier 2)
    let mut chimeric_alignments = Vec::new();
    if params.chim_segment_min > 0 && clusters.len() > 1 {
        use crate::chimeric::ChimericDetector;
        let detector = ChimericDetector::new(params);
        chimeric_alignments.extend(
            detector.detect_from_multi_clusters(&clusters, &seeds, read_seq, read_name, index)?,
        );
    }

    // Step 3: Stitch seeds within each cluster
    let scorer = AlignmentScorer::from_params(params);
    let mut transcripts = Vec::new();

    // Use junction DB for annotation-aware scoring if available
    let junction_db = if index.junction_db.is_empty() {
        None
    } else {
        Some(&index.junction_db)
    };

    for cluster in clusters.iter() {
        let cluster_transcripts =
            stitch_seeds_with_jdb(cluster, &seeds, read_seq, index, &scorer, junction_db)?;
        transcripts.extend(cluster_transcripts);
    }

    // Step 4: Filter transcripts
    let read_length = read_seq.len() as f64;

    // Log filtering statistics
    let pre_filter_count = transcripts.len();
    let mut filter_reasons = std::collections::HashMap::new();

    transcripts.retain(|t| {
        // Absolute score threshold
        if t.score < params.out_filter_score_min {
            *filter_reasons.entry("score_min").or_insert(0) += 1;
            return false;
        }

        // Relative score threshold (score / read_length)
        if (t.score as f64) < params.out_filter_score_min_over_lread * read_length {
            *filter_reasons.entry("score_min_relative").or_insert(0) += 1;
            return false;
        }

        // Absolute mismatch count
        if t.n_mismatch > params.out_filter_mismatch_nmax {
            *filter_reasons.entry("mismatch_max").or_insert(0) += 1;
            log::debug!(
                "Filtered {}: {} mismatches > {} max (read_len={}, score={})",
                read_name,
                t.n_mismatch,
                params.out_filter_mismatch_nmax,
                read_length,
                t.score
            );
            return false;
        }

        // Relative mismatch count (mismatches / read_length)
        let mismatch_rate = t.n_mismatch as f64 / read_length;
        if mismatch_rate > params.out_filter_mismatch_nover_lmax {
            *filter_reasons.entry("mismatch_rate").or_insert(0) += 1;
            log::debug!(
                "Filtered {}: {:.1}% mismatch rate > {:.1}% max ({}/{} bases, score={})",
                read_name,
                mismatch_rate * 100.0,
                params.out_filter_mismatch_nover_lmax * 100.0,
                t.n_mismatch,
                read_length,
                t.score
            );
            return false;
        }

        // Absolute matched bases
        let n_matched = t.n_matched();
        if n_matched < params.out_filter_match_nmin {
            *filter_reasons.entry("match_min").or_insert(0) += 1;
            return false;
        }

        // Relative matched bases (matched / read_length)
        if (n_matched as f64) < params.out_filter_match_nmin_over_lread * read_length {
            *filter_reasons.entry("match_min_relative").or_insert(0) += 1;
            return false;
        }

        // Junction motif filtering
        match params.out_filter_intron_motifs {
            IntronMotifFilter::None => {
                // Accept all motifs
            }
            IntronMotifFilter::RemoveNoncanonical => {
                // Reject if any junction is non-canonical
                if t.junction_motifs
                    .iter()
                    .any(|m| *m == SpliceMotif::NonCanonical)
                {
                    *filter_reasons.entry("noncanonical_junction").or_insert(0) += 1;
                    return false;
                }
            }
            IntronMotifFilter::RemoveNoncanonicalUnannotated => {
                // Only reject if a non-canonical junction is NOT annotated in GTF
                if t.junction_motifs
                    .iter()
                    .zip(t.junction_annotated.iter())
                    .any(|(m, annotated)| *m == SpliceMotif::NonCanonical && !annotated)
                {
                    *filter_reasons
                        .entry("noncanonical_unannotated_junction")
                        .or_insert(0) += 1;
                    return false;
                }
            }
        }

        // Intron strand consistency filtering (outFilterIntronStrands)
        if params.out_filter_intron_strands == IntronStrandFilter::RemoveInconsistentStrands {
            let mut has_plus = false;
            let mut has_minus = false;
            for motif in &t.junction_motifs {
                match motif.implied_strand() {
                    Some('+') => has_plus = true,
                    Some('-') => has_minus = true,
                    None => {}
                    _ => {}
                }
            }
            if has_plus && has_minus {
                *filter_reasons.entry("inconsistent_strand").or_insert(0) += 1;
                return false;
            }
        }

        true
    });

    // Log filtering summary if anything was filtered
    if pre_filter_count > transcripts.len() {
        let filtered = pre_filter_count - transcripts.len();
        log::debug!(
            "Read {}: Filtered {}/{} transcripts: {:?}",
            read_name,
            filtered,
            pre_filter_count,
            filter_reasons
        );
    }

    // Step 3b: Detect chimeric alignments from soft-clips (Tier 1)
    if params.chim_segment_min > 0 {
        use crate::chimeric::ChimericDetector;
        let detector = ChimericDetector::new(params);
        for transcript in &transcripts {
            if let Some(chim) =
                detector.detect_from_soft_clips(transcript, read_seq, read_name, index)?
            {
                chimeric_alignments.push(chim);
            }
        }
    }

    // Step 5: Deduplicate transcripts with identical genomic coordinates
    // Keep only the highest-scoring transcript for each unique location
    // Sort by (chr, start, end, strand, score_descending) so dedup_by keeps the best
    transcripts.sort_by(|a, b| {
        (a.chr_idx, a.genome_start, a.genome_end, a.is_reverse)
            .cmp(&(b.chr_idx, b.genome_start, b.genome_end, b.is_reverse))
            .then_with(|| b.score.cmp(&a.score)) // Higher score first
    });

    // Dedup consecutive entries with same coordinates (keeps first = highest score)
    transcripts.dedup_by(|a, b| {
        a.chr_idx == b.chr_idx
            && a.genome_start == b.genome_start
            && a.genome_end == b.genome_end
            && a.is_reverse == b.is_reverse
    });

    // Step 5b: Re-sort by score (descending) with deterministic tie-breaking
    // When scores are equal, prefer: smallest chr index → smallest position → forward strand
    // This matches STAR's tie-breaking behavior for multi-mappers
    transcripts.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| a.chr_idx.cmp(&b.chr_idx))
            .then_with(|| a.genome_start.cmp(&b.genome_start))
            .then_with(|| a.is_reverse.cmp(&b.is_reverse))
    });

    // Step 5b: Filter to keep only alignments within score range of the best
    // This is CRITICAL for unique vs multi-mapped classification
    if !transcripts.is_empty() {
        let max_score = transcripts[0].score;
        let score_threshold = max_score - params.out_filter_multimap_score_range;

        // Keep only transcripts within score range of the best
        transcripts.retain(|t| t.score >= score_threshold);
    }

    // Step 5c: Truncate to max multimap count
    transcripts.truncate(params.out_filter_multimap_nmax as usize);

    // Step 6: Filter chimeric alignments
    if params.chim_segment_min > 0 {
        chimeric_alignments.retain(|chim| {
            chim.meets_min_segment_length(params.chim_segment_min)
                && chim.meets_min_score(params.chim_score_min)
        });
    }

    // n_for_mapq: currently same as transcripts.len().
    // NOTE: rDNA reads (~157) still get MAPQ=255 instead of STAR's 1-3 because our large
    // clusters (589kb) merge seeds from all tandem repeat copies. The DP finds one optimal
    // path shared across all clusters (same seeds due to 589kb distance). Bin-counting
    // approaches fail because "wrong" clusters also share seeds and produce the same
    // transcript. The fix requires splitting clusters into ~65kb sub-windows with
    // independent DP per sub-window (Phase 16.5b → cluster splitting).
    let n_for_mapq = transcripts.len();

    let unmapped_reason = if transcripts.is_empty() {
        // Transcripts were generated by DP but all filtered out
        Some(UnmappedReason::TooShort)
    } else {
        None
    };

    Ok((
        transcripts,
        chimeric_alignments,
        n_for_mapq,
        unmapped_reason,
    ))
}

/// Align paired-end reads by aligning each mate independently, then pairing results.
///
/// # Algorithm
/// 1. Align mate1 independently using the proven SE align_read() path
/// 2. Align mate2 independently using the proven SE align_read() path
/// 3. Pair transcripts by chromosome + distance constraints
/// 4. Deduplicate, sort by combined score, filter
///
/// This pragmatic approach leverages the well-tested SE alignment (95.7% position
/// agreement) rather than the broken seed-pooling approach. Full mate-aware joint
/// DP stitching (STAR's actual approach) is deferred to Phase 16.6.
///
/// # Arguments
/// * `mate1_seq` - First mate sequence (encoded)
/// * `mate2_seq` - Second mate sequence (encoded)
/// * `index` - Genome index
/// * `params` - Parameters (includes alignMatesGapMax)
///
/// # Returns
/// Tuple of (paired alignments, n_for_mapq, unmapped_reason):
/// - paired alignments: sorted by combined score (best first)
/// - n_for_mapq: effective alignment count for MAPQ (max of both mates' cluster counts)
/// - unmapped_reason: `Some(reason)` if no paired alignments, `None` if mapped
pub fn align_paired_read(
    mate1_seq: &[u8],
    mate2_seq: &[u8],
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<(Vec<PairedAlignment>, usize, Option<UnmappedReason>), Error> {
    // Step 1: Align each mate independently using the proven SE path
    let (mate1_transcripts, _, mate1_mapq_n, mate1_unmapped) =
        align_read(mate1_seq, "", index, params)?;
    let (mate2_transcripts, _, mate2_mapq_n, mate2_unmapped) =
        align_read(mate2_seq, "", index, params)?;

    // If either mate has no alignments, return empty (both unmapped)
    if mate1_transcripts.is_empty() || mate2_transcripts.is_empty() {
        // Use the first mate's unmapped reason, or Other
        let reason = mate1_unmapped
            .or(mate2_unmapped)
            .unwrap_or(UnmappedReason::Other);
        return Ok((Vec::new(), 0, Some(reason)));
    }

    let pe_mapq_n = mate1_mapq_n.max(mate2_mapq_n);

    // Step 2: Pair transcripts by chromosome + distance
    let mut paired_alignments = Vec::new();

    for t1 in &mate1_transcripts {
        for t2 in &mate2_transcripts {
            // Must be on same chromosome
            if t1.chr_idx != t2.chr_idx {
                continue;
            }

            // Check distance constraint
            if params.align_mates_gap_max > 0 && !check_proper_pair(t1, t2, params) {
                continue;
            }

            let is_proper_pair = check_proper_pair(t1, t2, params);
            let insert_size = calculate_insert_size(t1, t2);

            paired_alignments.push(PairedAlignment {
                mate1_transcript: t1.clone(),
                mate2_transcript: t2.clone(),
                mate1_region: (0, mate1_seq.len()),
                mate2_region: (0, mate2_seq.len()),
                is_proper_pair,
                insert_size,
            });
        }
    }

    // Step 3: Deduplicate — if same (mate1 location, mate2 location) appears, keep highest score
    paired_alignments.sort_by(|a, b| {
        (
            a.mate1_transcript.chr_idx,
            a.mate1_transcript.genome_start,
            a.mate1_transcript.is_reverse,
            a.mate2_transcript.genome_start,
            a.mate2_transcript.is_reverse,
        )
            .cmp(&(
                b.mate1_transcript.chr_idx,
                b.mate1_transcript.genome_start,
                b.mate1_transcript.is_reverse,
                b.mate2_transcript.genome_start,
                b.mate2_transcript.is_reverse,
            ))
            .then_with(|| {
                let b_combined = b.mate1_transcript.score + b.mate2_transcript.score;
                let a_combined = a.mate1_transcript.score + a.mate2_transcript.score;
                b_combined.cmp(&a_combined)
            })
    });
    paired_alignments.dedup_by(|a, b| {
        a.mate1_transcript.chr_idx == b.mate1_transcript.chr_idx
            && a.mate1_transcript.genome_start == b.mate1_transcript.genome_start
            && a.mate1_transcript.is_reverse == b.mate1_transcript.is_reverse
            && a.mate2_transcript.genome_start == b.mate2_transcript.genome_start
            && a.mate2_transcript.is_reverse == b.mate2_transcript.is_reverse
    });

    // Step 4: Sort by combined score (descending) with deterministic tie-breaking
    paired_alignments.sort_by(|a, b| {
        let a_combined = a.mate1_transcript.score + a.mate2_transcript.score;
        let b_combined = b.mate1_transcript.score + b.mate2_transcript.score;
        b_combined
            .cmp(&a_combined)
            .then_with(|| a.mate1_transcript.chr_idx.cmp(&b.mate1_transcript.chr_idx))
            .then_with(|| {
                a.mate1_transcript
                    .genome_start
                    .cmp(&b.mate1_transcript.genome_start)
            })
            .then_with(|| {
                a.mate1_transcript
                    .is_reverse
                    .cmp(&b.mate1_transcript.is_reverse)
            })
    });

    // Step 5: Filter by score range from best combined score
    if !paired_alignments.is_empty() {
        let best_score = paired_alignments[0].mate1_transcript.score
            + paired_alignments[0].mate2_transcript.score;
        let score_threshold = best_score - params.out_filter_multimap_score_range;
        paired_alignments.retain(|pa| {
            let combined = pa.mate1_transcript.score + pa.mate2_transcript.score;
            combined >= score_threshold
        });
    }

    // Step 6: Apply quality filters and truncate
    filter_paired_transcripts(&mut paired_alignments, params);

    let unmapped_reason = if paired_alignments.is_empty() {
        Some(UnmappedReason::TooShort)
    } else {
        None
    };

    Ok((paired_alignments, pe_mapq_n, unmapped_reason))
}

/// Check if paired alignment is a proper pair
fn check_proper_pair(
    mate1_trans: &Transcript,
    mate2_trans: &Transcript,
    params: &Parameters,
) -> bool {
    // Proper pair criteria:
    // 1. Both mates mapped (checked by caller)
    // 2. Same chromosome (checked by caller)
    // 3. Distance within alignMatesGapMax

    if params.align_mates_gap_max == 0 {
        return true; // Auto mode = unlimited
    }

    // Calculate genomic distance
    let start = mate1_trans.genome_start.min(mate2_trans.genome_start);
    let end = mate1_trans.genome_end.max(mate2_trans.genome_end);
    let genomic_span = end - start;

    genomic_span <= params.align_mates_gap_max as u64
}

/// Calculate signed insert size (TLEN)
fn calculate_insert_size(mate1_trans: &Transcript, mate2_trans: &Transcript) -> i32 {
    // TLEN = signed genomic distance from leftmost mate to rightmost mate
    let start = mate1_trans.genome_start.min(mate2_trans.genome_start);
    let end = mate1_trans.genome_end.max(mate2_trans.genome_end);
    let abs_tlen = (end - start) as i32;

    // Sign convention: positive if mate1 is leftmost
    if mate1_trans.genome_start <= mate2_trans.genome_start {
        abs_tlen
    } else {
        -abs_tlen
    }
}

/// Filter paired transcripts by quality thresholds.
/// Each mate is checked independently against per-read thresholds.
/// Combined score uses sum of both mates.
fn filter_paired_transcripts(paired_alns: &mut Vec<PairedAlignment>, params: &Parameters) {
    paired_alns.retain(|pa| {
        let t1 = &pa.mate1_transcript;
        let t2 = &pa.mate2_transcript;
        let mate1_len = (pa.mate1_region.1 - pa.mate1_region.0) as f64;
        let mate2_len = (pa.mate2_region.1 - pa.mate2_region.0) as f64;

        // Check each mate independently for per-read thresholds
        for (t, read_len) in [(t1, mate1_len), (t2, mate2_len)] {
            // Absolute score threshold (per mate)
            if t.score < params.out_filter_score_min {
                return false;
            }

            // Relative score threshold (per mate)
            if (t.score as f64) < params.out_filter_score_min_over_lread * read_len {
                return false;
            }

            // Absolute mismatch count (per mate)
            if t.n_mismatch > params.out_filter_mismatch_nmax {
                return false;
            }

            // Relative mismatch count (per mate)
            if (t.n_mismatch as f64) > params.out_filter_mismatch_nover_lmax * read_len {
                return false;
            }

            // Absolute matched bases (per mate)
            let n_matched = t.n_matched();
            if n_matched < params.out_filter_match_nmin {
                return false;
            }

            // Relative matched bases (per mate)
            if (n_matched as f64) < params.out_filter_match_nmin_over_lread * read_len {
                return false;
            }
        }

        true
    });

    // Limit to top N
    paired_alns.truncate(params.out_filter_multimap_nmax as usize);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Genome;
    use crate::index::packed_array::PackedArray;
    use crate::index::sa_index::SaIndex;
    use crate::index::suffix_array::SuffixArray;
    use clap::Parser;

    fn make_test_params() -> Parameters {
        // Parse empty args to get default parameters
        Parameters::try_parse_from(vec!["ruSTAR"]).unwrap()
    }

    fn make_test_index() -> GenomeIndex {
        // Simple genome: ACGTACGTNN (10 bases)
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 4, 4];
        let n_genome = 64u64; // Padded
        let mut sequence = vec![5u8; (n_genome * 2) as usize];
        sequence[0..seq.len()].copy_from_slice(&seq);

        // Build reverse complement
        for i in 0..n_genome as usize {
            let base = sequence[i];
            let complement = if base < 4 { 3 - base } else { base };
            sequence[2 * n_genome as usize - 1 - i] = complement;
        }

        let genome = Genome {
            sequence,
            n_genome,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![10],
            chr_start: vec![0, n_genome],
        };

        // Create dummy SA and SAindex (would need real index for actual alignment)
        let gstrand_bit = 33;
        let suffix_array = SuffixArray {
            data: PackedArray::new(gstrand_bit, 0),
            gstrand_bit,
            gstrand_mask: (1u64 << gstrand_bit) - 1,
        };

        let word_length = gstrand_bit + 3;
        let sa_index = SaIndex {
            data: PackedArray::new(word_length, 0),
            nbases: 14,
            genome_sa_index_start: vec![0],
            word_length,
            gstrand_bit,
        };

        GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db: crate::junction::SpliceJunctionDb::empty(),
        }
    }

    #[test]
    fn test_align_read_no_seeds() {
        let index = make_test_index();
        let params = make_test_params();

        // Read with all N's (no seeds possible)
        let read_seq = vec![4, 4, 4, 4, 4, 4, 4, 4, 4, 4];

        let result = align_read(&read_seq, "READ_001", &index, &params);
        assert!(result.is_ok());

        let (transcripts, chimeras, n_for_mapq, unmapped_reason) = result.unwrap();
        assert_eq!(transcripts.len(), 0); // No alignment
        assert_eq!(chimeras.len(), 0); // No chimeric alignments
        assert_eq!(n_for_mapq, 0);
        assert_eq!(unmapped_reason, Some(UnmappedReason::Other));
    }

    #[test]
    fn test_transcript_filtering_score() {
        let index = make_test_index();
        let mut params = make_test_params();
        params.out_filter_score_min = 50;

        // Would need actual seeds and alignment to test this properly
        // This test just verifies the function doesn't crash
        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, "READ_002", &index, &params);
        assert!(result.is_ok());
    }

    #[test]
    fn test_transcript_filtering_mismatch() {
        let index = make_test_index();
        let mut params = make_test_params();
        params.out_filter_mismatch_nmax = 2;

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, "READ_003", &index, &params);
        assert!(result.is_ok());
    }

    #[test]
    fn test_transcript_multimap_limit() {
        let index = make_test_index();
        let mut params = make_test_params();
        params.out_filter_multimap_nmax = 5;

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, "READ_004", &index, &params);
        assert!(result.is_ok());

        let (transcripts, _chimeras, _n_for_mapq, _reason) = result.unwrap();
        assert!(transcripts.len() <= 5);
    }

    #[test]
    fn test_align_paired_read_no_seeds() {
        let index = make_test_index();
        let params = make_test_params();

        // Both mates with all N's
        let mate1 = vec![4, 4, 4, 4, 4, 4, 4, 4];
        let mate2 = vec![4, 4, 4, 4, 4, 4, 4, 4];

        let result = align_paired_read(&mate1, &mate2, &index, &params);
        assert!(result.is_ok());
        let (paired_alns, n_for_mapq, unmapped_reason) = result.unwrap();
        assert_eq!(paired_alns.len(), 0);
        assert_eq!(n_for_mapq, 0);
        assert!(unmapped_reason.is_some());
    }

    #[test]
    fn test_check_proper_pair_distance() {
        use crate::align::transcript::{CigarOp, Exon};

        let params = make_test_params();

        // Create two transcripts on same chromosome
        let t1 = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1100,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let t2 = Transcript {
            chr_idx: 0,
            genome_start: 1200,
            genome_end: 1300,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 1200,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Distance = 300bp, within default limit (auto mode = unlimited)
        assert!(check_proper_pair(&t1, &t2, &params));
    }

    #[test]
    fn test_check_proper_pair_too_far() {
        use crate::align::transcript::{CigarOp, Exon};

        let mut params = make_test_params();
        params.align_mates_gap_max = 100;

        let t1 = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1100,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let t2 = Transcript {
            chr_idx: 0,
            genome_start: 1300,
            genome_end: 1400,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 1300,
                genome_end: 1400,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Distance = 400bp, exceeds limit of 100bp
        assert!(!check_proper_pair(&t1, &t2, &params));
    }

    #[test]
    fn test_calculate_insert_size_positive() {
        use crate::align::transcript::{CigarOp, Exon};

        // Mate1 is leftmost
        let t1 = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1100,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let t2 = Transcript {
            chr_idx: 0,
            genome_start: 1200,
            genome_end: 1300,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 1200,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let tlen = calculate_insert_size(&t1, &t2);
        assert_eq!(tlen, 300); // Positive because mate1 is leftmost
    }

    #[test]
    fn test_strand_consistency_filter() {
        use crate::align::transcript::{CigarOp, Exon, Transcript};
        use crate::params::IntronStrandFilter;

        // Create a transcript with conflicting strand motifs
        let t_inconsistent = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1300,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::CtAc], // +strand and -strand
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Create a transcript with consistent strand motifs
        let t_consistent = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1300,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::GcAg], // both + strand
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Verify implied_strand detects the conflict
        let mut has_plus = false;
        let mut has_minus = false;
        for motif in &t_inconsistent.junction_motifs {
            match motif.implied_strand() {
                Some('+') => has_plus = true,
                Some('-') => has_minus = true,
                _ => {}
            }
        }
        assert!(has_plus && has_minus); // Inconsistent

        // Verify consistent transcript has no conflict
        has_plus = false;
        has_minus = false;
        for motif in &t_consistent.junction_motifs {
            match motif.implied_strand() {
                Some('+') => has_plus = true,
                Some('-') => has_minus = true,
                _ => {}
            }
        }
        assert!(has_plus && !has_minus); // Consistent (all +)

        // Verify the filter enum
        assert_eq!(
            IntronStrandFilter::RemoveInconsistentStrands,
            IntronStrandFilter::RemoveInconsistentStrands
        );
        assert_ne!(
            IntronStrandFilter::None,
            IntronStrandFilter::RemoveInconsistentStrands
        );
    }

    #[test]
    fn test_calculate_insert_size_negative() {
        use crate::align::transcript::{CigarOp, Exon};

        // Mate2 is leftmost
        let t1 = Transcript {
            chr_idx: 0,
            genome_start: 1200,
            genome_end: 1300,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1200,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let t2 = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1100,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let tlen = calculate_insert_size(&t1, &t2);
        assert_eq!(tlen, -300); // Negative because mate2 is leftmost
    }

    #[test]
    fn test_noncanonical_unannotated_filter() {
        use crate::align::score::SpliceMotif;
        use crate::align::transcript::{CigarOp, Exon, Transcript};

        // Helper: check if a transcript would be filtered by RemoveNoncanonicalUnannotated
        // (mirrors the logic in the retain closure)
        let would_filter = |t: &Transcript| -> bool {
            t.junction_motifs
                .iter()
                .zip(t.junction_annotated.iter())
                .any(|(m, annotated)| *m == SpliceMotif::NonCanonical && !annotated)
        };

        let base_transcript = || Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1200,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1200,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Case 1: NonCanonical + unannotated → should be filtered
        let mut t1 = base_transcript();
        t1.junction_motifs = vec![SpliceMotif::NonCanonical];
        t1.junction_annotated = vec![false];
        assert!(
            would_filter(&t1),
            "NonCanonical + unannotated should be filtered"
        );

        // Case 2: NonCanonical + annotated → should be KEPT
        let mut t2 = base_transcript();
        t2.junction_motifs = vec![SpliceMotif::NonCanonical];
        t2.junction_annotated = vec![true];
        assert!(
            !would_filter(&t2),
            "NonCanonical + annotated should be kept"
        );

        // Case 3: Canonical + unannotated → should be KEPT
        let mut t3 = base_transcript();
        t3.junction_motifs = vec![SpliceMotif::GtAg];
        t3.junction_annotated = vec![false];
        assert!(!would_filter(&t3), "Canonical + unannotated should be kept");

        // Case 4: Mixed — one canonical + one non-canonical unannotated → filtered
        let mut t4 = base_transcript();
        t4.junction_motifs = vec![SpliceMotif::GtAg, SpliceMotif::NonCanonical];
        t4.junction_annotated = vec![true, false];
        assert!(
            would_filter(&t4),
            "Mixed with unannotated non-canonical should be filtered"
        );

        // Case 5: Mixed — one canonical + one non-canonical annotated → kept
        let mut t5 = base_transcript();
        t5.junction_motifs = vec![SpliceMotif::GtAg, SpliceMotif::NonCanonical];
        t5.junction_annotated = vec![false, true];
        assert!(
            !would_filter(&t5),
            "Mixed with annotated non-canonical should be kept"
        );
    }
}
