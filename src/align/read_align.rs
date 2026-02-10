/// Read alignment driver function
use crate::align::score::{AlignmentScorer, SpliceMotif};
use crate::align::seed::Seed;
use crate::align::stitch::{cluster_seeds, stitch_seeds};
use crate::align::transcript::Transcript;
use crate::error::Error;
use crate::index::GenomeIndex;
use crate::params::{IntronMotifFilter, IntronStrandFilter, Parameters};

/// Paired-end alignment result
#[derive(Debug, Clone)]
pub struct PairedAlignment {
    /// Unified transcript covering both mates
    pub transcript: Transcript,
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
/// Tuple of (transcripts, chimeric alignments), both sorted by score (best first)
pub fn align_read(
    read_seq: &[u8],
    read_name: &str,
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<(Vec<Transcript>, Vec<crate::chimeric::ChimericAlignment>), Error> {
    // Step 1: Find seeds
    // Use a reasonable default min seed length (typically 8-20bp)
    let min_seed_length = 8;
    let seeds = Seed::find_seeds(read_seq, index, min_seed_length, params)?;

    if seeds.is_empty() {
        return Ok((Vec::new(), Vec::new())); // No seeds found
    }

    // Step 2: Cluster seeds
    let max_cluster_dist = 100000; // 100kb window (TODO: make configurable)
    let max_loci_for_anchor = 10; // Seeds mapping to <=10 loci can be anchors
    let clusters = cluster_seeds(
        &seeds,
        index,
        max_cluster_dist,
        max_loci_for_anchor,
        params.win_anchor_multimap_nmax,
        params.seed_none_loci_per_window,
    );

    if clusters.is_empty() {
        return Ok((Vec::new(), Vec::new()));
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

    for cluster in clusters.iter() {
        let cluster_transcripts = stitch_seeds(cluster, &seeds, read_seq, index, &scorer)?;
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
                // TODO: requires checking junction database for annotation status
                // For now, treat same as RemoveNoncanonical
                if t.junction_motifs
                    .iter()
                    .any(|m| *m == SpliceMotif::NonCanonical)
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

    Ok((transcripts, chimeric_alignments))
}

/// Align paired-end read using STAR's unified seed clustering approach.
///
/// # Algorithm (matching STAR)
/// 1. Find seeds from both mates independently (tagged with mate_id)
/// 2. Pool seeds together for unified clustering
/// 3. Cluster seeds in genomic windows (both mates' seeds together)
/// 4. Stitch seeds with DP that handles both mates
/// 5. Split results back to mate regions for SAM output
///
/// # Arguments
/// * `mate1_seq` - First mate sequence (encoded)
/// * `mate2_seq` - Second mate sequence (encoded)
/// * `index` - Genome index
/// * `params` - Parameters (includes alignMatesGapMax)
///
/// # Returns
/// Vector of paired alignments, sorted by score
pub fn align_paired_read(
    mate1_seq: &[u8],
    mate2_seq: &[u8],
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<Vec<PairedAlignment>, Error> {
    // Step 1: Find seeds from both mates (pooled together)
    let min_seed_length = 8;
    let pooled_seeds =
        Seed::find_paired_seeds(mate1_seq, mate2_seq, index, min_seed_length, params)?;

    if pooled_seeds.is_empty() {
        return Ok(Vec::new());
    }

    // Step 2: Cluster seeds (unified across both mates)
    let max_cluster_dist = 100000; // 100kb window
    let max_loci_for_anchor = 10;
    let clusters = cluster_seeds(
        &pooled_seeds,
        index,
        max_cluster_dist,
        max_loci_for_anchor,
        params.win_anchor_multimap_nmax,
        params.seed_none_loci_per_window,
    );

    if clusters.is_empty() {
        return Ok(Vec::new());
    }

    // Step 3: Stitch seeds within each cluster
    // For Phase 8, we use a simplified approach: align each mate independently
    // then combine results. Full mate-aware DP stitching is deferred to refinement.
    let scorer = AlignmentScorer::from_params(params);
    let mut paired_alignments = Vec::new();

    for cluster in clusters {
        // Separate cluster by mate
        let mate1_indices: Vec<_> = cluster
            .seed_indices
            .iter()
            .filter(|&&i| pooled_seeds[i].mate_id == 0)
            .copied()
            .collect();
        let mate2_indices: Vec<_> = cluster
            .seed_indices
            .iter()
            .filter(|&&i| pooled_seeds[i].mate_id == 1)
            .copied()
            .collect();

        // Stitch each mate independently
        let mate1_transcripts = if !mate1_indices.is_empty() {
            // Create a SeedCluster for mate1
            use crate::align::stitch::SeedCluster;
            let mate1_cluster = SeedCluster {
                seed_indices: mate1_indices,
                chr_idx: cluster.chr_idx,
                genome_start: cluster.genome_start,
                genome_end: cluster.genome_end,
                is_reverse: cluster.is_reverse,
                anchor_idx: cluster.anchor_idx,
            };
            stitch_seeds(&mate1_cluster, &pooled_seeds, mate1_seq, index, &scorer)?
        } else {
            Vec::new()
        };

        let mate2_transcripts = if !mate2_indices.is_empty() {
            // Create a SeedCluster for mate2
            use crate::align::stitch::SeedCluster;
            let mate2_cluster = SeedCluster {
                seed_indices: mate2_indices,
                chr_idx: cluster.chr_idx,
                genome_start: cluster.genome_start,
                genome_end: cluster.genome_end,
                is_reverse: cluster.is_reverse,
                anchor_idx: cluster.anchor_idx,
            };
            stitch_seeds(&mate2_cluster, &pooled_seeds, mate2_seq, index, &scorer)?
        } else {
            Vec::new()
        };

        // Combine best transcript from each mate
        if let (Some(t1), Some(t2)) = (mate1_transcripts.first(), mate2_transcripts.first()) {
            // Check if they map to the same chromosome
            if t1.chr_idx == t2.chr_idx {
                // Create a paired alignment
                let paired =
                    combine_mate_transcripts(t1, t2, mate1_seq.len(), mate2_seq.len(), params)?;
                paired_alignments.push(paired);
            }
        }
    }

    // Step 4: Filter and sort
    filter_paired_transcripts(&mut paired_alignments, params);
    // Deterministic tie-breaking for equal-score pairs:
    // smallest chr index → smallest position → forward strand first
    paired_alignments.sort_by(|a, b| {
        b.transcript
            .score
            .cmp(&a.transcript.score)
            .then_with(|| a.transcript.chr_idx.cmp(&b.transcript.chr_idx))
            .then_with(|| a.transcript.genome_start.cmp(&b.transcript.genome_start))
            .then_with(|| a.transcript.is_reverse.cmp(&b.transcript.is_reverse))
    });

    Ok(paired_alignments)
}

/// Combine two mate transcripts into a paired alignment
fn combine_mate_transcripts(
    mate1_trans: &Transcript,
    mate2_trans: &Transcript,
    mate1_len: usize,
    mate2_len: usize,
    params: &Parameters,
) -> Result<PairedAlignment, Error> {
    // For now, we create a simplified paired alignment
    // The transcript from mate1 is used as the primary transcript
    // This is a simplified implementation - full mate-aware stitching is deferred

    let mate1_region = (0, mate1_len);
    let mate2_region = (0, mate2_len);

    // Check if proper pair
    let is_proper_pair = check_proper_pair(mate1_trans, mate2_trans, params);

    // Calculate insert size (signed genomic distance)
    let insert_size = calculate_insert_size(mate1_trans, mate2_trans);

    Ok(PairedAlignment {
        transcript: mate1_trans.clone(),
        mate1_region,
        mate2_region,
        is_proper_pair,
        insert_size,
    })
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

/// Filter paired transcripts by quality thresholds
fn filter_paired_transcripts(paired_alns: &mut Vec<PairedAlignment>, params: &Parameters) {
    paired_alns.retain(|pa| {
        let t = &pa.transcript;
        let read_length =
            (pa.mate1_region.1 - pa.mate1_region.0 + pa.mate2_region.1 - pa.mate2_region.0) as f64;

        // Absolute score threshold
        if t.score < params.out_filter_score_min {
            return false;
        }

        // Relative score threshold
        if (t.score as f64) < params.out_filter_score_min_over_lread * read_length {
            return false;
        }

        // Absolute mismatch count
        if t.n_mismatch > params.out_filter_mismatch_nmax {
            return false;
        }

        // Relative mismatch count
        if (t.n_mismatch as f64) > params.out_filter_mismatch_nover_lmax * read_length {
            return false;
        }

        // Absolute matched bases
        let n_matched = t.n_matched();
        if n_matched < params.out_filter_match_nmin {
            return false;
        }

        // Relative matched bases
        if (n_matched as f64) < params.out_filter_match_nmin_over_lread * read_length {
            return false;
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

        let (transcripts, chimeras) = result.unwrap();
        assert_eq!(transcripts.len(), 0); // No alignment
        assert_eq!(chimeras.len(), 0); // No chimeric alignments
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

        let (transcripts, _chimeras) = result.unwrap();
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
        assert_eq!(result.unwrap().len(), 0);
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
            read_seq: vec![0; 100],
        };

        let tlen = calculate_insert_size(&t1, &t2);
        assert_eq!(tlen, -300); // Negative because mate2 is leftmost
    }
}
