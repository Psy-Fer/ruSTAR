/// Read alignment driver function
use crate::align::score::{AlignmentScorer, SpliceMotif};
use crate::align::seed::Seed;
use crate::align::stitch::{
    adjust_mate2_coords, cluster_seeds, finalize_transcript, find_mate_boundary,
    split_working_transcript, stitch_seeds_with_jdb_debug, stitch_seeds_working,
};
use crate::align::transcript::Transcript;
use crate::error::Error;
use crate::index::GenomeIndex;
use crate::params::{IntronMotifFilter, IntronStrandFilter, Parameters};
use crate::stats::UnmappedReason;

/// Result of aligning a single read: (transcripts, chimeric_alignments, n_for_mapq, unmapped_reason)
pub type AlignReadResult = (
    Vec<Transcript>,
    Vec<crate::chimeric::ChimericAlignment>,
    usize,
    Option<UnmappedReason>,
);

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
    /// Pre-finalization combined-read stitching score (wt.score before split_working_transcript).
    /// Used for multi-mapper score-range ranking — matches STAR's combined-read score approach,
    /// which is consistent across duplicate loci. Post-finalization per-mate scores diverge due
    /// to independent extension into soft-clip regions.
    pub combined_wt_score: i32,
    /// Pre-finalization combined-read coverage: sum of exon read spans in the joint WT before
    /// split_working_transcript. Mirrors STAR's nMatch on the joint combined transcript used in
    /// mappedFilter. For valid pairs ≈ len1+len2; for overlapping false-positives ≈ one mate len.
    pub combined_n_match: u32,
}

/// Result of paired-end alignment, covering all mapping outcomes.
#[derive(Debug, Clone)]
pub enum PairedAlignmentResult {
    /// Both mates mapped and paired successfully
    BothMapped(Box<PairedAlignment>),
    /// Only one mate mapped; rescue failed or was not attempted for the other
    HalfMapped {
        /// Transcript of the mapped mate
        mapped_transcript: Transcript,
        /// true = mate1 is the mapped mate, false = mate2 is mapped
        mate1_is_mapped: bool,
    },
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
) -> Result<AlignReadResult, Error> {
    let debug_read = !params.read_name_filter.is_empty() && read_name == params.read_name_filter;

    // Step 1: Find seeds (seedMapMin from params)
    let min_seed_length = params.seed_map_min;
    let seeds = Seed::find_seeds(
        read_seq,
        index,
        min_seed_length,
        params,
        if debug_read { read_name } else { "" },
    )?;

    if debug_read {
        let total_positions: usize = seeds.iter().map(|s| s.sa_end - s.sa_start).sum();
        eprintln!(
            "[DEBUG {}] Seeds: {} seeds, {} total SA positions, read_len={}",
            read_name,
            seeds.len(),
            total_positions,
            read_seq.len()
        );
        let n_lr = seeds.iter().filter(|s| !s.search_rc).count();
        let n_rl = seeds.iter().filter(|s| s.search_rc).count();
        eprintln!(
            "  {} seeds total: {} L→R (sparse), {} R→L (sparse)",
            seeds.len(),
            n_lr,
            n_rl
        );
        for (i, s) in seeds.iter().enumerate() {
            let n_loci = s.sa_end - s.sa_start;
            eprintln!(
                "  seed[{}]: read_pos={}, len={}, n_loci={}, search_rc={}, sa=[{},{})",
                i, s.read_pos, s.length, n_loci, s.search_rc, s.sa_start, s.sa_end
            );
        }
    }

    if seeds.is_empty() {
        if debug_read {
            eprintln!("[DEBUG {}] No seeds found — unmapped", read_name);
        }
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Step 2: Cluster seeds (STAR's bin-based windowing)
    // seed_per_window_nmax capacity eviction is handled inside cluster_seeds()
    let clusters = cluster_seeds(&seeds, index, params, read_seq.len());

    if debug_read {
        eprintln!(
            "[DEBUG {}] Clusters: {} clusters",
            read_name,
            clusters.len()
        );
        for (i, cluster) in clusters.iter().take(10).enumerate() {
            let chr_name = if cluster.chr_idx < index.genome.chr_name.len() {
                &index.genome.chr_name[cluster.chr_idx]
            } else {
                "unknown"
            };
            eprintln!(
                "  cluster[{}]: chr={}, is_reverse={}, seeds={}, anchor_bin={}",
                i,
                chr_name,
                cluster.is_reverse,
                cluster.alignments.len(),
                cluster.anchor_bin,
            );
            for (j, wa) in cluster.alignments.iter().enumerate() {
                let chr_pos = wa.genome_pos.saturating_sub(
                    if cluster.chr_idx < index.genome.chr_start.len() {
                        index.genome.chr_start[cluster.chr_idx]
                    } else {
                        0
                    },
                ) + 1; // 1-based
                eprintln!(
                    "    wa[{}]: read_pos={}, len={}, genome_pos={} ({}:{}), n_rep={}, is_anchor={}",
                    j,
                    wa.read_pos,
                    wa.length,
                    wa.genome_pos,
                    chr_name,
                    chr_pos,
                    wa.n_rep,
                    wa.is_anchor
                );
                if j >= 5 {
                    eprintln!(
                        "    ... ({} more WA entries)",
                        cluster.alignments.len() - j - 1
                    );
                    break;
                }
            }
        }
    }

    if clusters.is_empty() {
        if debug_read {
            eprintln!("[DEBUG {}] No clusters — unmapped", read_name);
        }
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Cap total clusters (alignWindowsPerReadNmax)
    let mut clusters = clusters;
    clusters.truncate(params.align_windows_per_read_nmax);

    // NOTE: STAR's winReadCoverageRelativeMin filter is long-reads-only
    // (#ifdef COMPILE_FOR_LONG_READS in stitchPieces.cpp). Standard STAR
    // does NOT filter clusters by seed coverage. Removed to match STAR.

    // Step 2b: Detect chimeric alignments from multi-cluster seeds (Tier 2)
    let mut chimeric_alignments = Vec::new();
    if params.chim_segment_min > 0 && clusters.len() > 1 {
        use crate::chimeric::ChimericDetector;
        let detector = ChimericDetector::new(params);
        chimeric_alignments
            .extend(detector.detect_from_multi_clusters(&clusters, read_seq, read_name, index)?);
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

    for (ci, cluster) in clusters.iter().enumerate() {
        let debug_name = if debug_read { read_name } else { "" };
        let cluster_transcripts = stitch_seeds_with_jdb_debug(
            cluster,
            read_seq,
            index,
            &scorer,
            junction_db,
            params.align_transcripts_per_window_nmax,
            debug_name,
        )?;
        if debug_read {
            eprintln!(
                "[DEBUG {}] Cluster[{}]: {} transcripts from DP",
                read_name,
                ci,
                cluster_transcripts.len()
            );
            for (ti, t) in cluster_transcripts.iter().enumerate().take(5) {
                let chr_name = if t.chr_idx < index.genome.chr_name.len() {
                    &index.genome.chr_name[t.chr_idx]
                } else {
                    "unknown"
                };
                let cigar_str: String = t.cigar.iter().map(|op| format!("{}", op)).collect();
                eprintln!(
                    "  transcript[{}]: chr={}:{}-{} ({}) score={} mm={} junctions={} cigar={}",
                    ti,
                    chr_name,
                    t.genome_start,
                    t.genome_end,
                    if t.is_reverse { "-" } else { "+" },
                    t.score,
                    t.n_mismatch,
                    t.n_junction,
                    cigar_str
                );
            }
        }
        transcripts.extend(cluster_transcripts);
    }

    // Step 4a: Deduplicate and score-range filter — BEFORE quality filters.
    // STAR order: multMapSelect (score-range) → mappedFilter (quality gates).
    // Doing quality filters first is wrong: it can remove the high-scoring primary,
    // leaving a lower-scoring secondary that then passes as the "best" alignment.

    // Deduplicate transcripts with identical genomic coordinates AND CIGAR.
    transcripts.sort_by(|a, b| {
        (
            a.chr_idx,
            a.genome_start,
            a.genome_end,
            a.is_reverse,
            &a.cigar,
        )
            .cmp(&(
                b.chr_idx,
                b.genome_start,
                b.genome_end,
                b.is_reverse,
                &b.cigar,
            ))
            .then_with(|| b.score.cmp(&a.score))
    });
    transcripts.dedup_by(|a, b| {
        a.chr_idx == b.chr_idx
            && a.genome_start == b.genome_start
            && a.genome_end == b.genome_end
            && a.is_reverse == b.is_reverse
            && a.cigar == b.cigar
    });

    // Sort by score descending (deterministic tie-breaking).
    transcripts.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| a.n_junction.cmp(&b.n_junction))
            .then_with(|| a.chr_idx.cmp(&b.chr_idx))
            .then_with(|| a.genome_start.cmp(&b.genome_start))
            .then_with(|| a.is_reverse.cmp(&b.is_reverse))
    });

    // Score-range filter: keep only alignments within outFilterMultimapScoreRange of the best.
    // (STAR's multMapSelect step — must run before quality filters.)
    if !transcripts.is_empty() {
        let max_score = transcripts[0].score;
        let score_threshold = max_score - params.out_filter_multimap_score_range;
        transcripts.retain(|t| t.score >= score_threshold);
    }

    // Multimap count check: too many loci → unmapped.
    if transcripts.len() > params.out_filter_multimap_nmax as usize {
        let n_loci = transcripts.len();
        transcripts.clear();
        return Ok((
            transcripts,
            chimeric_alignments,
            n_loci,
            Some(UnmappedReason::TooManyLoci),
        ));
    }

    // Step 4: Quality filters (STAR's mappedFilter — runs after score-range selection).
    // STAR uses (Lread-1) for relative thresholds and casts to integer
    // (ReadAlign_mappedFilter.cpp lines 8-9)
    let read_length = read_seq.len() as f64;
    let lread_m1 = (read_seq.len() as f64) - 1.0;

    // Log filtering statistics
    let pre_filter_count = transcripts.len();
    let mut filter_reasons = std::collections::HashMap::new();

    transcripts.retain(|t| {
        // Absolute score threshold
        if t.score < params.out_filter_score_min {
            *filter_reasons.entry("score_min").or_insert(0) += 1;
            return false;
        }

        // Relative score threshold: STAR casts to intScore (i32)
        if t.score < (params.out_filter_score_min_over_lread * lread_m1) as i32 {
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

        // Relative matched bases: STAR casts to uint (u32)
        if n_matched < (params.out_filter_match_nmin_over_lread * lread_m1) as u32 {
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
                if t.junction_motifs.contains(&SpliceMotif::NonCanonical) {
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
        // STAR's RemoveInconsistentStrands removes transcripts that have junctions
        // with mixed intron strand (some imply + strand, some imply - strand).
        // This handles chimeric/impossible transcripts spanning both strands.
        // Note: a reverse-strand read CAN have + strand motifs (antisense reads
        // from + strand genes in unstranded RNA-seq) — this is valid and STAR
        // keeps such reads. Only mixed-strand within one transcript is filtered.
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

    if debug_read {
        eprintln!(
            "[DEBUG {}] After quality filters: {}/{} transcripts remain (reasons: {:?})",
            read_name,
            transcripts.len(),
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

    // Note: STAR sometimes finds 2 equivalent indel placements in homopolymer runs
    // via its recursive stitcher's seed exploration (NH=2 instead of NH=1 for ~5 reads).
    // Generating equivalents post-hoc causes more harm than good (41 false NH=2 vs 5 fixed).
    // The root cause is jR scanning placing insertions at different positions — fixing that
    // would be a better approach than post-hoc enumeration.

    // Step 6: Filter chimeric alignments
    if params.chim_segment_min > 0 {
        chimeric_alignments.retain(|chim| {
            chim.meets_min_segment_length(params.chim_segment_min)
                && chim.meets_min_score(params.chim_score_min)
        });
    }

    // n_for_mapq = transcripts.len() after dedup and filtering.
    // Multi-transcript DP (Phase 16.10) produces multiple transcripts per window
    // for tandem repeats (e.g. rDNA), yielding correct NH → correct MAPQ.
    let n_for_mapq = transcripts.len();

    if debug_read {
        eprintln!(
            "[DEBUG {}] Final: {} transcripts, n_for_mapq={}",
            read_name,
            transcripts.len(),
            n_for_mapq
        );
        for (i, t) in transcripts.iter().enumerate() {
            let chr_name = if t.chr_idx < index.genome.chr_name.len() {
                &index.genome.chr_name[t.chr_idx]
            } else {
                "unknown"
            };
            let cigar_str: String = t.cigar.iter().map(|op| format!("{}", op)).collect();
            eprintln!(
                "  FINAL[{}]: chr={}:{}-{} ({}) score={} mm={} junctions={} cigar={}",
                i,
                chr_name,
                t.genome_start,
                t.genome_end,
                if t.is_reverse { "-" } else { "+" },
                t.score,
                t.n_mismatch,
                t.n_junction,
                cigar_str
            );
        }
    }

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

/// STAR's fragment spacer base (MARK_FRAG_SPACER_BASE = 11 in IncludeDefine.h:174).
/// Not a valid genome base (A=0,C=1,G=2,T=3,N=4), so no seed can span this byte.
const MATE_SPACER_BASE: u8 = 11;

/// Build the STAR-faithful combined PE read: [mate1_fwd][SPACER=11][RC(mate2)].
///
/// Seeds from mate1 (0..len1) map to the forward strand.
/// Seeds from RC(mate2) (len1+1..len1+1+len2) also map to the forward strand for
/// proper pairs (mate2 is reverse-strand, so RC(mate2) is forward-strand).
/// This allows both mates' seeds to cluster together in a single forward-strand window.
fn build_combined_read(mate1: &[u8], mate2: &[u8]) -> Vec<u8> {
    let mut v = Vec::with_capacity(mate1.len() + 1 + mate2.len());
    v.extend_from_slice(mate1);
    v.push(MATE_SPACER_BASE);
    // Append RC(mate2)
    for &b in mate2.iter().rev() {
        v.push(if b < 4 { 3 - b } else { b });
    }
    v
}

/// Assign mate_id to seeds found on the combined read.
///
/// Both L→R and R→L seeds: after find_seeds() converts R→L read_pos to combined-read
/// (L→R) space (seed.rs line 261-262), seed.read_pos is in original combined-read space for all.
/// combined = [mate1(0..len1)][SPACER][RC(mate2)(len1+1..)]
///   read_pos < len1  → mate_id=0 (mate1)
///   read_pos > len1  → mate_id=1 (mate2 RC portion)
fn assign_seed_mate_ids(seeds: &mut Vec<Seed>, len1: usize, _len2: usize) {
    for seed in seeds.iter_mut() {
        seed.mate_id = if seed.read_pos < len1 { 0 } else { 1 };
    }
    // Drop any seed that spans the spacer (read_pos == len1 in combined-read space).
    // This is a guard — the spacer byte=11 is not in the genome, so no seed should span it.
    seeds.retain(|s| s.read_pos != len1);
}

/// Align paired-end reads using STAR's joint combined-read approach (Phase 16.11).
///
/// # Algorithm (STAR-faithful)
/// 1. Build combined read [mate1_fwd][SPACER=11][RC(mate2)] — STAR's canonical PE format
/// 2. Find seeds on combined read; assign mate_id (0=mate1, 1=mate2) by read position
/// 3. Cluster seeds together (both mates in same forward-strand window for proper pairs)
/// 4. Stitch clusters with mate-boundary detection (canonSJ=-3 in STAR)
/// 5. Split joint transcripts at the mate boundary into two SAM records
/// 6. Fall back to independent SE alignment for any clusters not captured jointly
///
/// # Arguments
/// * `mate1_seq` - First mate sequence (encoded)
/// * `mate2_seq` - Second mate sequence (encoded)
/// * `index` - Genome index
/// * `params` - Parameters (includes alignMatesGapMax)
///
/// # Returns
/// Tuple of (paired alignment results, n_for_mapq, unmapped_reason)
pub fn align_paired_read(
    mate1_seq: &[u8],
    mate2_seq: &[u8],
    read_name: &str,
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<(Vec<PairedAlignmentResult>, usize, Option<UnmappedReason>), Error> {
    let len1 = mate1_seq.len();
    let len2 = mate2_seq.len();
    const SPACER_LEN: usize = 1;

    // STAR-faithful joint PE alignment via combined-read approach.
    // Combined = [mate1_fwd][SPACER=11][RC(mate2)].
    // Forward cluster: mate1_fwd and RC(mate2) both map to + strand → proper pair
    //   mate_id=0 = mate1 (read pos [0, len1)), mate_id=1 = mate2 ([len1+1, total_len))
    // Reverse cluster: mate2 and RC(mate1) map to + strand via RC combined read
    //   In RC combined: [mate2][SPACER_RC][RC(mate1)]; after stitch coord conversion:
    //   mate_id=0 = mate2 (read pos [0, len2)), mate_id=1 = mate1 ([len2+1, total_len))
    let combined = build_combined_read(mate1_seq, mate2_seq);
    let rc_mate2: Vec<u8> = mate2_seq
        .iter()
        .rev()
        .map(|&b| if b < 4 { 3 - b } else { b })
        .collect();
    let rc_mate1: Vec<u8> = mate1_seq
        .iter()
        .rev()
        .map(|&b| if b < 4 { 3 - b } else { b })
        .collect();

    let debug_pe = !params.read_name_filter.is_empty() && read_name == params.read_name_filter;

    let mut seeds = Seed::find_seeds(&combined, index, params.seed_map_min, params, "")?;
    assign_seed_mate_ids(&mut seeds, len1, len2);

    if debug_pe {
        eprintln!(
            "[DEBUG-PE] Combined read len={}, seeds={}",
            combined.len(),
            seeds.len()
        );
        for s in &seeds {
            eprintln!(
                "  seed: read_pos={}, len={}, sa_range=[{},{}), rc={}, mate_id={}",
                s.read_pos, s.length, s.sa_start, s.sa_end, s.search_rc, s.mate_id
            );
        }
    }

    let clusters = cluster_seeds(&seeds, index, params, combined.len());

    if debug_pe {
        eprintln!("[DEBUG-PE] Clusters: {}", clusters.len());
        for (i, c) in clusters.iter().enumerate() {
            eprintln!(
                "  cluster[{}]: is_reverse={}, chr={}, {} alignments",
                i,
                c.is_reverse,
                c.chr_idx,
                c.alignments.len()
            );
            for wa in &c.alignments {
                eprintln!(
                    "    WA: read_pos={}, len={}, genome_pos={}, sa_pos={}, anchor={}, mate_id={}",
                    wa.read_pos, wa.length, wa.genome_pos, wa.sa_pos, wa.is_anchor, wa.mate_id
                );
            }
        }
    }
    let scorer = AlignmentScorer::from_params(params);
    let junction_db = if index.junction_db.is_empty() {
        None
    } else {
        Some(&index.junction_db)
    };

    // STAR uses combined read length (len1+spacer+len2) as denominator for PE score threshold.
    // This is stricter than summing individual thresholds when they're computed separately.
    let combined_score_threshold =
        (params.out_filter_score_min_over_lread * (len1 + SPACER_LEN + len2 - 1) as f64) as i32;

    let mut joint_pairs: Vec<PairedAlignment> = Vec::new();
    let mut mate1_candidates: Vec<Transcript> = Vec::new();
    let mut mate2_candidates: Vec<Transcript> = Vec::new();

    for cluster in clusters.iter().take(params.align_windows_per_read_nmax) {
        let cluster_is_reverse = cluster.is_reverse;

        let debug_stitch_name: &str = if debug_pe { read_name } else { "" };
        let (wts, stitch_cluster, _stitch_is_reverse, _stitch_read) = stitch_seeds_working(
            cluster,
            &combined,
            index,
            &scorer,
            junction_db,
            params.align_transcripts_per_window_nmax,
            params.align_mates_gap_max as u64,
            debug_stitch_name,
        )?;

        for wt in wts {
            if wt.exons.is_empty() {
                continue;
            }

            if debug_pe {
                let has_m0 = wt.exons.iter().any(|e| e.mate_id == 0);
                let has_m1 = wt.exons.iter().any(|e| e.mate_id == 1);
                eprintln!(
                    "[DEBUG-WT] score={} n_exons={} has_m0={} has_m1={} cluster_rev={}",
                    wt.score,
                    wt.exons.len(),
                    has_m0,
                    has_m1,
                    cluster_is_reverse
                );
                for (i, e) in wt.exons.iter().enumerate() {
                    eprintln!(
                        "  exon[{}]: r={}-{} g={}-{} mate={}",
                        i, e.read_start, e.read_end, e.genome_start, e.genome_end, e.mate_id
                    );
                }
            }

            if !cluster_is_reverse {
                // --- Forward cluster ---
                // mate_id=0 = mate1 ([0, len1)), mate_id=1 = mate2 ([len1+1, total_len))
                match find_mate_boundary(&wt) {
                    Some(boundary_idx) => {
                        // STAR combined-read score filter (ReadAlign_mappedFilter.cpp line 8):
                        // trBest->maxScore < outFilterScoreMinOverLread * (Lread-1).
                        // Apply scoreGenomicLengthLog2scale penalty (stitchWindowAligns.cpp:262-265)
                        // before comparing to threshold. STAR stores Score=max(0,Score+penalty).
                        // Must check BEFORE split_working_transcript because split sets both
                        // wt1.score=wt.score and wt2.score=wt.score (approximations), and using
                        // wt1.score+wt2.score would double-count, inflating the effective score.
                        let g_span = wt
                            .exons
                            .last()
                            .unwrap()
                            .genome_end
                            .saturating_sub(wt.exons.first().unwrap().genome_start);
                        let length_penalty = scorer.genomic_length_penalty(g_span);
                        let adjusted_score = (wt.score + length_penalty).max(0);
                        // Log near-threshold failures: reads where the genomicLength penalty
                        // causes the score to fall below threshold. These are the ~35 regression
                        // reads introduced by Phase 16.31. Fires for any read pair, not just debug.
                        if wt.score >= combined_score_threshold
                            && adjusted_score < combined_score_threshold
                        {
                            eprintln!(
                                "[NEAR-THRESHOLD-FWD] {} raw={} penalty={} adj={} threshold={} n_exons={} gspan={}",
                                read_name,
                                wt.score,
                                length_penalty,
                                adjusted_score,
                                combined_score_threshold,
                                wt.exons.len(),
                                g_span
                            );
                        }
                        if adjusted_score < combined_score_threshold {
                            continue;
                        }
                        let combined_wt_score = adjusted_score;
                        if debug_pe {
                            eprintln!(
                                "[RUSTAR-FINALIZE-PREEXT] raw={} penalty={} adj={} n_exons={} n_mismatch={}",
                                wt.score,
                                length_penalty,
                                adjusted_score,
                                wt.exons.len(),
                                wt.n_mismatch
                            );
                            for (i, e) in wt.exons.iter().enumerate() {
                                eprintln!(
                                    "[RUSTAR-FINALIZE-EXON[{}]] rStart={} gStart={} len={} iFrag={}",
                                    i,
                                    e.read_start,
                                    e.genome_start,
                                    e.read_end - e.read_start,
                                    e.mate_id
                                );
                            }
                        }
                        // Pre-split combined-read coverage: mirrors STAR's nMatch on the joint
                        // combined transcript (used in mappedFilter). Sum of exon read spans
                        // before finalize_transcript inflates per-mate contributions.
                        let combined_n_match: u32 = wt
                            .exons
                            .iter()
                            .map(|e| (e.read_end - e.read_start) as u32)
                            .sum();

                        // Joint: split at mate1→mate2 boundary
                        let (wt1, wt2) =
                            match split_working_transcript(&wt, boundary_idx, len1, SPACER_LEN) {
                                Some(pair) => pair,
                                None => continue,
                            };
                        if wt1.exons.is_empty() || wt2.exons.is_empty() {
                            continue;
                        }

                        let t1 = match finalize_transcript(
                            &wt1,
                            mate1_seq,
                            index,
                            &scorer,
                            &stitch_cluster,
                            false,
                            false,
                        ) {
                            Some(t) => t,
                            None => continue,
                        };
                        let mut t2 = match finalize_transcript(
                            &wt2,
                            &rc_mate2,
                            index,
                            &scorer,
                            &stitch_cluster,
                            false,
                            true,
                        ) {
                            Some(t) => t,
                            None => continue,
                        };
                        t2.is_reverse = true;
                        t2.read_seq = mate2_seq.to_vec();

                        if t1.chr_idx != t2.chr_idx {
                            continue;
                        }

                        // STAR stitchWindowAligns.cpp:189-218: overlap consistency checks.
                        // When mates overlap in the genome, verify the overlap is consistent.
                        // STAR applies these AFTER backward extension (post-extension EX_G).
                        // Check 1 uses the estimated post-extension genome_start of mate1:
                        //   post_ext_start ≈ wt1.first.genome_start - wt1.first.read_start
                        // (left extension fills in exactly read_start bases when full extension occurs).
                        // This fixes short-insert pairs where raw genome_start > mate2.start due to
                        // the mate1 left soft clip, which post-extension would align correctly.
                        {
                            // STAR check2_lhs uses stitched exon end; use post-finalization t1.genome_end.
                            let left_end = t1.genome_end;
                            let right_start = wt2.exons.first().unwrap().genome_start;
                            if left_end > right_start {
                                // Check 1: post-extension leftMateStart > rightMateStart
                                let wt1_first = wt1.exons.first().unwrap();
                                let left_start_ext = wt1_first
                                    .genome_start
                                    .saturating_sub(wt1_first.read_start as u64);
                                if left_start_ext > right_start {
                                    continue;
                                }
                                // Check 2: leftMateEnd > rightMateEnd (estimated post right-extension)
                                let last_m2 = wt2.exons.last().unwrap();
                                let right_end =
                                    last_m2.genome_start + (len2 - last_m2.read_start) as u64;
                                if left_end > right_end {
                                    continue;
                                }
                            }
                        }

                        // Reject negative insert size (mate2 ends before mate1 starts)
                        if t2.genome_end <= t1.genome_start {
                            continue;
                        }

                        // D5: junction consistency in overlap region (FR: t1=left mate)
                        if !pe_junctions_consistent(&t1, &t2) {
                            continue;
                        }

                        let is_proper_pair = check_proper_pair(&t1, &t2, params);
                        let insert_size = calculate_insert_size(&t1, &t2);
                        joint_pairs.push(PairedAlignment {
                            mate1_transcript: t1,
                            mate2_transcript: t2,
                            mate1_region: (0, len1),
                            mate2_region: (0, len2),
                            is_proper_pair,
                            insert_size,
                            combined_wt_score,
                            combined_n_match,
                        });
                    }
                    None => {
                        // Single-mate: mate_id=0 = mate1, mate_id=1 = mate2
                        let mate = wt.exons.iter().find(|e| e.mate_id != 2).map(|e| e.mate_id);
                        match mate {
                            Some(0) => {
                                if let Some(t) = finalize_transcript(
                                    &wt,
                                    mate1_seq,
                                    index,
                                    &scorer,
                                    &stitch_cluster,
                                    false,
                                    false,
                                ) {
                                    mate1_candidates.push(t);
                                }
                            }
                            Some(1) => {
                                // Positions [len1+1, ..] → adjust to [0, len2)
                                if let Some(wt2) = adjust_mate2_coords(&wt, len1, SPACER_LEN)
                                    && let Some(mut t) = finalize_transcript(
                                        &wt2,
                                        &rc_mate2,
                                        index,
                                        &scorer,
                                        &stitch_cluster,
                                        false,
                                        false,
                                    )
                                {
                                    t.is_reverse = true;
                                    t.read_seq = mate2_seq.to_vec();
                                    mate2_candidates.push(t);
                                }
                            }
                            _ => {}
                        }
                    }
                }
            } else {
                // --- Reverse cluster ---
                // After stitch coord conversion, stitch_read = RC(combined) = [mate2][SPACER_RC][RC(mate1)].
                // Mate identity determined by read POSITION (not mate_id, which is inconsistent
                // for L→R vs R→L seeds after coordinate flip):
                //   read_end <= len2         → mate2 region  → finalize with mate2_seq
                //   read_start >= len2+1     → mate1 region  → finalize with rc_mate1 (adj by -(len2+1))
                let spacer_end = len2 + SPACER_LEN; // = len2 + 1
                let has_mate2_exons = wt
                    .exons
                    .iter()
                    .any(|e| e.read_start < spacer_end && e.read_end <= len2);
                let has_mate1_exons = wt.exons.iter().any(|e| e.read_start >= spacer_end);

                if has_mate2_exons && has_mate1_exons {
                    // Joint: mate2 at [0, len2), mate1 at [len2+1, ..) in stitch_read
                    // Boundary = index of first exon with read_start >= spacer_end (first mate1 exon)
                    let boundary_idx = wt
                        .exons
                        .iter()
                        .position(|e| e.read_start >= spacer_end)
                        .unwrap_or(wt.exons.len());
                    if boundary_idx == 0 || boundary_idx == wt.exons.len() {
                        continue;
                    }
                    // Validate mate2 exons: all must fit within [0, len2)
                    if wt.exons[..boundary_idx].iter().any(|e| e.read_end > len2) {
                        continue;
                    }
                    // Validate mate1 exons: all must start at >= spacer_end
                    if wt.exons[boundary_idx..]
                        .iter()
                        .any(|e| e.read_start < spacer_end)
                    {
                        continue;
                    }

                    // STAR combined-read score filter: check adjusted score before split.
                    // Apply scoreGenomicLengthLog2scale penalty (stitchWindowAligns.cpp:262-265).
                    // split_working_transcript sets wt1.score=wt.score and wt2.score=wt.score,
                    // so using wt1.score+wt2.score would double-count.
                    let g_span = wt
                        .exons
                        .last()
                        .unwrap()
                        .genome_end
                        .saturating_sub(wt.exons.first().unwrap().genome_start);
                    let length_penalty = scorer.genomic_length_penalty(g_span);
                    let adjusted_score = (wt.score + length_penalty).max(0);
                    if wt.score >= combined_score_threshold
                        && adjusted_score < combined_score_threshold
                    {
                        eprintln!(
                            "[NEAR-THRESHOLD-REV] {} raw={} penalty={} adj={} threshold={} n_exons={} gspan={}",
                            read_name,
                            wt.score,
                            length_penalty,
                            adjusted_score,
                            combined_score_threshold,
                            wt.exons.len(),
                            g_span
                        );
                    }
                    if adjusted_score < combined_score_threshold {
                        continue;
                    }
                    let combined_wt_score = adjusted_score;
                    if debug_pe {
                        eprintln!(
                            "[RUSTAR-FINALIZE-PREEXT-REV] raw={} penalty={} adj={} n_exons={} n_mismatch={}",
                            wt.score,
                            length_penalty,
                            adjusted_score,
                            wt.exons.len(),
                            wt.n_mismatch
                        );
                        for (i, e) in wt.exons.iter().enumerate() {
                            eprintln!(
                                "[RUSTAR-FINALIZE-EXON-REV[{}]] rStart={} gStart={} len={} iFrag={}",
                                i,
                                e.read_start,
                                e.genome_start,
                                e.read_end - e.read_start,
                                e.mate_id
                            );
                        }
                    }
                    // Pre-split combined-read coverage: mirrors STAR's nMatch on the joint
                    // combined transcript (used in mappedFilter). Stitch_read for reverse cluster
                    // is [mate2_fwd | SPACER_RC | RC(mate1_fwd)]; exon read spans cover same range.
                    let combined_n_match: u32 = wt
                        .exons
                        .iter()
                        .map(|e| (e.read_end - e.read_start) as u32)
                        .sum();

                    // Build wt2 (mate2 portion, no position adjustment needed)
                    let (wt2, wt1) =
                        match split_working_transcript(&wt, boundary_idx, len2, SPACER_LEN) {
                            Some(pair) => pair,
                            None => continue,
                        };
                    if wt2.exons.is_empty() || wt1.exons.is_empty() {
                        continue;
                    }

                    let t2 = match finalize_transcript(
                        &wt2,
                        mate2_seq,
                        index,
                        &scorer,
                        &stitch_cluster,
                        cluster_is_reverse,
                        false,
                    ) {
                        Some(t) => t,
                        None => continue,
                    };
                    let mut t1 = match finalize_transcript(
                        &wt1,
                        &rc_mate1,
                        index,
                        &scorer,
                        &stitch_cluster,
                        cluster_is_reverse,
                        true,
                    ) {
                        Some(t) => t,
                        None => continue,
                    };
                    t1.is_reverse = true;
                    t1.read_seq = mate1_seq.to_vec();

                    if t1.chr_idx != t2.chr_idx {
                        continue;
                    }

                    // STAR stitchWindowAligns.cpp:189-218: overlap consistency checks.
                    // In the reverse cluster: wt2 = mate2 content (left in genome for RF pairs),
                    // wt1 = mate1 content (right in genome).
                    // Applied pre-finalization using wt2/wt1 seed coords.
                    {
                        // STAR check2_lhs = EX_G + EX_L of left mate's last stitched exon.
                        // STAR's stitcher extends seeds to full match boundaries; ruSTAR's WT
                        // stops at the last seed's end. Use post-finalization t2.genome_end.
                        let left_end = t2.genome_end;
                        let right_start = wt1.exons.first().unwrap().genome_start;
                        if left_end > right_start {
                            // Check 1: estimated post-extension left mate start > right mate start → reject.
                            // Use the same post-extension estimate as the forward cluster (Phase 16.30):
                            // left_start_ext = first_exon.genome_start - first_exon.read_start.
                            // This accounts for short-insert overlapping pairs where mate2's first seed
                            // starts after mate1's first seed in genome space but, after left-extension,
                            // mate2's start is correctly at or before mate1's start.
                            let first_left = wt2.exons.first().unwrap();
                            let left_start_ext = first_left
                                .genome_start
                                .saturating_sub(first_left.read_start as u64);
                            if left_start_ext > right_start {
                                continue;
                            }
                            // Check 2: left mate ends after estimated right mate end → reject
                            let last_right = wt1.exons.last().unwrap();
                            let right_end = last_right.genome_start
                                + (len1 as u64).saturating_sub(last_right.read_start as u64);
                            if left_end > right_end {
                                continue;
                            }
                        }
                    }

                    // Reject negative insert size (mate1 ends before mate2 starts)
                    if t1.genome_end <= t2.genome_start {
                        continue;
                    }

                    // D5: junction consistency in overlap region (RF: t2=left mate)
                    if !pe_junctions_consistent(&t2, &t1) {
                        continue;
                    }

                    let is_proper_pair = check_proper_pair(&t1, &t2, params);
                    let insert_size = calculate_insert_size(&t1, &t2);
                    joint_pairs.push(PairedAlignment {
                        mate1_transcript: t1,
                        mate2_transcript: t2,
                        mate1_region: (0, len1),
                        mate2_region: (0, len2),
                        is_proper_pair,
                        insert_size,
                        combined_wt_score,
                        combined_n_match,
                    });
                } else if has_mate2_exons && !has_mate1_exons {
                    // Mate2-only: all exons in [0, len2) → finalize with mate2_seq
                    if wt.exons.iter().any(|e| e.read_end > len2) {
                        continue;
                    }
                    if let Some(t) = finalize_transcript(
                        &wt,
                        mate2_seq,
                        index,
                        &scorer,
                        &stitch_cluster,
                        cluster_is_reverse,
                        false,
                    ) {
                        mate2_candidates.push(t);
                    }
                } else if has_mate1_exons && !has_mate2_exons {
                    // Mate1-only: all exons in [len2+1, ..] → adjust to [0, len1), finalize with rc_mate1
                    if let Some(wt1) = adjust_mate2_coords(&wt, len2, SPACER_LEN)
                        && let Some(mut t) = finalize_transcript(
                            &wt1,
                            &rc_mate1,
                            index,
                            &scorer,
                            &stitch_cluster,
                            cluster_is_reverse,
                            false,
                        )
                    {
                        t.is_reverse = true;
                        t.read_seq = mate1_seq.to_vec();
                        mate1_candidates.push(t);
                    }
                }
            }
        }
    }

    if debug_pe {
        eprintln!(
            "[DEBUG-PE] After cluster loop: joint_pairs={}, mate1_cands={}, mate2_cands={}",
            joint_pairs.len(),
            mate1_candidates.len(),
            mate2_candidates.len()
        );
        for (i, pa) in joint_pairs.iter().enumerate() {
            eprintln!(
                "  pair[{}]: M1 chr={} pos={} rev={} score={} | M2 chr={} pos={} rev={} score={} | wt_score={}",
                i,
                pa.mate1_transcript.chr_idx,
                pa.mate1_transcript.genome_start,
                pa.mate1_transcript.is_reverse,
                pa.mate1_transcript.score,
                pa.mate2_transcript.chr_idx,
                pa.mate2_transcript.genome_start,
                pa.mate2_transcript.is_reverse,
                pa.mate2_transcript.score,
                pa.combined_wt_score
            );
        }
    }

    // --- Decision tree: dedup, sort, score-filter joint pairs, then fallbacks ---
    joint_pairs.sort_by(|a, b| {
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
    joint_pairs.dedup_by(|a, b| {
        a.mate1_transcript.chr_idx == b.mate1_transcript.chr_idx
            && a.mate1_transcript.genome_start == b.mate1_transcript.genome_start
            && a.mate1_transcript.is_reverse == b.mate1_transcript.is_reverse
            && a.mate2_transcript.genome_start == b.mate2_transcript.genome_start
            && a.mate2_transcript.is_reverse == b.mate2_transcript.is_reverse
    });
    joint_pairs.sort_by(|a, b| {
        b.combined_wt_score
            .cmp(&a.combined_wt_score)
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
    if !joint_pairs.is_empty() {
        // STAR's stitchWindowAligns score-range filter uses the combined WT score (post genomic
        // length penalty) to rank and filter multi-mappers. Using per-mate finalized scores
        // (t1.score + t2.score) fails because finalize_transcript's extendAlign can equalize
        // scores across spurious splice pairs: all paired transcripts get the same per-mate score
        // even when their combined WT scores differ by 4-5 points (e.g., a spurious 144kb intron
        // pair scores wt=288 vs correct unspliced wt=292, but both produce mate scores of 148+146).
        let best_score = joint_pairs
            .iter()
            .map(|pa| pa.combined_wt_score)
            .max()
            .unwrap_or(0);
        let score_threshold = best_score - params.out_filter_multimap_score_range;
        joint_pairs.retain(|pa| pa.combined_wt_score >= score_threshold);
    }
    filter_paired_transcripts(&mut joint_pairs, params);

    // STAR mappedFilter: nBestMate > outFilterMultimapNmax → unmapped(multi).
    // Applied after score-range filter and quality filter, same as STAR's mappedFilter ordering.
    if joint_pairs.len() > params.out_filter_multimap_nmax as usize {
        return Ok((Vec::new(), 0, Some(UnmappedReason::TooManyLoci)));
    }

    // 1. BothMapped from joint combined-read path
    if !joint_pairs.is_empty() {
        let pe_mapq_n = joint_pairs.len().max(1);
        let results = joint_pairs
            .into_iter()
            .map(|pa| PairedAlignmentResult::BothMapped(Box::new(pa)))
            .collect();
        return Ok((results, pe_mapq_n, None));
    }

    // 2. STAR: mono-mate transcript score (≈len_mate) < 0.66*(len1+len2+1) → "unmapped: too short"
    // STAR's mappedFilter uses combined read Lread as denominator, so single-mate candidates
    // always fail the threshold. STAR has no rescue step — return unmapped.
    Ok((Vec::new(), 0, Some(UnmappedReason::TooShort)))
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
/// STAR's mappedFilter (ReadAlign_mappedFilter.cpp) applies ALL quality checks to the combined
/// read as a single unit (Lread = len1+1+len2). There are no per-mate quality filters for PE.
fn filter_paired_transcripts(paired_alns: &mut Vec<PairedAlignment>, params: &Parameters) {
    paired_alns.retain(|pa| {
        let t1 = &pa.mate1_transcript;
        let t2 = &pa.mate2_transcript;
        let mate1_len = (pa.mate1_region.1 - pa.mate1_region.0) as f64;
        let mate2_len = (pa.mate2_region.1 - pa.mate2_region.0) as f64;
        // Lread-1 for the combined read: (len1+1+len2) - 1 = len1 + len2
        let combined_lread_m1 = mate1_len + mate2_len;

        let combined_nm = t1.n_mismatch + t2.n_mismatch;

        // STAR's mappedFilter applies nMatch against Lread-1 on the JOINT combined transcript
        // (ReadAlign_mappedFilter.cpp line 21). For overlapping pairs (both mates at same genome
        // region), STAR's joint transcript covers only one mate's worth → nMatch ≈ 137 < 198.
        // ruSTAR extends each mate independently after split, inflating per-mate n_matched().
        //
        // Fix: use the pre-split combined-read coverage captured before split_working_transcript,
        // which directly mirrors STAR's joint transcript nMatch without extension inflation.
        let combined_match = pa.combined_n_match;

        // Score: trBest->maxScore < outFilterScoreMin || < outFilterScoreMinOverLread*(Lread-1)
        // STAR checks the combined WT score (ReadAlign_mappedFilter.cpp), not the sum of per-mate
        // finalized scores. Using t1.score + t2.score incorrectly double-applies the genomic
        // length penalty (once per mate), rejecting valid overlapping short-insert pairs where
        // per-mate spans are small (penalty ≈ −2 each → sum penalty −4 vs combined penalty −2).
        // combined_wt_score was already verified ≥ combined_score_threshold before the split.
        let combined_score = pa.combined_wt_score;
        if combined_score < params.out_filter_score_min
            || combined_score < (params.out_filter_score_min_over_lread * combined_lread_m1) as i32
        {
            return false;
        }

        // Match bases: trBest->nMatch < outFilterMatchNmin || < outFilterMatchNminOverLread*(Lread-1)
        if combined_match < params.out_filter_match_nmin
            || combined_match < (params.out_filter_match_nmin_over_lread * combined_lread_m1) as u32
        {
            return false;
        }

        // Mismatches: nMM > outFilterMismatchNmaxTotal || nMM/rLength > outFilterMismatchNoverLmax
        // outFilterMismatchNmaxTotal = min(outFilterMismatchNmax, outFilterMismatchNoverReadLmax*(len1+len2))
        if combined_nm > params.out_filter_mismatch_nmax
            || (combined_nm as f64)
                > params.out_filter_mismatch_nover_lmax * (mate1_len + mate2_len)
        {
            return false;
        }

        true
    });

    // STAR's mappedFilter: if nTr > outFilterMultimapNmax → unmapped=multi (clear all).
    // Unlike SE which truncates for output, PE must REJECT the entire pair when over the limit.
    if paired_alns.len() > params.out_filter_multimap_nmax as usize {
        paired_alns.clear();
    }
}

/// Extract splice junctions from a Transcript's CIGAR as (donor, acceptor) pairs.
/// Junction coords are in genomic space (0-based). Junctions only exist where CigarOp::RefSkip is.
fn extract_junctions_from_cigar(t: &Transcript) -> Vec<(u64, u64)> {
    let mut junctions = Vec::new();
    let mut genome_pos = t.genome_start;
    for op in &t.cigar {
        use crate::align::transcript::CigarOp;
        match op {
            CigarOp::Match(n) | CigarOp::Equal(n) | CigarOp::Diff(n) | CigarOp::Del(n) => {
                genome_pos += *n as u64;
            }
            CigarOp::RefSkip(n) => {
                let donor = genome_pos;
                let acceptor = genome_pos + *n as u64;
                junctions.push((donor, acceptor));
                genome_pos = acceptor;
            }
            CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) => {}
        }
    }
    junctions
}

/// D5: Check junction consistency in the overlapping region of paired-end mates.
/// When mates overlap in the genome, every splice junction in the overlap from the
/// left mate must appear in the right mate, and vice versa.
/// Implements STAR stitchWindowAligns.cpp check after overlap detection.
///
/// `left` is the mate with the lower genome_start, `right` is the other mate.
pub(crate) fn pe_junctions_consistent(left: &Transcript, right: &Transcript) -> bool {
    // Overlapping region: [overlap_start, overlap_end)
    let overlap_start = left.genome_start.max(right.genome_start);
    let overlap_end = left.genome_end.min(right.genome_end);
    if overlap_start >= overlap_end {
        return true; // No overlap — nothing to check
    }

    let left_juncs = extract_junctions_from_cigar(left);
    let right_juncs = extract_junctions_from_cigar(right);

    // Every junction from left that falls within the overlap must be in right too
    for (donor, acceptor) in &left_juncs {
        if *donor >= overlap_start
            && *acceptor <= overlap_end
            && !right_juncs.iter().any(|(d, a)| d == donor && a == acceptor)
        {
            return false;
        }
    }
    // Every junction from right that falls within the overlap must be in left too
    for (donor, acceptor) in &right_juncs {
        if *donor >= overlap_start
            && *acceptor <= overlap_end
            && !left_juncs.iter().any(|(d, a)| d == donor && a == acceptor)
        {
            return false;
        }
    }
    true
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

        let result = align_paired_read(&mate1, &mate2, "test", &index, &params);
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

        // Create a transcript with conflicting strand motifs (mixed + and - within one transcript)
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

        // Create a transcript with consistent strand motifs (all + strand)
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

        // Note: STAR's RemoveInconsistentStrands filters transcripts where
        // junctions have MIXED implied strand (both + and - within one transcript).
        // It does NOT compare junction strand vs alignment strand — a reverse-strand
        // read at a GT/AG junction (antisense of + strand gene) is valid and kept.

        // Verify mixed-strand is detected
        let mut has_plus = false;
        let mut has_minus = false;
        for motif in &t_inconsistent.junction_motifs {
            match motif.implied_strand() {
                Some('+') => has_plus = true,
                Some('-') => has_minus = true,
                _ => {}
            }
        }
        assert!(has_plus && has_minus); // Inconsistent (mixed strands)

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

        // Also verify: single CT/AC on forward-strand or single GT/AG on reverse-strand
        // are NOT filtered (STAR keeps these — antisense reads are valid)
        let ctac_only = vec![SpliceMotif::CtAc];
        let (mut hp, mut hm) = (false, false);
        for m in &ctac_only {
            match m.implied_strand() {
                Some('+') => hp = true,
                Some('-') => hm = true,
                _ => {}
            }
        }
        assert!(!hp && hm); // Only minus → NOT mixed → NOT filtered

        // Verify the filter enum
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

    #[test]
    fn test_align_paired_both_unmapped() {
        // Both mates are all N's → both unmapped → empty Vec
        let index = make_test_index();
        let params = make_test_params();

        let mate1 = vec![4, 4, 4, 4, 4, 4, 4, 4];
        let mate2 = vec![4, 4, 4, 4, 4, 4, 4, 4];

        let (results, n_for_mapq, unmapped_reason) =
            align_paired_read(&mate1, &mate2, "test", &index, &params).unwrap();
        assert!(results.is_empty(), "Both unmapped should return empty Vec");
        assert_eq!(n_for_mapq, 0);
        assert!(unmapped_reason.is_some(), "Should have unmapped reason");
    }

    #[test]
    fn test_paired_alignment_result_enum_variants() {
        use crate::align::transcript::{CigarOp, Exon};

        let transcript = Transcript {
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

        // Test BothMapped variant
        let both = PairedAlignmentResult::BothMapped(Box::new(PairedAlignment {
            mate1_transcript: transcript.clone(),
            mate2_transcript: transcript.clone(),
            mate1_region: (0, 100),
            mate2_region: (0, 100),
            is_proper_pair: true,
            insert_size: 200,
            combined_wt_score: 0,
            combined_n_match: 200,
        }));
        assert!(matches!(both, PairedAlignmentResult::BothMapped(_)));

        // Test HalfMapped variant
        let half = PairedAlignmentResult::HalfMapped {
            mapped_transcript: transcript,
            mate1_is_mapped: true,
        };
        assert!(matches!(half, PairedAlignmentResult::HalfMapped { .. }));

        // Verify mate1_is_mapped
        if let PairedAlignmentResult::HalfMapped {
            mate1_is_mapped, ..
        } = half
        {
            assert!(mate1_is_mapped);
        }
    }
}
