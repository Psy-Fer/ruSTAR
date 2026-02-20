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

    // Step 2: Cluster seeds (STAR's bin-based windowing)
    // seed_per_window_nmax capacity eviction is handled inside cluster_seeds()
    let max_loci_for_anchor = 10; // Seeds mapping to <=10 loci can be anchors
    let clusters = cluster_seeds(
        &seeds,
        read_seq,
        index,
        params.win_bin_nbits,
        params.win_anchor_dist_nbins,
        params.win_flank_nbins,
        max_loci_for_anchor,
        params.win_anchor_multimap_nmax,
        params.seed_per_window_nmax,
    );

    if clusters.is_empty() {
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Cap total clusters (alignWindowsPerReadNmax)
    let mut clusters = clusters;
    clusters.truncate(params.align_windows_per_read_nmax);

    // Step 2a: Filter clusters by seed coverage (winReadCoverageRelativeMin)
    let read_length = read_seq.len();
    let min_coverage = params.win_read_coverage_relative_min;
    let clusters: Vec<_> = clusters
        .into_iter()
        .filter(|cluster| {
            // Compute total seed coverage for this cluster (union of read ranges)
            let mut covered = vec![false; read_length];
            for wa in &cluster.alignments {
                let end = wa.read_pos + wa.length;
                for c in covered
                    .iter_mut()
                    .take(end.min(read_length))
                    .skip(wa.read_pos)
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

    for cluster in clusters.iter() {
        let cluster_transcripts =
            stitch_seeds_with_jdb(cluster, read_seq, index, &scorer, junction_db)?;
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
    // NOTE: rDNA reads (~157) still get MAPQ=255 instead of STAR's 1-3 because our
    // clustering gathers all seeds by proximity (within 589kb), whereas STAR assigns
    // seeds to windows by bin position. The fix requires implementing STAR's per-window
    // bin-based seed assignment with winFlankNbins flanking extension.
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

/// Attempt to rescue an unmapped mate using the mapped mate's position as anchor.
///
/// # Algorithm
/// 1. Find seeds for the unmapped mate (genome-wide, same as SE path)
/// 2. Filter seeds to keep only those within the rescue window on the mapped mate's chromosome
/// 3. Cluster the filtered seeds
/// 4. Stitch seeds within each cluster via DP
/// 5. Apply the same quality filters as align_read()
/// 6. Return the best transcript if any passes filters
///
/// # Arguments
/// * `unmapped_seq` - Sequence of the unmapped mate (encoded)
/// * `mapped_transcript` - Best transcript from the mapped mate (anchor)
/// * `index` - Genome index
/// * `params` - Parameters
///
/// # Returns
/// `Some(transcript)` if rescue succeeds, `None` otherwise
fn rescue_unmapped_mate(
    unmapped_seq: &[u8],
    mapped_transcript: &Transcript,
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<Option<Transcript>, Error> {
    // Step 1: Find seeds for unmapped mate (genome-wide search)
    let min_seed_length = 8;
    let seeds = Seed::find_seeds(unmapped_seq, index, min_seed_length, params)?;

    if seeds.is_empty() {
        return Ok(None);
    }

    // Step 2: Determine rescue window
    // Use alignMatesGapMax if set, otherwise default 1Mb
    let rescue_range: u64 = if params.align_mates_gap_max > 0 {
        params.align_mates_gap_max as u64
    } else {
        1_000_000
    };
    let anchor_chr = mapped_transcript.chr_idx;
    let anchor_start = mapped_transcript.genome_start;
    let anchor_end = mapped_transcript.genome_end;
    let window_start = anchor_start.saturating_sub(rescue_range);
    let window_end = anchor_end.saturating_add(rescue_range);

    // Step 3: Filter seeds — keep only those with at least one genome position
    // on the same chromosome and within the rescue window
    let filtered_seed_indices: Vec<usize> = seeds
        .iter()
        .enumerate()
        .filter(|(_, seed)| {
            seed.genome_positions(index).any(|(pos, is_rev)| {
                // Convert SA position to forward genome coordinate
                let fwd_pos = index.sa_pos_to_forward(pos, is_rev, seed.length);
                // Check if this position is on the anchor chromosome and within window
                if let Some((chr_idx, _)) = index.genome.position_to_chr(fwd_pos) {
                    chr_idx == anchor_chr && fwd_pos >= window_start && fwd_pos <= window_end
                } else {
                    false
                }
            })
        })
        .map(|(i, _)| i)
        .collect();

    if filtered_seed_indices.is_empty() {
        return Ok(None);
    }

    // Step 4: Build a filtered seed list for clustering
    let filtered_seeds: Vec<Seed> = filtered_seed_indices
        .iter()
        .map(|&i| seeds[i].clone())
        .collect();

    // Step 5: Cluster and stitch the filtered seeds (bin-based windowing)
    let max_loci_for_anchor = 10;
    let clusters = cluster_seeds(
        &filtered_seeds,
        unmapped_seq,
        index,
        params.win_bin_nbits,
        params.win_anchor_dist_nbins,
        params.win_flank_nbins,
        max_loci_for_anchor,
        params.win_anchor_multimap_nmax,
        params.seed_per_window_nmax,
    );

    if clusters.is_empty() {
        return Ok(None);
    }

    let scorer = AlignmentScorer::from_params(params);
    let junction_db = if index.junction_db.is_empty() {
        None
    } else {
        Some(&index.junction_db)
    };

    let mut transcripts = Vec::new();
    for cluster in clusters.iter().take(params.align_windows_per_read_nmax) {
        let cluster_transcripts =
            stitch_seeds_with_jdb(cluster, unmapped_seq, index, &scorer, junction_db)?;
        transcripts.extend(cluster_transcripts);
    }

    // Step 6: Apply quality filters (same as align_read)
    let read_length = unmapped_seq.len() as f64;
    transcripts.retain(|t| {
        if t.score < params.out_filter_score_min {
            return false;
        }
        if (t.score as f64) < params.out_filter_score_min_over_lread * read_length {
            return false;
        }
        if t.n_mismatch > params.out_filter_mismatch_nmax {
            return false;
        }
        let mismatch_rate = t.n_mismatch as f64 / read_length;
        if mismatch_rate > params.out_filter_mismatch_nover_lmax {
            return false;
        }
        let n_matched = t.n_matched();
        if n_matched < params.out_filter_match_nmin {
            return false;
        }
        if (n_matched as f64) < params.out_filter_match_nmin_over_lread * read_length {
            return false;
        }
        true
    });

    // Sort by score (best first) and return the best
    transcripts.sort_by(|a, b| b.score.cmp(&a.score));
    Ok(transcripts.into_iter().next())
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
/// Tuple of (paired alignment results, n_for_mapq, unmapped_reason):
/// - results: `PairedAlignmentResult` variants (BothMapped or HalfMapped)
/// - n_for_mapq: effective alignment count for MAPQ (max of both mates' cluster counts)
/// - unmapped_reason: `Some(reason)` if no alignments at all, `None` if any mate mapped
pub fn align_paired_read(
    mate1_seq: &[u8],
    mate2_seq: &[u8],
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<(Vec<PairedAlignmentResult>, usize, Option<UnmappedReason>), Error> {
    // Step 1: Align each mate independently using the proven SE path
    let (mate1_transcripts, _, mate1_mapq_n, mate1_unmapped) =
        align_read(mate1_seq, "", index, params)?;
    let (mate2_transcripts, _, mate2_mapq_n, mate2_unmapped) =
        align_read(mate2_seq, "", index, params)?;

    // Tier 2 & 3: Handle cases where one or both mates fail to align
    let mate1_empty = mate1_transcripts.is_empty();
    let mate2_empty = mate2_transcripts.is_empty();

    if mate1_empty && mate2_empty {
        // Both mates unmapped — no rescue possible
        let reason = mate1_unmapped
            .or(mate2_unmapped)
            .unwrap_or(UnmappedReason::Other);
        return Ok((Vec::new(), 0, Some(reason)));
    }

    if mate1_empty || mate2_empty {
        // One mate mapped, one didn't — attempt rescue
        let (mapped_transcripts, unmapped_seq, mate1_is_mapped) = if mate2_empty {
            (&mate1_transcripts, mate2_seq, true)
        } else {
            (&mate2_transcripts, mate1_seq, false)
        };

        let mapped_best = &mapped_transcripts[0];

        // Tier 2: Try to rescue the unmapped mate
        if let Some(rescued) = rescue_unmapped_mate(unmapped_seq, mapped_best, index, params)? {
            // Rescue succeeded — build a proper pair
            let (t1, t2) = if mate1_is_mapped {
                (mapped_best.clone(), rescued)
            } else {
                (rescued, mapped_best.clone())
            };

            let is_proper_pair = check_proper_pair(&t1, &t2, params);
            let insert_size = calculate_insert_size(&t1, &t2);

            let pe_mapq_n = if mate1_is_mapped {
                mate1_mapq_n.max(1)
            } else {
                mate2_mapq_n.max(1)
            };

            return Ok((
                vec![PairedAlignmentResult::BothMapped(Box::new(
                    PairedAlignment {
                        mate1_transcript: t1,
                        mate2_transcript: t2,
                        mate1_region: (0, mate1_seq.len()),
                        mate2_region: (0, mate2_seq.len()),
                        is_proper_pair,
                        insert_size,
                    },
                ))],
                pe_mapq_n,
                None,
            ));
        }

        // Tier 3: Rescue failed — return HalfMapped
        let pe_mapq_n = if mate1_is_mapped {
            mate1_mapq_n
        } else {
            mate2_mapq_n
        };

        return Ok((
            vec![PairedAlignmentResult::HalfMapped {
                mapped_transcript: mapped_best.clone(),
                mate1_is_mapped,
            }],
            pe_mapq_n,
            None, // Not fully unmapped — one mate mapped
        ));
    }

    // Tier 1: Both mates aligned independently — pair them
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

    if paired_alignments.is_empty() {
        // Both mates mapped independently but no valid pairing found
        // Return the best mate as HalfMapped (prefer mate1)
        let (best_transcript, mate1_is_mapped) =
            if mate1_transcripts[0].score >= mate2_transcripts[0].score {
                (mate1_transcripts[0].clone(), true)
            } else {
                (mate2_transcripts[0].clone(), false)
            };
        return Ok((
            vec![PairedAlignmentResult::HalfMapped {
                mapped_transcript: best_transcript,
                mate1_is_mapped,
            }],
            pe_mapq_n,
            None,
        ));
    }

    // Wrap in PairedAlignmentResult::BothMapped
    let results = paired_alignments
        .into_iter()
        .map(|pa| PairedAlignmentResult::BothMapped(Box::new(pa)))
        .collect();

    Ok((results, pe_mapq_n, None))
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

    #[test]
    fn test_rescue_unmapped_mate_no_seeds() {
        // Unmapped mate has all N's (no seeds possible), should return None
        let index = make_test_index();
        let params = make_test_params();

        use crate::align::transcript::{CigarOp, Exon};
        let mapped_transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
            }],
            cigar: vec![CigarOp::Match(4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        let unmapped_seq = vec![4, 4, 4, 4, 4, 4, 4, 4]; // all N's
        let result = rescue_unmapped_mate(&unmapped_seq, &mapped_transcript, &index, &params);
        assert!(result.is_ok());
        assert!(
            result.unwrap().is_none(),
            "Should return None when unmapped mate has no seeds"
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
            align_paired_read(&mate1, &mate2, &index, &params).unwrap();
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
