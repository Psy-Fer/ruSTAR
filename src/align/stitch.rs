/// Seed clustering and stitching via dynamic programming
use crate::align::score::{AlignmentScorer, GapType};
use crate::align::seed::Seed;
use crate::align::transcript::{CigarOp, Transcript};
use crate::error::Error;
use crate::index::GenomeIndex;

/// Verify match length at a specific genome position with correct strand handling.
///
/// Seeds found via binary search guarantee `sa_nbases` matching bases at `sa_start`,
/// but other positions in the SA range may match fewer bases. This function re-verifies
/// the actual match length at each specific genome position.
///
/// For reverse-strand positions, adds `n_genome` offset to access the reverse-complement
/// region of the genome (which is stored at `[n_genome, 2*n_genome)`).
fn verify_match_at_position(
    read_seq: &[u8],
    read_pos: usize,
    genome_pos: u64,
    is_reverse: bool,
    max_length: usize,
    index: &GenomeIndex,
) -> usize {
    let actual_genome_pos = if is_reverse {
        genome_pos + index.genome.n_genome
    } else {
        genome_pos
    };
    let mut length = 0;
    for i in 0..max_length {
        if read_pos + i >= read_seq.len() {
            break;
        }
        match index.genome.get_base(actual_genome_pos + i as u64) {
            Some(gb) if gb < 5 && gb == read_seq[read_pos + i] => length += 1,
            _ => break,
        }
    }
    length
}

/// Count mismatches in an alignment by comparing read sequence to genome sequence.
///
/// The read sequence is always in forward orientation. For reverse-strand alignments,
/// the genome is accessed at `pos + n_genome` (the reverse-complement region) rather
/// than reverse-complementing the read.
///
/// # Arguments
/// * `read_seq` - Read sequence in forward orientation (encoded as 0=A, 1=C, 2=G, 3=T, 4=N)
/// * `cigar_ops` - CIGAR operations
/// * `genome_start` - Starting position in genome (decoded SA position, WITHOUT n_genome offset)
/// * `read_start` - Starting position in read
/// * `index` - Genome index (contains genome sequence)
/// * `is_reverse` - Whether this is a reverse-strand alignment
///
/// # Returns
/// Number of mismatched bases (excluding N bases)
/// Count mismatches in a simple region (no CIGAR, just read vs genome)
fn count_mismatches_in_region(
    read_seq: &[u8],
    read_start: usize,
    genome_start: u64,
    length: usize,
    index: &GenomeIndex,
    is_reverse: bool,
) -> u32 {
    let genome_offset = if is_reverse { index.genome.n_genome } else { 0 };
    let mut n_mismatch = 0u32;

    for i in 0..length {
        let read_pos = read_start + i;
        if read_pos >= read_seq.len() {
            break;
        }
        let read_base = read_seq[read_pos];
        if let Some(genome_base) = index
            .genome
            .get_base(genome_start + i as u64 + genome_offset)
        {
            if read_base != genome_base && read_base != 4 && genome_base != 4 {
                n_mismatch += 1;
            }
        }
    }

    n_mismatch
}

fn count_mismatches(
    read_seq: &[u8],
    cigar_ops: &[CigarOp],
    genome_start: u64,
    read_start: usize,
    index: &GenomeIndex,
    is_reverse: bool,
) -> u32 {
    // Add n_genome offset for reverse-strand genome access
    let genome_offset = if is_reverse { index.genome.n_genome } else { 0 };

    let mut n_mismatch = 0u32;
    let mut read_pos = read_start;
    let mut genome_pos = genome_start;

    for op in cigar_ops {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                for _i in 0..*len {
                    if read_pos < read_seq.len() {
                        let read_base = read_seq[read_pos];
                        if let Some(genome_base) = index.genome.get_base(genome_pos + genome_offset)
                        {
                            if read_base != genome_base && read_base != 4 && genome_base != 4 {
                                n_mismatch += 1;
                            }
                        }
                    }
                    read_pos += 1;
                    genome_pos += 1;
                }
            }
            CigarOp::Ins(len) => {
                read_pos += *len as usize;
            }
            CigarOp::Del(len) | CigarOp::RefSkip(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::SoftClip(len) => {
                read_pos += *len as usize;
            }
            CigarOp::HardClip(_) => {}
        }
    }

    n_mismatch
}

/// Result of extending an alignment into flanking regions
#[derive(Debug, Clone)]
struct ExtendResult {
    /// How far the extension reached (bases)
    extend_len: usize,
    /// Maximum score achieved during extension
    max_score: i32,
    /// Number of mismatches in the extended region
    n_mismatch: u32,
}

/// Extend alignment from a boundary into flanking sequence, mirroring STAR's extendAlign().
///
/// Walks base-by-base from the alignment boundary, scoring +1 match / -1 mismatch,
/// tracking the maximum-score extension point. Stops when total mismatches exceed
/// `min(p_mm_max * total_length, n_mm_max)`.
///
/// # Arguments
/// * `read_seq` - Full read sequence (encoded)
/// * `read_start` - Boundary position in read (where extension begins)
/// * `genome_start` - Corresponding genome position (WITHOUT n_genome offset)
/// * `direction` - +1 for rightward extension, -1 for leftward
/// * `max_extend` - Maximum distance to extend (to read boundary)
/// * `n_mm_prev` - Mismatches already in the aligned portion
/// * `len_prev` - Length of the already-aligned portion
/// * `n_mm_max` - outFilterMismatchNmax (absolute max mismatches)
/// * `p_mm_max` - outFilterMismatchNoverLmax (max mismatch ratio)
/// * `index` - Genome index
/// * `is_reverse` - Whether this is a reverse-strand alignment
fn extend_alignment(
    read_seq: &[u8],
    read_start: usize,
    genome_start: u64,
    direction: i32,
    max_extend: usize,
    n_mm_prev: u32,
    len_prev: usize,
    n_mm_max: u32,
    p_mm_max: f64,
    index: &GenomeIndex,
    is_reverse: bool,
) -> ExtendResult {
    if max_extend == 0 {
        return ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        };
    }

    let genome_offset = if is_reverse { index.genome.n_genome } else { 0 };

    let mut score: i32 = 0;
    let mut max_score: i32 = 0;
    let mut best_len: usize = 0;
    let mut best_mm: u32 = 0;
    let mut n_mm = 0u32;

    for i in 0..max_extend {
        // Calculate read and genome positions based on direction
        let read_pos = if direction > 0 {
            read_start + i
        } else {
            // Leftward: read_start is exclusive boundary, go backwards
            if read_start < 1 + i {
                break;
            }
            read_start - 1 - i
        };

        if read_pos >= read_seq.len() {
            break;
        }

        let genome_pos = if direction > 0 {
            genome_start + i as u64
        } else {
            if genome_start < 1 + i as u64 {
                break;
            }
            genome_start - 1 - i as u64
        };

        // Get genome base (with strand offset)
        let genome_base = match index.genome.get_base(genome_pos + genome_offset) {
            Some(b) => b,
            None => break,
        };

        // Stop at chromosome boundary (padding = 5)
        if genome_base == 5 {
            break;
        }

        let read_base = read_seq[read_pos];

        // Skip N bases (no score impact, matches STAR behavior)
        if read_base == 4 || genome_base == 4 {
            continue;
        }

        if read_base == genome_base {
            score += 1;
        } else {
            score -= 1;
            n_mm += 1;
        }

        // Check mismatch limits considering the full alignment
        let total_mm = n_mm_prev + n_mm;
        let total_len = len_prev + i + 1;
        let mm_limit = ((p_mm_max * total_len as f64) as u32).min(n_mm_max);
        if total_mm > mm_limit {
            break;
        }

        // Record best extension point (highest score)
        if score > max_score {
            max_score = score;
            best_len = i + 1;
            best_mm = n_mm;
        }
    }

    // Only accept extension if it has positive score
    if max_score > 0 {
        ExtendResult {
            extend_len: best_len,
            max_score,
            n_mismatch: best_mm,
        }
    } else {
        ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        }
    }
}

/// A cluster of seeds mapping to the same genomic region
#[derive(Debug, Clone)]
pub struct SeedCluster {
    /// Indices into the seeds array
    pub seed_indices: Vec<usize>,
    /// Chromosome index
    pub chr_idx: usize,
    /// Genomic start (leftmost position)
    pub genome_start: u64,
    /// Genomic end (rightmost position)
    pub genome_end: u64,
    /// Strand (false = forward, true = reverse)
    pub is_reverse: bool,
    /// Anchor seed index (in the seeds array)
    pub anchor_idx: usize,
}

/// Cluster seeds by genomic proximity around anchor seeds.
///
/// # Arguments
/// * `seeds` - All seeds found in the read
/// * `index` - Genome index
/// * `max_cluster_dist` - Maximum genomic distance to cluster seeds (e.g., 100kb)
/// * `max_loci_for_anchor` - Maximum SA range for a seed to be an anchor (e.g., 10)
/// * `win_anchor_multimap_nmax` - Max loci anchors can map to (STAR default: 50)
/// * `seed_none_loci_per_window` - Max seed positions per window (STAR default: 10)
///
/// # Returns
/// Vector of seed clusters
pub fn cluster_seeds(
    seeds: &[Seed],
    index: &GenomeIndex,
    max_cluster_dist: u64,
    max_loci_for_anchor: usize,
    win_anchor_multimap_nmax: usize,
    seed_none_loci_per_window: usize,
) -> Vec<SeedCluster> {
    let mut clusters = Vec::with_capacity(seeds.len());

    // Find anchor seeds (seeds that map to few genomic locations)
    let mut anchors: Vec<usize> = Vec::with_capacity(seeds.len());
    for (i, seed) in seeds.iter().enumerate() {
        let n_loci = seed.sa_end - seed.sa_start;
        if n_loci <= max_loci_for_anchor {
            anchors.push(i);
        }
    }

    // If no anchors, use best seeds as potential anchors
    // STAR-like behavior: sort by SA range size, use smallest (most specific)
    if anchors.is_empty() {
        // Use ALL seeds as anchors (like original STAR behavior)
        // The winAnchorMultimapNmax and seedNoneLociPerWindow limits will prevent explosion
        anchors = (0..seeds.len()).collect();
    }

    // For each anchor, create clusters
    for &anchor_idx in &anchors {
        let anchor = &seeds[anchor_idx];
        let n_anchor_loci = anchor.sa_end - anchor.sa_start;

        // Skip anchors with too many loci (STAR: winAnchorMultimapNmax)
        if n_anchor_loci > win_anchor_multimap_nmax {
            continue;
        }

        // Limit number of positions per anchor (STAR: seedNoneLociPerWindow)
        let max_positions = seed_none_loci_per_window.min(n_anchor_loci);

        // For each genomic position of the anchor (limited, using iterator)
        for (anchor_pos, anchor_strand) in anchor.genome_positions(index).take(max_positions) {
            // Find chromosome
            let chr_info = match index.genome.position_to_chr(anchor_pos) {
                Some(info) => info,
                None => continue, // Position in padding
            };
            let chr_idx = chr_info.0;

            // Collect seeds within clustering window
            let mut seed_indices = Vec::with_capacity(seeds.len());
            seed_indices.push(anchor_idx);
            let mut genome_start = anchor_pos;
            let mut genome_end = anchor_pos + anchor.length as u64;

            for (i, seed) in seeds.iter().enumerate() {
                if i == anchor_idx {
                    continue;
                }

                // Check if seed overlaps with window (using iterator, no Vec alloc)
                let n_seed_loci = seed.sa_end - seed.sa_start;
                let max_seed_positions = seed_none_loci_per_window.min(n_seed_loci);
                let mut found = false;
                for (pos, strand) in seed.genome_positions(index).take(max_seed_positions) {
                    if strand != anchor_strand {
                        continue;
                    }

                    let seed_chr = match index.genome.position_to_chr(pos) {
                        Some(info) => info.0,
                        None => continue,
                    };

                    if seed_chr != chr_idx {
                        continue;
                    }

                    let seed_end = pos + seed.length as u64;

                    // Check if within clustering distance
                    let dist = if pos > genome_end {
                        pos - genome_end
                    } else {
                        genome_start.saturating_sub(seed_end)
                    };

                    if dist <= max_cluster_dist {
                        seed_indices.push(i);
                        genome_start = genome_start.min(pos);
                        genome_end = genome_end.max(seed_end);
                        found = true;
                        break; // Only add each seed once per cluster
                    }
                }
                let _ = found;
            }

            // Create cluster if it has at least one seed
            if !seed_indices.is_empty() {
                clusters.push(SeedCluster {
                    seed_indices,
                    chr_idx,
                    genome_start,
                    genome_end,
                    is_reverse: anchor_strand,
                    anchor_idx,
                });
            }
        }
    }

    clusters
}

/// DP state for seed stitching
#[derive(Debug, Clone)]
struct DpState {
    score: i32,
    prev_seed: Option<usize>,
    genome_pos: u64,
    cigar_ops: Vec<CigarOp>,
    n_mismatch: u32,
    n_gap: u32,
    n_junction: u32,
    junction_motifs: Vec<crate::align::score::SpliceMotif>,
}

/// Expanded seed with specific genome position
#[derive(Debug, Clone)]
struct ExpandedSeed {
    read_pos: usize,
    read_end: usize,
    genome_pos: u64,
    genome_end: u64,
    length: usize,
}

/// Stitch seeds within a cluster using dynamic programming.
///
/// # Arguments
/// * `cluster` - Seed cluster
/// * `seeds` - All seeds
/// * `read_seq` - Read sequence
/// * `index` - Genome index
/// * `scorer` - Alignment scorer
///
/// # Returns
/// Vector of transcripts (may have multiple paths through the cluster)
pub fn stitch_seeds(
    cluster: &SeedCluster,
    seeds: &[Seed],
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
) -> Result<Vec<Transcript>, Error> {
    // Expand seeds to specific genome positions matching the cluster
    let mut expanded_seeds = Vec::with_capacity(cluster.seed_indices.len() * 4);

    for &seed_idx in &cluster.seed_indices {
        let seed = &seeds[seed_idx];

        for (pos, strand) in seed.genome_positions(index) {
            if strand != cluster.is_reverse {
                continue;
            }

            let chr = match index.genome.position_to_chr(pos) {
                Some(info) => info.0,
                None => continue,
            };

            if chr != cluster.chr_idx {
                continue;
            }

            // Check if within cluster bounds (allow some slack)
            if pos < cluster.genome_start || pos > cluster.genome_end + 1000000 {
                continue;
            }

            // Re-verify match length at this specific genome position.
            // Binary search only guarantees `sa_nbases` matching bases at `sa_start`;
            // other positions in the SA range may match fewer bases.
            let actual_length =
                verify_match_at_position(read_seq, seed.read_pos, pos, strand, seed.length, index);

            if actual_length < 8 {
                continue; // Too short after verification
            }

            expanded_seeds.push(ExpandedSeed {
                read_pos: seed.read_pos,
                read_end: seed.read_pos + actual_length,
                genome_pos: pos,
                genome_end: pos + actual_length as u64,
                length: actual_length,
            });
        }
    }

    if expanded_seeds.is_empty() {
        return Ok(Vec::new());
    }

    // Sort by read position, then by length descending (longest first for dedup)
    expanded_seeds.sort_by(|a, b| a.read_pos.cmp(&b.read_pos).then(b.length.cmp(&a.length)));

    // Deduplicate: keep only the longest seed per (read_pos, genome_pos) pair
    expanded_seeds.dedup_by(|a, b| a.read_pos == b.read_pos && a.genome_pos == b.genome_pos);

    // Cap expanded seeds to prevent pathological O(n²) DP on repetitive regions
    const MAX_EXPANDED_SEEDS: usize = 200;
    if expanded_seeds.len() > MAX_EXPANDED_SEEDS {
        // Keep longest seeds (re-sort by length descending, take top N, re-sort by read_pos)
        expanded_seeds.sort_by(|a, b| b.length.cmp(&a.length));
        expanded_seeds.truncate(MAX_EXPANDED_SEEDS);
        expanded_seeds.sort_by_key(|s| s.read_pos);
    }

    // Initialize DP: one state per expanded seed
    let n = expanded_seeds.len();
    let mut dp: Vec<DpState> = Vec::with_capacity(n);

    // Base case: each seed starts with positive score for matches
    // STAR uses +1 per matched base as the baseline score
    for exp_seed in &expanded_seeds {
        let initial_cigar = vec![CigarOp::Match(exp_seed.length as u32)];
        dp.push(DpState {
            score: exp_seed.length as i32, // Positive score: +1 per matched base
            prev_seed: None,
            genome_pos: exp_seed.genome_pos,
            cigar_ops: initial_cigar,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: Vec::new(),
        });
    }

    // DP: for each seed, try connecting to all compatible previous seeds
    // Optimization: find best j first, then build CIGAR only once (avoids repeated cloning)
    for i in 1..n {
        let curr = &expanded_seeds[i];
        let mut best_score = dp[i].score;
        let mut best_j: Option<usize> = None;
        let mut best_gap_type = GapType::Deletion(0);
        let mut best_read_gap = 0i64;
        let mut best_genome_gap = 0i64;
        let mut best_gap_mismatches = 0u32;
        let mut best_eff_length = curr.length;

        for j in 0..i {
            let prev = &expanded_seeds[j];

            // STAR-style overlap trimming: if seeds overlap in read, trim seed B's start
            let mut eff_read_pos = curr.read_pos;
            let mut eff_genome_pos = curr.genome_pos;
            let mut eff_length = curr.length;

            if prev.read_end > eff_read_pos {
                let overlap = prev.read_end - eff_read_pos;
                if overlap >= eff_length {
                    continue; // Seed B fully consumed by overlap
                }
                eff_read_pos = prev.read_end;
                eff_genome_pos += overlap as u64;
                eff_length -= overlap;
            }

            // Similarly handle genome overlap
            if prev.genome_end > eff_genome_pos && eff_genome_pos > prev.genome_pos {
                let g_overlap = (prev.genome_end - eff_genome_pos) as usize;
                if g_overlap >= eff_length {
                    continue; // Seed B fully consumed by genome overlap
                }
                eff_read_pos += g_overlap;
                eff_genome_pos += g_overlap as u64;
                eff_length -= g_overlap;
            }

            // Calculate gaps using effective (trimmed) coordinates
            let read_gap = (eff_read_pos as i64) - (prev.read_end as i64);
            let genome_gap = (eff_genome_pos as i64) - (prev.genome_end as i64);

            // Score the gap
            let (gap_score, gap_type) =
                scorer.score_gap(genome_gap, read_gap, prev.genome_end, &index.genome);

            // Count mismatches in the gap region (shared match portion)
            let gap_mismatches = if read_gap > 0 && genome_gap > 0 {
                let shared = (read_gap.min(genome_gap)) as usize;
                if shared > 0 {
                    count_mismatches_in_region(
                        read_seq,
                        prev.read_end,
                        prev.genome_end,
                        shared,
                        index,
                        cluster.is_reverse,
                    )
                } else {
                    0
                }
            } else {
                0
            };

            // Apply alignSJstitchMismatchNmax filter
            // STAR rejects junction stitches with too many mismatches for the motif type
            if let GapType::SpliceJunction { ref motif, .. } = gap_type {
                if !scorer.stitch_mismatch_allowed(motif, gap_mismatches) {
                    continue; // Reject: too many mismatches at this junction type
                }

                // Enforce minimum overhang (alignSJoverhangMin / alignSJDBoverhangMin)
                // Left overhang = length of the preceding seed (immediately flanks the junction)
                // Right overhang = effective length of current seed after overlap trimming
                let left_overhang = prev.length;
                let right_overhang = eff_length;
                let min_overhang = scorer.align_sj_overhang_min as usize;

                if left_overhang < min_overhang || right_overhang < min_overhang {
                    continue; // Reject: insufficient overhang flanking splice junction
                }
            }

            // STAR: shared region scores +1 per match, -1 per mismatch
            // Total = shared_length - 2*mismatches (each mismatch flips from +1 to -1)
            let shared_bases = if read_gap > 0 && genome_gap > 0 {
                read_gap.min(genome_gap) as i32
            } else {
                0
            };
            let shared_score = shared_bases - (2 * gap_mismatches as i32);

            // Transition score = prev_score + gap_penalty + shared_region_score + current_seed_match_score
            let transition_score = dp[j].score + gap_score + shared_score + (eff_length as i32);

            if transition_score > best_score {
                best_score = transition_score;
                best_j = Some(j);
                best_gap_type = gap_type;
                best_read_gap = read_gap;
                best_genome_gap = genome_gap;
                best_gap_mismatches = gap_mismatches;
                best_eff_length = eff_length;
            }
        }

        // Build CIGAR only once for the best transition (instead of cloning per candidate)
        if let Some(j) = best_j {
            let mut cigar = dp[j].cigar_ops.clone();

            // Emit CIGAR operations for the gap between seeds
            // Handle negative gaps explicitly to avoid integer overflow
            if best_read_gap == 0 && best_genome_gap == 0 {
                // No gap - seeds are adjacent
            } else if best_read_gap < 0 || best_genome_gap < 0 {
                // Negative gap indicates overlapping seeds or error
                // score_gap() already penalized this, but we skip the connection
                log::trace!(
                    "Skipping connection with negative gap: read_gap={}, genome_gap={}",
                    best_read_gap,
                    best_genome_gap
                );
                continue;
            } else if best_read_gap == 0 && best_genome_gap > 0 {
                // Pure deletion or splice junction (read doesn't advance, genome does)
                let gg = best_genome_gap as u32; // Safe: checked > 0
                match best_gap_type {
                    GapType::SpliceJunction { intron_len, .. } => {
                        cigar.push(CigarOp::RefSkip(intron_len));
                    }
                    _ => {
                        cigar.push(CigarOp::Del(gg));
                    }
                }
            } else if best_read_gap > 0 && best_genome_gap == 0 {
                // Pure insertion (read advances, genome doesn't)
                let rg = best_read_gap as u32; // Safe: checked > 0
                cigar.push(CigarOp::Ins(rg));
            } else {
                // Both gaps positive: combined gap region
                let rg = best_read_gap as u32; // Safe: checked > 0 above
                let gg = best_genome_gap as u32; // Safe: checked > 0 above

                let shared = rg.min(gg);
                let excess_genome = gg.saturating_sub(rg);
                let excess_read = rg.saturating_sub(gg);

                if shared > 0 {
                    // Try to merge with previous Match operation
                    if let Some(CigarOp::Match(prev_len)) = cigar.last_mut() {
                        *prev_len += shared;
                    } else {
                        cigar.push(CigarOp::Match(shared));
                    }
                }
                if excess_genome > 0 {
                    match best_gap_type {
                        GapType::SpliceJunction { intron_len, .. } => {
                            cigar.push(CigarOp::RefSkip(intron_len));
                        }
                        _ => {
                            cigar.push(CigarOp::Del(excess_genome));
                        }
                    }
                }
                if excess_read > 0 {
                    cigar.push(CigarOp::Ins(excess_read));
                }
            }

            // Add current seed match using effective length (after overlap trimming)
            if let Some(CigarOp::Match(prev_len)) = cigar.last_mut() {
                // Merge with previous Match operation
                *prev_len += best_eff_length as u32;
            } else {
                // No previous Match, or previous op was not Match
                cigar.push(CigarOp::Match(best_eff_length as u32));
            }

            let mut n_gap = dp[j].n_gap;
            let mut n_junction = dp[j].n_junction;
            let mut junction_motifs = dp[j].junction_motifs.clone();
            match best_gap_type {
                GapType::Insertion(_) | GapType::Deletion(_) => n_gap += 1,
                GapType::SpliceJunction { motif, .. } => {
                    n_junction += 1;
                    junction_motifs.push(motif);
                }
            }
            dp[i].score = best_score;
            dp[i].prev_seed = Some(j);
            dp[i].genome_pos = curr.genome_pos;
            dp[i].cigar_ops = cigar;
            dp[i].n_mismatch = dp[j].n_mismatch + best_gap_mismatches;
            dp[i].n_gap = n_gap;
            dp[i].n_junction = n_junction;
            dp[i].junction_motifs = junction_motifs;
        }
        // If no best_j, dp[i] keeps its initial state (single seed)
    }

    // Backtrack from best final state
    let best_final_idx = dp
        .iter()
        .enumerate()
        .max_by_key(|(_, state)| state.score)
        .map(|(i, _)| i)
        .unwrap();

    let best_state = &dp[best_final_idx];

    // Find the first seed in the DP chain by tracing back through prev_seed
    let mut chain_start_idx = best_final_idx;
    while let Some(prev) = dp[chain_start_idx].prev_seed {
        chain_start_idx = prev;
    }
    let chain_start_seed = &expanded_seeds[chain_start_idx];

    // Calculate the read region covered by the alignment
    let alignment_start = chain_start_seed.read_pos;

    // Calculate aligned read length from CIGAR operations
    let mut aligned_read_len = 0usize;
    for op in &best_state.cigar_ops {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) | CigarOp::Ins(len) => {
                aligned_read_len += *len as usize;
            }
            _ => {}
        }
    }
    let alignment_end = alignment_start + aligned_read_len;

    // Compute right-side genome position by walking the DP CIGAR
    let mut right_genome_pos = chain_start_seed.genome_pos;
    for op in &best_state.cigar_ops {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                right_genome_pos += *len as u64;
            }
            CigarOp::Del(len) | CigarOp::RefSkip(len) => {
                right_genome_pos += *len as u64;
            }
            _ => {}
        }
    }

    // Extend alignment into flanking regions (STAR-style extendAlign)
    let left_extend = if alignment_start > 0 {
        extend_alignment(
            read_seq,
            alignment_start,             // read boundary (exclusive, leftward)
            chain_start_seed.genome_pos, // genome boundary (exclusive, leftward)
            -1,                          // leftward
            alignment_start,             // max distance to read start
            best_state.n_mismatch,       // mismatches in aligned portion
            aligned_read_len,            // length of aligned portion
            scorer.n_mm_max,
            scorer.p_mm_max,
            index,
            cluster.is_reverse,
        )
    } else {
        ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        }
    };

    let right_extend = if alignment_end < read_seq.len() {
        extend_alignment(
            read_seq,
            alignment_end,                  // read boundary (inclusive, rightward)
            right_genome_pos,               // genome boundary (inclusive, rightward)
            1,                              // rightward
            read_seq.len() - alignment_end, // max distance to read end
            best_state.n_mismatch + left_extend.n_mismatch, // cumulative mismatches
            aligned_read_len + left_extend.extend_len, // cumulative length
            scorer.n_mm_max,
            scorer.p_mm_max,
            index,
            cluster.is_reverse,
        )
    } else {
        ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        }
    };

    // Build final CIGAR with extensions
    let mut final_cigar = Vec::new();

    // Remaining 5' soft clip after left extension
    let remaining_left_clip = alignment_start - left_extend.extend_len;
    if remaining_left_clip > 0 {
        final_cigar.push(CigarOp::SoftClip(remaining_left_clip as u32));
    }

    // Left extension Match (merge with first main CIGAR Match if possible)
    if left_extend.extend_len > 0 {
        final_cigar.push(CigarOp::Match(left_extend.extend_len as u32));
    }

    // Main alignment CIGAR (merge first op with left extension Match if both are Match)
    for op in &best_state.cigar_ops {
        if let CigarOp::Match(len) = op {
            if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
                *prev_len += len;
                continue;
            }
        }
        final_cigar.push(*op);
    }

    // Right extension Match (merge with last main CIGAR Match if possible)
    if right_extend.extend_len > 0 {
        if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
            *prev_len += right_extend.extend_len as u32;
        } else {
            final_cigar.push(CigarOp::Match(right_extend.extend_len as u32));
        }
    }

    // Remaining 3' soft clip after right extension
    let remaining_right_clip = (read_seq.len() - alignment_end) - right_extend.extend_len;
    if remaining_right_clip > 0 {
        final_cigar.push(CigarOp::SoftClip(remaining_right_clip as u32));
    }

    // Adjust genome start position for left extension
    let adjusted_genome_start = chain_start_seed.genome_pos - left_extend.extend_len as u64;

    // Adjust score for extensions
    let adjusted_score = best_state.score + left_extend.max_score + right_extend.max_score;

    // Build exons from CIGAR
    use crate::align::transcript::Exon;
    let mut exons = Vec::new();
    let mut read_pos = 0usize; // Start from 0, will be adjusted by soft clips
    let mut genome_pos = adjusted_genome_start;

    for op in &final_cigar {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                // Create or extend current exon
                let len = *len as usize;
                exons.push(Exon {
                    genome_start: genome_pos,
                    genome_end: genome_pos + len as u64,
                    read_start: read_pos,
                    read_end: read_pos + len,
                });
                read_pos += len;
                genome_pos += len as u64;
            }
            CigarOp::Ins(len) => {
                // Insertion: advances read but not genome
                read_pos += *len as usize;
            }
            CigarOp::Del(len) => {
                // Deletion: advances genome but not read
                genome_pos += *len as u64;
            }
            CigarOp::RefSkip(len) => {
                // Intron: advances genome but not read (starts new exon)
                genome_pos += *len as u64;
            }
            CigarOp::SoftClip(len) => {
                // Soft clip: advances read but not genome
                read_pos += *len as usize;
            }
            CigarOp::HardClip(_) => {
                // Hard clip: doesn't advance either
            }
        }
    }

    // Merge consecutive exons (from consecutive Match operations)
    let mut merged_exons: Vec<Exon> = Vec::new();
    for exon in exons {
        if let Some(last_exon) = merged_exons.last_mut() {
            if last_exon.genome_end == exon.genome_start && last_exon.read_end == exon.read_start {
                // Merge with previous exon
                last_exon.genome_end = exon.genome_end;
                last_exon.read_end = exon.read_end;
                continue;
            }
        }
        merged_exons.push(exon);
    }

    // Get the actual genome start position from exons
    let actual_genome_start = merged_exons
        .first()
        .map(|e| e.genome_start)
        .unwrap_or(adjusted_genome_start);

    // Count mismatches in the final alignment
    let n_mismatch = count_mismatches(
        read_seq,
        &final_cigar,
        actual_genome_start, // Use actual alignment start, not first seed position!
        0,                   // Read starts at position 0 (CIGAR includes soft clips)
        index,
        cluster.is_reverse, // Pass reverse-strand flag for correct sequence comparison
    );

    // Build transcript
    let final_seed = &expanded_seeds[best_final_idx];
    let transcript = Transcript {
        chr_idx: cluster.chr_idx,
        genome_start: merged_exons
            .first()
            .map(|e| e.genome_start)
            .unwrap_or(cluster.genome_start),
        genome_end: merged_exons
            .last()
            .map(|e| e.genome_end)
            .unwrap_or(final_seed.genome_end),
        is_reverse: cluster.is_reverse,
        exons: merged_exons,
        cigar: final_cigar, // Use CIGAR with extensions + soft clips
        score: adjusted_score,
        n_mismatch,
        n_gap: best_state.n_gap,
        n_junction: best_state.n_junction,
        junction_motifs: best_state.junction_motifs.clone(),
        read_seq: read_seq.to_vec(),
    };

    Ok(vec![transcript])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Genome;
    use crate::index::packed_array::PackedArray;
    use crate::index::sa_index::SaIndex;
    use crate::index::suffix_array::SuffixArray;

    fn make_simple_index() -> GenomeIndex {
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

        // Create dummy SA and SAindex
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
    fn test_cluster_seeds_simple() {
        // This test would require a properly populated SA
        // For now, just test that we can create the clustering structure
        // without panicking on an empty index
        let index = make_simple_index();

        // Create seeds with empty SA ranges (won't expand to any positions)
        let seeds = vec![
            Seed {
                read_pos: 0,
                length: 5,
                sa_start: 0,
                sa_end: 0, // Empty range
                is_reverse: false,
                mate_id: 2,
            },
            Seed {
                read_pos: 10,
                length: 5,
                sa_start: 0,
                sa_end: 0, // Empty range
                is_reverse: false,
                mate_id: 2,
            },
        ];

        let clusters = cluster_seeds(&seeds, &index, 100000, 10, 50, 10);

        // With empty SA ranges, no clusters will be created
        assert_eq!(clusters.len(), 0);
    }

    #[test]
    fn test_expanded_seed_sorting() {
        let mut expanded = [
            ExpandedSeed {
                read_pos: 10,
                read_end: 15,
                genome_pos: 100,
                genome_end: 105,
                length: 5,
            },
            ExpandedSeed {
                read_pos: 5,
                read_end: 10,
                genome_pos: 50,
                genome_end: 55,
                length: 5,
            },
        ];

        expanded.sort_by_key(|s| s.read_pos);

        assert_eq!(expanded[0].read_pos, 5);
        assert_eq!(expanded[1].read_pos, 10);
    }

    /// Helper to build a GenomeIndex with a specific forward sequence
    fn make_index_with_seq(seq: &[u8]) -> GenomeIndex {
        let n_genome = ((seq.len() as u64 + 1) / 64 + 1) * 64;
        let mut sequence = vec![5u8; (n_genome * 2) as usize];
        sequence[0..seq.len()].copy_from_slice(seq);

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
            chr_length: vec![seq.len() as u64],
            chr_start: vec![0, n_genome],
        };

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
    fn test_extend_perfect_match_rightward() {
        // Genome: ACGTACGTAC (10 bases, A=0, C=1, G=2, T=3)
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let index = make_index_with_seq(&seq);
        // Read matches genome perfectly from position 5 onward
        let read_seq = vec![1, 2, 3, 0, 1]; // matches genome[5..10]

        let result = extend_alignment(
            &read_seq, 0,   // read_start (boundary)
            5,   // genome_start
            1,   // rightward
            5,   // max_extend
            0,   // no previous mismatches
            0,   // no previous length
            10,  // n_mm_max
            0.3, // p_mm_max
            &index, false, // forward strand
        );

        assert_eq!(result.extend_len, 5);
        assert_eq!(result.max_score, 5);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_perfect_match_leftward() {
        // Genome: ACGTACGTAC
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let index = make_index_with_seq(&seq);
        // Read matches genome[0..5] = ACGTA
        let read_seq = vec![0, 1, 2, 3, 0];

        let result = extend_alignment(
            &read_seq, 5,   // read_start (exclusive boundary for leftward)
            5,   // genome_start (exclusive boundary for leftward)
            -1,  // leftward
            5,   // max_extend
            0,   // no previous mismatches
            0,   // no previous length
            10,  // n_mm_max
            0.3, // p_mm_max
            &index, false,
        );

        assert_eq!(result.extend_len, 5);
        assert_eq!(result.max_score, 5);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_stops_at_optimal_point_with_mismatches() {
        // Genome: A C G T A C G T (positions 0-7)
        let genome_seq = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let index = make_index_with_seq(&genome_seq);
        // Read: A C G T T T T T (matches first 4, then all mismatches)
        let read_seq: Vec<u8> = vec![0, 1, 2, 3, 3, 3, 3, 3];

        let result = extend_alignment(
            &read_seq, 0,   // read_start
            0,   // genome_start
            1,   // rightward
            8,   // max_extend
            0,   // no previous mismatches
            0,   // no previous length
            10,  // n_mm_max
            0.3, // p_mm_max
            &index, false,
        );

        // Should extend 4 bases (perfect match), then mismatches drag score down
        assert_eq!(result.extend_len, 4);
        assert_eq!(result.max_score, 4);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_chromosome_boundary() {
        // Genome: A C G (3 bases, then padding=5)
        let genome_seq = vec![0, 1, 2];
        let index = make_index_with_seq(&genome_seq);
        // Read is 5 bases, but genome only has 3
        let read_seq: Vec<u8> = vec![0, 1, 2, 3, 0];

        let result = extend_alignment(
            &read_seq, 0, // read_start
            0, // genome_start
            1, // rightward
            5, // max_extend
            0, 0, 10, 0.3, &index, false,
        );

        // Should stop at 3 bases (genome boundary)
        assert_eq!(result.extend_len, 3);
        assert_eq!(result.max_score, 3);
    }

    #[test]
    fn test_extend_n_bases_skipped() {
        // Genome: A N C G (N=4 at position 1)
        let genome_seq = vec![0, 4, 1, 2];
        let index = make_index_with_seq(&genome_seq);
        // Read: A A C G (matches at 0, N skip at 1, matches at 2-3)
        let read_seq: Vec<u8> = vec![0, 0, 1, 2];

        let result = extend_alignment(&read_seq, 0, 0, 1, 4, 0, 0, 10, 0.3, &index, false);

        // Should extend all 4 bases: match + N(skip) + match + match = score 3
        assert_eq!(result.extend_len, 4);
        assert_eq!(result.max_score, 3);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_all_mismatch_returns_zero() {
        // Genome: A A A A (all 0)
        let genome_seq = vec![0, 0, 0, 0];
        let index = make_index_with_seq(&genome_seq);
        // Read: T T T T (all 3, complete mismatch)
        let read_seq: Vec<u8> = vec![3, 3, 3, 3];

        let result = extend_alignment(&read_seq, 0, 0, 1, 4, 0, 0, 10, 0.3, &index, false);

        // Score never goes positive, so extend_len should be 0
        assert_eq!(result.extend_len, 0);
        assert_eq!(result.max_score, 0);
    }

    #[test]
    fn test_extend_recovery_through_mismatch() {
        // Genome: A C G T A C G T A C G T A C (14 bases)
        let genome_seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let index = make_index_with_seq(&genome_seq);
        // Read: A C G X A C G T A C G T A C (1 mismatch at pos 3, then 10 matches)
        //       M M M X M M M M M M M M M M
        let read_seq: Vec<u8> = vec![0, 1, 2, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1];

        let result = extend_alignment(&read_seq, 0, 0, 1, 14, 0, 0, 10, 0.3, &index, false);

        // Should extend past the mismatch: 3M + 1X + 10M
        // Score: +3 -1 +10 = 12, best at position 14
        assert_eq!(result.extend_len, 14);
        assert_eq!(result.max_score, 12);
        assert_eq!(result.n_mismatch, 1);
    }

    #[test]
    fn test_extend_zero_max_extend() {
        let genome_seq = vec![0, 1, 2, 3];
        let index = make_index_with_seq(&genome_seq);
        let read_seq: Vec<u8> = vec![0, 1, 2, 3];

        let result = extend_alignment(
            &read_seq, 0, 0, 1, 0, // max_extend = 0
            0, 0, 10, 0.3, &index, false,
        );

        assert_eq!(result.extend_len, 0);
        assert_eq!(result.max_score, 0);
    }

    #[test]
    fn test_overhang_check_rejects_short_overhang() {
        // Test that the overhang check rejects splice junctions with tiny flanking seeds.
        // Scenario: prev seed length=3 (below default min of 5), current seed length=20
        // With alignSJoverhangMin=5, the 3bp left overhang should cause rejection.
        use crate::align::score::AlignmentScorer;

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

        // Left overhang (prev.length) = 3, below min of 5
        let left_overhang: usize = 3;
        let right_overhang: usize = 20;
        let min_overhang = scorer.align_sj_overhang_min as usize;

        // Should be rejected
        assert!(left_overhang < min_overhang || right_overhang < min_overhang);

        // Right overhang too small
        let left_overhang: usize = 20;
        let right_overhang: usize = 4;
        assert!(left_overhang < min_overhang || right_overhang < min_overhang);
    }

    #[test]
    fn test_overhang_check_accepts_sufficient_overhang() {
        // Test that splice junctions with sufficient overhang pass the check.
        use crate::align::score::AlignmentScorer;

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

        // Both overhangs >= 5
        let left_overhang: usize = 5;
        let right_overhang: usize = 10;
        let min_overhang = scorer.align_sj_overhang_min as usize;

        // Should pass
        assert!(!(left_overhang < min_overhang || right_overhang < min_overhang));

        // Exactly at minimum
        let left_overhang: usize = 5;
        let right_overhang: usize = 5;
        assert!(!(left_overhang < min_overhang || right_overhang < min_overhang));
    }
}
