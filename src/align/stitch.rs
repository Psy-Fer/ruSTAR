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
        });
    }

    // DP: for each seed, try connecting to all compatible previous seeds
    // Optimization: find best j first, then build CIGAR only once (avoids repeated cloning)
    for i in 1..n {
        let curr = &expanded_seeds[i];
        let mut best_score = dp[i].score;
        let mut best_j: Option<usize> = None;
        let mut best_gap_score = 0i32;
        let mut best_gap_type = GapType::Deletion(0);
        let mut best_read_gap = 0i64;
        let mut best_genome_gap = 0i64;

        for j in 0..i {
            let prev = &expanded_seeds[j];

            // Check compatibility: no overlap in read, consistent genome order
            if prev.read_end > curr.read_pos {
                continue; // Overlapping in read
            }

            if prev.genome_end > curr.genome_pos && curr.genome_pos > prev.genome_pos {
                continue; // Overlapping in genome (not a clean gap)
            }

            // Calculate gaps
            let read_gap = (curr.read_pos - prev.read_end) as i64;
            let genome_gap = (curr.genome_pos - prev.genome_end) as i64;

            // Score the gap
            let (gap_score, gap_type) =
                scorer.score_gap(genome_gap, read_gap, prev.genome_end, &index.genome);

            // Transition score = prev_score + gap_penalty + current_seed_match_score
            // Current seed contributes +1 per matched base (STAR convention)
            let transition_score = dp[j].score + gap_score + (curr.length as i32);

            if transition_score > best_score {
                best_score = transition_score;
                best_j = Some(j);
                best_gap_score = gap_score;
                best_gap_type = gap_type;
                best_read_gap = read_gap;
                best_genome_gap = genome_gap;
            }
        }

        // Build CIGAR only once for the best transition (instead of cloning per candidate)
        if let Some(j) = best_j {
            let mut cigar = dp[j].cigar_ops.clone();

            // Emit CIGAR operations for the gap between seeds
            // Handle negative gaps explicitly to avoid integer overflow
            let has_gap;

            if best_read_gap == 0 && best_genome_gap == 0 {
                // No gap - seeds are adjacent
                has_gap = false;
            } else if best_read_gap < 0 || best_genome_gap < 0 {
                // Negative gap indicates overlapping seeds or error
                // score_gap() already penalized this, but we skip the connection
                log::warn!(
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
                has_gap = true;
            } else if best_read_gap > 0 && best_genome_gap == 0 {
                // Pure insertion (read advances, genome doesn't)
                let rg = best_read_gap as u32; // Safe: checked > 0
                cigar.push(CigarOp::Ins(rg));
                has_gap = true;
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
                has_gap = true;
            }

            // Add current seed match (always try to merge with previous Match op)
            if let Some(CigarOp::Match(prev_len)) = cigar.last_mut() {
                // Merge with previous Match operation
                *prev_len += curr.length as u32;
            } else {
                // No previous Match, or previous op was not Match
                cigar.push(CigarOp::Match(curr.length as u32));
            }

            let mut n_gap = dp[j].n_gap;
            let mut n_junction = dp[j].n_junction;
            match best_gap_type {
                GapType::Insertion(_) | GapType::Deletion(_) => n_gap += 1,
                GapType::SpliceJunction { .. } => n_junction += 1,
            }
            let _ = best_gap_score;

            dp[i].score = best_score;
            dp[i].prev_seed = Some(j);
            dp[i].genome_pos = curr.genome_pos;
            dp[i].cigar_ops = cigar;
            dp[i].n_mismatch = 0;
            dp[i].n_gap = n_gap;
            dp[i].n_junction = n_junction;
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

    let mut final_cigar = Vec::new();

    // Add 5' soft clip if alignment doesn't start at read position 0
    if alignment_start > 0 {
        final_cigar.push(CigarOp::SoftClip(alignment_start as u32));
    }

    // Add the main alignment CIGAR
    final_cigar.extend(best_state.cigar_ops.clone());

    // Add 3' soft clip if alignment doesn't end at read end
    if alignment_end < read_seq.len() {
        let clip_len = read_seq.len() - alignment_end;
        final_cigar.push(CigarOp::SoftClip(clip_len as u32));
    }

    // Build exons from CIGAR
    use crate::align::transcript::Exon;
    let mut exons = Vec::new();
    let mut read_pos = 0usize; // Start from 0, will be adjusted by soft clips
    let mut genome_pos = chain_start_seed.genome_pos;

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
        .unwrap_or(chain_start_seed.genome_pos);

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
        cigar: final_cigar, // Use CIGAR with soft clips
        score: best_state.score,
        n_mismatch,
        n_gap: best_state.n_gap,
        n_junction: best_state.n_junction,
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
}
