/// Seed clustering and stitching via dynamic programming
use crate::align::score::{AlignmentScorer, GapType};
use crate::align::seed::Seed;
use crate::align::transcript::{CigarOp, Transcript};
use crate::error::Error;
use crate::index::GenomeIndex;

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
///
/// # Returns
/// Vector of seed clusters
pub fn cluster_seeds(
    seeds: &[Seed],
    index: &GenomeIndex,
    max_cluster_dist: u64,
    max_loci_for_anchor: usize,
) -> Vec<SeedCluster> {
    let mut clusters = Vec::new();

    // Find anchor seeds (seeds that map to few genomic locations)
    let mut anchors: Vec<usize> = Vec::new();
    for (i, seed) in seeds.iter().enumerate() {
        let n_loci = seed.sa_end - seed.sa_start;
        if n_loci <= max_loci_for_anchor {
            anchors.push(i);
        }
    }

    // If no anchors, use all seeds as potential anchors
    if anchors.is_empty() {
        anchors = (0..seeds.len()).collect();
    }

    // For each anchor, create clusters
    for &anchor_idx in &anchors {
        let anchor = &seeds[anchor_idx];
        let anchor_positions = anchor.get_genome_positions(index);

        // For each genomic position of the anchor
        for (anchor_pos, anchor_strand) in anchor_positions {
            // Find chromosome
            let chr_info = match index.genome.position_to_chr(anchor_pos) {
                Some(info) => info,
                None => continue, // Position in padding
            };
            let chr_idx = chr_info.0;

            // Collect seeds within clustering window
            let mut seed_indices = vec![anchor_idx];
            let mut genome_start = anchor_pos;
            let mut genome_end = anchor_pos + anchor.length as u64;

            for (i, seed) in seeds.iter().enumerate() {
                if i == anchor_idx {
                    continue;
                }

                // Check if seed overlaps with window
                let seed_positions = seed.get_genome_positions(index);
                for (pos, strand) in seed_positions {
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
                        break; // Only add each seed once per cluster
                    }
                }
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
    let mut expanded_seeds = Vec::new();

    for &seed_idx in &cluster.seed_indices {
        let seed = &seeds[seed_idx];
        let positions = seed.get_genome_positions(index);

        for (pos, strand) in positions {
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

            // Check if within cluster bounds
            let seed_end = pos + seed.length as u64;
            if pos >= cluster.genome_start && seed_end <= cluster.genome_end + 1000000 {
                // Allow some slack
                expanded_seeds.push(ExpandedSeed {
                    read_pos: seed.read_pos,
                    read_end: seed.read_pos + seed.length,
                    genome_pos: pos,
                    genome_end: seed_end,
                    length: seed.length,
                });
            }
        }
    }

    if expanded_seeds.is_empty() {
        return Ok(Vec::new());
    }

    // Sort by read position
    expanded_seeds.sort_by_key(|s| s.read_pos);

    // Initialize DP: one state per expanded seed
    let n = expanded_seeds.len();
    let mut dp: Vec<DpState> = Vec::with_capacity(n);

    // Base case: first seed
    for exp_seed in &expanded_seeds {
        let initial_cigar = vec![CigarOp::Match(exp_seed.length as u32)];
        dp.push(DpState {
            score: 0, // Initial seed has no gap penalty
            prev_seed: None,
            genome_pos: exp_seed.genome_pos,
            cigar_ops: initial_cigar,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
        });
    }

    // DP: for each seed, try connecting to all compatible previous seeds
    for i in 1..n {
        let curr = &expanded_seeds[i];
        let mut best_score = dp[i].score;
        let mut best_prev = dp[i].prev_seed;
        let mut best_cigar = dp[i].cigar_ops.clone();
        let mut best_n_gap = dp[i].n_gap;
        let mut best_n_junction = dp[i].n_junction;

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

            let transition_score = dp[j].score + gap_score;

            if transition_score > best_score {
                best_score = transition_score;
                best_prev = Some(j);

                // Build CIGAR: previous CIGAR + gap + current match
                let mut cigar = dp[j].cigar_ops.clone();

                // Add gap operation
                match gap_type {
                    GapType::Insertion(len) => {
                        cigar.push(CigarOp::Ins(len));
                    }
                    GapType::Deletion(len) => {
                        if len > 0 {
                            cigar.push(CigarOp::Del(len));
                        }
                    }
                    GapType::SpliceJunction { intron_len, .. } => {
                        cigar.push(CigarOp::RefSkip(intron_len));
                    }
                }

                // Add current seed match
                cigar.push(CigarOp::Match(curr.length as u32));

                best_cigar = cigar;

                // Update gap and junction counts
                best_n_gap = dp[j].n_gap;
                best_n_junction = dp[j].n_junction;

                match gap_type {
                    GapType::Insertion(_) | GapType::Deletion(_) => best_n_gap += 1,
                    GapType::SpliceJunction { .. } => best_n_junction += 1,
                }
            }
        }

        // Update DP state
        dp[i].score = best_score;
        dp[i].prev_seed = best_prev;
        dp[i].genome_pos = curr.genome_pos;
        dp[i].cigar_ops = best_cigar;
        dp[i].n_mismatch = 0; // TODO: count mismatches in gaps
        dp[i].n_gap = best_n_gap;
        dp[i].n_junction = best_n_junction;
    }

    // Backtrack from best final state
    let best_final_idx = dp
        .iter()
        .enumerate()
        .max_by_key(|(_, state)| state.score)
        .map(|(i, _)| i)
        .unwrap();

    let best_state = &dp[best_final_idx];
    let final_seed = &expanded_seeds[best_final_idx];

    // Build transcript
    let transcript = Transcript {
        chr_idx: cluster.chr_idx,
        genome_start: cluster.genome_start,
        genome_end: final_seed.genome_end,
        is_reverse: cluster.is_reverse,
        exons: vec![], // TODO: build exons from CIGAR
        cigar: best_state.cigar_ops.clone(),
        score: best_state.score,
        n_mismatch: best_state.n_mismatch,
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
            },
            Seed {
                read_pos: 10,
                length: 5,
                sa_start: 0,
                sa_end: 0, // Empty range
                is_reverse: false,
            },
        ];

        let clusters = cluster_seeds(&seeds, &index, 100000, 10);

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
