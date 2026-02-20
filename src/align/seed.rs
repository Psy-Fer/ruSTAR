use crate::error::Error;
use crate::index::GenomeIndex;
use crate::io::fastq::complement_base;
use crate::params::Parameters;

/// A seed represents an exact match between a read position and genome location(s).
#[derive(Debug, Clone)]
pub struct Seed {
    /// Position in the read where this seed starts
    pub read_pos: usize,

    /// Length of the exact match
    pub length: usize,

    /// Range in the suffix array [start, end) where this k-mer appears
    pub sa_start: usize,
    pub sa_end: usize,

    /// Whether this seed is on the reverse strand of the read
    pub is_reverse: bool,

    /// Whether this seed was found via R→L (reverse-complement) search.
    /// When true, genome_positions() converts coordinates back to forward orientation.
    pub search_rc: bool,

    /// Mate identifier for paired-end reads
    /// 0 = mate1, 1 = mate2, 2 = single-end (default)
    pub mate_id: u8,
}

impl Seed {
    /// Find all seeds for a read sequence using MMP (Maximal Mappable Prefix) search.
    ///
    /// For each position in the read, performs binary search on the suffix array
    /// to find the longest exact match.
    ///
    /// # Arguments
    /// * `read_seq` - Read sequence (encoded as 0=A, 1=C, 2=G, 3=T)
    /// * `index` - Genome index with SA and SAindex
    /// * `min_seed_length` - Minimum seed length to report (typically 8-20)
    /// * `params` - Parameters including seedMultimapNmax
    ///
    /// # Returns
    /// Vector of seeds found in the read
    pub fn find_seeds(
        read_seq: &[u8],
        index: &GenomeIndex,
        min_seed_length: usize,
        params: &Parameters,
    ) -> Result<Vec<Seed>, Error> {
        let mut seeds = Vec::new();
        let read_len = read_seq.len();

        // Search L→R (forward direction on read)
        for read_pos in 0..read_len {
            let result =
                find_seed_at_position(read_seq, read_pos, index, min_seed_length, false, params)?;
            if let Some(seed) = result.seed {
                seeds.push(seed);
                // Cap total seeds per read (STAR: seedPerReadNmax)
                if seeds.len() >= params.seed_per_read_nmax {
                    return Ok(seeds);
                }
            }
        }

        // Search R→L: reverse-complement the read and run L→R search on it
        let rc_read = reverse_complement_read(read_seq);
        for rc_pos in 0..rc_read.len() {
            let result =
                find_seed_at_position(&rc_read, rc_pos, index, min_seed_length, false, params)?;
            if let Some(mut seed) = result.seed {
                // Convert RC read position back to original read coordinates
                seed.read_pos = read_len - rc_pos - seed.length;
                seed.search_rc = true;
                seeds.push(seed);
                if seeds.len() >= params.seed_per_read_nmax {
                    break;
                }
            }
        }

        Ok(seeds)
    }

    /// Find all seeds for paired-end reads using unified seed pooling.
    ///
    /// This implements STAR's hybrid approach: seeds from both mates are found
    /// independently, tagged with their mate origin, then pooled together for
    /// unified clustering.
    ///
    /// # Arguments
    /// * `mate1_seq` - First mate sequence (encoded)
    /// * `mate2_seq` - Second mate sequence (encoded)
    /// * `index` - Genome index with SA and SAindex
    /// * `min_seed_length` - Minimum seed length to report
    /// * `params` - Parameters including seedMultimapNmax
    ///
    /// # Returns
    /// Vector of seeds from both mates, tagged with mate_id (0 or 1)
    pub fn find_paired_seeds(
        mate1_seq: &[u8],
        mate2_seq: &[u8],
        index: &GenomeIndex,
        min_seed_length: usize,
        params: &Parameters,
    ) -> Result<Vec<Seed>, Error> {
        // Find seeds from mate1 (tag with mate_id = 0)
        let mut seeds = Self::find_seeds(mate1_seq, index, min_seed_length, params)?;
        for seed in &mut seeds {
            seed.mate_id = 0;
        }

        // Find seeds from mate2 (tag with mate_id = 1)
        // IMPORTANT: read_pos is relative to mate2 start (will be adjusted during stitching)
        let mut seeds2 = Self::find_seeds(mate2_seq, index, min_seed_length, params)?;
        for seed in &mut seeds2 {
            seed.mate_id = 1;
        }

        // Pool seeds together
        seeds.extend(seeds2);

        Ok(seeds)
    }

    /// Get all genome positions for this seed.
    ///
    /// Expands the SA range to actual genome positions.
    pub fn get_genome_positions(&self, index: &GenomeIndex) -> Vec<(u64, bool)> {
        self.genome_positions(index).collect()
    }

    /// Iterate over genome positions for this seed without allocating.
    ///
    /// Returns an iterator that lazily decodes SA entries.
    /// For R→L seeds (search_rc == true), converts positions back:
    /// (pos, is_rev) → (n_genome - pos - length, !is_rev)
    /// Positions where the conversion would underflow are filtered out.
    pub fn genome_positions<'a>(
        &'a self,
        index: &'a GenomeIndex,
    ) -> impl Iterator<Item = (u64, bool)> + 'a {
        let search_rc = self.search_rc;
        let length = self.length as u64;
        let n_genome = index.genome.n_genome;
        (self.sa_start..self.sa_end).filter_map(move |sa_idx| {
            let sa_entry = index.suffix_array.get(sa_idx);
            let (pos, is_rev) = index.suffix_array.decode(sa_entry);
            if search_rc {
                if pos + length <= n_genome {
                    Some((n_genome - pos - length, !is_rev))
                } else {
                    None // Position would span past genome boundary
                }
            } else {
                Some((pos, is_rev))
            }
        })
    }
}

/// Reverse-complement an encoded read sequence.
///
/// Reverses the order and complements each base (A↔T, C↔G).
fn reverse_complement_read(read_seq: &[u8]) -> Vec<u8> {
    read_seq.iter().rev().map(|&b| complement_base(b)).collect()
}

/// Result of an MMP (Maximal Mappable Prefix) search at a single position.
/// Always provides the advance length for Lmapped tracking, even when no
/// seed is stored (matching STAR's behavior).
struct MmpResult {
    /// The seed to store, if it passed all filters (multimap, min length)
    seed: Option<Seed>,
    /// MMP length to advance by (>= 1). Used for Lmapped tracking regardless
    /// of whether a seed was stored.
    advance: usize,
}

/// Search one direction using STAR's sparse starting positions with Lmapped tracking.
///
/// Divides read into Nstart evenly-spaced starting positions. From each start,
/// does successive MMP searches forward, advancing past found seeds (Lmapped).
/// Produces ~20-40 seeds per 150bp read (vs ~100+ with dense every-position search)
/// while guaranteeing no gap > seedMapMin (default 5bp).
///
/// NOTE: Bug-fixed but dormant — our DP stitcher needs dense (every-position) seeds.
/// Sparse seeds cause false splices (4.3% vs 2.2%) because mismatch-position seeds
/// map to spurious locations without enough neighboring seeds to vote them down.
/// Activate when DP is adapted for sparse seeds (extension-based gap filling).
#[allow(dead_code)]
fn search_direction_sparse(
    read_seq: &[u8],
    original_read_len: usize,
    index: &GenomeIndex,
    min_seed_length: usize,
    params: &Parameters,
    effective_start_lmax: usize,
    is_rc: bool,
    seeds: &mut Vec<Seed>,
) -> Result<(), Error> {
    let read_len = read_seq.len();

    // Calculate Nstart and Lstart (STAR formula: Nstart = ceil(read_len / effective_start_lmax))
    let nstart = read_len.div_ceil(effective_start_lmax);
    let lstart = read_len / nstart.max(1);

    for istart in 0..nstart {
        let start_pos = (istart * lstart).min(read_len);
        let mut pos = start_pos;

        // From this starting position, search forward with Lmapped tracking.
        // Always search at least once per start position, then use seedMapMin
        // for continuation (STAR default: 5).
        loop {
            if pos >= read_len {
                break;
            }

            let result =
                find_seed_at_position(read_seq, pos, index, min_seed_length, false, params)?;

            if let Some(mut seed) = result.seed {
                // Apply seedSearchLmax cap
                if params.seed_search_lmax > 0 && seed.length > params.seed_search_lmax {
                    seed.length = params.seed_search_lmax;
                }

                seed.search_rc = is_rc;

                // Convert RC read_pos back to original read coordinates
                if is_rc {
                    seed.read_pos = original_read_len - seed.read_pos - seed.length;
                }

                seeds.push(seed);

                if seeds.len() >= params.seed_per_read_nmax {
                    return Ok(());
                }
            }

            pos += result.advance; // Always advance by MMP length (matches STAR)

            // Check if enough remaining read to justify continuing
            if pos + params.seed_map_min >= read_len {
                break;
            }
        }
    }

    Ok(())
}

/// Find a seed starting at a specific position in the read.
///
/// Returns an MmpResult that always provides the MMP advance length for Lmapped
/// tracking, even when no seed is stored. This matches STAR's behavior where
/// `maxMappableLength2strands()` always returns the MMP length, and `Lmapped += L`
/// always advances — regardless of whether the seed passes filters.
fn find_seed_at_position(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    min_seed_length: usize,
    is_reverse: bool,
    params: &Parameters,
) -> Result<MmpResult, Error> {
    if read_pos >= read_seq.len() {
        return Ok(MmpResult {
            seed: None,
            advance: 1,
        });
    }

    // Extract k-mer for SAindex lookup
    let sa_nbases = index.sa_index.nbases as usize;
    let remaining = read_seq.len() - read_pos;

    if remaining < min_seed_length {
        return Ok(MmpResult {
            seed: None,
            advance: 1,
        });
    }

    // Use SAindex to get initial SA range
    let lookup_len = remaining.min(sa_nbases);
    let mut kmer_idx = 0u64;
    let mut has_n = false;

    for i in 0..lookup_len {
        let base = read_seq[read_pos + i];
        if base >= 4 {
            has_n = true;
            break;
        }
        kmer_idx = (kmer_idx << 2) | (base as u64);
    }

    if has_n {
        return Ok(MmpResult {
            seed: None,
            advance: 1,
        }); // N base — no match possible
    }

    // Look up in SAindex
    let (sa_pos, is_present) = index.sa_index.lookup(kmer_idx, lookup_len as u32);

    if !is_present {
        return Ok(MmpResult {
            seed: None,
            advance: 1,
        }); // K-mer not in index
    }

    // Binary search to find exact range
    let (sa_start, sa_end) =
        binary_search_sa(read_seq, read_pos, index, sa_pos as usize, lookup_len)?;

    if sa_start >= sa_end {
        return Ok(MmpResult {
            seed: None,
            advance: 1,
        }); // SA range empty
    }

    // Find maximum mappable prefix length and narrow SA range (STAR's maxMappableLength).
    // Binary searches within the SA range, extending match while narrowing both boundaries.
    let (match_length, narrowed_start, narrowed_end) =
        max_mappable_length(read_seq, read_pos, index, sa_start, sa_end, lookup_len);
    let advance = match_length.max(1);

    // Check seedMultimapNmax: filter seeds that map to too many loci
    // Key fix: still advance by MMP length even when seed is not stored
    // Uses narrowed range (accurate loci count, not overestimated k-mer range)
    let n_loci = narrowed_end - narrowed_start;
    if n_loci > params.seed_multimap_nmax {
        return Ok(MmpResult {
            seed: None,
            advance,
        });
    }

    if match_length >= min_seed_length {
        Ok(MmpResult {
            seed: Some(Seed {
                read_pos,
                length: match_length,
                sa_start: narrowed_start,
                sa_end: narrowed_end,
                is_reverse,
                search_rc: false,
                mate_id: 2, // Single-end default
            }),
            advance,
        })
    } else {
        Ok(MmpResult {
            seed: None,
            advance,
        })
    }
}

/// Binary search the SA around a starting position to find exact range.
fn binary_search_sa(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    _hint_pos: usize,
    min_match: usize,
) -> Result<(usize, usize), Error> {
    let sa_len = index.suffix_array.len();

    // Find lower bound (first position where suffix >= query)
    let mut left = 0;
    let mut right = sa_len;

    while left < right {
        let mid = (left + right) / 2;
        let cmp = compare_suffix_to_query(read_seq, read_pos, index, mid, min_match)?;

        if cmp < 0 {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    let sa_start = left;

    // Find upper bound (first position where suffix > query)
    // Reuse left starting from sa_start
    right = sa_len;

    while left < right {
        let mid = (left + right) / 2;
        let cmp = compare_suffix_to_query(read_seq, read_pos, index, mid, min_match)?;

        if cmp <= 0 {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    let sa_end = left;

    Ok((sa_start, sa_end))
}

/// Compare a genome suffix (from SA) to a read query sequence.
///
/// Returns: -1 if suffix < query, 0 if equal (up to min_match), 1 if suffix > query
fn compare_suffix_to_query(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    sa_idx: usize,
    min_match: usize,
) -> Result<i32, Error> {
    let sa_entry = index.suffix_array.get(sa_idx);
    let (genome_pos, is_reverse) = index.suffix_array.decode(sa_entry);

    let genome_start = if is_reverse {
        genome_pos as usize + index.genome.n_genome as usize
    } else {
        genome_pos as usize
    };

    let max_cmp = min_match.min(read_seq.len() - read_pos);

    for i in 0..max_cmp {
        let read_base = read_seq[read_pos + i];
        let genome_idx = genome_start + i;

        if genome_idx >= index.genome.sequence.len() {
            return Ok(-1); // Genome suffix ended
        }

        let genome_base = index.genome.sequence[genome_idx];

        if genome_base >= 5 {
            return Ok(-1); // Hit padding
        }

        match genome_base.cmp(&read_base) {
            std::cmp::Ordering::Less => return Ok(-1),
            std::cmp::Ordering::Greater => return Ok(1),
            std::cmp::Ordering::Equal => continue,
        }
    }

    Ok(0) // Equal up to min_match
}

/// Overflow-safe median of two unsigned integers.
/// Equivalent to STAR's medianUint2: a/2 + b/2 + (a%2 + b%2)/2
fn median_uint2(a: usize, b: usize) -> usize {
    a / 2 + b / 2 + (a % 2 + b % 2) / 2
}

/// Compare read to genome at a specific SA position, starting from offset l_start.
/// Returns (total_match_length, is_read_greater_at_mismatch).
/// Ports STAR's compareSeqToGenome (SuffixArrayFuns.cpp).
///
/// Starts comparing from offset `l_start` (bases 0..l_start are assumed to match).
/// Walks forward until a mismatch, end of read, or genome padding.
fn compare_seq_to_genome(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    sa_idx: usize,
    l_start: usize,
) -> (usize, bool) {
    let sa_entry = index.suffix_array.get(sa_idx);
    let (genome_pos, is_reverse) = index.suffix_array.decode(sa_entry);

    let genome_start = if is_reverse {
        genome_pos as usize + index.genome.n_genome as usize
    } else {
        genome_pos as usize
    };

    let remaining = read_seq.len() - read_pos;
    let mut match_len = l_start;

    for i in l_start..remaining {
        let genome_idx = genome_start + i;

        if genome_idx >= index.genome.sequence.len() {
            // Past end of genome array — treat like padding (STAR: comp_res > 0)
            return (match_len, true);
        }

        let genome_base = index.genome.sequence[genome_idx];

        if genome_base >= 5 {
            // Padding character — STAR returns comp_res > 0 (read > genome)
            return (match_len, true);
        }

        let read_base = read_seq[read_pos + i];

        if read_base != genome_base {
            return (match_len, read_base > genome_base);
        }

        match_len += 1;
    }

    // Matched all remaining bases — STAR returns comp_res < 0 (genome >= read)
    (match_len, false)
}

/// Find maximum mappable prefix length within SA range [sa_start, sa_end).
/// Binary searches the range while extending match length, then narrows to
/// all positions matching the maximum length.
/// Returns (match_length, narrowed_sa_start, narrowed_sa_end_exclusive).
/// Ports STAR's maxMappableLength (SuffixArrayFuns.cpp).
fn max_mappable_length(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    sa_start: usize,
    sa_end: usize,
    l_initial: usize,
) -> (usize, usize, usize) {
    let remaining = read_seq.len() - read_pos;

    // Single element: just compare
    if sa_start + 1 >= sa_end {
        let (l, _) = compare_seq_to_genome(read_seq, read_pos, index, sa_start, l_initial);
        return (l, sa_start, sa_start + 1);
    }

    // Convert to inclusive range (STAR convention internally)
    let mut i1 = sa_start;
    let mut i2 = sa_end - 1;

    let (mut l1, _) = compare_seq_to_genome(read_seq, read_pos, index, i1, l_initial);
    let (mut l2, _) = compare_seq_to_genome(read_seq, read_pos, index, i2, l_initial);

    let mut l = l1.min(l2);
    let mut l3 = l;
    let mut i3 = i1;

    // Track history for find_mult_range
    let (mut i1a, mut l1a) = (i1, l1);
    let (mut i1b, mut l1b) = (i1, l1);
    let (mut i2a, mut l2a) = (i2, l2);
    let (mut i2b, mut l2b) = (i2, l2);

    // Binary search within SA range
    while i1 + 1 < i2 {
        i3 = median_uint2(i1, i2);
        let comp3;
        (l3, comp3) = compare_seq_to_genome(read_seq, read_pos, index, i3, l);

        if l3 == remaining {
            break; // Perfect match found
        }

        if comp3 {
            // read > genome at mismatch: move left boundary up
            i1a = i1b;
            l1a = l1b;
            i1b = i1;
            l1b = l1;
            i1 = i3;
            l1 = l3;
        } else {
            // read <= genome at mismatch: move right boundary down
            i2a = i2b;
            l2a = l2b;
            i2b = i2;
            l2b = l2;
            i2 = i3;
            l2 = l3;
        }

        l = l1.min(l2);
    }

    // Pick the best match length
    if l3 < remaining {
        if l1 > l2 {
            l3 = l1;
            i3 = i1;
        } else {
            l3 = l2;
            i3 = i2;
        }
    }

    // Find narrowed range using find_mult_range
    let narrowed_start = find_mult_range(
        read_seq, read_pos, index, remaining, i3, l3, i1, l1, i1a, l1a, i1b, l1b,
    );
    let narrowed_end = find_mult_range(
        read_seq, read_pos, index, remaining, i3, l3, i2, l2, i2a, l2a, i2b, l2b,
    );

    // Convert back to exclusive end
    (l3, narrowed_start, narrowed_end + 1)
}

/// Binary search to find the SA boundary where match length transitions
/// from >= l3 to < l3. Used to narrow the SA range to only positions
/// matching the maximum prefix length.
/// Ports STAR's funFindMultRange.
#[allow(clippy::too_many_arguments)]
fn find_mult_range(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    _remaining: usize,
    mut i3: usize,
    l3: usize,
    mut i1: usize,
    l1: usize,
    i1a: usize,
    l1a: usize,
    i1b: usize,
    l1b: usize,
) -> usize {
    if l1 >= l3 {
        return i1; // Boundary already at target length
    }

    // Use closest cached intermediate with L >= l3 to narrow search space
    if l1b >= l3 {
        i3 = i1b;
    } else if l1a >= l3 {
        i3 = i1a;
    }

    // Binary search between i3 (L >= l3) and i1 (L < l3)
    // Start comparison from the minimum known matching length
    let mut l_start = l1;

    while (i3 + 1 < i1) || (i1 + 1 < i3) {
        let i4 = median_uint2(i3, i1);
        let (l4, _) = compare_seq_to_genome(read_seq, read_pos, index, i4, l_start);

        if l4 >= l3 {
            i3 = i4;
        } else {
            i1 = i4;
            l_start = l4; // Update comparison start to shorter match
        }
    }

    i3
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Parameters;
    use clap::Parser;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_test_index(sequence: &str) -> GenomeIndex {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "{}", sequence).unwrap();

        let dir = tempfile::tempdir().unwrap();

        let args = vec![
            "ruSTAR",
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            file.path().to_str().unwrap(),
            "--genomeDir",
            dir.path().to_str().unwrap(),
            "--genomeChrBinNbits",
            "2",
            "--genomeSAindexNbases",
            "2",
        ];

        let params = Parameters::parse_from(args);
        GenomeIndex::build(&params).unwrap()
    }

    fn encode_sequence(seq: &str) -> Vec<u8> {
        seq.bytes()
            .map(|b| match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            })
            .collect()
    }

    #[test]
    fn find_exact_match() {
        let index = make_test_index("ACGTACGT");
        let read = encode_sequence("ACGT");

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params).unwrap();

        // Should find at least one seed
        assert!(!seeds.is_empty());

        // First seed should be at position 0 with length 4
        assert_eq!(seeds[0].read_pos, 0);
        assert_eq!(seeds[0].length, 4);
    }

    #[test]
    fn min_seed_length_filter() {
        let index = make_test_index("AAAAAAAA");
        let read = encode_sequence("AAA");

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        // With min_seed_length=4, should find nothing (read is only 3bp)
        let seeds = Seed::find_seeds(&read, &index, 4, &params).unwrap();
        assert!(seeds.is_empty());

        // With min_seed_length=2, should find seeds
        let seeds = Seed::find_seeds(&read, &index, 2, &params).unwrap();
        assert!(!seeds.is_empty());
    }

    #[test]
    fn no_match() {
        let index = make_test_index("ACAC");
        let read = encode_sequence("GGGG");

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 2, &params).unwrap();

        // No seeds should be found (GGGG not in ACAC or its reverse complement GTGT)
        assert!(seeds.is_empty());
    }

    #[test]
    fn get_genome_positions() {
        let index = make_test_index("ACGTACGT");
        let read = encode_sequence("ACGT");

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params).unwrap();
        assert!(!seeds.is_empty());

        // Get positions for first seed
        let positions = seeds[0].get_genome_positions(&index);
        assert!(!positions.is_empty());

        // Should have at least one valid position
        for (pos, _is_reverse) in positions {
            assert!(pos < index.genome.n_genome);
        }
    }

    #[test]
    fn test_single_end_mate_id() {
        let index = make_test_index("ACGTACGT");
        let read = encode_sequence("ACGT");

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params).unwrap();
        assert!(!seeds.is_empty());

        // Single-end seeds should have mate_id = 2
        for seed in seeds {
            assert_eq!(seed.mate_id, 2);
        }
    }

    #[test]
    fn test_find_paired_seeds() {
        let index = make_test_index("ACGTACGTTTGGCCAA");
        let mate1 = encode_sequence("ACGT");
        let mate2 = encode_sequence("TTGG");

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_paired_seeds(&mate1, &mate2, &index, 4, &params).unwrap();

        // Should have seeds from both mates
        let mate1_seeds: Vec<_> = seeds.iter().filter(|s| s.mate_id == 0).collect();
        let mate2_seeds: Vec<_> = seeds.iter().filter(|s| s.mate_id == 1).collect();

        assert!(!mate1_seeds.is_empty(), "Should have mate1 seeds");
        assert!(!mate2_seeds.is_empty(), "Should have mate2 seeds");

        // Verify mate1 seeds have correct read positions
        for seed in mate1_seeds {
            assert!(seed.read_pos < mate1.len());
        }

        // Verify mate2 seeds have correct read positions (relative to mate2)
        for seed in mate2_seeds {
            assert!(seed.read_pos < mate2.len());
        }
    }

    #[test]
    fn test_paired_seeds_pooling() {
        let index = make_test_index("ACGTACGT");
        let mate1 = encode_sequence("ACGT");
        let mate2 = encode_sequence("ACGT");

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_paired_seeds(&mate1, &mate2, &index, 4, &params).unwrap();

        // Should have roughly double the seeds (one set from each mate)
        let mate1_count = seeds.iter().filter(|s| s.mate_id == 0).count();
        let mate2_count = seeds.iter().filter(|s| s.mate_id == 1).count();

        assert!(mate1_count > 0);
        assert!(mate2_count > 0);
        assert_eq!(seeds.len(), mate1_count + mate2_count);
    }

    #[test]
    fn test_reverse_complement_read() {
        // ACGT → RC = ACGT (palindrome)
        let read = encode_sequence("ACGT");
        let rc = reverse_complement_read(&read);
        assert_eq!(rc, encode_sequence("ACGT"));

        // AACC → RC = GGTT
        let read2 = encode_sequence("AACC");
        let rc2 = reverse_complement_read(&read2);
        assert_eq!(rc2, encode_sequence("GGTT"));

        // Single base
        let read3 = encode_sequence("A");
        let rc3 = reverse_complement_read(&read3);
        assert_eq!(rc3, encode_sequence("T"));

        // N bases preserved
        let read4 = vec![0, 4, 1]; // A, N, C
        let rc4 = reverse_complement_read(&read4);
        assert_eq!(rc4, vec![2, 4, 3]); // G, N, T
    }

    #[test]
    fn test_rl_seeds_found() {
        // Genome has ACGTACGT. The RC of that is ACGTACGT (palindrome),
        // so L→R already finds everything. Use an asymmetric sequence instead.
        // Genome: AACCGGTT — RC genome half has AACCGGTT too.
        // Read: CCGG — L→R finds it at pos 2 in genome.
        // RC of read: CCGG — R→L also finds it.
        // So with this palindromic example, R→L seeds duplicate L→R.
        // Instead use: Genome = AACCTTGG, Read = CCAAGGTT (= RC of AACCTTGG)
        // The read itself won't match L→R in forward genome, but its RC (AACCTTGG) will.
        let index = make_test_index("AACCTTGG");
        // Read is RC of genome: CCAAGGTT
        let read = encode_sequence("CCAAGGTT");

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params).unwrap();

        // Should have R→L seeds (search_rc == true)
        let rc_seeds: Vec<_> = seeds.iter().filter(|s| s.search_rc).collect();
        let lr_seeds: Vec<_> = seeds.iter().filter(|s| !s.search_rc).collect();

        // R→L search should find seeds because RC(read) = AACCTTGG matches genome
        assert!(
            !rc_seeds.is_empty(),
            "R→L search should find seeds (RC of read matches genome). All seeds: {:?}",
            seeds
        );

        // Verify R→L seeds have valid read positions
        for seed in &rc_seeds {
            assert!(
                seed.read_pos + seed.length <= read.len(),
                "R→L seed read_pos {} + length {} exceeds read length {}",
                seed.read_pos,
                seed.length,
                read.len()
            );
        }

        // Total seeds should be more than just L→R
        assert!(
            seeds.len() > lr_seeds.len(),
            "Total seeds ({}) should exceed L→R seeds ({})",
            seeds.len(),
            lr_seeds.len()
        );
    }

    #[test]
    fn test_shared_seed_cap() {
        // Test that combined L→R + R→L respects seedPerReadNmax
        let index = make_test_index("ACGTACGTACGTACGT");
        let read = encode_sequence("ACGTACGT");

        let args = vec![
            "ruSTAR",
            "--runMode",
            "alignReads",
            "--seedPerReadNmax",
            "3",
        ];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params).unwrap();
        assert!(
            seeds.len() <= 3,
            "Total seeds ({}) should respect seedPerReadNmax=3",
            seeds.len()
        );
    }

    #[test]
    fn test_sparse_nstart_calculation() {
        // Verify Nstart = ceil(read_len / effective_start_lmax)
        // 150 / 50 = 3 exact → Nstart=3, Lstart=50
        assert_eq!(150_usize.div_ceil(50), 3);
        assert_eq!(150 / 3_usize.max(1), 50);

        // 151 / 50 = 3.02 → Nstart=4, Lstart=37
        assert_eq!(151_usize.div_ceil(50), 4);
        assert_eq!(151 / 4_usize.max(1), 37);

        // 30 / 50 = 0.6 → Nstart=1, Lstart=30
        assert_eq!(30_usize.div_ceil(50), 1);
        assert_eq!(30 / 1_usize.max(1), 30);

        // Edge: 50 / 50 = 1 → Nstart=1, Lstart=50
        assert_eq!(50_usize.div_ceil(50), 1);
        assert_eq!(50 / 1_usize.max(1), 50);

        // 100 / 50 = 2 → Nstart=2, Lstart=50
        assert_eq!(100_usize.div_ceil(50), 2);
        assert_eq!(100 / 2_usize.max(1), 50);
    }

    #[test]
    fn test_sparse_rc_read_pos_conversion() {
        // Genome: AACCTTGG, read is RC of genome: CCAAGGTT
        // RC(read) = AACCTTGG matches forward genome at pos 0
        // When searching R→L (is_rc=true), seeds found in the RC read have
        // positions relative to the RC read. They must be converted back to
        // original read coordinates: read_pos = original_read_len - rc_pos - length
        let index = make_test_index("AACCTTGG");
        let read = encode_sequence("CCAAGGTT");

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params).unwrap();

        for seed in &seeds {
            // All seeds (L→R and R→L) must have valid read positions
            assert!(
                seed.read_pos + seed.length <= read.len(),
                "Seed read_pos {} + length {} exceeds read len {} (search_rc={})",
                seed.read_pos,
                seed.length,
                read.len(),
                seed.search_rc
            );
        }

        // Should have R→L seeds
        let rc_seeds: Vec<_> = seeds.iter().filter(|s| s.search_rc).collect();
        assert!(
            !rc_seeds.is_empty(),
            "Should have R→L seeds with valid read positions"
        );
    }

    #[test]
    fn test_sparse_fewer_seeds_than_dense() {
        // With a longer genome and read, sparse search should produce fewer seeds
        // than the old dense (every-position) search
        let genome_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let index = make_test_index(genome_seq);

        // Use a read that's long enough for multiple start positions
        let read = encode_sequence("ACGTACGTACGTACGTACGTACGT"); // 24bp

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let sparse_seeds = Seed::find_seeds(&read, &index, 4, &params).unwrap();

        // Count how many seeds dense would produce (every position that has a match)
        let mut dense_count = 0;
        for read_pos in 0..read.len() {
            let result = find_seed_at_position(&read, read_pos, &index, 4, false, &params).unwrap();
            if result.seed.is_some() {
                dense_count += 1;
            }
        }
        // Also count R→L dense seeds
        let rc_read = reverse_complement_read(&read);
        for rc_pos in 0..rc_read.len() {
            let result =
                find_seed_at_position(&rc_read, rc_pos, &index, 4, false, &params).unwrap();
            if result.seed.is_some() {
                dense_count += 1;
            }
        }

        assert!(
            sparse_seeds.len() <= dense_count,
            "Sparse ({}) should produce <= dense ({}) seeds",
            sparse_seeds.len(),
            dense_count
        );
    }

    #[test]
    fn test_rc_seed_genome_positions() {
        // Genome: AACCTTGG, read RC = CCAAGGTT
        // RC(read) = AACCTTGG matches forward genome
        let index = make_test_index("AACCTTGG");
        let read = encode_sequence("CCAAGGTT");

        let args = vec!["ruSTAR", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params).unwrap();

        for seed in &seeds {
            if seed.search_rc {
                let positions: Vec<_> = seed.genome_positions(&index).collect();
                assert!(
                    !positions.is_empty(),
                    "RC seed should have genome positions"
                );

                for (pos, _is_rev) in &positions {
                    // Converted positions should be valid (within genome)
                    assert!(
                        *pos < index.genome.n_genome,
                        "Converted position {} should be < n_genome {}",
                        pos,
                        index.genome.n_genome
                    );
                }
            }
        }
    }
}
