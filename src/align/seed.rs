use crate::error::Error;
use crate::index::GenomeIndex;

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
    ///
    /// # Returns
    /// Vector of seeds found in the read
    pub fn find_seeds(
        read_seq: &[u8],
        index: &GenomeIndex,
        min_seed_length: usize,
    ) -> Result<Vec<Seed>, Error> {
        let mut seeds = Vec::new();

        // Search forward strand
        for read_pos in 0..read_seq.len() {
            if let Some(seed) = find_seed_at_position(
                read_seq,
                read_pos,
                index,
                min_seed_length,
                false,
            )? {
                seeds.push(seed);
            }
        }

        Ok(seeds)
    }

    /// Get all genome positions for this seed.
    ///
    /// Expands the SA range to actual genome positions.
    pub fn get_genome_positions(&self, index: &GenomeIndex) -> Vec<(u64, bool)> {
        let mut positions = Vec::new();

        for sa_idx in self.sa_start..self.sa_end {
            let sa_entry = index.suffix_array.get(sa_idx);
            let (pos, is_reverse) = index.suffix_array.decode(sa_entry);
            positions.push((pos, is_reverse));
        }

        positions
    }
}

/// Find a seed starting at a specific position in the read.
fn find_seed_at_position(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    min_seed_length: usize,
    is_reverse: bool,
) -> Result<Option<Seed>, Error> {
    if read_pos >= read_seq.len() {
        return Ok(None);
    }

    // Extract k-mer for SAindex lookup
    let sa_nbases = index.sa_index.nbases as usize;
    let remaining = read_seq.len() - read_pos;

    if remaining < min_seed_length {
        return Ok(None);
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
        return Ok(None); // Skip k-mers containing N
    }

    // Look up in SAindex
    let (sa_pos, is_present) = index.sa_index.lookup(kmer_idx, lookup_len as u32);

    if !is_present {
        return Ok(None);
    }

    // Binary search to find exact range
    let (sa_start, sa_end) = binary_search_sa(
        read_seq,
        read_pos,
        index,
        sa_pos as usize,
        lookup_len,
    )?;

    if sa_start >= sa_end {
        return Ok(None);
    }

    // Extend match as far as possible
    let match_length = extend_match(read_seq, read_pos, index, sa_start)?;

    if match_length >= min_seed_length {
        Ok(Some(Seed {
            read_pos,
            length: match_length,
            sa_start,
            sa_end,
            is_reverse,
        }))
    } else {
        Ok(None)
    }
}

/// Binary search the SA around a starting position to find exact range.
fn binary_search_sa(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    hint_pos: usize,
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
    let mut left = sa_start;
    let mut right = sa_len;

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

/// Extend a match as far as possible.
fn extend_match(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    sa_idx: usize,
) -> Result<usize, Error> {
    let sa_entry = index.suffix_array.get(sa_idx);
    let (genome_pos, is_reverse) = index.suffix_array.decode(sa_entry);

    let genome_start = if is_reverse {
        genome_pos as usize + index.genome.n_genome as usize
    } else {
        genome_pos as usize
    };

    let mut length = 0;
    let max_len = read_seq.len() - read_pos;

    for i in 0..max_len {
        let genome_idx = genome_start + i;

        if genome_idx >= index.genome.sequence.len() {
            break;
        }

        let genome_base = index.genome.sequence[genome_idx];

        if genome_base >= 5 || genome_base != read_seq[read_pos + i] {
            break;
        }

        length += 1;
    }

    Ok(length)
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

        let seeds = Seed::find_seeds(&read, &index, 4).unwrap();

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

        // With min_seed_length=4, should find nothing (read is only 3bp)
        let seeds = Seed::find_seeds(&read, &index, 4).unwrap();
        assert!(seeds.is_empty());

        // With min_seed_length=2, should find seeds
        let seeds = Seed::find_seeds(&read, &index, 2).unwrap();
        assert!(!seeds.is_empty());
    }

    #[test]
    fn no_match() {
        let index = make_test_index("ACAC");
        let read = encode_sequence("GGGG");

        let seeds = Seed::find_seeds(&read, &index, 2).unwrap();

        // No seeds should be found (GGGG not in ACAC or its reverse complement GTGT)
        assert!(seeds.is_empty());
    }

    #[test]
    fn get_genome_positions() {
        let index = make_test_index("ACGTACGT");
        let read = encode_sequence("ACGT");

        let seeds = Seed::find_seeds(&read, &index, 4).unwrap();
        assert!(!seeds.is_empty());

        // Get positions for first seed
        let positions = seeds[0].get_genome_positions(&index);
        assert!(!positions.is_empty());

        // Should have at least one valid position
        for (pos, _is_reverse) in positions {
            assert!(pos < index.genome.n_genome);
        }
    }
}
