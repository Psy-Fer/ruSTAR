use crate::error::Error;
use crate::genome::Genome;
use crate::index::packed_array::PackedArray;
use crate::index::suffix_array::SuffixArray;

/// SA index for fast k-mer lookup during binary search.
///
/// For each k-mer prefix (up to genomeSAindexNbases length), stores the
/// SA range where suffixes starting with that prefix can be found.
#[derive(Clone)]
pub struct SaIndex {
    /// Length of indexed k-mers (typically 14)
    pub nbases: u32,

    /// Cumulative counts for each k-mer level
    /// Length = nbases + 1
    /// genomeSAindexStart[k] = total number of k-mers of length <= k
    pub genome_sa_index_start: Vec<u64>,

    /// Packed array storing SA index entries
    /// Each entry is (gstrand_bit + 3) bits:
    /// - Bits 0..gstrand_bit: SA position (or max value if prefix absent)
    /// - Bit gstrand_bit+1: "contains N" flag
    /// - Bit gstrand_bit+2: "prefix absent" flag
    pub data: PackedArray,

    /// Word length for SAindex entries (gstrand_bit + 3)
    pub word_length: u32,

    /// Strand bit position (from SA)
    pub gstrand_bit: u32,
}

impl SaIndex {
    /// Calculate the total number of k-mer prefixes to index.
    ///
    /// Returns sum of 4^1 + 4^2 + ... + 4^nbases = (4^(nbases+1) - 4) / 3
    pub fn calculate_num_indices(nbases: u32) -> u64 {
        if nbases == 0 {
            return 0;
        }
        // Sum of geometric series: (4^(n+1) - 4) / 3
        let power = 4u64.pow(nbases + 1);
        (power - 4) / 3
    }

    /// Build SA index from genome and sorted suffix array.
    ///
    /// # Arguments
    /// * `genome` - The genome sequence
    /// * `sa` - The sorted suffix array
    /// * `nbases` - K-mer length to index (genomeSAindexNbases, typically 14)
    pub fn build(genome: &Genome, sa: &SuffixArray, nbases: u32) -> Result<Self, Error> {
        let gstrand_bit = sa.gstrand_bit;
        let word_length = gstrand_bit + 3; // +2 for flags, +1 for STAR's formula

        // Calculate genomeSAindexStart array
        let mut genome_sa_index_start = vec![0u64; (nbases + 1) as usize];
        genome_sa_index_start[0] = 0;
        for k in 1..=nbases {
            genome_sa_index_start[k as usize] =
                genome_sa_index_start[(k - 1) as usize] + 4u64.pow(k);
        }

        let num_indices = Self::calculate_num_indices(nbases);

        log::info!(
            "Building SA index: nbases={}, num_indices={}, word_length={}",
            nbases,
            num_indices,
            word_length
        );

        // Initialize packed array with "absent" markers
        let mut data = PackedArray::new(word_length, num_indices as usize);
        let absent_marker = (1u64 << (gstrand_bit + 2)) | ((1u64 << gstrand_bit) - 1);

        for i in 0..num_indices as usize {
            data.write(i, absent_marker);
        }

        // Iterate through SA and record first occurrence of each k-mer
        for sa_idx in 0..sa.len() {
            let sa_entry = sa.get(sa_idx);
            let (pos, is_reverse) = sa.decode(sa_entry);

            // Adjust for strand
            let genome_pos = if is_reverse {
                pos as usize + genome.n_genome as usize
            } else {
                pos as usize
            };

            // Extract k-mers of all lengths up to nbases
            for k in 1..=nbases {
                if genome_pos + (k as usize) > genome.sequence.len() {
                    break;
                }

                // Build k-mer index
                let mut kmer_idx = 0u64;
                let mut has_n = false;

                for offset in 0..k {
                    let base = genome.sequence[genome_pos + offset as usize];
                    if base >= 4 {
                        // N or padding
                        has_n = true;
                        break;
                    }
                    kmer_idx = (kmer_idx << 2) | (base as u64);
                }

                if has_n {
                    continue; // Skip k-mers containing N
                }

                // Calculate index in SAindex array
                let sai_pos = genome_sa_index_start[(k - 1) as usize] + kmer_idx;

                // Check if this k-mer hasn't been seen yet
                let current_entry = data.read(sai_pos as usize);
                let is_absent = (current_entry >> (gstrand_bit + 2)) & 1 != 0;

                if is_absent {
                    // Record first occurrence
                    let entry = sa_idx as u64;

                    // Set "contains N" flag if needed (already checked, so clear)
                    // Set "prefix absent" flag to 0 (present)

                    data.write(sai_pos as usize, entry);
                }
            }
        }

        Ok(SaIndex {
            nbases,
            genome_sa_index_start,
            data,
            word_length,
            gstrand_bit,
        })
    }

    /// Get SA range for a k-mer lookup.
    ///
    /// Returns (start_sa_index, is_present) where is_present indicates
    /// whether this k-mer exists in the genome.
    pub fn lookup(&self, kmer_idx: u64, k: u32) -> (u64, bool) {
        if k == 0 || k > self.nbases {
            return (0, false);
        }

        let sai_pos = self.genome_sa_index_start[(k - 1) as usize] + kmer_idx;
        if sai_pos >= self.data.len() as u64 {
            return (0, false);
        }

        let entry = self.data.read(sai_pos as usize);

        // Check "prefix absent" flag (bit gstrand_bit + 2)
        let is_absent = (entry >> (self.gstrand_bit + 2)) & 1 != 0;

        if is_absent {
            return (0, false);
        }

        // Extract SA position
        let sa_pos = entry & ((1u64 << self.gstrand_bit) - 1);
        (sa_pos, true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Parameters;
    use clap::Parser;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_test_index(sequence: &str, bin_nbits: u32, sa_nbases: u32) -> SaIndex {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">test").unwrap();
        writeln!(file, "{}", sequence).unwrap();

        let bin_nbits_str = bin_nbits.to_string();
        let sa_nbases_str = sa_nbases.to_string();
        let args = vec![
            "ruSTAR",
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            file.path().to_str().unwrap(),
            "--genomeChrBinNbits",
            &bin_nbits_str,
            "--genomeSAindexNbases",
            &sa_nbases_str,
        ];

        let params = Parameters::parse_from(args);
        let genome = Genome::from_fasta(&params).unwrap();
        let sa = SuffixArray::build(&genome).unwrap();
        SaIndex::build(&genome, &sa, sa_nbases).unwrap()
    }

    #[test]
    fn calculate_num_indices() {
        // For nbases=1: 4^1 = 4
        assert_eq!(SaIndex::calculate_num_indices(1), 4);

        // For nbases=2: 4^1 + 4^2 = 4 + 16 = 20
        assert_eq!(SaIndex::calculate_num_indices(2), 20);

        // For nbases=3: 4 + 16 + 64 = 84
        assert_eq!(SaIndex::calculate_num_indices(3), 84);

        // For nbases=14 (STAR default): ~357 million
        let n14 = SaIndex::calculate_num_indices(14);
        assert_eq!(n14, 357_913_940);
    }

    #[test]
    fn build_simple_index() {
        let sai = make_test_index("ACGT", 2, 2);

        assert_eq!(sai.nbases, 2);
        assert_eq!(sai.genome_sa_index_start, vec![0, 4, 20]);
        assert!(!sai.data.is_empty());
    }

    #[test]
    fn lookup_present_kmer() {
        let sai = make_test_index("AAAA", 2, 2);

        // K-mer "AA" (00 in 2-bit) should be present
        let kmer_idx = 0b00_00; // AA
        let (sa_pos, is_present) = sai.lookup(kmer_idx, 2);

        assert!(is_present, "K-mer 'AA' should be present in 'AAAA'");
        assert!(sa_pos < sai.data.len() as u64);
    }

    #[test]
    fn lookup_absent_kmer() {
        let sai = make_test_index("ACAC", 2, 2);

        // K-mer "GG" (22 in 2-bit) should be absent from "ACAC" (and its revcomp "GTGT")
        let kmer_idx = 0b10_10; // GG
        let (_sa_pos, is_present) = sai.lookup(kmer_idx, 2);

        assert!(!is_present, "K-mer 'GG' should be absent from 'ACAC'");
    }

    #[test]
    fn genome_sa_index_start_progression() {
        let sai = make_test_index("ACGT", 2, 4);

        // Should be [0, 4, 20, 84, 340]
        assert_eq!(sai.genome_sa_index_start[0], 0);
        assert_eq!(sai.genome_sa_index_start[1], 4); // 4^1
        assert_eq!(sai.genome_sa_index_start[2], 20); // 4 + 4^2
        assert_eq!(sai.genome_sa_index_start[3], 84); // 4 + 16 + 4^3
        assert_eq!(sai.genome_sa_index_start[4], 340); // 4 + 16 + 64 + 4^4
    }
}
