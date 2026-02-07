/// Scoring functions for alignment gaps and splice junctions
use crate::genome::Genome;
use crate::params::Parameters;

/// Alignment scorer with user-defined penalties
#[derive(Debug, Clone)]
pub struct AlignmentScorer {
    /// Canonical splice junction penalty (GT-AG)
    pub score_gap: i32,
    /// Non-canonical splice junction penalty
    pub score_gap_noncan: i32,
    /// GC-AG splice junction penalty
    pub score_gap_gcag: i32,
    /// AT-AC splice junction penalty
    pub score_gap_atac: i32,
    /// Deletion open penalty
    pub score_del_open: i32,
    /// Deletion extension penalty (per base)
    pub score_del_base: i32,
    /// Insertion open penalty
    pub score_ins_open: i32,
    /// Insertion extension penalty (per base)
    pub score_ins_base: i32,
    /// Minimum intron length (gaps >= this are treated as splice junctions)
    pub align_intron_min: u32,
    /// Bonus for annotated splice junctions (from GTF)
    pub sjdb_score: i32,
}

impl AlignmentScorer {
    /// Create scorer from parameters
    pub fn from_params(params: &Parameters) -> Self {
        Self {
            score_gap: params.score_gap,
            score_gap_noncan: params.score_gap_noncan,
            score_gap_gcag: params.score_gap_gcag,
            score_gap_atac: params.score_gap_atac,
            score_del_open: params.score_del_open,
            score_del_base: params.score_del_base,
            score_ins_open: params.score_ins_open,
            score_ins_base: params.score_ins_base,
            align_intron_min: params.align_intron_min,
            sjdb_score: params.sjdb_score,
        }
    }

    /// Apply annotation bonus to junction score
    ///
    /// # Arguments
    /// * `base_score` - Base score from motif
    /// * `annotated` - Whether the junction is annotated in GTF
    ///
    /// # Returns
    /// Adjusted score with annotation bonus applied
    pub fn score_annotated_junction(&self, base_score: i32, annotated: bool) -> i32 {
        if annotated {
            base_score + self.sjdb_score
        } else {
            base_score
        }
    }

    /// Score a gap between two aligned regions
    ///
    /// # Arguments
    /// - `genome_gap`: Gap in genome coordinates (can be negative for insertions)
    /// - `read_gap`: Gap in read coordinates
    /// - `genome_pos`: Starting genomic position (for motif detection)
    /// - `genome`: Genome reference (for motif detection)
    ///
    /// # Returns
    /// `(score, gap_type)`
    pub fn score_gap(
        &self,
        genome_gap: i64,
        read_gap: i64,
        genome_pos: u64,
        genome: &Genome,
    ) -> (i32, GapType) {
        match (genome_gap, read_gap) {
            // Insertion: read advances but genome doesn't
            (0, rg) if rg > 0 => {
                let len = rg as u32;
                let score = self.score_ins_open + self.score_ins_base * len as i32;
                (score, GapType::Insertion(len))
            }
            // Deletion or splice junction: genome advances but read doesn't
            (gg, 0) if gg > 0 => {
                let len = gg as u32;
                if len >= self.align_intron_min {
                    // Splice junction
                    let motif = self.detect_splice_motif(genome_pos, len, genome);
                    let score = self.score_splice_junction(&motif);
                    (
                        score,
                        GapType::SpliceJunction {
                            intron_len: len,
                            motif,
                        },
                    )
                } else {
                    // Deletion
                    let score = self.score_del_open + self.score_del_base * len as i32;
                    (score, GapType::Deletion(len))
                }
            }
            // Both advance: not a simple gap (mismatches/matches)
            _ => (0, GapType::Deletion(0)),
        }
    }

    /// Detect splice junction motif
    ///
    /// # Arguments
    /// - `donor_pos`: Position of the donor site (first base after exon)
    /// - `intron_len`: Length of the intron
    /// - `genome`: Genome reference
    ///
    /// # Returns
    /// The detected splice motif
    pub fn detect_splice_motif(
        &self,
        donor_pos: u64,
        intron_len: u32,
        genome: &Genome,
    ) -> SpliceMotif {
        // Read 2bp donor and 2bp acceptor
        // Donor: donor_pos, donor_pos+1
        // Acceptor: donor_pos+intron_len-2, donor_pos+intron_len-1
        let d1 = genome.get_base(donor_pos);
        let d2 = genome.get_base(donor_pos + 1);
        let a1 = genome.get_base(donor_pos + intron_len as u64 - 2);
        let a2 = genome.get_base(donor_pos + intron_len as u64 - 1);

        // Check if all bases are valid
        match (d1, d2, a1, a2) {
            (Some(d1), Some(d2), Some(a1), Some(a2)) => {
                // GT-AG: (2,3,0,2)
                if d1 == 2 && d2 == 3 && a1 == 0 && a2 == 2 {
                    return SpliceMotif::GtAg;
                }
                // GC-AG: (2,1,0,2)
                if d1 == 2 && d2 == 1 && a1 == 0 && a2 == 2 {
                    return SpliceMotif::GcAg;
                }
                // AT-AC: (0,3,0,1)
                if d1 == 0 && d2 == 3 && a1 == 0 && a2 == 1 {
                    return SpliceMotif::AtAc;
                }
                SpliceMotif::NonCanonical
            }
            _ => SpliceMotif::NonCanonical,
        }
    }

    /// Score a splice junction based on motif
    fn score_splice_junction(&self, motif: &SpliceMotif) -> i32 {
        match motif {
            SpliceMotif::GtAg => self.score_gap,
            SpliceMotif::GcAg => self.score_gap_gcag,
            SpliceMotif::AtAc => self.score_gap_atac,
            SpliceMotif::NonCanonical => self.score_gap_noncan,
        }
    }
}

/// Splice junction motif types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpliceMotif {
    /// GT-AG (canonical)
    GtAg,
    /// GC-AG (semi-canonical)
    GcAg,
    /// AT-AC (semi-canonical)
    AtAc,
    /// Non-canonical
    NonCanonical,
}

/// Gap type between aligned regions
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GapType {
    /// Insertion in read
    Insertion(u32),
    /// Deletion in read
    Deletion(u32),
    /// Splice junction
    SpliceJunction { intron_len: u32, motif: SpliceMotif },
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_genome(seq: &[u8]) -> Genome {
        // Create simple genome with one chromosome
        let n_genome = ((seq.len() as u64 + 1) / 64 + 1) * 64; // Pad to 64-byte boundary
        let mut sequence = vec![5u8; (n_genome * 2) as usize];

        // Copy forward sequence
        sequence[0..seq.len()].copy_from_slice(seq);

        // Build reverse complement
        for i in 0..n_genome as usize {
            let base = sequence[i];
            let complement = if base < 4 { 3 - base } else { base };
            sequence[2 * n_genome as usize - 1 - i] = complement;
        }

        Genome {
            sequence,
            n_genome,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![seq.len() as u64],
            chr_start: vec![0, n_genome],
        }
    }

    #[test]
    fn test_detect_gtag_motif() {
        // Sequence layout (0-based):
        // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
        // A  A [G  T  C  C  C  C  C  C  C  C  A  G] A  A
        //       ^donor                    ^acceptor
        // Intron from position 2 to 14 (exclusive), length = 12
        // Donor: positions 2,3 (GT)
        // Acceptor: positions 12,13 (AG)
        // Bases: A=0, C=1, G=2, T=3
        let seq = vec![
            0, 0, // AA (positions 0-1)
            2, 3, // GT (positions 2-3, donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's (positions 4-11, intron body)
            0, 2, // AG (positions 12-13, acceptor)
            0, 0, // AA (positions 14-15)
        ];
        let genome = make_test_genome(&seq);

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
        };

        // Intron from position 2, length 12 (spans positions 2-13 inclusive)
        let motif = scorer.detect_splice_motif(2, 12, &genome);
        assert_eq!(motif, SpliceMotif::GtAg);

        let score = scorer.score_splice_junction(&motif);
        assert_eq!(score, 0); // Canonical
    }

    #[test]
    fn test_detect_gcag_motif() {
        // GC-AG motif: (2,1,0,2)
        let seq = vec![
            0, 0, // AA
            2, 1, // GC (donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's
            0, 2, // AG (acceptor)
            0, 0, // AA
        ];
        let genome = make_test_genome(&seq);

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
        };

        let motif = scorer.detect_splice_motif(2, 12, &genome);
        assert_eq!(motif, SpliceMotif::GcAg);

        let score = scorer.score_splice_junction(&motif);
        assert_eq!(score, -4);
    }

    #[test]
    fn test_detect_atac_motif() {
        // AT-AC motif: (0,3,0,1)
        let seq = vec![
            0, 0, // AA
            0, 3, // AT (donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's
            0, 1, // AC (acceptor)
            0, 0, // AA
        ];
        let genome = make_test_genome(&seq);

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
        };

        let motif = scorer.detect_splice_motif(2, 12, &genome);
        assert_eq!(motif, SpliceMotif::AtAc);

        let score = scorer.score_splice_junction(&motif);
        assert_eq!(score, -8);
    }

    #[test]
    fn test_detect_noncanonical_motif() {
        // Some random motif
        let seq = vec![
            0, 0, 1, 1, // AA CC (donor)
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // 10 C's
            3, 3, 0, 0, // TT AA (acceptor)
        ];
        let genome = make_test_genome(&seq);

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
        };

        let motif = scorer.detect_splice_motif(2, 12, &genome);
        assert_eq!(motif, SpliceMotif::NonCanonical);

        let score = scorer.score_splice_junction(&motif);
        assert_eq!(score, -8);
    }

    #[test]
    fn test_score_gap_insertion() {
        let genome = make_test_genome(&[0, 1, 2, 3]);
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
        };

        let (score, gap_type) = scorer.score_gap(0, 5, 0, &genome);
        assert_eq!(score, -2 + (-2 * 5)); // open + 5*base
        assert_eq!(gap_type, GapType::Insertion(5));
    }

    #[test]
    fn test_score_gap_deletion() {
        let genome = make_test_genome(&[0, 1, 2, 3]);
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
        };

        // Small gap (< align_intron_min) is deletion
        let (score, gap_type) = scorer.score_gap(10, 0, 0, &genome);
        assert_eq!(score, -2 + (-2 * 10));
        assert_eq!(gap_type, GapType::Deletion(10));
    }

    #[test]
    fn test_score_gap_splice_junction() {
        // Create genome with GT-AG motif
        let seq = vec![
            0, 0, 2, 3, // AA GT (donor)
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // 22 C's
            0, 2, 0, 0, // AG AA (acceptor)
        ];
        let genome = make_test_genome(&seq);

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
        };

        // Gap starting at position 2 (GT), length 26 (>= 21) is splice junction
        let (score, gap_type) = scorer.score_gap(26, 0, 2, &genome);
        assert_eq!(score, 0); // Canonical GT-AG
        assert!(matches!(
            gap_type,
            GapType::SpliceJunction {
                intron_len: 26,
                motif: SpliceMotif::GtAg
            }
        ));
    }

    #[test]
    fn test_annotated_junction_bonus() {
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
        };

        // Annotated junction should get bonus
        let annotated_score = scorer.score_annotated_junction(0, true);
        assert_eq!(annotated_score, 2);

        // Novel junction should not get bonus
        let novel_score = scorer.score_annotated_junction(0, false);
        assert_eq!(novel_score, 0);

        // Bonus applies to any base score
        let annotated_noncanon = scorer.score_annotated_junction(-8, true);
        assert_eq!(annotated_noncanon, -6); // -8 + 2
    }
}
