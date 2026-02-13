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
    /// Max mismatches for stitching SJs: [non-canonical, GT-AG, GC-AG, AT-AC]
    /// -1 means unlimited
    pub align_sj_stitch_mismatch_nmax: [i32; 4],
    /// Max absolute number of mismatches for alignment extension (outFilterMismatchNmax)
    pub n_mm_max: u32,
    /// Max ratio of mismatches to total alignment length (outFilterMismatchNoverLmax)
    pub p_mm_max: f64,
    /// Minimum overhang for splice junctions (alignSJoverhangMin, default 5)
    pub align_sj_overhang_min: u32,
    /// Minimum overhang for annotated splice junctions (alignSJDBoverhangMin, default 3)
    pub align_sjdb_overhang_min: u32,
    /// Maximum intron length (alignIntronMax, 0 = use default 589824)
    pub align_intron_max: u32,
    /// Extra score log-scaled with genomic length: scale * log2(genomicLength)
    pub score_genomic_length_log2_scale: f64,
}

impl AlignmentScorer {
    /// Create a minimal scorer for motif detection only (used in junction recording)
    pub fn from_params_minimal() -> Self {
        Self {
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
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
        }
    }

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
            align_sj_stitch_mismatch_nmax: [
                params.align_sj_stitch_mismatch_nmax[0],
                params.align_sj_stitch_mismatch_nmax[1],
                params.align_sj_stitch_mismatch_nmax[2],
                params.align_sj_stitch_mismatch_nmax[3],
            ],
            n_mm_max: params.out_filter_mismatch_nmax,
            p_mm_max: params.out_filter_mismatch_nover_lmax,
            align_sj_overhang_min: params.align_sj_overhang_min,
            align_sjdb_overhang_min: params.align_sjdb_overhang_min,
            align_intron_max: if params.align_intron_max == 0 {
                params.win_bin_window_dist() as u32
            } else {
                params.align_intron_max
            },
            score_genomic_length_log2_scale: params.score_genomic_length_log2_scale,
        }
    }

    /// Compute genomic length penalty for a transcript.
    /// STAR formula: ceil(log2(genomicLength) * scale - 0.5), clamped so score >= 0.
    pub fn genomic_length_penalty(&self, genomic_span: u64) -> i32 {
        if self.score_genomic_length_log2_scale == 0.0 || genomic_span == 0 {
            return 0;
        }
        ((genomic_span as f64).log2() * self.score_genomic_length_log2_scale - 0.5).ceil() as i32
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

    /// Check if the number of mismatches is allowed for this junction motif type
    pub fn stitch_mismatch_allowed(&self, motif: &SpliceMotif, n_mismatch: u32) -> bool {
        let idx = match motif {
            SpliceMotif::NonCanonical => 0,
            SpliceMotif::GtAg | SpliceMotif::CtAc => 1,
            SpliceMotif::GcAg | SpliceMotif::CtGc => 2,
            SpliceMotif::AtAc | SpliceMotif::GtAt => 3,
        };
        let max_mm = self.align_sj_stitch_mismatch_nmax[idx];
        max_mm < 0 || n_mismatch <= max_mm as u32
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
        self.score_gap_with_strand(genome_gap, read_gap, genome_pos, genome, false, 0)
    }

    /// Score a gap between consecutive seeds, with strand-aware motif detection.
    ///
    /// For reverse-strand reads, `genome_pos` is a raw SA position in the RC genome.
    /// The donor position must be converted to forward genome coordinates for motif
    /// detection: `forward_donor = n_genome - rc_donor - intron_len`.
    pub fn score_gap_with_strand(
        &self,
        genome_gap: i64,
        read_gap: i64,
        genome_pos: u64,
        genome: &Genome,
        is_reverse: bool,
        n_genome: u64,
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
                if len >= self.align_intron_min && len <= self.align_intron_max {
                    // Splice junction — detect motif on forward genome
                    let donor = if is_reverse {
                        n_genome - genome_pos - len as u64
                    } else {
                        genome_pos
                    };
                    let motif = self.detect_splice_motif(donor, len, genome);
                    let score = self.score_splice_junction(&motif);
                    (
                        score,
                        GapType::SpliceJunction {
                            intron_len: len,
                            motif,
                        },
                    )
                } else {
                    // Deletion (too short for intron, or exceeds max intron length)
                    let score = self.score_del_open + self.score_del_base * len as i32;
                    (score, GapType::Deletion(len))
                }
            }
            // Both advance: combined gap (handled by CIGAR builder in stitch.rs)
            // Score the net indel portion
            (gg, rg) if gg > 0 && rg > 0 => {
                let excess = gg - rg;
                if excess > 0 {
                    let del_len = excess as u32;
                    if del_len >= self.align_intron_min && del_len <= self.align_intron_max {
                        let rc_donor = genome_pos + rg as u64;
                        let donor = if is_reverse {
                            n_genome - rc_donor - del_len as u64
                        } else {
                            rc_donor
                        };
                        let motif = self.detect_splice_motif(donor, del_len, genome);
                        let score = self.score_splice_junction(&motif);
                        (
                            score,
                            GapType::SpliceJunction {
                                intron_len: del_len,
                                motif,
                            },
                        )
                    } else {
                        let score = self.score_del_open + self.score_del_base * del_len as i32;
                        (score, GapType::Deletion(del_len))
                    }
                } else if excess < 0 {
                    let ins_len = (-excess) as u32;
                    let score = self.score_ins_open + self.score_ins_base * ins_len as i32;
                    (score, GapType::Insertion(ins_len))
                } else {
                    // Equal gaps: no net indel
                    (0, GapType::Deletion(0))
                }
            }
            // Other cases (negative gaps, etc.)
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
        // Read 2bp donor and 2bp acceptor from the FORWARD genome
        // Donor: donor_pos, donor_pos+1
        // Acceptor: donor_pos+intron_len-2, donor_pos+intron_len-1
        // Always read forward strand — motif pattern determines the strand
        let d1 = genome.get_base(donor_pos);
        let d2 = genome.get_base(donor_pos + 1);
        let a1 = genome.get_base(donor_pos + intron_len as u64 - 2);
        let a2 = genome.get_base(donor_pos + intron_len as u64 - 1);

        // Check if all bases are valid
        // A=0, C=1, G=2, T=3
        match (d1, d2, a1, a2) {
            (Some(d1), Some(d2), Some(a1), Some(a2)) => {
                // Forward-strand motifs
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
                // Reverse-strand motifs (reverse complement on forward genome)
                // CT-AC: (1,3,0,1) — reverse complement of GT-AG
                if d1 == 1 && d2 == 3 && a1 == 0 && a2 == 1 {
                    return SpliceMotif::CtAc;
                }
                // CT-GC: (1,3,2,1) — reverse complement of GC-AG
                if d1 == 1 && d2 == 3 && a1 == 2 && a2 == 1 {
                    return SpliceMotif::CtGc;
                }
                // GT-AT: (2,3,0,3) — reverse complement of AT-AC
                if d1 == 2 && d2 == 3 && a1 == 0 && a2 == 3 {
                    return SpliceMotif::GtAt;
                }
                SpliceMotif::NonCanonical
            }
            _ => SpliceMotif::NonCanonical,
        }
    }

    /// Score a splice junction based on motif
    fn score_splice_junction(&self, motif: &SpliceMotif) -> i32 {
        match motif {
            SpliceMotif::GtAg | SpliceMotif::CtAc => self.score_gap,
            SpliceMotif::GcAg | SpliceMotif::CtGc => self.score_gap_gcag,
            SpliceMotif::AtAc | SpliceMotif::GtAt => self.score_gap_atac,
            SpliceMotif::NonCanonical => self.score_gap_noncan,
        }
    }
}

/// Splice junction motif types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpliceMotif {
    /// GT-AG (canonical, + strand)
    GtAg,
    /// CT-AC (canonical, - strand; reverse complement of GT-AG)
    CtAc,
    /// GC-AG (semi-canonical, + strand)
    GcAg,
    /// CT-GC (semi-canonical, - strand; reverse complement of GC-AG)
    CtGc,
    /// AT-AC (semi-canonical, + strand)
    AtAc,
    /// GT-AT (semi-canonical, - strand; reverse complement of AT-AC)
    GtAt,
    /// Non-canonical
    NonCanonical,
}

impl SpliceMotif {
    /// Get the implied transcript strand from this motif.
    /// Forward-strand motifs (GT-AG, GC-AG, AT-AC) → Some('+')
    /// Reverse-strand motifs (CT-AC, CT-GC, GT-AT) → Some('-')
    /// Non-canonical → None (no strand information)
    pub fn implied_strand(&self) -> Option<char> {
        match self {
            SpliceMotif::GtAg | SpliceMotif::GcAg | SpliceMotif::AtAc => Some('+'),
            SpliceMotif::CtAc | SpliceMotif::CtGc | SpliceMotif::GtAt => Some('-'),
            SpliceMotif::NonCanonical => None,
        }
    }

    /// Get the motif filter category index for outSJfilter* parameters.
    /// 0 = non-canonical, 1 = GT/AG or CT/AC, 2 = GC/AG or CT/GC, 3 = AT/AC or GT/AT
    pub fn filter_category(&self) -> usize {
        match self {
            SpliceMotif::NonCanonical => 0,
            SpliceMotif::GtAg | SpliceMotif::CtAc => 1,
            SpliceMotif::GcAg | SpliceMotif::CtGc => 2,
            SpliceMotif::AtAc | SpliceMotif::GtAt => 3,
        }
    }

    /// Get the filter category from an encoded motif value (0-6 as in SJ.out.tab).
    /// 0→0 (non-canonical), 1|2→1 (GT/AG family), 3|4→2 (GC/AG family), 5|6→3 (AT/AC family)
    pub fn filter_category_from_encoded(encoded: u8) -> usize {
        match encoded {
            1 | 2 => 1,
            3 | 4 => 2,
            5 | 6 => 3,
            _ => 0, // 0 or any unknown = non-canonical
        }
    }
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
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
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
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
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
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
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
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
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
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
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
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
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
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
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
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
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

    #[test]
    fn test_detect_reverse_complement_motifs() {
        // Test all 3 reverse-complement motifs on the forward genome
        // These appear at minus-strand gene splice sites

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
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
        };

        // CT-AC motif: (1,3,0,1) — reverse complement of GT-AG
        let seq_ctac = vec![
            0, 0, // AA
            1, 3, // CT (donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's
            0, 1, // AC (acceptor)
            0, 0, // AA
        ];
        let genome_ctac = make_test_genome(&seq_ctac);
        let motif = scorer.detect_splice_motif(2, 12, &genome_ctac);
        assert_eq!(motif, SpliceMotif::CtAc);
        // Should score same as canonical GT-AG
        assert_eq!(scorer.score_splice_junction(&motif), 0);

        // CT-GC motif: (1,3,2,1) — reverse complement of GC-AG
        let seq_ctgc = vec![
            0, 0, // AA
            1, 3, // CT (donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's
            2, 1, // GC (acceptor)
            0, 0, // AA
        ];
        let genome_ctgc = make_test_genome(&seq_ctgc);
        let motif = scorer.detect_splice_motif(2, 12, &genome_ctgc);
        assert_eq!(motif, SpliceMotif::CtGc);
        assert_eq!(scorer.score_splice_junction(&motif), -4);

        // GT-AT motif: (2,3,0,3) — reverse complement of AT-AC
        let seq_gtat = vec![
            0, 0, // AA
            2, 3, // GT (donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's
            0, 3, // AT (acceptor)
            0, 0, // AA
        ];
        let genome_gtat = make_test_genome(&seq_gtat);
        let motif = scorer.detect_splice_motif(2, 12, &genome_gtat);
        assert_eq!(motif, SpliceMotif::GtAt);
        assert_eq!(scorer.score_splice_junction(&motif), -8);
    }

    #[test]
    fn test_align_intron_max_default() {
        // When align_intron_max is 0, from_params should resolve to 589824
        use clap::Parser;
        let params = crate::params::Parameters::try_parse_from(vec!["ruSTAR"]).unwrap();
        assert_eq!(params.align_intron_max, 0);
        let scorer = AlignmentScorer::from_params(&params);
        assert_eq!(scorer.align_intron_max, 589_824);
    }

    #[test]
    fn test_align_intron_max_custom() {
        // Custom alignIntronMax should be passed through directly
        use clap::Parser;
        let params =
            crate::params::Parameters::try_parse_from(vec!["ruSTAR", "--alignIntronMax", "100000"])
                .unwrap();
        assert_eq!(params.align_intron_max, 100_000);
        let scorer = AlignmentScorer::from_params(&params);
        assert_eq!(scorer.align_intron_max, 100_000);
    }

    #[test]
    fn test_gap_at_intron_max_is_splice_junction() {
        // A gap exactly at alignIntronMax should still be a splice junction
        // Create genome large enough with GT-AG motif at boundaries
        let mut seq = vec![0u8; 600_000]; // ~600kb genome
        // Place GT at position 100
        seq[100] = 2; // G
        seq[101] = 3; // T
        // Place AG at position 100 + 589824 - 2, 100 + 589824 - 1
        let acceptor_pos = 100 + 589_824 - 2;
        if acceptor_pos + 1 < seq.len() {
            seq[acceptor_pos] = 0; // A
            seq[acceptor_pos + 1] = 2; // G
        }
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
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
        };

        // Gap of exactly 589824 starting at position 100 should be splice junction
        let (score, gap_type) = scorer.score_gap(589_824, 0, 100, &genome);
        assert_eq!(score, 0); // GT-AG canonical
        assert!(matches!(
            gap_type,
            GapType::SpliceJunction {
                intron_len: 589_824,
                ..
            }
        ));
    }

    #[test]
    fn test_implied_strand() {
        // Forward-strand motifs
        assert_eq!(SpliceMotif::GtAg.implied_strand(), Some('+'));
        assert_eq!(SpliceMotif::GcAg.implied_strand(), Some('+'));
        assert_eq!(SpliceMotif::AtAc.implied_strand(), Some('+'));
        // Reverse-strand motifs
        assert_eq!(SpliceMotif::CtAc.implied_strand(), Some('-'));
        assert_eq!(SpliceMotif::CtGc.implied_strand(), Some('-'));
        assert_eq!(SpliceMotif::GtAt.implied_strand(), Some('-'));
        // Non-canonical
        assert_eq!(SpliceMotif::NonCanonical.implied_strand(), None);
    }

    #[test]
    fn test_filter_category() {
        assert_eq!(SpliceMotif::NonCanonical.filter_category(), 0);
        assert_eq!(SpliceMotif::GtAg.filter_category(), 1);
        assert_eq!(SpliceMotif::CtAc.filter_category(), 1);
        assert_eq!(SpliceMotif::GcAg.filter_category(), 2);
        assert_eq!(SpliceMotif::CtGc.filter_category(), 2);
        assert_eq!(SpliceMotif::AtAc.filter_category(), 3);
        assert_eq!(SpliceMotif::GtAt.filter_category(), 3);
    }

    #[test]
    fn test_filter_category_from_encoded() {
        assert_eq!(SpliceMotif::filter_category_from_encoded(0), 0); // non-canonical
        assert_eq!(SpliceMotif::filter_category_from_encoded(1), 1); // GT/AG
        assert_eq!(SpliceMotif::filter_category_from_encoded(2), 1); // CT/AC
        assert_eq!(SpliceMotif::filter_category_from_encoded(3), 2); // GC/AG
        assert_eq!(SpliceMotif::filter_category_from_encoded(4), 2); // CT/GC
        assert_eq!(SpliceMotif::filter_category_from_encoded(5), 3); // AT/AC
        assert_eq!(SpliceMotif::filter_category_from_encoded(6), 3); // GT/AT
        assert_eq!(SpliceMotif::filter_category_from_encoded(7), 0); // unknown → non-canonical
    }

    #[test]
    fn test_gap_exceeding_intron_max_is_deletion() {
        // A gap exceeding alignIntronMax should be treated as deletion
        let seq = vec![0u8; 100]; // Small genome, gap length doesn't need real sequence
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
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 1000, // Small max for testing
            score_genomic_length_log2_scale: -0.25,
        };

        // Gap of 1001 (> 1000 max) should be deletion, not splice junction
        let (_score, gap_type) = scorer.score_gap(1001, 0, 0, &genome);
        assert!(matches!(gap_type, GapType::Deletion(1001)));

        // Gap of 1000 (== max) should still be splice junction
        let (_score, gap_type) = scorer.score_gap(1000, 0, 0, &genome);
        assert!(matches!(
            gap_type,
            GapType::SpliceJunction {
                intron_len: 1000,
                ..
            }
        ));

        // Gap of 21 (== min) should be splice junction
        let (_score, gap_type) = scorer.score_gap(21, 0, 0, &genome);
        assert!(matches!(
            gap_type,
            GapType::SpliceJunction { intron_len: 21, .. }
        ));

        // Gap of 20 (< min) should be deletion
        let (_score, gap_type) = scorer.score_gap(20, 0, 0, &genome);
        assert!(matches!(gap_type, GapType::Deletion(20)));
    }
}
