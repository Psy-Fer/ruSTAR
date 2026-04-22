//! STAR-faithful splice junction insertion into the genome index.
//!
//! Ports `source/sjdbPrepare.cpp` + `source/sjdbBuildIndex.cpp` from STAR.
//! At `genomeGenerate` time, STAR extracts the flanking `sjdbOverhang`
//! bases on each side of every GTF-derived splice junction, concatenates
//! them into a `Gsj` buffer, appends that buffer to the `Genome` binary,
//! and extends the suffix array to index the new bases. This module
//! provides the same machinery for ruSTAR so that the generated
//! `Genome` / `SA` / `SAindex` / `sjdbInfo.txt` / `sjdbList.out.tab` files
//! match STAR's byte-for-byte.
//!
//! The module is orchestrated from `index::GenomeIndex::build` after the
//! base suffix array has been built.

/// STAR's `sjdbMotif` encoding (see `source/sjdbPrepare.cpp:37-50`).
///
/// `0` = non-canonical / unknown.
/// Values 1-6 correspond to: GT/AG(+), CT/AC(-), GC/AG(+), CT/GC(-),
/// AT/AC(+), GT/AT(-), where each pair is the dinucleotide at the donor /
/// acceptor side of the intron (base-encoded A=0, C=1, G=2, T=3).
pub const MOTIF_NON_CANONICAL: u8 = 0;
pub const MOTIF_GT_AG: u8 = 1;
pub const MOTIF_CT_AC: u8 = 2;
pub const MOTIF_GC_AG: u8 = 3;
pub const MOTIF_CT_GC: u8 = 4;
pub const MOTIF_AT_AC: u8 = 5;
pub const MOTIF_GT_AT: u8 = 6;

/// Classify the splice motif at an intron's donor/acceptor boundary.
///
/// `s` is the 0-based genome-absolute position of the intron's FIRST base
/// (the base immediately after the donor exon's last base). `e` is the
/// 0-based position of the intron's LAST base (immediately before the
/// acceptor exon's first base). Both must be valid indices into `genome`.
///
/// Returns one of the `MOTIF_*` constants.
pub fn classify_motif(genome: &[u8], s: u64, e: u64) -> u8 {
    let si = s as usize;
    let ei = e as usize;
    if ei + 1 > genome.len() || si + 1 >= genome.len() {
        return MOTIF_NON_CANONICAL;
    }
    let b0 = genome[si]; // intron first base
    let b1 = genome[si + 1];
    let b2 = genome[ei - 1];
    let b3 = genome[ei];

    match (b0, b1, b2, b3) {
        (2, 3, 0, 2) => MOTIF_GT_AG,
        (1, 3, 0, 1) => MOTIF_CT_AC,
        (2, 1, 0, 2) => MOTIF_GC_AG,
        (1, 3, 2, 1) => MOTIF_CT_GC,
        (0, 3, 0, 1) => MOTIF_AT_AC,
        (2, 3, 0, 3) => MOTIF_GT_AT,
        _ => MOTIF_NON_CANONICAL,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Build a helper genome buffer. Positions beyond the written slice
    // default to the placeholder padding byte (5) so motif checks at edges
    // are exercised safely.
    fn make_genome(bases: &[(usize, u8)]) -> Vec<u8> {
        let max_idx = bases.iter().map(|(i, _)| *i).max().unwrap_or(0);
        let mut g = vec![5u8; max_idx + 10];
        for (i, b) in bases {
            g[*i] = *b;
        }
        g
    }

    #[test]
    fn classify_gt_ag_forward() {
        // S..E intron = G T ... A G → motif 1
        let g = make_genome(&[(100, 2), (101, 3), (199, 0), (200, 2)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_GT_AG);
    }

    #[test]
    fn classify_ct_ac_reverse() {
        // C T ... A C → motif 2
        let g = make_genome(&[(100, 1), (101, 3), (199, 0), (200, 1)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_CT_AC);
    }

    #[test]
    fn classify_gc_ag() {
        let g = make_genome(&[(100, 2), (101, 1), (199, 0), (200, 2)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_GC_AG);
    }

    #[test]
    fn classify_ct_gc() {
        let g = make_genome(&[(100, 1), (101, 3), (199, 2), (200, 1)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_CT_GC);
    }

    #[test]
    fn classify_at_ac() {
        let g = make_genome(&[(100, 0), (101, 3), (199, 0), (200, 1)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_AT_AC);
    }

    #[test]
    fn classify_gt_at() {
        let g = make_genome(&[(100, 2), (101, 3), (199, 0), (200, 3)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_GT_AT);
    }

    #[test]
    fn classify_non_canonical() {
        // Non-matching combination → 0
        let g = make_genome(&[(100, 0), (101, 0), (199, 0), (200, 0)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_NON_CANONICAL);
    }

    #[test]
    fn classify_out_of_bounds() {
        let g = vec![5u8; 10];
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_NON_CANONICAL);
    }
}
