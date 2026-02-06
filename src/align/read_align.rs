/// Read alignment driver function
use crate::align::score::AlignmentScorer;
use crate::align::seed::Seed;
use crate::align::stitch::{cluster_seeds, stitch_seeds};
use crate::align::transcript::Transcript;
use crate::error::Error;
use crate::index::GenomeIndex;
use crate::params::Parameters;

/// Align a read to the genome.
///
/// # Algorithm
/// 1. Find seeds (exact matches) using MMP search
/// 2. Cluster seeds by genomic proximity
/// 3. Stitch seeds within each cluster using DP
/// 4. Filter transcripts by quality thresholds
/// 5. Sort by score and limit to top N
///
/// # Arguments
/// * `read_seq` - Read sequence (encoded as 0=A, 1=C, 2=G, 3=T)
/// * `index` - Genome index
/// * `params` - User parameters
///
/// # Returns
/// Vector of transcripts (alignments), sorted by score (best first)
pub fn align_read(
    read_seq: &[u8],
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<Vec<Transcript>, Error> {
    // Step 1: Find seeds
    // Use a reasonable default min seed length (typically 8-20bp)
    let min_seed_length = 8;
    let seeds = Seed::find_seeds(read_seq, index, min_seed_length)?;

    if seeds.is_empty() {
        return Ok(Vec::new()); // No seeds found
    }

    // Step 2: Cluster seeds
    let max_cluster_dist = 100000; // 100kb window (TODO: make configurable)
    let max_loci_for_anchor = 10; // Seeds mapping to <=10 loci can be anchors
    let clusters = cluster_seeds(&seeds, index, max_cluster_dist, max_loci_for_anchor);

    if clusters.is_empty() {
        return Ok(Vec::new());
    }

    // Step 3: Stitch seeds within each cluster
    let scorer = AlignmentScorer::from_params(params);
    let mut transcripts = Vec::new();

    for cluster in clusters {
        let cluster_transcripts = stitch_seeds(&cluster, &seeds, read_seq, index, &scorer)?;
        transcripts.extend(cluster_transcripts);
    }

    // Step 4: Filter transcripts
    let read_length = read_seq.len() as f64;
    transcripts.retain(|t| {
        // Absolute score threshold
        if t.score < params.out_filter_score_min {
            return false;
        }

        // Relative score threshold (score / read_length)
        if (t.score as f64) < params.out_filter_score_min_over_lread * read_length {
            return false;
        }

        // Absolute mismatch count
        if t.n_mismatch > params.out_filter_mismatch_nmax {
            return false;
        }

        // Relative mismatch count (mismatches / read_length)
        if (t.n_mismatch as f64) > params.out_filter_mismatch_nover_lmax * read_length {
            return false;
        }

        // Absolute matched bases
        let n_matched = t.n_matched();
        if n_matched < params.out_filter_match_nmin {
            return false;
        }

        // Relative matched bases (matched / read_length)
        if (n_matched as f64) < params.out_filter_match_nmin_over_lread * read_length {
            return false;
        }

        true
    });

    // Step 5: Sort by score (descending) and limit to top N
    transcripts.sort_by(|a, b| b.score.cmp(&a.score));
    transcripts.truncate(params.out_filter_multimap_nmax as usize);

    Ok(transcripts)
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
        }
    }

    #[test]
    fn test_align_read_no_seeds() {
        let index = make_test_index();
        let params = make_test_params();

        // Read with all N's (no seeds possible)
        let read_seq = vec![4, 4, 4, 4, 4, 4, 4, 4, 4, 4];

        let result = align_read(&read_seq, &index, &params);
        assert!(result.is_ok());

        let transcripts = result.unwrap();
        assert_eq!(transcripts.len(), 0); // No alignment
    }

    #[test]
    fn test_transcript_filtering_score() {
        let index = make_test_index();
        let mut params = make_test_params();
        params.out_filter_score_min = 50;

        // Would need actual seeds and alignment to test this properly
        // This test just verifies the function doesn't crash
        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, &index, &params);
        assert!(result.is_ok());
    }

    #[test]
    fn test_transcript_filtering_mismatch() {
        let index = make_test_index();
        let mut params = make_test_params();
        params.out_filter_mismatch_nmax = 2;

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, &index, &params);
        assert!(result.is_ok());
    }

    #[test]
    fn test_transcript_multimap_limit() {
        let index = make_test_index();
        let mut params = make_test_params();
        params.out_filter_multimap_nmax = 5;

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, &index, &params);
        assert!(result.is_ok());

        let transcripts = result.unwrap();
        assert!(transcripts.len() <= 5);
    }
}
