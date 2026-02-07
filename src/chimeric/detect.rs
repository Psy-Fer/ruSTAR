// Chimeric alignment detection algorithms

use crate::align::SeedCluster;
use crate::align::seed::Seed;
use crate::align::stitch::stitch_seeds;
use crate::align::transcript::Transcript;
use crate::chimeric::score::{calculate_repeat_length, classify_junction_type};
use crate::chimeric::segment::{ChimericAlignment, ChimericSegment};
use crate::error::Error;
use crate::index::GenomeIndex;
use crate::params::Parameters;

/// Chimeric alignment detector
pub struct ChimericDetector<'a> {
    params: &'a Parameters,
}

impl<'a> ChimericDetector<'a> {
    /// Create a new chimeric detector
    pub fn new(params: &'a Parameters) -> Self {
        Self { params }
    }

    /// Detect chimeric alignments from soft-clipped transcripts (Tier 1)
    ///
    /// Triggers when:
    /// - Transcript has >20% soft-clipped bases
    /// - Soft-clipped segment >= chimSegmentMin
    pub fn detect_from_soft_clips(
        &self,
        transcript: &Transcript,
        read_seq: &[u8],
        _read_name: &str,
        _index: &GenomeIndex,
    ) -> Result<Option<ChimericAlignment>, Error> {
        let (left_clip, right_clip) = transcript.count_soft_clips();
        let total_clip = left_clip + right_clip;
        let read_len = read_seq.len() as u32;

        // Check if soft-clipping exceeds threshold
        if total_clip as f64 / read_len as f64 <= 0.20 {
            return Ok(None);
        }

        // Check if clipped segment meets minimum length
        let min_seg = self.params.chim_segment_min;
        if left_clip < min_seg && right_clip < min_seg {
            return Ok(None);
        }

        // For now, we mark this as a potential chimera but don't re-map the clipped region
        // (Tier 3 re-mapping is deferred to Phase 12.2)
        // Instead, we just log that a chimeric candidate was detected but not processed
        Ok(None)
    }

    /// Detect chimeric alignments from multi-cluster seeds (Tier 2)
    ///
    /// Triggers when:
    /// - Seeds cluster on different chromosomes
    /// - Seeds cluster on different strands (same chromosome)
    /// - Seeds cluster with large genomic distance (>1Mb, same chr/strand)
    pub fn detect_from_multi_clusters(
        &self,
        clusters: &[SeedCluster],
        seeds: &[Seed],
        read_seq: &[u8],
        read_name: &str,
        index: &GenomeIndex,
    ) -> Result<Vec<ChimericAlignment>, Error> {
        let mut chimeras = Vec::new();

        // Find cluster pairs with chimeric signatures
        for i in 0..clusters.len() {
            for j in (i + 1)..clusters.len() {
                if self.is_chimeric_signature(&clusters[i], &clusters[j]) {
                    // Try to build chimeric alignment from these clusters
                    if let Some(chim) = self.build_chimeric_from_clusters(
                        &clusters[i],
                        &clusters[j],
                        seeds,
                        read_seq,
                        read_name,
                        index,
                    )? {
                        chimeras.push(chim);
                    }
                }
            }
        }

        Ok(chimeras)
    }

    /// Check if two clusters represent a chimeric signature
    fn is_chimeric_signature(&self, c1: &SeedCluster, c2: &SeedCluster) -> bool {
        // Different chromosomes
        if c1.chr_idx != c2.chr_idx {
            return true;
        }

        // Different strands (same chromosome)
        if c1.is_reverse != c2.is_reverse {
            return true;
        }

        // Large genomic distance (same chr/strand)
        let distance = genomic_distance(c1, c2);
        if distance > 1_000_000 {
            return true;
        }

        false
    }

    /// Build chimeric alignment from two clusters
    fn build_chimeric_from_clusters(
        &self,
        cluster1: &SeedCluster,
        cluster2: &SeedCluster,
        seeds: &[Seed],
        read_seq: &[u8],
        read_name: &str,
        index: &GenomeIndex,
    ) -> Result<Option<ChimericAlignment>, Error> {
        if cluster1.seed_indices.is_empty() || cluster2.seed_indices.is_empty() {
            return Ok(None);
        }

        // Stitch each cluster independently using existing stitch_seeds
        use crate::align::score::AlignmentScorer;
        let scorer = AlignmentScorer::from_params(self.params);

        let transcripts1 = stitch_seeds(cluster1, seeds, read_seq, index, &scorer)?;
        let transcripts2 = stitch_seeds(cluster2, seeds, read_seq, index, &scorer)?;

        if transcripts1.is_empty() || transcripts2.is_empty() {
            return Ok(None);
        }

        // Take best transcript from each cluster
        let t1 = &transcripts1[0];
        let t2 = &transcripts2[0];

        // Determine donor/acceptor based on read position
        let (donor_t, acceptor_t) = if t1.exons[0].read_start < t2.exons[0].read_start {
            (t1, t2)
        } else {
            (t2, t1)
        };

        // Convert transcripts to chimeric segments
        let donor = transcript_to_segment(donor_t)?;
        let acceptor = transcript_to_segment(acceptor_t)?;

        // Check minimum segment lengths
        if !donor.meets_min_length(self.params.chim_segment_min)
            || !acceptor.meets_min_length(self.params.chim_segment_min)
        {
            return Ok(None);
        }

        // Classify junction type
        let junction_type = classify_junction_type(
            &index.genome,
            donor.chr_idx,
            donor.genome_end,
            donor.is_reverse,
            acceptor.chr_idx,
            acceptor.genome_start,
            acceptor.is_reverse,
        );

        // Calculate repeat lengths
        let (repeat_len_donor, repeat_len_acceptor) = calculate_repeat_length(
            &index.genome,
            donor.chr_idx,
            donor.genome_end,
            acceptor.chr_idx,
            acceptor.genome_start,
            20, // max check distance
        );

        // Create chimeric alignment
        let chim = ChimericAlignment::new(
            donor,
            acceptor,
            junction_type,
            repeat_len_donor,
            repeat_len_acceptor,
            read_seq.to_vec(),
            read_name.to_string(),
        );

        Ok(Some(chim))
    }
}

/// Calculate genomic distance between two clusters
fn genomic_distance(c1: &SeedCluster, c2: &SeedCluster) -> u64 {
    if c1.chr_idx != c2.chr_idx {
        return u64::MAX;
    }

    if c1.genome_end < c2.genome_start {
        c2.genome_start - c1.genome_end
    } else if c2.genome_end < c1.genome_start {
        c1.genome_start - c2.genome_end
    } else {
        0 // overlapping
    }
}

/// Convert a transcript to a chimeric segment
fn transcript_to_segment(transcript: &Transcript) -> Result<ChimericSegment, Error> {
    if transcript.exons.is_empty() {
        return Err(Error::Alignment(
            "Cannot convert empty transcript to segment".to_string(),
        ));
    }

    // Get overall bounds
    let read_start = transcript.exons[0].read_start;
    let read_end = transcript.exons.last().unwrap().read_end;

    Ok(ChimericSegment::new(
        transcript.chr_idx,
        transcript.genome_start,
        transcript.genome_end,
        transcript.is_reverse,
        read_start,
        read_end,
        transcript.cigar.clone(),
        transcript.score,
        transcript.n_mismatch,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[test]
    fn test_genomic_distance_same_chr() {
        let c1 = SeedCluster {
            seed_indices: vec![0],
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            anchor_idx: 0,
        };
        let c2 = SeedCluster {
            seed_indices: vec![1],
            chr_idx: 0,
            genome_start: 1200,
            genome_end: 1300,
            is_reverse: false,
            anchor_idx: 1,
        };

        assert_eq!(genomic_distance(&c1, &c2), 100);
        assert_eq!(genomic_distance(&c2, &c1), 100);
    }

    #[test]
    fn test_genomic_distance_overlapping() {
        let c1 = SeedCluster {
            seed_indices: vec![0],
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1200,
            is_reverse: false,
            anchor_idx: 0,
        };
        let c2 = SeedCluster {
            seed_indices: vec![1],
            chr_idx: 0,
            genome_start: 1100,
            genome_end: 1300,
            is_reverse: false,
            anchor_idx: 1,
        };

        assert_eq!(genomic_distance(&c1, &c2), 0);
    }

    #[test]
    fn test_genomic_distance_different_chr() {
        let c1 = SeedCluster {
            seed_indices: vec![0],
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            anchor_idx: 0,
        };
        let c2 = SeedCluster {
            seed_indices: vec![1],
            chr_idx: 1,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            anchor_idx: 1,
        };

        assert_eq!(genomic_distance(&c1, &c2), u64::MAX);
    }

    #[test]
    fn test_is_chimeric_signature_different_chr() {
        let params = Parameters::try_parse_from(vec!["ruSTAR"]).unwrap();
        let detector = ChimericDetector::new(&params);

        let c1 = SeedCluster {
            seed_indices: vec![0],
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            anchor_idx: 0,
        };
        let c2 = SeedCluster {
            seed_indices: vec![1],
            chr_idx: 1,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            anchor_idx: 1,
        };

        assert!(detector.is_chimeric_signature(&c1, &c2));
    }

    #[test]
    fn test_is_chimeric_signature_strand_break() {
        let params = Parameters::try_parse_from(vec!["ruSTAR"]).unwrap();
        let detector = ChimericDetector::new(&params);

        let c1 = SeedCluster {
            seed_indices: vec![0],
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            anchor_idx: 0,
        };
        let c2 = SeedCluster {
            seed_indices: vec![1],
            chr_idx: 0,
            genome_start: 1200,
            genome_end: 1300,
            is_reverse: true,
            anchor_idx: 1,
        };

        assert!(detector.is_chimeric_signature(&c1, &c2));
    }

    #[test]
    fn test_is_chimeric_signature_large_distance() {
        let params = Parameters::try_parse_from(vec!["ruSTAR"]).unwrap();
        let detector = ChimericDetector::new(&params);

        let c1 = SeedCluster {
            seed_indices: vec![0],
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            anchor_idx: 0,
        };
        let c2 = SeedCluster {
            seed_indices: vec![1],
            chr_idx: 0,
            genome_start: 2_000_000,
            genome_end: 2_000_100,
            is_reverse: false,
            anchor_idx: 1,
        };

        assert!(detector.is_chimeric_signature(&c1, &c2));
    }

    #[test]
    fn test_is_chimeric_signature_close_same_strand() {
        let params = Parameters::try_parse_from(vec!["ruSTAR"]).unwrap();
        let detector = ChimericDetector::new(&params);

        let c1 = SeedCluster {
            seed_indices: vec![0],
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            anchor_idx: 0,
        };
        let c2 = SeedCluster {
            seed_indices: vec![1],
            chr_idx: 0,
            genome_start: 1200,
            genome_end: 1300,
            is_reverse: false,
            anchor_idx: 1,
        };

        assert!(!detector.is_chimeric_signature(&c1, &c2));
    }
}
