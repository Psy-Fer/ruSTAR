/// MAPQ (mapping quality) calculation
///
/// Following STAR's approach:
/// - Unique mappers (n=1): use outSAMmapqUnique (default 255)
/// - Multi-mappers (n>1): -10*log10(1 - 1/n), capped at 255
/// - Unmapped (n=0): 0

/// Calculate MAPQ score based on number of alignments
///
/// # Arguments
/// * `n_alignments` - Number of valid alignments for this read
/// * `mapq_unique` - MAPQ value for unique mappers (typically 255)
///
/// # Returns
/// MAPQ score (0-255)
pub fn calculate_mapq(n_alignments: usize, mapq_unique: u8) -> u8 {
    match n_alignments {
        0 => 0,           // Unmapped
        1 => mapq_unique, // Unique mapper
        n => {
            // Multiple alignments: -10*log10(1 - 1/n)
            // This gives the probability that the chosen alignment is correct
            let p_correct = 1.0 - 1.0 / (n as f64);
            let mapq = -10.0 * p_correct.log10();
            mapq.round().min(255.0) as u8
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mapq_unmapped() {
        assert_eq!(calculate_mapq(0, 255), 0);
    }

    #[test]
    fn test_mapq_unique() {
        assert_eq!(calculate_mapq(1, 255), 255);
        assert_eq!(calculate_mapq(1, 60), 60);
    }

    #[test]
    fn test_mapq_multi() {
        // 2 alignments: -10*log10(0.5) ≈ 3.01
        assert_eq!(calculate_mapq(2, 255), 3);

        // 10 alignments: -10*log10(0.9) ≈ 0.46
        assert_eq!(calculate_mapq(10, 255), 0);

        // 100 alignments: -10*log10(0.99) ≈ 0.04
        assert_eq!(calculate_mapq(100, 255), 0);
    }

    #[test]
    fn test_mapq_capped() {
        // Even with very high alignment quality, cap at 255
        assert_eq!(calculate_mapq(1, 255), 255);
    }

    #[test]
    fn test_mapq_formula() {
        // Verify the formula for a few known values
        // n=2: p=0.5, -10*log10(0.5) = 3.01
        let mapq = calculate_mapq(2, 255);
        assert!(mapq >= 3 && mapq <= 3);

        // n=3: p=0.667, -10*log10(0.667) = 1.76
        let mapq = calculate_mapq(3, 255);
        assert!(mapq >= 2 && mapq <= 2);

        // n=5: p=0.8, -10*log10(0.8) = 0.97
        let mapq = calculate_mapq(5, 255);
        assert!(mapq >= 1 && mapq <= 1);
    }
}
