/// Alignment statistics tracking and reporting
use log::info;

/// Tracks alignment statistics for a read mapping run
#[derive(Default, Debug)]
pub struct AlignmentStats {
    /// Total number of reads processed
    pub total_reads: u64,
    /// Reads that mapped uniquely (exactly 1 locus)
    pub uniquely_mapped: u64,
    /// Reads that mapped to multiple loci (2-N loci)
    pub multi_mapped: u64,
    /// Reads that did not map
    pub unmapped: u64,
    /// Reads that mapped to too many loci (exceeds outFilterMultimapNmax)
    pub too_many_loci: u64,
}

impl AlignmentStats {
    /// Create new statistics tracker
    pub fn new() -> Self {
        Self::default()
    }

    /// Record an alignment result
    ///
    /// # Arguments
    /// * `n_alignments` - Number of valid alignments found for the read
    /// * `max_multimaps` - Maximum allowed multi-map count (outFilterMultimapNmax)
    pub fn record_alignment(&mut self, n_alignments: usize, max_multimaps: usize) {
        self.total_reads += 1;
        match n_alignments {
            0 => self.unmapped += 1,
            1 => self.uniquely_mapped += 1,
            n if n <= max_multimaps => self.multi_mapped += 1,
            _ => self.too_many_loci += 1,
        }
    }

    /// Print summary statistics to log
    pub fn print_summary(&self) {
        if self.total_reads == 0 {
            info!("No reads processed");
            return;
        }

        info!("=== Alignment Summary ===");
        info!("Number of input reads: {}", self.total_reads);
        info!(
            "Uniquely mapped reads: {} ({:.2}%)",
            self.uniquely_mapped,
            100.0 * self.uniquely_mapped as f64 / self.total_reads as f64
        );
        info!(
            "Multi-mapped reads: {} ({:.2}%)",
            self.multi_mapped,
            100.0 * self.multi_mapped as f64 / self.total_reads as f64
        );
        info!(
            "Unmapped reads: {} ({:.2}%)",
            self.unmapped,
            100.0 * self.unmapped as f64 / self.total_reads as f64
        );
        if self.too_many_loci > 0 {
            info!(
                "Reads with too many loci: {} ({:.2}%)",
                self.too_many_loci,
                100.0 * self.too_many_loci as f64 / self.total_reads as f64
            );
        }

        // Calculate mapped percentage
        let mapped = self.uniquely_mapped + self.multi_mapped;
        info!(
            "Total mapped: {} ({:.2}%)",
            mapped,
            100.0 * mapped as f64 / self.total_reads as f64
        );
    }

    /// Get percentage of uniquely mapped reads
    pub fn unique_percent(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            100.0 * self.uniquely_mapped as f64 / self.total_reads as f64
        }
    }

    /// Get percentage of multi-mapped reads
    pub fn multi_percent(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            100.0 * self.multi_mapped as f64 / self.total_reads as f64
        }
    }

    /// Get percentage of unmapped reads
    pub fn unmapped_percent(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            100.0 * self.unmapped as f64 / self.total_reads as f64
        }
    }

    /// Get total mapped reads (unique + multi)
    pub fn total_mapped(&self) -> u64 {
        self.uniquely_mapped + self.multi_mapped
    }

    /// Get percentage of mapped reads
    pub fn mapped_percent(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            100.0 * self.total_mapped() as f64 / self.total_reads as f64
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stats_default() {
        let stats = AlignmentStats::default();
        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.uniquely_mapped, 0);
        assert_eq!(stats.multi_mapped, 0);
        assert_eq!(stats.unmapped, 0);
        assert_eq!(stats.too_many_loci, 0);
    }

    #[test]
    fn test_record_unique() {
        let mut stats = AlignmentStats::new();
        stats.record_alignment(1, 10);
        assert_eq!(stats.total_reads, 1);
        assert_eq!(stats.uniquely_mapped, 1);
        assert_eq!(stats.multi_mapped, 0);
        assert_eq!(stats.unmapped, 0);
    }

    #[test]
    fn test_record_unmapped() {
        let mut stats = AlignmentStats::new();
        stats.record_alignment(0, 10);
        assert_eq!(stats.total_reads, 1);
        assert_eq!(stats.uniquely_mapped, 0);
        assert_eq!(stats.multi_mapped, 0);
        assert_eq!(stats.unmapped, 1);
    }

    #[test]
    fn test_record_multi() {
        let mut stats = AlignmentStats::new();
        stats.record_alignment(5, 10);
        assert_eq!(stats.total_reads, 1);
        assert_eq!(stats.uniquely_mapped, 0);
        assert_eq!(stats.multi_mapped, 1);
        assert_eq!(stats.unmapped, 0);
    }

    #[test]
    fn test_record_too_many() {
        let mut stats = AlignmentStats::new();
        stats.record_alignment(15, 10);
        assert_eq!(stats.total_reads, 1);
        assert_eq!(stats.uniquely_mapped, 0);
        assert_eq!(stats.multi_mapped, 0);
        assert_eq!(stats.unmapped, 0);
        assert_eq!(stats.too_many_loci, 1);
    }

    #[test]
    fn test_multiple_reads() {
        let mut stats = AlignmentStats::new();
        stats.record_alignment(1, 10); // unique
        stats.record_alignment(0, 10); // unmapped
        stats.record_alignment(5, 10); // multi
        stats.record_alignment(1, 10); // unique
        stats.record_alignment(15, 10); // too many

        assert_eq!(stats.total_reads, 5);
        assert_eq!(stats.uniquely_mapped, 2);
        assert_eq!(stats.multi_mapped, 1);
        assert_eq!(stats.unmapped, 1);
        assert_eq!(stats.too_many_loci, 1);
    }

    #[test]
    fn test_percentages() {
        let mut stats = AlignmentStats::new();
        stats.record_alignment(1, 10); // unique
        stats.record_alignment(0, 10); // unmapped
        stats.record_alignment(5, 10); // multi
        stats.record_alignment(1, 10); // unique

        assert_eq!(stats.total_reads, 4);
        assert!((stats.unique_percent() - 50.0).abs() < 0.01);
        assert!((stats.multi_percent() - 25.0).abs() < 0.01);
        assert!((stats.unmapped_percent() - 25.0).abs() < 0.01);
        assert!((stats.mapped_percent() - 75.0).abs() < 0.01);
    }

    #[test]
    fn test_empty_stats() {
        let stats = AlignmentStats::new();
        assert_eq!(stats.unique_percent(), 0.0);
        assert_eq!(stats.multi_percent(), 0.0);
        assert_eq!(stats.unmapped_percent(), 0.0);
        assert_eq!(stats.mapped_percent(), 0.0);
        assert_eq!(stats.total_mapped(), 0);
    }
}
