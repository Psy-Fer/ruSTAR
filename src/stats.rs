/// Alignment statistics tracking and reporting
use log::info;
use std::sync::atomic::{AtomicU64, Ordering};

/// Tracks alignment statistics for a read mapping run
/// Thread-safe using atomic counters
#[derive(Debug)]
pub struct AlignmentStats {
    /// Total number of reads processed
    pub total_reads: AtomicU64,
    /// Reads that mapped uniquely (exactly 1 locus)
    pub uniquely_mapped: AtomicU64,
    /// Reads that mapped to multiple loci (2-N loci)
    pub multi_mapped: AtomicU64,
    /// Reads that did not map
    pub unmapped: AtomicU64,
    /// Reads that mapped to too many loci (exceeds outFilterMultimapNmax)
    pub too_many_loci: AtomicU64,
}

impl Default for AlignmentStats {
    fn default() -> Self {
        Self {
            total_reads: AtomicU64::new(0),
            uniquely_mapped: AtomicU64::new(0),
            multi_mapped: AtomicU64::new(0),
            unmapped: AtomicU64::new(0),
            too_many_loci: AtomicU64::new(0),
        }
    }
}

impl AlignmentStats {
    /// Create new statistics tracker
    pub fn new() -> Self {
        Self::default()
    }

    /// Record an alignment result (thread-safe)
    ///
    /// # Arguments
    /// * `n_alignments` - Number of valid alignments found for the read
    /// * `max_multimaps` - Maximum allowed multi-map count (outFilterMultimapNmax)
    pub fn record_alignment(&self, n_alignments: usize, max_multimaps: usize) {
        self.total_reads.fetch_add(1, Ordering::Relaxed);
        match n_alignments {
            0 => {
                self.unmapped.fetch_add(1, Ordering::Relaxed);
            }
            1 => {
                self.uniquely_mapped.fetch_add(1, Ordering::Relaxed);
            }
            n if n <= max_multimaps => {
                self.multi_mapped.fetch_add(1, Ordering::Relaxed);
            }
            _ => {
                self.too_many_loci.fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    /// Print summary statistics to log
    pub fn print_summary(&self) {
        // Load all atomics once at start
        let total_reads = self.total_reads.load(Ordering::Relaxed);
        let uniquely_mapped = self.uniquely_mapped.load(Ordering::Relaxed);
        let multi_mapped = self.multi_mapped.load(Ordering::Relaxed);
        let unmapped = self.unmapped.load(Ordering::Relaxed);
        let too_many_loci = self.too_many_loci.load(Ordering::Relaxed);

        if total_reads == 0 {
            info!("No reads processed");
            return;
        }

        info!("=== Alignment Summary ===");
        info!("Number of input reads: {}", total_reads);
        info!(
            "Uniquely mapped reads: {} ({:.2}%)",
            uniquely_mapped,
            100.0 * uniquely_mapped as f64 / total_reads as f64
        );
        info!(
            "Multi-mapped reads: {} ({:.2}%)",
            multi_mapped,
            100.0 * multi_mapped as f64 / total_reads as f64
        );
        info!(
            "Unmapped reads: {} ({:.2}%)",
            unmapped,
            100.0 * unmapped as f64 / total_reads as f64
        );
        if too_many_loci > 0 {
            info!(
                "Reads with too many loci: {} ({:.2}%)",
                too_many_loci,
                100.0 * too_many_loci as f64 / total_reads as f64
            );
        }

        // Calculate mapped percentage
        let mapped = uniquely_mapped + multi_mapped;
        info!(
            "Total mapped: {} ({:.2}%)",
            mapped,
            100.0 * mapped as f64 / total_reads as f64
        );
    }

    /// Get percentage of uniquely mapped reads
    pub fn unique_percent(&self) -> f64 {
        let total_reads = self.total_reads.load(Ordering::Relaxed);
        let uniquely_mapped = self.uniquely_mapped.load(Ordering::Relaxed);
        if total_reads == 0 {
            0.0
        } else {
            100.0 * uniquely_mapped as f64 / total_reads as f64
        }
    }

    /// Get percentage of multi-mapped reads
    pub fn multi_percent(&self) -> f64 {
        let total_reads = self.total_reads.load(Ordering::Relaxed);
        let multi_mapped = self.multi_mapped.load(Ordering::Relaxed);
        if total_reads == 0 {
            0.0
        } else {
            100.0 * multi_mapped as f64 / total_reads as f64
        }
    }

    /// Get percentage of unmapped reads
    pub fn unmapped_percent(&self) -> f64 {
        let total_reads = self.total_reads.load(Ordering::Relaxed);
        let unmapped = self.unmapped.load(Ordering::Relaxed);
        if total_reads == 0 {
            0.0
        } else {
            100.0 * unmapped as f64 / total_reads as f64
        }
    }

    /// Get total mapped reads (unique + multi)
    pub fn total_mapped(&self) -> u64 {
        let uniquely_mapped = self.uniquely_mapped.load(Ordering::Relaxed);
        let multi_mapped = self.multi_mapped.load(Ordering::Relaxed);
        uniquely_mapped + multi_mapped
    }

    /// Get percentage of mapped reads
    pub fn mapped_percent(&self) -> f64 {
        let total_reads = self.total_reads.load(Ordering::Relaxed);
        if total_reads == 0 {
            0.0
        } else {
            100.0 * self.total_mapped() as f64 / total_reads as f64
        }
    }

    /// Get total number of reads processed
    pub fn total_reads(&self) -> u64 {
        self.total_reads.load(Ordering::Relaxed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stats_default() {
        let stats = AlignmentStats::new();
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 0);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.too_many_loci.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_record_unique() {
        let stats = AlignmentStats::new();
        stats.record_alignment(1, 10);
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 1);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_record_unmapped() {
        let stats = AlignmentStats::new();
        stats.record_alignment(0, 10);
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_record_multi() {
        let stats = AlignmentStats::new();
        stats.record_alignment(5, 10);
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 1);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_record_too_many() {
        let stats = AlignmentStats::new();
        stats.record_alignment(15, 10);
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.too_many_loci.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_multiple_reads() {
        let stats = AlignmentStats::new();
        stats.record_alignment(1, 10); // unique
        stats.record_alignment(0, 10); // unmapped
        stats.record_alignment(5, 10); // multi
        stats.record_alignment(1, 10); // unique
        stats.record_alignment(15, 10); // too many

        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 5);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 2);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 1);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 1);
        assert_eq!(stats.too_many_loci.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_percentages() {
        let stats = AlignmentStats::new();
        stats.record_alignment(1, 10); // unique
        stats.record_alignment(0, 10); // unmapped
        stats.record_alignment(5, 10); // multi
        stats.record_alignment(1, 10); // unique

        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 4);
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
