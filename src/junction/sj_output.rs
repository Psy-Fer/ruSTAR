/// Splice junction output (SJ.out.tab file generation)
///
/// Format (9 tab-separated columns):
/// 1. chromosome
/// 2. intron start (1-based)
/// 3. intron end (1-based)
/// 4. strand (0=undefined, 1=+, 2=-)
/// 5. motif (0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG, 4=CT/GC, 5=AT/AC, 6=GT/AT)
/// 6. annotated (0=no, 1=yes)
/// 7. unique-mapping reads
/// 8. multi-mapping reads
/// 9. maximum overhang
use crate::align::score::SpliceMotif;
use crate::error::Error;
use crate::genome::Genome;
use dashmap::DashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};

/// Key for junction statistics
#[derive(Hash, Eq, PartialEq, Clone, Debug)]
pub(crate) struct SjKey {
    pub chr_idx: usize,
    pub intron_start: u64,
    pub intron_end: u64,
    pub strand: u8,
    pub motif: u8, // Encoded motif value
}

/// Counts for a single junction
#[derive(Debug)]
pub(crate) struct SjCounts {
    pub unique_count: AtomicU32,
    pub multi_count: AtomicU32,
    pub max_overhang: AtomicU32,
    pub annotated: bool,
}

/// Thread-safe junction statistics accumulator
pub struct SpliceJunctionStats {
    /// Thread-safe map for parallel accumulation
    junctions: DashMap<SjKey, SjCounts>,
}

impl Clone for SpliceJunctionStats {
    fn clone(&self) -> Self {
        let new_map = DashMap::new();
        for entry in self.junctions.iter() {
            let key = entry.key().clone();
            let counts = entry.value();
            new_map.insert(
                key,
                SjCounts {
                    unique_count: AtomicU32::new(counts.unique_count.load(Ordering::Relaxed)),
                    multi_count: AtomicU32::new(counts.multi_count.load(Ordering::Relaxed)),
                    max_overhang: AtomicU32::new(counts.max_overhang.load(Ordering::Relaxed)),
                    annotated: counts.annotated,
                },
            );
        }
        Self { junctions: new_map }
    }
}

impl SpliceJunctionStats {
    /// Create new empty statistics accumulator
    pub fn new() -> Self {
        Self {
            junctions: DashMap::new(),
        }
    }

    /// Record a junction from alignment (thread-safe)
    ///
    /// # Arguments
    /// * `chr_idx` - Chromosome index
    /// * `start` - Intron start (1-based)
    /// * `end` - Intron end (1-based)
    /// * `strand` - Strand (0=unknown, 1=+, 2=-)
    /// * `motif` - Splice motif
    /// * `is_unique` - True if from unique mapping, false if multi-mapping
    /// * `overhang` - Minimum overhang on either side
    /// * `annotated` - True if junction is annotated in GTF
    pub fn record_junction(
        &self,
        chr_idx: usize,
        start: u64,
        end: u64,
        strand: u8,
        motif: SpliceMotif,
        is_unique: bool,
        overhang: u32,
        annotated: bool,
    ) {
        let key = SjKey {
            chr_idx,
            intron_start: start,
            intron_end: end,
            strand,
            motif: encode_motif(motif, strand),
        };

        // Get or create entry
        self.junctions
            .entry(key)
            .or_insert_with(|| SjCounts {
                unique_count: AtomicU32::new(0),
                multi_count: AtomicU32::new(0),
                max_overhang: AtomicU32::new(0),
                annotated,
            })
            .value()
            .record(is_unique, overhang);
    }

    /// Write SJ.out.tab file
    pub fn write_output(&self, output_path: &Path, genome: &Genome) -> Result<(), Error> {
        let file = File::create(output_path).map_err(|e| Error::io(e, output_path))?;
        let mut writer = BufWriter::new(file);

        // Collect all junctions
        let mut junctions: Vec<_> = self
            .junctions
            .iter()
            .map(|entry| {
                let key = entry.key();
                let counts = entry.value();
                (
                    key.chr_idx,
                    key.intron_start,
                    key.intron_end,
                    key.strand,
                    key.motif,
                    counts.annotated,
                    counts.unique_count.load(Ordering::Relaxed),
                    counts.multi_count.load(Ordering::Relaxed),
                    counts.max_overhang.load(Ordering::Relaxed),
                )
            })
            .collect();

        // Sort by chromosome, start, end
        junctions.sort_by(|a, b| {
            a.0.cmp(&b.0) // chr_idx
                .then(a.1.cmp(&b.1)) // start
                .then(a.2.cmp(&b.2)) // end
        });

        // Write each junction
        for (chr_idx, start, end, strand, motif, annotated, unique, multi, max_overhang) in
            junctions
        {
            let chr_name = genome
                .chr_name
                .get(chr_idx)
                .ok_or_else(|| Error::Index("Invalid chromosome index in junction".to_string()))?;

            // Convert from global genome coordinates to per-chromosome coordinates
            let chr_start = genome.chr_start[chr_idx];
            let chr_pos_start = start - chr_start;
            let chr_pos_end = end - chr_start;

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                chr_name,
                chr_pos_start,
                chr_pos_end,
                strand,
                motif,
                if annotated { 1 } else { 0 },
                unique,
                multi,
                max_overhang
            )
            .map_err(|e| Error::io(e, output_path))?;
        }

        writer.flush().map_err(|e| Error::io(e, output_path))?;

        log::info!(
            "Wrote {} junctions to {}",
            self.junctions.len(),
            output_path.display()
        );

        Ok(())
    }

    /// Get the number of unique junctions tracked
    pub fn len(&self) -> usize {
        self.junctions.len()
    }

    /// Check if any junctions have been recorded
    pub fn is_empty(&self) -> bool {
        self.junctions.is_empty()
    }

    /// Iterate over all junctions (for two-pass mode filtering)
    pub(crate) fn iter(
        &self,
    ) -> impl Iterator<Item = dashmap::mapref::multiple::RefMulti<'_, SjKey, SjCounts>> {
        self.junctions.iter()
    }
}

impl Default for SpliceJunctionStats {
    fn default() -> Self {
        Self::new()
    }
}

impl SjCounts {
    /// Record a junction occurrence (thread-safe)
    fn record(&self, is_unique: bool, overhang: u32) {
        if is_unique {
            self.unique_count.fetch_add(1, Ordering::Relaxed);
        } else {
            self.multi_count.fetch_add(1, Ordering::Relaxed);
        }

        // Update max overhang using compare-and-swap
        let mut current = self.max_overhang.load(Ordering::Relaxed);
        while overhang > current {
            match self.max_overhang.compare_exchange_weak(
                current,
                overhang,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(x) => current = x,
            }
        }
    }
}

/// Encode splice motif for SJ.out.tab
///
/// STAR convention:
/// 0 = non-canonical
/// 1 = GT/AG (canonical, + strand)
/// 2 = CT/AC (reverse complement of GT/AG, - strand)
/// 3 = GC/AG (+ strand)
/// 4 = CT/GC (reverse complement of GC/AG, - strand)
/// 5 = AT/AC (+ strand)
/// 6 = GT/AT (rare)
fn encode_motif(motif: SpliceMotif, strand: u8) -> u8 {
    match motif {
        SpliceMotif::GtAg => {
            if strand == 2 {
                2 // CT/AC (reverse complement)
            } else {
                1 // GT/AG
            }
        }
        SpliceMotif::CtAc => 2, // CT/AC (detected directly on forward genome)
        SpliceMotif::GcAg => {
            if strand == 2 {
                4 // CT/GC (reverse complement)
            } else {
                3 // GC/AG
            }
        }
        SpliceMotif::CtGc => 4, // CT/GC (detected directly on forward genome)
        SpliceMotif::AtAc => 5, // AT/AC
        SpliceMotif::GtAt => 6, // GT/AT (reverse complement of AT/AC)
        SpliceMotif::NonCanonical => 0, // Non-canonical
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sj_stats_new() {
        let stats = SpliceJunctionStats::new();
        assert_eq!(stats.len(), 0);
        assert!(stats.is_empty());
    }

    #[test]
    fn test_record_junction_unique() {
        let stats = SpliceJunctionStats::new();

        stats.record_junction(0, 100, 200, 1, SpliceMotif::GtAg, true, 10, false);

        assert_eq!(stats.len(), 1);

        let entry = stats.junctions.iter().next().unwrap();
        assert_eq!(entry.value().unique_count.load(Ordering::Relaxed), 1);
        assert_eq!(entry.value().multi_count.load(Ordering::Relaxed), 0);
        assert_eq!(entry.value().max_overhang.load(Ordering::Relaxed), 10);
        assert!(!entry.value().annotated);
    }

    #[test]
    fn test_record_junction_multi() {
        let stats = SpliceJunctionStats::new();

        stats.record_junction(0, 100, 200, 1, SpliceMotif::GtAg, false, 15, true);

        assert_eq!(stats.len(), 1);

        let entry = stats.junctions.iter().next().unwrap();
        assert_eq!(entry.value().unique_count.load(Ordering::Relaxed), 0);
        assert_eq!(entry.value().multi_count.load(Ordering::Relaxed), 1);
        assert_eq!(entry.value().max_overhang.load(Ordering::Relaxed), 15);
        assert!(entry.value().annotated);
    }

    #[test]
    fn test_record_junction_multiple_times() {
        let stats = SpliceJunctionStats::new();

        // Record same junction multiple times
        stats.record_junction(0, 100, 200, 1, SpliceMotif::GtAg, true, 10, false);
        stats.record_junction(0, 100, 200, 1, SpliceMotif::GtAg, true, 15, false);
        stats.record_junction(0, 100, 200, 1, SpliceMotif::GtAg, false, 12, false);

        assert_eq!(stats.len(), 1); // Same junction

        let entry = stats.junctions.iter().next().unwrap();
        assert_eq!(entry.value().unique_count.load(Ordering::Relaxed), 2);
        assert_eq!(entry.value().multi_count.load(Ordering::Relaxed), 1);
        assert_eq!(entry.value().max_overhang.load(Ordering::Relaxed), 15);
    }

    #[test]
    fn test_max_overhang_update() {
        let stats = SpliceJunctionStats::new();

        // Record with increasing overhangs
        stats.record_junction(0, 100, 200, 1, SpliceMotif::GtAg, true, 10, false);
        stats.record_junction(0, 100, 200, 1, SpliceMotif::GtAg, true, 20, false);
        stats.record_junction(0, 100, 200, 1, SpliceMotif::GtAg, true, 15, false); // Lower, should not update

        let entry = stats.junctions.iter().next().unwrap();
        assert_eq!(entry.value().max_overhang.load(Ordering::Relaxed), 20);
    }

    #[test]
    fn test_encode_motif_gtag() {
        assert_eq!(encode_motif(SpliceMotif::GtAg, 1), 1); // + strand
        assert_eq!(encode_motif(SpliceMotif::GtAg, 2), 2); // - strand (reverse complement)
    }

    #[test]
    fn test_encode_motif_gcag() {
        assert_eq!(encode_motif(SpliceMotif::GcAg, 1), 3); // + strand
        assert_eq!(encode_motif(SpliceMotif::GcAg, 2), 4); // - strand (reverse complement)
    }

    #[test]
    fn test_encode_motif_atac() {
        assert_eq!(encode_motif(SpliceMotif::AtAc, 1), 5);
        assert_eq!(encode_motif(SpliceMotif::AtAc, 2), 5); // No reverse complement encoding
    }

    #[test]
    fn test_encode_motif_noncanonical() {
        assert_eq!(encode_motif(SpliceMotif::NonCanonical, 1), 0);
        assert_eq!(encode_motif(SpliceMotif::NonCanonical, 2), 0);
    }

    #[test]
    fn test_sj_key_equality() {
        let key1 = SjKey {
            chr_idx: 0,
            intron_start: 100,
            intron_end: 200,
            strand: 1,
            motif: 1,
        };
        let key2 = SjKey {
            chr_idx: 0,
            intron_start: 100,
            intron_end: 200,
            strand: 1,
            motif: 1,
        };
        let key3 = SjKey {
            chr_idx: 0,
            intron_start: 100,
            intron_end: 200,
            strand: 2,
            motif: 2,
        };

        assert_eq!(key1, key2);
        assert_ne!(key1, key3); // Different strand and motif
    }

    #[test]
    fn test_write_output() {
        use tempfile::NamedTempFile;

        let stats = SpliceJunctionStats::new();
        stats.record_junction(0, 100, 200, 1, SpliceMotif::GtAg, true, 10, false);
        stats.record_junction(0, 300, 400, 2, SpliceMotif::GcAg, false, 15, true);

        let genome = Genome {
            sequence: vec![0; 1000],
            n_genome: 1000,
            n_chr_real: 1,
            chr_start: vec![0, 1000],
            chr_length: vec![1000],
            chr_name: vec!["chr1".to_string()],
        };

        let output_file = NamedTempFile::new().unwrap();
        stats.write_output(output_file.path(), &genome).unwrap();

        // Read and verify output
        let content = std::fs::read_to_string(output_file.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(lines.len(), 2);

        // First junction (sorted by position)
        let fields1: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(fields1[0], "chr1"); // chr
        assert_eq!(fields1[1], "100"); // start
        assert_eq!(fields1[2], "200"); // end
        assert_eq!(fields1[3], "1"); // strand
        assert_eq!(fields1[4], "1"); // motif (GT/AG)
        assert_eq!(fields1[5], "0"); // not annotated
        assert_eq!(fields1[6], "1"); // unique count
        assert_eq!(fields1[7], "0"); // multi count
        assert_eq!(fields1[8], "10"); // max overhang

        // Second junction
        let fields2: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(fields2[0], "chr1");
        assert_eq!(fields2[1], "300");
        assert_eq!(fields2[2], "400");
        assert_eq!(fields2[3], "2"); // - strand
        assert_eq!(fields2[4], "4"); // motif (CT/GC, reverse of GC/AG)
        assert_eq!(fields2[5], "1"); // annotated
        assert_eq!(fields2[6], "0"); // unique count
        assert_eq!(fields2[7], "1"); // multi count
        assert_eq!(fields2[8], "15"); // max overhang
    }
}
