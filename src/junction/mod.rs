/// Splice junction annotation and tracking
///
/// This module handles:
/// - GTF file parsing for gene/transcript/exon annotations
/// - Building a junction database from annotated exons
/// - Junction lookup during alignment (annotated vs novel)
/// - Junction statistics collection for SJ.out.tab output
mod gtf;
mod sj_output;

pub use sj_output::SpliceJunctionStats;

use crate::error::Error;
use crate::genome::Genome;
use std::collections::HashMap;
use std::path::Path;

/// Key for junction lookup: (chr_idx, intron_start, intron_end, strand)
#[derive(Hash, Eq, PartialEq, Clone, Debug)]
struct JunctionKey {
    chr_idx: usize,
    intron_start: u64,
    intron_end: u64,
    strand: u8, // 0=unknown, 1=+, 2=-
}

/// Information about a splice junction
#[derive(Debug, Clone)]
struct JunctionInfo {
    annotated: bool,
    // Future: gene_id, transcript_ids for provenance tracking
}

/// Splice junction database built from GTF annotations
pub struct SpliceJunctionDb {
    /// Map: (chr_idx, intron_start, intron_end, strand) → annotated
    junctions: HashMap<JunctionKey, JunctionInfo>,
}

impl SpliceJunctionDb {
    /// Create empty database (for no-GTF mode)
    pub fn empty() -> Self {
        Self {
            junctions: HashMap::new(),
        }
    }

    /// Build junction database from GTF file
    pub fn from_gtf(gtf_path: &Path, genome: &Genome) -> Result<Self, Error> {
        log::info!("Loading GTF annotations from: {}", gtf_path.display());

        // Parse GTF and extract exon features
        let exons = gtf::parse_gtf(gtf_path)?;
        log::debug!("Parsed {} exon features from GTF", exons.len());

        // Extract junctions from consecutive exons
        let junctions_vec = gtf::extract_junctions_from_exons(exons, genome)?;
        log::info!(
            "Extracted {} annotated junctions from GTF",
            junctions_vec.len()
        );

        // Build HashMap for fast lookup
        let mut junctions = HashMap::new();
        for (chr_idx, intron_start, intron_end, strand) in junctions_vec {
            let key = JunctionKey {
                chr_idx,
                intron_start,
                intron_end,
                strand,
            };
            junctions.insert(key, JunctionInfo { annotated: true });
        }

        Ok(Self { junctions })
    }

    /// Check if a junction is annotated in the GTF
    ///
    /// # Arguments
    /// * `chr_idx` - Chromosome index
    /// * `start` - Intron start position (last exon base + 1)
    /// * `end` - Intron end position (first exon base of next exon - 1)
    /// * `strand` - Strand (0=unknown, 1=+, 2=-)
    ///
    /// # Returns
    /// `true` if junction is annotated, `false` otherwise
    pub fn is_annotated(&self, chr_idx: usize, start: u64, end: u64, strand: u8) -> bool {
        let key = JunctionKey {
            chr_idx,
            intron_start: start,
            intron_end: end,
            strand,
        };
        self.junctions.get(&key).is_some_and(|info| info.annotated)
    }

    /// Get the number of annotated junctions in the database
    pub fn len(&self) -> usize {
        self.junctions.len()
    }

    /// Check if the database is empty
    pub fn is_empty(&self) -> bool {
        self.junctions.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_junction_db_empty() {
        let db = SpliceJunctionDb::empty();
        assert_eq!(db.len(), 0);
        assert!(db.is_empty());
        assert!(!db.is_annotated(0, 100, 200, 1));
    }

    #[test]
    fn test_junction_key_equality() {
        let key1 = JunctionKey {
            chr_idx: 0,
            intron_start: 100,
            intron_end: 200,
            strand: 1,
        };
        let key2 = JunctionKey {
            chr_idx: 0,
            intron_start: 100,
            intron_end: 200,
            strand: 1,
        };
        let key3 = JunctionKey {
            chr_idx: 0,
            intron_start: 100,
            intron_end: 200,
            strand: 2,
        };

        assert_eq!(key1, key2);
        assert_ne!(key1, key3); // Different strand
    }

    #[test]
    fn test_junction_lookup() {
        let mut db = SpliceJunctionDb::empty();

        // Manually insert a junction
        db.junctions.insert(
            JunctionKey {
                chr_idx: 0,
                intron_start: 100,
                intron_end: 200,
                strand: 1,
            },
            JunctionInfo { annotated: true },
        );

        // Should find annotated junction
        assert!(db.is_annotated(0, 100, 200, 1));

        // Should not find with different strand
        assert!(!db.is_annotated(0, 100, 200, 2));

        // Should not find with different coordinates
        assert!(!db.is_annotated(0, 101, 200, 1));
        assert!(!db.is_annotated(0, 100, 201, 1));
    }

    #[test]
    fn test_junction_strand_specific() {
        let mut db = SpliceJunctionDb::empty();

        // Add same junction coordinates but different strands
        db.junctions.insert(
            JunctionKey {
                chr_idx: 0,
                intron_start: 100,
                intron_end: 200,
                strand: 1,
            },
            JunctionInfo { annotated: true },
        );
        db.junctions.insert(
            JunctionKey {
                chr_idx: 0,
                intron_start: 100,
                intron_end: 200,
                strand: 2,
            },
            JunctionInfo { annotated: true },
        );

        assert_eq!(db.len(), 2);
        assert!(db.is_annotated(0, 100, 200, 1));
        assert!(db.is_annotated(0, 100, 200, 2));
        assert!(!db.is_annotated(0, 100, 200, 0)); // Unknown strand
    }
}
