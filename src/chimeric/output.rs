// Chimeric.out.junction file writer

use crate::chimeric::segment::ChimericAlignment;
use crate::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// Writer for Chimeric.out.junction file
pub struct ChimericJunctionWriter {
    writer: BufWriter<File>,
}

impl ChimericJunctionWriter {
    /// Create a new chimeric junction writer
    ///
    /// Creates file: {prefix}Chimeric.out.junction
    pub fn new(prefix: &str) -> Result<Self, Error> {
        let mut path = PathBuf::from(prefix);
        path.push("Chimeric.out.junction");

        let file = File::create(&path).map_err(|e| Error::io(e, &path))?;

        let writer = BufWriter::new(file);
        Ok(Self { writer })
    }

    /// Write a chimeric alignment to the file
    ///
    /// Format: 14 tab-separated columns
    /// 1. Donor chromosome
    /// 2. Donor breakpoint (1-based)
    /// 3. Donor strand (+/-)
    /// 4. Acceptor chromosome
    /// 5. Acceptor breakpoint (1-based)
    /// 6. Acceptor strand (+/-)
    /// 7. Junction type (0-6)
    /// 8. Repeat length donor
    /// 9. Repeat length acceptor
    /// 10. Read name
    /// 11. First segment start (1-based)
    /// 12. First segment CIGAR
    /// 13. Second segment start (1-based)
    /// 14. Second segment CIGAR
    pub fn write_alignment(
        &mut self,
        alignment: &ChimericAlignment,
        chr_names: &[String],
        read_name: &str,
    ) -> Result<(), Error> {
        // Get chromosome names
        let donor_chr = &chr_names[alignment.donor.chr_idx];
        let acceptor_chr = &chr_names[alignment.acceptor.chr_idx];

        // Get breakpoints (1-based)
        let donor_bp = alignment.donor_breakpoint();
        let acceptor_bp = alignment.acceptor_breakpoint();

        // Get strand symbols
        let donor_strand = alignment.donor_strand();
        let acceptor_strand = alignment.acceptor_strand();

        // Get junction type
        let junction_type = alignment.junction_type;

        // Get repeat lengths
        let repeat_donor = alignment.repeat_len_donor;
        let repeat_acceptor = alignment.repeat_len_acceptor;

        // Get segment start positions (1-based)
        let donor_start = alignment.donor.genome_start + 1;
        let acceptor_start = alignment.acceptor.genome_start + 1;

        // Convert CIGAR to string
        let donor_cigar = cigar_to_string(&alignment.donor.cigar);
        let acceptor_cigar = cigar_to_string(&alignment.acceptor.cigar);

        // Write line
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            donor_chr,
            donor_bp,
            donor_strand,
            acceptor_chr,
            acceptor_bp,
            acceptor_strand,
            junction_type,
            repeat_donor,
            repeat_acceptor,
            read_name,
            donor_start,
            donor_cigar,
            acceptor_start,
            acceptor_cigar,
        )
        .map_err(|e| Error::Chimeric(format!("Failed to write chimeric junction: {}", e)))?;

        Ok(())
    }

    /// Flush buffered data to disk
    pub fn flush(&mut self) -> Result<(), Error> {
        self.writer
            .flush()
            .map_err(|e| Error::Chimeric(format!("Failed to flush chimeric junction file: {}", e)))
    }
}

/// Convert CIGAR operations to CIGAR string
fn cigar_to_string(cigar: &[crate::align::transcript::CigarOp]) -> String {
    use crate::align::transcript::CigarOp;

    let mut result = String::new();
    for op in cigar {
        match op {
            CigarOp::Match(len) => result.push_str(&format!("{}M", len)),
            CigarOp::Equal(len) => result.push_str(&format!("{}=", len)),
            CigarOp::Diff(len) => result.push_str(&format!("{}X", len)),
            CigarOp::Ins(len) => result.push_str(&format!("{}I", len)),
            CigarOp::Del(len) => result.push_str(&format!("{}D", len)),
            CigarOp::RefSkip(len) => result.push_str(&format!("{}N", len)),
            CigarOp::SoftClip(len) => result.push_str(&format!("{}S", len)),
            CigarOp::HardClip(len) => result.push_str(&format!("{}H", len)),
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::CigarOp;
    use crate::chimeric::segment::{ChimericAlignment, ChimericSegment};
    use std::io::Read;
    use tempfile::tempdir;

    #[test]
    fn test_cigar_to_string() {
        let cigar = vec![
            CigarOp::Match(50),
            CigarOp::Ins(2),
            CigarOp::Del(3),
            CigarOp::RefSkip(1000),
            CigarOp::SoftClip(5),
        ];

        assert_eq!(cigar_to_string(&cigar), "50M2I3D1000N5S");
    }

    #[test]
    fn test_chimeric_junction_writer_creation() {
        let dir = tempdir().unwrap();
        let prefix = dir.path().to_str().unwrap();

        let writer = ChimericJunctionWriter::new(prefix);
        assert!(writer.is_ok());

        let mut path = PathBuf::from(prefix);
        path.push("Chimeric.out.junction");
        assert!(path.exists());
    }

    #[test]
    fn test_write_inter_chromosomal() {
        let dir = tempdir().unwrap();
        let prefix = dir.path().to_str().unwrap();

        let mut writer = ChimericJunctionWriter::new(prefix).unwrap();

        // Create mock chimeric alignment (chr9 -> chr22, BCR-ABL fusion)
        let donor = ChimericSegment::new(
            0,
            133738300,
            133738363,
            false,
            0,
            63,
            vec![CigarOp::Match(63)],
            100,
            2,
        );

        let acceptor = ChimericSegment::new(
            1,
            23632600,
            23632637,
            false,
            63,
            100,
            vec![CigarOp::Match(37)],
            80,
            1,
        );

        let alignment = ChimericAlignment::new(
            donor,
            acceptor,
            1, // GT/AG
            0,
            0,
            vec![0; 100],
            "READ_001".to_string(),
        );

        let chr_names = vec!["chr9".to_string(), "chr22".to_string()];

        writer
            .write_alignment(&alignment, &chr_names, "READ_001")
            .unwrap();
        writer.flush().unwrap();

        // Read file and verify
        let mut path = PathBuf::from(prefix);
        path.push("Chimeric.out.junction");

        let mut content = String::new();
        File::open(&path)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();

        let line = content.trim();
        let fields: Vec<&str> = line.split('\t').collect();

        assert_eq!(fields.len(), 14);
        assert_eq!(fields[0], "chr9"); // donor chr
        assert_eq!(fields[1], "133738363"); // donor breakpoint
        assert_eq!(fields[2], "+"); // donor strand
        assert_eq!(fields[3], "chr22"); // acceptor chr
        assert_eq!(fields[4], "23632601"); // acceptor breakpoint
        assert_eq!(fields[5], "+"); // acceptor strand
        assert_eq!(fields[6], "1"); // junction type
        assert_eq!(fields[7], "0"); // repeat donor
        assert_eq!(fields[8], "0"); // repeat acceptor
        assert_eq!(fields[9], "READ_001"); // read name
        assert_eq!(fields[10], "133738301"); // donor start (1-based)
        assert_eq!(fields[11], "63M"); // donor CIGAR
        assert_eq!(fields[12], "23632601"); // acceptor start (1-based)
        assert_eq!(fields[13], "37M"); // acceptor CIGAR
    }

    #[test]
    fn test_write_strand_break() {
        let dir = tempdir().unwrap();
        let prefix = dir.path().to_str().unwrap();

        let mut writer = ChimericJunctionWriter::new(prefix).unwrap();

        // Create mock chimeric alignment (same chr, opposite strands)
        let donor = ChimericSegment::new(
            0,
            1000,
            1050,
            false, // forward
            0,
            50,
            vec![CigarOp::Match(50)],
            100,
            1,
        );

        let acceptor = ChimericSegment::new(
            0,
            2000,
            2050,
            true, // reverse
            50,
            100,
            vec![CigarOp::Match(50)],
            100,
            1,
        );

        let alignment = ChimericAlignment::new(
            donor,
            acceptor,
            0, // non-canonical
            0,
            0,
            vec![0; 100],
            "READ_002".to_string(),
        );

        let chr_names = vec!["chr1".to_string()];

        writer
            .write_alignment(&alignment, &chr_names, "READ_002")
            .unwrap();
        writer.flush().unwrap();

        // Read file and verify
        let mut path = PathBuf::from(prefix);
        path.push("Chimeric.out.junction");

        let mut content = String::new();
        File::open(&path)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();

        let line = content.trim();
        let fields: Vec<&str> = line.split('\t').collect();

        assert_eq!(fields.len(), 14);
        assert_eq!(fields[0], "chr1"); // donor chr
        assert_eq!(fields[2], "+"); // donor strand
        assert_eq!(fields[3], "chr1"); // acceptor chr
        assert_eq!(fields[5], "-"); // acceptor strand (reverse)
        assert_eq!(fields[6], "0"); // junction type (non-canonical)
    }
}
