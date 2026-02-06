use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::Error;

/// A single chromosome from a FASTA file.
#[derive(Debug, Clone)]
pub struct Chromosome {
    pub name: String,
    pub sequence: Vec<u8>,
}

/// Parse FASTA files and extract chromosomes.
///
/// Matches STAR's behavior:
/// - Base encoding: A=0, C=1, G=2, T=3, N/other=4
/// - Control characters (ASCII < 32) are skipped
/// - Case-insensitive (a == A)
pub fn parse_fasta_files<P: AsRef<Path>>(paths: &[P]) -> Result<Vec<Chromosome>, Error> {
    let mut chromosomes = Vec::new();

    for path in paths {
        let path = path.as_ref();
        let file = File::open(path).map_err(|e| Error::io(e, path))?;
        let reader = BufReader::new(file);

        let mut current_name: Option<String> = None;
        let mut current_seq: Vec<u8> = Vec::new();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result.map_err(|e| Error::io(e, path))?;

            if line.is_empty() {
                continue;
            }

            if let Some(stripped) = line.strip_prefix('>') {
                // New chromosome header
                if let Some(name) = current_name.take() {
                    // Save the previous chromosome
                    chromosomes.push(Chromosome {
                        name,
                        sequence: current_seq,
                    });
                    current_seq = Vec::new();
                }

                // Extract chromosome name (everything after '>' up to first whitespace)
                let name = stripped
                    .split_whitespace()
                    .next()
                    .ok_or_else(|| {
                        Error::Fasta(format!(
                            "empty chromosome name at {}:{}",
                            path.display(),
                            line_num + 1
                        ))
                    })?
                    .to_string();

                current_name = Some(name);
            } else {
                // Sequence line
                if current_name.is_none() {
                    return Err(Error::Fasta(format!(
                        "sequence data before first header at {}:{}",
                        path.display(),
                        line_num + 1
                    )));
                }

                // Convert bases to numeric encoding, skipping control characters
                for &byte in line.as_bytes() {
                    if byte < 32 {
                        // Control character — skip
                        continue;
                    }

                    let encoded = match byte {
                        b'A' | b'a' => 0,
                        b'C' | b'c' => 1,
                        b'G' | b'g' => 2,
                        b'T' | b't' => 3,
                        _ => 4, // N and all other non-control characters
                    };

                    current_seq.push(encoded);
                }
            }
        }

        // Save the last chromosome
        if let Some(name) = current_name {
            chromosomes.push(Chromosome {
                name,
                sequence: current_seq,
            });
        }
    }

    if chromosomes.is_empty() {
        return Err(Error::Fasta(
            "no chromosomes found in FASTA files".to_string(),
        ));
    }

    Ok(chromosomes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn parse_single_chromosome() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, "NNNN").unwrap();

        let chroms = parse_fasta_files(&[file.path()]).unwrap();
        assert_eq!(chroms.len(), 1);
        assert_eq!(chroms[0].name, "chr1");
        assert_eq!(chroms[0].sequence, vec![0, 1, 2, 3, 4, 4, 4, 4]);
    }

    #[test]
    fn parse_multiple_chromosomes() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1 some comment").unwrap();
        writeln!(file, "ACG").unwrap();
        writeln!(file, ">chr2").unwrap();
        writeln!(file, "TGA").unwrap();

        let chroms = parse_fasta_files(&[file.path()]).unwrap();
        assert_eq!(chroms.len(), 2);
        assert_eq!(chroms[0].name, "chr1");
        assert_eq!(chroms[0].sequence, vec![0, 1, 2]);
        assert_eq!(chroms[1].name, "chr2");
        assert_eq!(chroms[1].sequence, vec![3, 2, 0]);
    }

    #[test]
    fn case_insensitive() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">test").unwrap();
        writeln!(file, "AaCcGgTt").unwrap();

        let chroms = parse_fasta_files(&[file.path()]).unwrap();
        assert_eq!(chroms[0].sequence, vec![0, 0, 1, 1, 2, 2, 3, 3]);
    }

    #[test]
    fn multiple_files() {
        let mut file1 = NamedTempFile::new().unwrap();
        writeln!(file1, ">chr1").unwrap();
        writeln!(file1, "AC").unwrap();

        let mut file2 = NamedTempFile::new().unwrap();
        writeln!(file2, ">chr2").unwrap();
        writeln!(file2, "GT").unwrap();

        let chroms = parse_fasta_files(&[file1.path(), file2.path()]).unwrap();
        assert_eq!(chroms.len(), 2);
        assert_eq!(chroms[0].name, "chr1");
        assert_eq!(chroms[1].name, "chr2");
    }

    #[test]
    fn empty_file_error() {
        let file = NamedTempFile::new().unwrap();
        let result = parse_fasta_files(&[file.path()]);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("no chromosomes"));
    }

    #[test]
    fn sequence_before_header_error() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "ACGT").unwrap();

        let result = parse_fasta_files(&[file.path()]);
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("sequence data before first header")
        );
    }
}
