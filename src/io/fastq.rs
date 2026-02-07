/// FASTQ reader with base encoding and decompression support
use crate::error::Error;
use flate2::read::GzDecoder;
use noodles::fastq;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::{Command, Stdio};

/// A read from a FASTQ file with encoded bases
#[derive(Debug, Clone)]
pub struct EncodedRead {
    /// Read identifier
    pub name: String,
    /// Base sequence encoded as 0=A, 1=C, 2=G, 3=T, 4=N
    pub sequence: Vec<u8>,
    /// Quality scores (raw FASTQ quality values)
    pub quality: Vec<u8>,
}

/// FASTQ reader that handles decompression and base encoding
pub struct FastqReader {
    inner: fastq::Reader<Box<dyn BufRead + Send>>,
}

impl FastqReader {
    /// Open a FASTQ file (plain or gzip compressed)
    ///
    /// # Arguments
    /// * `path` - Path to FASTQ file
    /// * `decompress_cmd` - Optional decompression command (e.g., "zcat" for .gz files)
    ///
    /// # Returns
    /// A FastqReader that iterates over encoded reads
    pub fn open(path: &Path, decompress_cmd: Option<&str>) -> Result<Self, Error> {
        let reader: Box<dyn BufRead + Send> = if let Some(cmd) = decompress_cmd {
            // Use external decompression command
            Self::open_with_command(path, cmd)?
        } else {
            // Auto-detect compression by file extension
            let path_str = path.to_string_lossy();
            let is_gzipped = path_str.ends_with(".gz") || path_str.ends_with(".gzip");

            let file = File::open(path).map_err(|e| Error::io(e, path))?;

            if is_gzipped {
                // Gzipped file
                Box::new(BufReader::new(GzDecoder::new(file)))
            } else {
                // Plain text FASTQ
                Box::new(BufReader::new(file))
            }
        };

        let fastq_reader = fastq::Reader::new(reader);

        Ok(Self {
            inner: fastq_reader,
        })
    }

    /// Open FASTQ file using external decompression command
    fn open_with_command(path: &Path, cmd: &str) -> Result<Box<dyn BufRead + Send>, Error> {
        let mut child = Command::new(cmd)
            .arg(path)
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|e| Error::io(e, path))?;

        let stdout = child.stdout.take().ok_or_else(|| {
            Error::from(std::io::Error::other(
                "failed to capture stdout from decompression command",
            ))
        })?;

        Ok(Box::new(BufReader::new(stdout)))
    }

    /// Get next read with encoded bases
    pub fn next_encoded(&mut self) -> Result<Option<EncodedRead>, Error> {
        match self.inner.records().next() {
            Some(Ok(record)) => {
                let name = std::str::from_utf8(record.name())
                    .map_err(|e| {
                        Error::from(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!("invalid UTF-8 in read name: {}", e),
                        ))
                    })?
                    .to_string();

                let sequence = record.sequence().iter().map(|&b| encode_base(b)).collect();

                let quality = record.quality_scores().to_vec();

                Ok(Some(EncodedRead {
                    name,
                    sequence,
                    quality,
                }))
            }
            Some(Err(e)) => Err(Error::from(e)),
            None => Ok(None),
        }
    }

    /// Read a batch of encoded reads for parallel processing
    ///
    /// # Arguments
    /// * `batch_size` - Maximum number of reads to return
    ///
    /// # Returns
    /// Vector of encoded reads (may be shorter than batch_size at end of file)
    pub fn read_batch(&mut self, batch_size: usize) -> Result<Vec<EncodedRead>, Error> {
        let mut batch = Vec::with_capacity(batch_size);
        for _ in 0..batch_size {
            match self.next_encoded()? {
                Some(read) => batch.push(read),
                None => break,
            }
        }
        Ok(batch)
    }
}

/// Convert FASTQ base character to genome encoding
///
/// # Arguments
/// * `base` - ASCII base character (A, C, G, T, N, or lowercase variants)
///
/// # Returns
/// Encoded base: 0=A, 1=C, 2=G, 3=T, 4=N (or any ambiguous base)
pub fn encode_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4, // N or any ambiguous base (R, Y, S, W, K, M, etc.)
    }
}

/// Decode genome encoding to ASCII base character
///
/// # Arguments
/// * `encoded` - Encoded base (0-4)
///
/// # Returns
/// ASCII base character (A, C, G, T, or N)
pub fn decode_base(encoded: u8) -> u8 {
    match encoded {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

/// Apply read clipping from 5' and 3' ends
///
/// # Arguments
/// * `seq` - Original sequence
/// * `qual` - Original quality scores
/// * `clip5p` - Number of bases to clip from 5' end
/// * `clip3p` - Number of bases to clip from 3' end
///
/// # Returns
/// Tuple of (clipped_sequence, clipped_quality)
pub fn clip_read(seq: &[u8], qual: &[u8], clip5p: usize, clip3p: usize) -> (Vec<u8>, Vec<u8>) {
    let len = seq.len();

    // Handle edge cases
    if clip5p + clip3p >= len {
        // Clipping removes entire read
        return (Vec::new(), Vec::new());
    }

    let start = clip5p;
    let end = len - clip3p;

    (seq[start..end].to_vec(), qual[start..end].to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'A'), 0);
        assert_eq!(encode_base(b'a'), 0);
        assert_eq!(encode_base(b'C'), 1);
        assert_eq!(encode_base(b'c'), 1);
        assert_eq!(encode_base(b'G'), 2);
        assert_eq!(encode_base(b'g'), 2);
        assert_eq!(encode_base(b'T'), 3);
        assert_eq!(encode_base(b't'), 3);
        assert_eq!(encode_base(b'N'), 4);
        assert_eq!(encode_base(b'n'), 4);
        // Ambiguous bases
        assert_eq!(encode_base(b'R'), 4);
        assert_eq!(encode_base(b'Y'), 4);
        assert_eq!(encode_base(b'S'), 4);
    }

    #[test]
    fn test_decode_base() {
        assert_eq!(decode_base(0), b'A');
        assert_eq!(decode_base(1), b'C');
        assert_eq!(decode_base(2), b'G');
        assert_eq!(decode_base(3), b'T');
        assert_eq!(decode_base(4), b'N');
        assert_eq!(decode_base(5), b'N'); // Invalid -> N
    }

    #[test]
    fn test_clip_read_none() {
        let seq = vec![0, 1, 2, 3, 0]; // ACGTA
        let qual = vec![30, 30, 30, 30, 30];

        let (clipped_seq, clipped_qual) = clip_read(&seq, &qual, 0, 0);
        assert_eq!(clipped_seq, seq);
        assert_eq!(clipped_qual, qual);
    }

    #[test]
    fn test_clip_read_5p() {
        let seq = vec![0, 1, 2, 3, 0]; // ACGTA
        let qual = vec![30, 30, 30, 30, 30];

        let (clipped_seq, clipped_qual) = clip_read(&seq, &qual, 2, 0);
        assert_eq!(clipped_seq, vec![2, 3, 0]); // GTA
        assert_eq!(clipped_qual, vec![30, 30, 30]);
    }

    #[test]
    fn test_clip_read_3p() {
        let seq = vec![0, 1, 2, 3, 0]; // ACGTA
        let qual = vec![30, 30, 30, 30, 30];

        let (clipped_seq, clipped_qual) = clip_read(&seq, &qual, 0, 2);
        assert_eq!(clipped_seq, vec![0, 1, 2]); // ACG
        assert_eq!(clipped_qual, vec![30, 30, 30]);
    }

    #[test]
    fn test_clip_read_both() {
        let seq = vec![0, 1, 2, 3, 0]; // ACGTA
        let qual = vec![30, 30, 30, 30, 30];

        let (clipped_seq, clipped_qual) = clip_read(&seq, &qual, 1, 1);
        assert_eq!(clipped_seq, vec![1, 2, 3]); // CGT
        assert_eq!(clipped_qual, vec![30, 30, 30]);
    }

    #[test]
    fn test_clip_read_entire() {
        let seq = vec![0, 1, 2, 3, 0]; // ACGTA
        let qual = vec![30, 30, 30, 30, 30];

        let (clipped_seq, clipped_qual) = clip_read(&seq, &qual, 3, 3);
        assert_eq!(clipped_seq, Vec::<u8>::new());
        assert_eq!(clipped_qual, Vec::<u8>::new());
    }

    #[test]
    fn test_fastq_reader_plain() {
        let mut tmpfile = NamedTempFile::new().unwrap();
        writeln!(tmpfile, "@read1").unwrap();
        writeln!(tmpfile, "ACGTN").unwrap();
        writeln!(tmpfile, "+").unwrap();
        writeln!(tmpfile, "IIIII").unwrap();
        writeln!(tmpfile, "@read2").unwrap();
        writeln!(tmpfile, "TGCA").unwrap();
        writeln!(tmpfile, "+").unwrap();
        writeln!(tmpfile, "HHHH").unwrap();
        tmpfile.flush().unwrap();

        let mut reader = FastqReader::open(tmpfile.path(), None).unwrap();

        let read1 = reader.next_encoded().unwrap().unwrap();
        assert_eq!(read1.name, "read1");
        assert_eq!(read1.sequence, vec![0, 1, 2, 3, 4]); // ACGTN
        assert_eq!(read1.quality.len(), 5);

        let read2 = reader.next_encoded().unwrap().unwrap();
        assert_eq!(read2.name, "read2");
        assert_eq!(read2.sequence, vec![3, 2, 1, 0]); // TGCA

        let read3 = reader.next_encoded().unwrap();
        assert!(read3.is_none());
    }

    #[test]
    fn test_fastq_reader_gzip() {
        use flate2::Compression;
        use flate2::write::GzEncoder;

        let tmpfile = tempfile::Builder::new()
            .suffix(".fastq.gz")
            .tempfile()
            .unwrap();
        let mut encoder = GzEncoder::new(tmpfile.as_file(), Compression::default());
        writeln!(encoder, "@read1").unwrap();
        writeln!(encoder, "ACGT").unwrap();
        writeln!(encoder, "+").unwrap();
        writeln!(encoder, "IIII").unwrap();
        encoder.finish().unwrap();

        let mut reader = FastqReader::open(tmpfile.path(), None).unwrap();

        let read1 = reader.next_encoded().unwrap().unwrap();
        assert_eq!(read1.name, "read1");
        assert_eq!(read1.sequence, vec![0, 1, 2, 3]); // ACGT
        assert_eq!(read1.quality.len(), 4);
    }
}
