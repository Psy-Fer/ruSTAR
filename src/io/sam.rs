/// SAM/BAM output writer with noodles
use crate::align::transcript::{CigarOp, Transcript};
use crate::error::Error;
use crate::genome::Genome;
use crate::io::fastq::decode_base;
use crate::mapq::calculate_mapq;
use crate::params::Parameters;
use noodles::sam;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};
use noodles::sam::header::record::value::{Map, map::Program};
use std::fs::File;
use std::io::BufWriter;
use std::num::NonZeroUsize;
use std::path::Path;

/// Buffer for SAM records built by parallel threads
#[derive(Default)]
pub struct BufferedSamRecords {
    pub records: Vec<RecordBuf>,
}

impl BufferedSamRecords {
    /// Create new buffer with capacity
    pub fn new() -> Self {
        Self {
            records: Vec::with_capacity(10000),
        }
    }

    /// Add a record to the buffer
    pub fn push(&mut self, record: RecordBuf) {
        self.records.push(record);
    }
}

/// SAM file writer
pub struct SamWriter {
    writer: sam::io::Writer<BufWriter<File>>,
    header: sam::Header,
}

impl SamWriter {
    /// Create a new SAM writer with header from genome index
    ///
    /// # Arguments
    /// * `output_path` - Path to output SAM file
    /// * `genome` - Genome index with chromosome information
    /// * `params` - Parameters (for @PG header)
    pub fn create(output_path: &Path, genome: &Genome, params: &Parameters) -> Result<Self, Error> {
        let file = File::create(output_path)?;
        let buf_writer = BufWriter::new(file);

        let header = build_sam_header(genome, params)?;
        let mut writer = sam::io::Writer::new(buf_writer);

        writer.write_header(&header)?;

        Ok(Self { writer, header })
    }

    /// Write alignment record(s) for a read
    ///
    /// # Arguments
    /// * `read_name` - Read identifier
    /// * `read_seq` - Read sequence (encoded)
    /// * `read_qual` - Quality scores
    /// * `transcripts` - Alignment transcripts (1 or more for multi-mappers)
    /// * `genome` - Genome index
    /// * `params` - Parameters
    pub fn write_alignment(
        &mut self,
        read_name: &str,
        read_seq: &[u8],
        read_qual: &[u8],
        transcripts: &[Transcript],
        genome: &Genome,
        params: &Parameters,
    ) -> Result<(), Error> {
        if transcripts.is_empty() {
            return Ok(());
        }

        let n_alignments = transcripts.len();
        let mapq = calculate_mapq(n_alignments, params.out_sam_mapq_unique);

        for (hit_index, transcript) in transcripts.iter().enumerate() {
            let record = transcript_to_record(
                transcript,
                read_name,
                read_seq,
                read_qual,
                genome,
                mapq,
                n_alignments,
                hit_index + 1, // 1-based
            )?;

            self.writer.write_alignment_record(&self.header, &record)?;
        }

        Ok(())
    }

    /// Write unmapped read
    ///
    /// # Arguments
    /// * `read_name` - Read identifier
    /// * `read_seq` - Read sequence (encoded)
    /// * `read_qual` - Quality scores
    pub fn write_unmapped(
        &mut self,
        read_name: &str,
        read_seq: &[u8],
        read_qual: &[u8],
    ) -> Result<(), Error> {
        let record = Self::build_unmapped_record(read_name, read_seq, read_qual)?;
        self.writer.write_alignment_record(&self.header, &record)?;
        Ok(())
    }

    /// Write batch of buffered records (for parallel processing)
    ///
    /// # Arguments
    /// * `batch` - Slice of records to write
    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        for record in batch {
            self.writer.write_alignment_record(&self.header, record)?;
        }
        Ok(())
    }

    /// Build unmapped record (without writing)
    ///
    /// # Arguments
    /// * `read_name` - Read identifier
    /// * `read_seq` - Read sequence (encoded)
    /// * `read_qual` - Quality scores
    pub fn build_unmapped_record(
        read_name: &str,
        read_seq: &[u8],
        read_qual: &[u8],
    ) -> Result<RecordBuf, Error> {
        let mut record = RecordBuf::default();

        // Name
        record.name_mut().replace(read_name.into());

        // FLAGS: 0x4 (unmapped)
        let flags = sam::alignment::record::Flags::UNMAPPED;
        *record.flags_mut() = flags;

        // Sequence (decode from genome encoding)
        let seq_bytes: Vec<u8> = read_seq.iter().map(|&b| decode_base(b)).collect();
        *record.sequence_mut() = Sequence::from(seq_bytes);

        // Quality scores
        *record.quality_scores_mut() = QualityScores::from(read_qual.to_vec());

        Ok(record)
    }

    /// Build alignment records (without writing) for a read
    ///
    /// # Arguments
    /// * `read_name` - Read identifier
    /// * `read_seq` - Read sequence (encoded)
    /// * `read_qual` - Quality scores
    /// * `transcripts` - Alignment transcripts (1 or more for multi-mappers)
    /// * `genome` - Genome index
    /// * `params` - Parameters
    pub fn build_alignment_records(
        read_name: &str,
        read_seq: &[u8],
        read_qual: &[u8],
        transcripts: &[Transcript],
        genome: &Genome,
        params: &Parameters,
    ) -> Result<Vec<RecordBuf>, Error> {
        if transcripts.is_empty() {
            return Ok(Vec::new());
        }

        let n_alignments = transcripts.len();
        let mapq = calculate_mapq(n_alignments, params.out_sam_mapq_unique);

        let mut records = Vec::with_capacity(n_alignments);
        for (hit_index, transcript) in transcripts.iter().enumerate() {
            let record = transcript_to_record(
                transcript,
                read_name,
                read_seq,
                read_qual,
                genome,
                mapq,
                n_alignments,
                hit_index + 1, // 1-based
            )?;
            records.push(record);
        }

        Ok(records)
    }
}

/// Build SAM header from genome
fn build_sam_header(genome: &Genome, _params: &Parameters) -> Result<sam::Header, Error> {
    let mut builder = sam::Header::builder();

    // @HD line (default version and unsorted)
    builder = builder.set_header(Default::default());

    // @SQ lines for each chromosome
    for i in 0..genome.n_chr_real {
        let name = &genome.chr_name[i];
        let length = genome.chr_length[i] as usize;

        let length_nz = NonZeroUsize::new(length)
            .ok_or_else(|| Error::Index(format!("chromosome {} has zero length", name)))?;

        builder = builder.add_reference_sequence(
            name.as_str(),
            Map::<sam::header::record::value::map::ReferenceSequence>::new(length_nz),
        );
    }

    // @PG line
    builder = builder.add_program("ruSTAR", Map::<Program>::default());

    Ok(builder.build())
}

/// Convert Transcript to SAM record
fn transcript_to_record(
    transcript: &Transcript,
    read_name: &str,
    read_seq: &[u8],
    read_qual: &[u8],
    genome: &Genome,
    mapq: u8,
    _n_alignments: usize,
    _hit_index: usize,
) -> Result<RecordBuf, Error> {
    let mut record = RecordBuf::default();

    // Name
    record.name_mut().replace(read_name.into());

    // FLAGS
    let mut flags = sam::alignment::record::Flags::empty();
    if transcript.is_reverse {
        flags |= sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
    }
    *record.flags_mut() = flags;

    // RNAME (reference sequence name)
    if transcript.chr_idx >= genome.n_chr_real {
        return Err(Error::Alignment(format!(
            "invalid chromosome index {} (max {})",
            transcript.chr_idx,
            genome.n_chr_real - 1
        )));
    }
    *record.reference_sequence_id_mut() = Some(transcript.chr_idx);

    // POS (1-based)
    let pos = (transcript.genome_start + 1) as usize;
    *record.alignment_start_mut() = Some(
        pos.try_into()
            .map_err(|e| Error::Alignment(format!("invalid alignment position {}: {}", pos, e)))?,
    );

    // MAPQ
    *record.mapping_quality_mut() = MappingQuality::new(mapq);

    // CIGAR
    let cigar = convert_cigar(&transcript.cigar)?;
    *record.cigar_mut() = cigar;

    // Sequence (decode from genome encoding)
    let seq_bytes: Vec<u8> = read_seq.iter().map(|&b| decode_base(b)).collect();
    *record.sequence_mut() = Sequence::from(seq_bytes);

    // Quality scores
    *record.quality_scores_mut() = QualityScores::from(read_qual.to_vec());

    // Optional tags
    // TODO: Add SAM tags (AS, NM, NH, HI, nM, jM, jI) in Phase 6 refinement
    // For now, minimal tags are sufficient for basic alignment output

    Ok(record)
}

/// Convert ruSTAR CigarOp to noodles Cigar
fn convert_cigar(ops: &[CigarOp]) -> Result<sam::alignment::record_buf::Cigar, Error> {
    use sam::alignment::record::cigar::op::Kind;

    let mut cigar = sam::alignment::record_buf::Cigar::default();

    for op in ops {
        let kind = match op {
            CigarOp::Match(_) => Kind::Match,
            CigarOp::Equal(_) => Kind::SequenceMatch,
            CigarOp::Diff(_) => Kind::SequenceMismatch,
            CigarOp::Ins(_) => Kind::Insertion,
            CigarOp::Del(_) => Kind::Deletion,
            CigarOp::RefSkip(_) => Kind::Skip,
            CigarOp::SoftClip(_) => Kind::SoftClip,
            CigarOp::HardClip(_) => Kind::HardClip,
        };

        let len = op.len() as usize;
        let noodles_op = sam::alignment::record::cigar::Op::new(kind, len);
        cigar.as_mut().push(noodles_op);
    }

    Ok(cigar)
}

// TODO: SAM optional tags will be added in Phase 6 refinement
// The noodles API for tags requires careful handling of lifetimes
// For now, we'll focus on getting the core alignment working

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Genome;
    use clap::Parser;
    use tempfile::NamedTempFile;

    fn make_test_genome() -> Genome {
        Genome {
            sequence: vec![0, 1, 2, 3, 0, 1, 2, 3], // ACGTACGT
            n_genome: 8,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![8],
            chr_start: vec![0, 8],
        }
    }

    #[test]
    fn test_convert_cigar() {
        let ops = vec![
            CigarOp::Match(50),
            CigarOp::Ins(3),
            CigarOp::Del(2),
            CigarOp::RefSkip(100),
            CigarOp::Match(50),
        ];

        let cigar = convert_cigar(&ops).unwrap();
        assert_eq!(cigar.as_ref().len(), 5);
    }

    #[test]
    fn test_build_sam_header() {
        let genome = make_test_genome();
        let params = Parameters::parse_from(vec!["ruSTAR", "--readFilesIn", "test.fq"]);

        let header = build_sam_header(&genome, &params).unwrap();

        // Check that we have reference sequences
        assert_eq!(header.reference_sequences().len(), 1);

        // Check that we have a program line (just check header is valid)
        assert_eq!(header.reference_sequences().len(), 1);
    }

    #[test]
    fn test_sam_writer_creation() {
        let genome = make_test_genome();
        let params = Parameters::parse_from(vec!["ruSTAR", "--readFilesIn", "test.fq"]);

        let tmpfile = NamedTempFile::new().unwrap();
        let writer = SamWriter::create(tmpfile.path(), &genome, &params);
        assert!(writer.is_ok());
    }

    #[test]
    fn test_write_unmapped() {
        let genome = make_test_genome();
        let params = Parameters::parse_from(vec!["ruSTAR", "--readFilesIn", "test.fq"]);

        let tmpfile = NamedTempFile::new().unwrap();
        let mut writer = SamWriter::create(tmpfile.path(), &genome, &params).unwrap();

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let read_qual = vec![30, 30, 30, 30];

        let result = writer.write_unmapped("read1", &read_seq, &read_qual);
        assert!(result.is_ok());
    }

    #[test]
    fn test_transcript_to_record() {
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 10,
            genome_end: 60,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(50)],
            score: 100,
            n_mismatch: 2,
            n_gap: 0,
            n_junction: 0,
            read_seq: vec![0, 1, 2, 3],
        };

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
        );
        assert!(record.is_ok());

        let record = record.unwrap();
        assert_eq!(
            record.name().map(|n| n.to_string()),
            Some("read1".to_string())
        );
        assert_eq!(record.reference_sequence_id(), Some(0));
        assert_eq!(record.alignment_start().map(|p| usize::from(p)), Some(11)); // 1-based
    }

    // TODO: Add tests for SAM tags when they are implemented
}
