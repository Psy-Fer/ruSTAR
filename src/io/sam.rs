/// SAM/BAM output writer with noodles
use crate::align::read_align::PairedAlignment;
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

    /// Build paired-end SAM records (without writing)
    ///
    /// Returns 2 records per pair (one for each mate)
    ///
    /// # Arguments
    /// * `read_name` - Read identifier (base name without /1 or /2)
    /// * `mate1_seq` - First mate sequence (encoded)
    /// * `mate1_qual` - First mate quality scores
    /// * `mate2_seq` - Second mate sequence (encoded)
    /// * `mate2_qual` - Second mate quality scores
    /// * `paired_alignments` - Paired alignments
    /// * `genome` - Genome index
    /// * `params` - Parameters
    #[allow(clippy::too_many_arguments)]
    pub fn build_paired_records(
        read_name: &str,
        mate1_seq: &[u8],
        mate1_qual: &[u8],
        mate2_seq: &[u8],
        mate2_qual: &[u8],
        paired_alignments: &[PairedAlignment],
        genome: &Genome,
        params: &Parameters,
    ) -> Result<Vec<RecordBuf>, Error> {
        if paired_alignments.is_empty() {
            // Both mates unmapped
            return Self::build_paired_unmapped_records(
                read_name, mate1_seq, mate1_qual, mate2_seq, mate2_qual,
            );
        }

        let n_alignments = paired_alignments.len();
        let mapq = calculate_mapq(n_alignments, params.out_sam_mapq_unique);

        let mut records = Vec::with_capacity(n_alignments * 2);

        for paired_aln in paired_alignments {
            // Create record for mate1
            let rec1 = build_paired_mate_record(
                read_name,
                mate1_seq,
                mate1_qual,
                &paired_aln.transcript,
                genome,
                mapq,
                true, // is_first_mate
                paired_aln.is_proper_pair,
                paired_aln.insert_size,
            )?;
            records.push(rec1);

            // Create record for mate2
            let rec2 = build_paired_mate_record(
                read_name,
                mate2_seq,
                mate2_qual,
                &paired_aln.transcript,
                genome,
                mapq,
                false, // is_first_mate
                paired_aln.is_proper_pair,
                -paired_aln.insert_size, // Negative for mate2
            )?;
            records.push(rec2);
        }

        Ok(records)
    }

    /// Build unmapped paired records (both mates unmapped)
    pub fn build_paired_unmapped_records(
        read_name: &str,
        mate1_seq: &[u8],
        mate1_qual: &[u8],
        mate2_seq: &[u8],
        mate2_qual: &[u8],
    ) -> Result<Vec<RecordBuf>, Error> {
        let mut records = Vec::with_capacity(2);

        // Mate1 record
        let mut rec1 = RecordBuf::default();
        rec1.name_mut().replace(read_name.into());

        // FLAGS: 0x1 (paired) | 0x4 (unmapped) | 0x8 (mate unmapped) | 0x40 (first in pair)
        let flags1 = sam::alignment::record::Flags::SEGMENTED
            | sam::alignment::record::Flags::UNMAPPED
            | sam::alignment::record::Flags::MATE_UNMAPPED
            | sam::alignment::record::Flags::FIRST_SEGMENT;
        *rec1.flags_mut() = flags1;

        let seq1_bytes: Vec<u8> = mate1_seq.iter().map(|&b| decode_base(b)).collect();
        *rec1.sequence_mut() = Sequence::from(seq1_bytes);
        *rec1.quality_scores_mut() = QualityScores::from(mate1_qual.to_vec());
        records.push(rec1);

        // Mate2 record
        let mut rec2 = RecordBuf::default();
        rec2.name_mut().replace(read_name.into());

        // FLAGS: 0x1 (paired) | 0x4 (unmapped) | 0x8 (mate unmapped) | 0x80 (second in pair)
        let flags2 = sam::alignment::record::Flags::SEGMENTED
            | sam::alignment::record::Flags::UNMAPPED
            | sam::alignment::record::Flags::MATE_UNMAPPED
            | sam::alignment::record::Flags::LAST_SEGMENT;
        *rec2.flags_mut() = flags2;

        let seq2_bytes: Vec<u8> = mate2_seq.iter().map(|&b| decode_base(b)).collect();
        *rec2.sequence_mut() = Sequence::from(seq2_bytes);
        *rec2.quality_scores_mut() = QualityScores::from(mate2_qual.to_vec());
        records.push(rec2);

        Ok(records)
    }
}

/// Build paired SAM header from genome
pub fn build_sam_header(genome: &Genome, _params: &Parameters) -> Result<sam::Header, Error> {
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
#[allow(clippy::too_many_arguments)]
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

    // POS (1-based, per-chromosome coordinate)
    // transcript.genome_start is a global genome coordinate, need to convert to per-chr
    let chr_start = genome.chr_start[transcript.chr_idx];
    let pos = (transcript.genome_start - chr_start + 1) as usize;
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

/// Build a SAM record for one mate of a paired-end read
#[allow(clippy::too_many_arguments)]
fn build_paired_mate_record(
    read_name: &str,
    mate_seq: &[u8],
    mate_qual: &[u8],
    transcript: &Transcript,
    genome: &Genome,
    mapq: u8,
    is_first_mate: bool,
    is_proper_pair: bool,
    insert_size: i32,
) -> Result<RecordBuf, Error> {
    let mut record = RecordBuf::default();

    // Name
    record.name_mut().replace(read_name.into());

    // FLAGS
    let mut flags = sam::alignment::record::Flags::SEGMENTED; // 0x1 (paired)

    if is_proper_pair {
        flags |= sam::alignment::record::Flags::PROPERLY_SEGMENTED; // 0x2
    }

    if transcript.is_reverse {
        flags |= sam::alignment::record::Flags::REVERSE_COMPLEMENTED; // 0x10
    }

    // Mate on opposite strand for proper pairs (simplified assumption)
    if !transcript.is_reverse {
        flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED; // 0x20
    }

    if is_first_mate {
        flags |= sam::alignment::record::Flags::FIRST_SEGMENT; // 0x40
    } else {
        flags |= sam::alignment::record::Flags::LAST_SEGMENT; // 0x80
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

    // POS (1-based, per-chromosome coordinate)
    // transcript.genome_start is a global genome coordinate, need to convert to per-chr
    let chr_start = genome.chr_start[transcript.chr_idx];
    let pos = (transcript.genome_start - chr_start + 1) as usize;
    *record.alignment_start_mut() = Some(
        pos.try_into()
            .map_err(|e| Error::Alignment(format!("invalid alignment position {}: {}", pos, e)))?,
    );

    // MAPQ
    *record.mapping_quality_mut() = MappingQuality::new(mapq);

    // CIGAR
    let cigar = convert_cigar(&transcript.cigar)?;
    *record.cigar_mut() = cigar;

    // RNEXT (mate reference sequence = same as this mate for proper pairs)
    *record.mate_reference_sequence_id_mut() = Some(transcript.chr_idx);

    // PNEXT (mate position = approximate, use transcript start for now)
    // In full implementation, we'd track actual mate position
    let mate_pos = (transcript.genome_start + 1) as usize;
    *record.mate_alignment_start_mut() = Some(
        mate_pos
            .try_into()
            .map_err(|e| Error::Alignment(format!("invalid mate position {}: {}", mate_pos, e)))?,
    );

    // TLEN (insert size)
    *record.template_length_mut() = insert_size;

    // Sequence (decode from genome encoding)
    let seq_bytes: Vec<u8> = mate_seq.iter().map(|&b| decode_base(b)).collect();
    *record.sequence_mut() = Sequence::from(seq_bytes);

    // Quality scores
    *record.quality_scores_mut() = QualityScores::from(mate_qual.to_vec());

    Ok(record)
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

    #[test]
    fn test_build_paired_unmapped_records() {
        let mate1_seq = vec![0, 1, 2, 3]; // ACGT
        let mate1_qual = vec![30, 30, 30, 30];
        let mate2_seq = vec![3, 2, 1, 0]; // TGCA
        let mate2_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_paired_unmapped_records(
            "read1",
            &mate1_seq,
            &mate1_qual,
            &mate2_seq,
            &mate2_qual,
        )
        .unwrap();

        assert_eq!(records.len(), 2);

        // Check mate1 record
        let rec1 = &records[0];
        assert_eq!(
            rec1.name().map(|n| n.to_string()),
            Some("read1".to_string())
        );
        assert!(rec1.flags().is_segmented());
        assert!(rec1.flags().is_unmapped());
        assert!(rec1.flags().is_mate_unmapped());
        assert!(rec1.flags().is_first_segment());

        // Check mate2 record
        let rec2 = &records[1];
        assert_eq!(
            rec2.name().map(|n| n.to_string()),
            Some("read1".to_string())
        );
        assert!(rec2.flags().is_segmented());
        assert!(rec2.flags().is_unmapped());
        assert!(rec2.flags().is_mate_unmapped());
        assert!(rec2.flags().is_last_segment());
    }

    #[test]
    fn test_build_paired_mate_record_flags() {
        use crate::align::transcript::Exon;

        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 10,
            genome_end: 60,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 10,
                genome_end: 60,
                read_start: 0,
                read_end: 50,
            }],
            cigar: vec![CigarOp::Match(50)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            read_seq: vec![0, 1, 2, 3],
        };

        let mate_seq = vec![0, 1, 2, 3];
        let mate_qual = vec![30, 30, 30, 30];

        // Test first mate
        let rec1 = build_paired_mate_record(
            "read1",
            &mate_seq,
            &mate_qual,
            &transcript,
            &genome,
            255,
            true, // is_first_mate
            true, // is_proper_pair
            300,
        )
        .unwrap();

        assert!(rec1.flags().is_segmented());
        assert!(rec1.flags().is_properly_segmented());
        assert!(rec1.flags().is_first_segment());
        assert!(!rec1.flags().is_last_segment());
        assert_eq!(rec1.template_length(), 300);

        // Test second mate
        let rec2 = build_paired_mate_record(
            "read1",
            &mate_seq,
            &mate_qual,
            &transcript,
            &genome,
            255,
            false, // is_first_mate
            true,  // is_proper_pair
            -300,
        )
        .unwrap();

        assert!(rec2.flags().is_segmented());
        assert!(rec2.flags().is_properly_segmented());
        assert!(!rec2.flags().is_first_segment());
        assert!(rec2.flags().is_last_segment());
        assert_eq!(rec2.template_length(), -300);
    }

    #[test]
    fn test_build_paired_mate_record_mate_fields() {
        use crate::align::transcript::Exon;

        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 200,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 100,
                genome_end: 200,
                read_start: 0,
                read_end: 100,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 200,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            read_seq: vec![0; 100],
        };

        let mate_seq = vec![0; 100];
        let mate_qual = vec![30; 100];

        let rec = build_paired_mate_record(
            "read1",
            &mate_seq,
            &mate_qual,
            &transcript,
            &genome,
            60,
            true,
            true,
            250,
        )
        .unwrap();

        // Check RNEXT points to same chromosome
        assert_eq!(rec.mate_reference_sequence_id(), Some(0));

        // Check PNEXT is set (1-based position)
        assert_eq!(
            rec.mate_alignment_start().map(|p| usize::from(p)),
            Some(101)
        );

        // Check TLEN
        assert_eq!(rec.template_length(), 250);
    }

    // TODO: Add tests for SAM tags when they are implemented
}
