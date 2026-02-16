/// SAM/BAM output writer with noodles
use crate::align::read_align::PairedAlignment;
use crate::align::transcript::{CigarOp, Transcript};
use crate::error::Error;
use crate::genome::Genome;
use crate::io::fastq::{complement_base, decode_base};
use crate::junction::encode_motif;
use crate::mapq::calculate_mapq;
use crate::params::Parameters;
use bstr::BString;
use noodles::sam;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::data::field::value::Array;
use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};
use noodles::sam::header::record::value::{Map, map::Program};
use std::collections::HashSet;
use std::fmt::Write as FmtWrite;
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
        n_for_mapq: usize,
    ) -> Result<(), Error> {
        if transcripts.is_empty() {
            return Ok(());
        }

        let n_alignments = transcripts.len();
        let max_output = if params.out_sam_mult_nmax < 0 {
            n_alignments
        } else {
            (params.out_sam_mult_nmax as usize).min(n_alignments)
        };
        let effective_n = n_alignments.max(n_for_mapq);
        let mapq = calculate_mapq(effective_n, params.out_sam_mapq_unique);
        let mut attrs = params.sam_attribute_set();
        if params.out_sam_strand_field != "intronMotif" {
            attrs.remove("XS");
        }

        for (hit_index, transcript) in transcripts.iter().take(max_output).enumerate() {
            let record = transcript_to_record(
                transcript,
                read_name,
                read_seq,
                read_qual,
                genome,
                mapq,
                max_output,    // NH = number of reported alignments
                hit_index + 1, // 1-based
                &attrs,
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
        n_for_mapq: usize,
    ) -> Result<Vec<RecordBuf>, Error> {
        if transcripts.is_empty() {
            return Ok(Vec::new());
        }

        let n_alignments = transcripts.len();
        let max_output = if params.out_sam_mult_nmax < 0 {
            n_alignments
        } else {
            (params.out_sam_mult_nmax as usize).min(n_alignments)
        };
        let effective_n = n_alignments.max(n_for_mapq);
        let mapq = calculate_mapq(effective_n, params.out_sam_mapq_unique);
        let mut attrs = params.sam_attribute_set();
        if params.out_sam_strand_field != "intronMotif" {
            attrs.remove("XS");
        }

        let mut records = Vec::with_capacity(max_output);
        for (hit_index, transcript) in transcripts.iter().take(max_output).enumerate() {
            let record = transcript_to_record(
                transcript,
                read_name,
                read_seq,
                read_qual,
                genome,
                mapq,
                max_output,    // NH = number of reported alignments
                hit_index + 1, // 1-based
                &attrs,
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
        n_for_mapq: usize,
    ) -> Result<Vec<RecordBuf>, Error> {
        if paired_alignments.is_empty() {
            // Both mates unmapped
            return Self::build_paired_unmapped_records(
                read_name, mate1_seq, mate1_qual, mate2_seq, mate2_qual,
            );
        }

        let n_alignments = paired_alignments.len();
        let max_output = if params.out_sam_mult_nmax < 0 {
            n_alignments
        } else {
            (params.out_sam_mult_nmax as usize).min(n_alignments)
        };
        let effective_n = n_alignments.max(n_for_mapq);
        let mapq = calculate_mapq(effective_n, params.out_sam_mapq_unique);
        let mut attrs = params.sam_attribute_set();
        if params.out_sam_strand_field != "intronMotif" {
            attrs.remove("XS");
        }

        let mut records = Vec::with_capacity(max_output * 2);

        for (pair_idx, paired_aln) in paired_alignments.iter().take(max_output).enumerate() {
            let hit_index = pair_idx + 1; // 1-based

            // Create record for mate1 (this=mate1, mate=mate2)
            let rec1 = build_paired_mate_record(
                read_name,
                mate1_seq,
                mate1_qual,
                &paired_aln.mate1_transcript,
                &paired_aln.mate2_transcript,
                genome,
                mapq,
                true, // is_first_mate
                paired_aln.is_proper_pair,
                paired_aln.insert_size,
                max_output, // NH = number of reported alignments
                hit_index,
                &attrs,
            )?;
            records.push(rec1);

            // Create record for mate2 (this=mate2, mate=mate1)
            let rec2 = build_paired_mate_record(
                read_name,
                mate2_seq,
                mate2_qual,
                &paired_aln.mate2_transcript,
                &paired_aln.mate1_transcript,
                genome,
                mapq,
                false, // is_first_mate
                paired_aln.is_proper_pair,
                -paired_aln.insert_size, // Negative for mate2
                max_output,              // NH = number of reported alignments
                hit_index,
                &attrs,
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
    n_alignments: usize,
    hit_index: usize,
    attrs: &HashSet<String>,
) -> Result<RecordBuf, Error> {
    let mut record = RecordBuf::default();

    // Name
    record.name_mut().replace(read_name.into());

    // FLAGS
    let mut flags = sam::alignment::record::Flags::empty();
    if transcript.is_reverse {
        flags |= sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
    }
    if hit_index > 1 {
        flags |= sam::alignment::record::Flags::SECONDARY;
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

    // Sequence and quality scores
    // Per SAM spec: when FLAG & 16 (reverse strand), SEQ is the reverse complement
    // of the original read, and QUAL is reversed.
    if transcript.is_reverse {
        // Reverse complement the sequence
        let seq_bytes: Vec<u8> = read_seq
            .iter()
            .rev()
            .map(|&b| decode_base(complement_base(b)))
            .collect();
        *record.sequence_mut() = Sequence::from(seq_bytes);

        // Reverse the quality scores
        let mut qual = read_qual.to_vec();
        qual.reverse();
        *record.quality_scores_mut() = QualityScores::from(qual);
    } else {
        let seq_bytes: Vec<u8> = read_seq.iter().map(|&b| decode_base(b)).collect();
        *record.sequence_mut() = Sequence::from(seq_bytes);
        *record.quality_scores_mut() = QualityScores::from(read_qual.to_vec());
    }

    // Optional tags: gated by --outSAMattributes
    let data = record.data_mut();
    if attrs.contains("NH") {
        data.insert(Tag::ALIGNMENT_HIT_COUNT, Value::from(n_alignments as i32));
    }
    if attrs.contains("HI") {
        data.insert(Tag::HIT_INDEX, Value::from(hit_index as i32));
    }
    if attrs.contains("AS") {
        data.insert(Tag::ALIGNMENT_SCORE, Value::from(transcript.score));
    }
    if attrs.contains("NM") {
        data.insert(
            Tag::EDIT_DISTANCE,
            Value::from(compute_edit_distance(transcript)),
        );
    }
    if attrs.contains("nM") {
        data.insert(
            Tag::new(b'n', b'M'),
            Value::from(transcript.n_mismatch as i32),
        );
    }
    if attrs.contains("XS") {
        if let Some(xs_strand) = derive_xs_strand(transcript) {
            data.insert(Tag::new(b'X', b'S'), Value::Character(xs_strand as u8));
        }
    }
    if attrs.contains("jM") {
        if let Some(jm) = build_jm_tag(transcript) {
            data.insert(Tag::new(b'j', b'M'), jm);
        }
    }
    if attrs.contains("jI") {
        if let Some(ji) = build_ji_tag(transcript, chr_start) {
            data.insert(Tag::new(b'j', b'I'), ji);
        }
    }
    if attrs.contains("MD") {
        let md = build_md_tag(transcript, read_seq, genome, transcript.is_reverse);
        data.insert(Tag::new(b'M', b'D'), Value::String(BString::from(md)));
    }

    Ok(record)
}

/// Compute edit distance (NM tag): mismatches + inserted bases + deleted bases
fn compute_edit_distance(transcript: &Transcript) -> i32 {
    let indel_bases: u32 = transcript
        .cigar
        .iter()
        .filter_map(|op| match op {
            CigarOp::Ins(n) | CigarOp::Del(n) => Some(*n),
            _ => None,
        })
        .sum();
    (transcript.n_mismatch + indel_bases) as i32
}

/// Derive XS strand tag from transcript junction motifs.
/// Returns Some('+') or Some('-') if all junctions agree on strand.
/// Returns None if no junctions, all non-canonical, or conflicting strands.
fn derive_xs_strand(transcript: &Transcript) -> Option<char> {
    let mut strand: Option<char> = None;
    for motif in &transcript.junction_motifs {
        if let Some(s) = motif.implied_strand() {
            match strand {
                None => strand = Some(s),
                Some(prev) if prev != s => return None,
                _ => {}
            }
        }
    }
    strand
}

/// Build jM tag: array of junction motif codes (one per intron/RefSkip in CIGAR).
///
/// Encoding: 0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG, 4=CT/GC, 5=AT/AC, 6=GT/AT.
/// Add +20 if junction is annotated in GTF.
fn build_jm_tag(transcript: &Transcript) -> Option<Value> {
    if transcript.junction_motifs.is_empty() {
        return None;
    }
    let motifs: Vec<i8> = transcript
        .junction_motifs
        .iter()
        .zip(
            transcript
                .junction_annotated
                .iter()
                .chain(std::iter::repeat(&false)),
        )
        .map(|(motif, &annotated)| {
            let code = encode_motif(*motif) as i8;
            if annotated { code + 20 } else { code }
        })
        .collect();
    Some(Value::Array(Array::Int8(motifs)))
}

/// Build jI tag: array of intron start/end coordinates (1-based, per-chromosome).
///
/// Format: [start1, end1, start2, end2, ...] where start is first intronic base
/// and end is last intronic base (both 1-based, inclusive).
fn build_ji_tag(transcript: &Transcript, chr_start: u64) -> Option<Value> {
    if transcript.n_junction == 0 {
        return None;
    }
    let mut coords: Vec<i32> = Vec::new();
    let mut genome_pos = transcript.genome_start;
    for op in &transcript.cigar {
        match op {
            CigarOp::RefSkip(n) => {
                let intron_start = (genome_pos - chr_start + 1) as i32; // 1-based
                let intron_end = (genome_pos + *n as u64 - chr_start) as i32; // 1-based inclusive
                coords.push(intron_start);
                coords.push(intron_end);
                genome_pos += *n as u64;
            }
            CigarOp::Match(n) | CigarOp::Equal(n) | CigarOp::Diff(n) | CigarOp::Del(n) => {
                genome_pos += *n as u64;
            }
            _ => {} // Ins, SoftClip, HardClip don't consume reference
        }
    }
    Some(Value::Array(Array::Int32(coords)))
}

/// Build MD tag string: matches/mismatches/deletions relative to reference.
///
/// Format: "10A5^AC6" = 10 match, A mismatch, 5 match, 2bp deletion (AC), 6 match.
/// The MD tag describes the reference sequence for positions that differ from the read.
fn build_md_tag(
    transcript: &Transcript,
    read_seq: &[u8],
    genome: &Genome,
    is_reverse: bool,
) -> String {
    // Build the SAM-order sequence (RC for reverse strand)
    let sam_seq: Vec<u8> = if is_reverse {
        read_seq.iter().rev().map(|&b| complement_base(b)).collect()
    } else {
        read_seq.to_vec()
    };

    let mut md = String::new();
    let mut match_count: u32 = 0;
    let mut genome_pos = transcript.genome_start;
    let mut read_pos: usize = 0;

    for op in &transcript.cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Equal(n) | CigarOp::Diff(n) => {
                for _ in 0..*n {
                    let ref_base = genome.get_base(genome_pos).unwrap_or(4);
                    let read_base = if read_pos < sam_seq.len() {
                        sam_seq[read_pos]
                    } else {
                        4
                    };
                    if ref_base == read_base {
                        match_count += 1;
                    } else {
                        write!(md, "{}", match_count).unwrap();
                        match_count = 0;
                        md.push(decode_base(ref_base) as char);
                    }
                    genome_pos += 1;
                    read_pos += 1;
                }
            }
            CigarOp::Del(n) => {
                write!(md, "{}", match_count).unwrap();
                match_count = 0;
                md.push('^');
                for _ in 0..*n {
                    let ref_base = genome.get_base(genome_pos).unwrap_or(4);
                    md.push(decode_base(ref_base) as char);
                    genome_pos += 1;
                }
            }
            CigarOp::Ins(n) => {
                read_pos += *n as usize;
            }
            CigarOp::RefSkip(n) => {
                genome_pos += *n as u64;
            }
            CigarOp::SoftClip(n) => {
                read_pos += *n as usize;
            }
            CigarOp::HardClip(_) => {}
        }
    }
    // Emit trailing match count
    write!(md, "{}", match_count).unwrap();
    md
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
    mate_transcript: &Transcript,
    genome: &Genome,
    mapq: u8,
    is_first_mate: bool,
    is_proper_pair: bool,
    insert_size: i32,
    n_alignments: usize,
    hit_index: usize,
    attrs: &HashSet<String>,
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

    // Mate reverse flag from the actual mate's alignment strand
    if mate_transcript.is_reverse {
        flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED; // 0x20
    }

    if is_first_mate {
        flags |= sam::alignment::record::Flags::FIRST_SEGMENT; // 0x40
    } else {
        flags |= sam::alignment::record::Flags::LAST_SEGMENT; // 0x80
    }

    if hit_index > 1 {
        flags |= sam::alignment::record::Flags::SECONDARY; // 0x100
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

    // RNEXT (mate reference sequence from the mate's actual alignment)
    *record.mate_reference_sequence_id_mut() = Some(mate_transcript.chr_idx);

    // PNEXT (mate position from the mate's actual alignment, per-chromosome coords)
    let mate_chr_start = genome.chr_start[mate_transcript.chr_idx];
    let mate_pos = (mate_transcript.genome_start - mate_chr_start + 1) as usize;
    *record.mate_alignment_start_mut() = Some(
        mate_pos
            .try_into()
            .map_err(|e| Error::Alignment(format!("invalid mate position {}: {}", mate_pos, e)))?,
    );

    // TLEN (insert size)
    *record.template_length_mut() = insert_size;

    // Sequence and quality scores (reverse complement for reverse strand)
    if transcript.is_reverse {
        let seq_bytes: Vec<u8> = mate_seq
            .iter()
            .rev()
            .map(|&b| decode_base(complement_base(b)))
            .collect();
        *record.sequence_mut() = Sequence::from(seq_bytes);

        let mut qual = mate_qual.to_vec();
        qual.reverse();
        *record.quality_scores_mut() = QualityScores::from(qual);
    } else {
        let seq_bytes: Vec<u8> = mate_seq.iter().map(|&b| decode_base(b)).collect();
        *record.sequence_mut() = Sequence::from(seq_bytes);
        *record.quality_scores_mut() = QualityScores::from(mate_qual.to_vec());
    }

    // Optional tags: gated by --outSAMattributes
    let data = record.data_mut();
    if attrs.contains("NH") {
        data.insert(Tag::ALIGNMENT_HIT_COUNT, Value::from(n_alignments as i32));
    }
    if attrs.contains("HI") {
        data.insert(Tag::HIT_INDEX, Value::from(hit_index as i32));
    }
    if attrs.contains("AS") {
        data.insert(Tag::ALIGNMENT_SCORE, Value::from(transcript.score));
    }
    if attrs.contains("NM") {
        data.insert(
            Tag::EDIT_DISTANCE,
            Value::from(compute_edit_distance(transcript)),
        );
    }
    if attrs.contains("nM") {
        data.insert(
            Tag::new(b'n', b'M'),
            Value::from(transcript.n_mismatch as i32),
        );
    }
    if attrs.contains("XS") {
        if let Some(xs_strand) = derive_xs_strand(transcript) {
            data.insert(Tag::new(b'X', b'S'), Value::Character(xs_strand as u8));
        }
    }
    if attrs.contains("jM") {
        if let Some(jm) = build_jm_tag(transcript) {
            data.insert(Tag::new(b'j', b'M'), jm);
        }
    }
    if attrs.contains("jI") {
        if let Some(ji) = build_ji_tag(transcript, chr_start) {
            data.insert(Tag::new(b'j', b'I'), ji);
        }
    }
    if attrs.contains("MD") {
        let md = build_md_tag(transcript, mate_seq, genome, transcript.is_reverse);
        data.insert(Tag::new(b'M', b'D'), Value::String(BString::from(md)));
    }

    Ok(record)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::score::SpliceMotif;
    use crate::genome::Genome;
    use clap::Parser;
    use tempfile::NamedTempFile;

    /// Build an attribute set with all tags enabled (for tests that don't care about filtering)
    fn all_attrs() -> HashSet<String> {
        ["NH", "HI", "AS", "NM", "nM", "XS", "jM", "jI", "MD"]
            .iter()
            .map(|s| s.to_string())
            .collect()
    }

    /// Build the standard attribute set (NH, HI, AS, NM, nM)
    fn standard_attrs() -> HashSet<String> {
        ["NH", "HI", "AS", "NM", "nM"]
            .iter()
            .map(|s| s.to_string())
            .collect()
    }

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
            junction_motifs: vec![],
            junction_annotated: vec![],
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
            &standard_attrs(),
        );
        assert!(record.is_ok());

        let record = record.unwrap();
        assert_eq!(
            record.name().map(|n| n.to_string()),
            Some("read1".to_string())
        );
        assert_eq!(record.reference_sequence_id(), Some(0));
        assert_eq!(record.alignment_start().map(|p| usize::from(p)), Some(11)); // 1-based
        // hit_index=1, so NOT secondary
        assert!(!record.flags().is_secondary());
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

        let mate1_transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
            }],
            cigar: vec![CigarOp::Match(4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        let mate2_transcript = Transcript {
            chr_idx: 0,
            genome_start: 4,
            genome_end: 7,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 4,
                genome_end: 7,
                read_start: 0,
                read_end: 3,
            }],
            cigar: vec![CigarOp::Match(3)],
            score: 90,
            n_mismatch: 1,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2],
        };

        let mate_seq = vec![0, 1, 2, 3];
        let mate_qual = vec![30, 30, 30, 30];

        // Test first mate (forward), mate2 is reverse → 0x20 should be set
        let rec1 = build_paired_mate_record(
            "read1",
            &mate_seq,
            &mate_qual,
            &mate1_transcript,
            &mate2_transcript,
            &genome,
            255,
            true, // is_first_mate
            true, // is_proper_pair
            300,
            1, // n_alignments
            1, // hit_index
            &standard_attrs(),
        )
        .unwrap();

        assert!(rec1.flags().is_segmented());
        assert!(rec1.flags().is_properly_segmented());
        assert!(rec1.flags().is_first_segment());
        assert!(!rec1.flags().is_last_segment());
        assert!(!rec1.flags().is_reverse_complemented()); // mate1 is forward
        assert!(rec1.flags().is_mate_reverse_complemented()); // mate2 is reverse
        assert!(!rec1.flags().is_secondary());
        assert_eq!(rec1.template_length(), 300);

        // Test second mate (reverse), mate1 is forward → 0x20 should NOT be set
        let rec2 = build_paired_mate_record(
            "read1",
            &mate_seq,
            &mate_qual,
            &mate2_transcript,
            &mate1_transcript,
            &genome,
            255,
            false, // is_first_mate
            true,  // is_proper_pair
            -300,
            1, // n_alignments
            1, // hit_index
            &standard_attrs(),
        )
        .unwrap();

        assert!(rec2.flags().is_segmented());
        assert!(rec2.flags().is_properly_segmented());
        assert!(!rec2.flags().is_first_segment());
        assert!(rec2.flags().is_last_segment());
        assert!(rec2.flags().is_reverse_complemented()); // mate2 is reverse
        assert!(!rec2.flags().is_mate_reverse_complemented()); // mate1 is forward
        assert!(!rec2.flags().is_secondary());
        assert_eq!(rec2.template_length(), -300);
    }

    #[test]
    fn test_build_paired_mate_record_mate_fields() {
        use crate::align::transcript::Exon;

        let genome = make_test_genome();

        // Mate1 at position 0 (chr_start=0, so per-chr pos = 1)
        let this_transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
            }],
            cigar: vec![CigarOp::Match(4)],
            score: 200,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 4],
        };

        // Mate2 at position 4 (chr_start=0, so per-chr pos = 5)
        let mate_transcript = Transcript {
            chr_idx: 0,
            genome_start: 4,
            genome_end: 7,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 4,
                genome_end: 7,
                read_start: 0,
                read_end: 3,
            }],
            cigar: vec![CigarOp::Match(3)],
            score: 150,
            n_mismatch: 1,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 3],
        };

        let mate_seq = vec![0; 4];
        let mate_qual = vec![30; 4];

        let rec = build_paired_mate_record(
            "read1",
            &mate_seq,
            &mate_qual,
            &this_transcript,
            &mate_transcript,
            &genome,
            60,
            true,
            true,
            250,
            1, // n_alignments
            1, // hit_index
            &standard_attrs(),
        )
        .unwrap();

        // RNEXT = mate's chr_idx
        assert_eq!(rec.mate_reference_sequence_id(), Some(0));

        // PNEXT = mate's per-chr position (genome_start=4, chr_start=0 → pos=5)
        assert_eq!(rec.mate_alignment_start().map(|p| usize::from(p)), Some(5));

        // Check TLEN
        assert_eq!(rec.template_length(), 250);

        // Tags should reflect THIS transcript (mate1), not the mate
        let data = rec.data();
        assert_eq!(
            data.get(&Tag::ALIGNMENT_SCORE),
            Some(&Value::from(200_i32)),
            "AS should be from this transcript (200), not mate (150)"
        );
        assert_eq!(
            data.get(&Tag::EDIT_DISTANCE),
            Some(&Value::from(0_i32)),
            "NM should be from this transcript (0), not mate (1)"
        );
    }

    #[test]
    fn test_tags_nh_hi_as_nm() {
        let genome = make_test_genome();

        // Transcript with 2 mismatches and a 3bp deletion → NM = 2 + 3 = 5
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 60,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(20), CigarOp::Del(3), CigarOp::Match(30)],
            score: 100,
            n_mismatch: 2,
            n_gap: 1,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            3, // n_alignments
            2, // hit_index
            &standard_attrs(),
        )
        .unwrap();

        // hit_index=2 → secondary
        assert!(record.flags().is_secondary());

        let data = record.data();
        assert_eq!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT),
            Some(&Value::from(3_i32)),
            "NH tag should be 3"
        );
        assert_eq!(
            data.get(&Tag::HIT_INDEX),
            Some(&Value::from(2_i32)),
            "HI tag should be 2"
        );
        assert_eq!(
            data.get(&Tag::ALIGNMENT_SCORE),
            Some(&Value::from(100_i32)),
            "AS tag should be 100"
        );
        assert_eq!(
            data.get(&Tag::EDIT_DISTANCE),
            Some(&Value::from(5_i32)),
            "NM tag should be 5 (2 mismatches + 3 deleted bases)"
        );
        assert_eq!(
            data.get(&Tag::new(b'n', b'M')),
            Some(&Value::from(2_i32)),
            "nM tag should be 2 (mismatches only, no indels)"
        );
    }

    #[test]
    fn test_edit_distance_computation() {
        // Pure match: NM = n_mismatch only
        let t1 = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 50,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(50)],
            score: 100,
            n_mismatch: 3,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };
        assert_eq!(compute_edit_distance(&t1), 3);

        // Match + Ins + Del: NM = n_mismatch + ins_len + del_len
        let t2 = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 60,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(20),
                CigarOp::Ins(5),
                CigarOp::Match(10),
                CigarOp::Del(7),
                CigarOp::Match(20),
            ],
            score: 80,
            n_mismatch: 1,
            n_gap: 2,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };
        assert_eq!(compute_edit_distance(&t2), 13); // 1 + 5 + 7

        // RefSkip (splice junction) should NOT count toward NM
        let t3 = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(25),
                CigarOp::RefSkip(1000),
                CigarOp::Match(25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };
        assert_eq!(compute_edit_distance(&t3), 0);

        // Soft clips should NOT count toward NM
        let t4 = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 40,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::SoftClip(10),
                CigarOp::Match(40),
                CigarOp::SoftClip(10),
            ],
            score: 40,
            n_mismatch: 2,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };
        assert_eq!(compute_edit_distance(&t4), 2);
    }

    #[test]
    fn test_transcript_to_record_has_tags() {
        // Verify the existing test_transcript_to_record scenario also has tags
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
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1, // unique mapper
            1,
            &standard_attrs(),
        )
        .unwrap();

        let data = record.data();
        // Unique mapper: NH=1, HI=1
        assert_eq!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT),
            Some(&Value::from(1_i32))
        );
        assert_eq!(data.get(&Tag::HIT_INDEX), Some(&Value::from(1_i32)));
        assert_eq!(data.get(&Tag::ALIGNMENT_SCORE), Some(&Value::from(100_i32)));
        // NM = 2 mismatches, no indels
        assert_eq!(data.get(&Tag::EDIT_DISTANCE), Some(&Value::from(2_i32)));
        // nM = 2 mismatches only (same as NM when no indels)
        assert_eq!(data.get(&Tag::new(b'n', b'M')), Some(&Value::from(2_i32)));
        // XS not in standard attrs
        assert_eq!(data.get(&Tag::new(b'X', b'S')), None);
    }

    #[test]
    fn test_secondary_flag() {
        let genome = make_test_genome();
        let params = Parameters::parse_from(vec!["ruSTAR", "--readFilesIn", "test.fq"]);

        let transcripts = vec![
            Transcript {
                chr_idx: 0,
                genome_start: 0,
                genome_end: 50,
                is_reverse: false,
                exons: vec![],
                cigar: vec![CigarOp::Match(50)],
                score: 100,
                n_mismatch: 0,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
                read_seq: vec![0; 4],
            },
            Transcript {
                chr_idx: 0,
                genome_start: 2,
                genome_end: 52,
                is_reverse: false,
                exons: vec![],
                cigar: vec![CigarOp::Match(50)],
                score: 98,
                n_mismatch: 1,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
                read_seq: vec![0; 4],
            },
            Transcript {
                chr_idx: 0,
                genome_start: 4,
                genome_end: 54,
                is_reverse: true,
                exons: vec![],
                cigar: vec![CigarOp::Match(50)],
                score: 96,
                n_mismatch: 2,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
                read_seq: vec![0; 4],
            },
        ];

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_alignment_records(
            "read1",
            &read_seq,
            &read_qual,
            &transcripts,
            &genome,
            &params,
            1,
        )
        .unwrap();

        assert_eq!(records.len(), 3);

        // Record 0 (HI=1): NOT secondary
        assert!(!records[0].flags().is_secondary());
        // Record 1 (HI=2): IS secondary
        assert!(records[1].flags().is_secondary());
        // Record 2 (HI=3): IS secondary + reverse complemented
        assert!(records[2].flags().is_secondary());
        assert!(records[2].flags().is_reverse_complemented());
    }

    #[test]
    fn test_xs_tag_spliced() {
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(25),
                CigarOp::RefSkip(100),
                CigarOp::Match(25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![false],
            read_seq: vec![0; 4],
        };

        let read_seq = vec![0, 1, 2, 3];
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
            &all_attrs(),
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            Some(&Value::Character(b'+')),
            "XS should be '+' for GT/AG motif"
        );
    }

    #[test]
    fn test_xs_tag_unspliced() {
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 50,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(50)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 4],
        };

        let read_seq = vec![0, 1, 2, 3];
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
            &all_attrs(),
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            None,
            "XS should NOT be present for unspliced reads"
        );
    }

    #[test]
    fn test_xs_tag_reverse_strand() {
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(25),
                CigarOp::RefSkip(100),
                CigarOp::Match(25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::CtAc],
            junction_annotated: vec![false],
            read_seq: vec![0; 4],
        };

        let read_seq = vec![0, 1, 2, 3];
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
            &all_attrs(),
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            Some(&Value::Character(b'-')),
            "XS should be '-' for CT/AC motif"
        );
    }

    #[test]
    fn test_xs_tag_conflicting_motifs() {
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 300,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(25),
                CigarOp::RefSkip(100),
                CigarOp::Match(25),
                CigarOp::RefSkip(100),
                CigarOp::Match(25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::CtAc], // +strand and -strand
            junction_annotated: vec![false, false],
            read_seq: vec![0; 4],
        };

        let read_seq = vec![0, 1, 2, 3];
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
            &all_attrs(),
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            None,
            "XS should NOT be present when junction motifs conflict on strand"
        );
    }

    #[test]
    fn test_xs_not_emitted_when_disabled() {
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(25),
                CigarOp::RefSkip(100),
                CigarOp::Match(25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![false],
            read_seq: vec![0; 4],
        };

        let read_seq = vec![0, 1, 2, 3];
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
            &standard_attrs(), // XS not in standard attrs
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            None,
            "XS should NOT be present when not in attribute set"
        );
    }

    #[test]
    fn test_out_sam_mult_nmax() {
        let genome = make_test_genome();
        let params = Parameters::parse_from(vec![
            "ruSTAR",
            "--readFilesIn",
            "test.fq",
            "--outSAMmultNmax",
            "3",
        ]);

        let transcripts: Vec<Transcript> = (0..5)
            .map(|i| Transcript {
                chr_idx: 0,
                genome_start: i as u64,
                genome_end: (i + 50) as u64,
                is_reverse: false,
                exons: vec![],
                cigar: vec![CigarOp::Match(50)],
                score: 100 - i,
                n_mismatch: 0,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
                read_seq: vec![0; 4],
            })
            .collect();

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_alignment_records(
            "read1",
            &read_seq,
            &read_qual,
            &transcripts,
            &genome,
            &params,
            1,
        )
        .unwrap();

        // Only 3 records output despite 5 transcripts
        assert_eq!(records.len(), 3);

        // NH should be 3 (number of reported alignments)
        for rec in &records {
            let data = rec.data();
            assert_eq!(
                data.get(&Tag::ALIGNMENT_HIT_COUNT),
                Some(&Value::from(3_i32)),
                "NH should be 3 (capped by outSAMmultNmax)"
            );
        }

        // HI should be 1, 2, 3
        assert_eq!(
            records[0].data().get(&Tag::HIT_INDEX),
            Some(&Value::from(1_i32))
        );
        assert_eq!(
            records[1].data().get(&Tag::HIT_INDEX),
            Some(&Value::from(2_i32))
        );
        assert_eq!(
            records[2].data().get(&Tag::HIT_INDEX),
            Some(&Value::from(3_i32))
        );

        // First is primary, rest are secondary
        assert!(!records[0].flags().is_secondary());
        assert!(records[1].flags().is_secondary());
        assert!(records[2].flags().is_secondary());
    }

    #[test]
    fn test_build_jm_tag_basic() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(25),
                CigarOp::RefSkip(100),
                CigarOp::Match(25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![false],
            read_seq: vec![],
        };

        let jm = build_jm_tag(&transcript);
        assert!(jm.is_some());
        // GT/AG = motif code 1, not annotated → 1
        assert_eq!(jm.unwrap(), Value::Array(Array::Int8(vec![1])));
    }

    #[test]
    fn test_build_jm_tag_annotated() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(25),
                CigarOp::RefSkip(100),
                CigarOp::Match(25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![true],
            read_seq: vec![],
        };

        let jm = build_jm_tag(&transcript);
        assert!(jm.is_some());
        // GT/AG = motif code 1, annotated → 1 + 20 = 21
        assert_eq!(jm.unwrap(), Value::Array(Array::Int8(vec![21])));
    }

    #[test]
    fn test_build_jm_tag_empty() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 50,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(50)],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        assert!(build_jm_tag(&transcript).is_none());
    }

    #[test]
    fn test_build_jm_tag_multiple_junctions() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 400,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(25),
                CigarOp::RefSkip(100),
                CigarOp::Match(25),
                CigarOp::RefSkip(100),
                CigarOp::Match(25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::CtAc],
            junction_annotated: vec![true, false],
            read_seq: vec![],
        };

        let jm = build_jm_tag(&transcript);
        assert!(jm.is_some());
        // GT/AG annotated=21, CT/AC not annotated=2
        assert_eq!(jm.unwrap(), Value::Array(Array::Int8(vec![21, 2])));
    }

    #[test]
    fn test_build_ji_tag_basic() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 325,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(25),
                CigarOp::RefSkip(200),
                CigarOp::Match(25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![false],
            read_seq: vec![],
        };

        // chr_start=0, genome_start=100, intron starts at 125, ends at 324
        let ji = build_ji_tag(&transcript, 0);
        assert!(ji.is_some());
        // Intron start: 100+25 - 0 + 1 = 126, Intron end: 100+25+200 - 0 = 325
        assert_eq!(ji.unwrap(), Value::Array(Array::Int32(vec![126, 325])));
    }

    #[test]
    fn test_build_ji_tag_empty() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 50,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(50)],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        assert!(build_ji_tag(&transcript, 0).is_none());
    }

    #[test]
    fn test_build_md_tag_perfect_match() {
        // Genome: ACGTACGT (A=0,C=1,G=2,T=3)
        let genome = make_test_genome();
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(4)],
            score: 4,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        // Read exactly matches genome[0..4] = ACGT
        let read_seq = vec![0, 1, 2, 3];
        let md = build_md_tag(&transcript, &read_seq, &genome, false);
        assert_eq!(md, "4");
    }

    #[test]
    fn test_build_md_tag_mismatches() {
        // Genome: ACGTACGT
        let genome = make_test_genome();
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(4)],
            score: 2,
            n_mismatch: 2,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 0, 2, 0], // A,A,G,A vs genome A,C,G,T
        };

        // Position 1: read=A, ref=C → mismatch (C in MD)
        // Position 3: read=A, ref=T → mismatch (T in MD)
        let read_seq = vec![0, 0, 2, 0]; // AAGA
        let md = build_md_tag(&transcript, &read_seq, &genome, false);
        assert_eq!(md, "1C1T0");
    }

    #[test]
    fn test_build_md_tag_deletion() {
        // Genome: ACGTACGT
        let genome = make_test_genome();
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 6,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(2), CigarOp::Del(2), CigarOp::Match(2)],
            score: 4,
            n_mismatch: 0,
            n_gap: 1,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 0, 1], // AC + AC (genome AC^GT AC)
        };

        // Read: A,C,[del G,T],A,C
        let read_seq = vec![0, 1, 0, 1];
        let md = build_md_tag(&transcript, &read_seq, &genome, false);
        assert_eq!(md, "2^GT2");
    }

    #[test]
    fn test_build_md_tag_insertion() {
        // Genome: ACGTACGT
        let genome = make_test_genome();
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(2), CigarOp::Ins(2), CigarOp::Match(2)],
            score: 4,
            n_mismatch: 0,
            n_gap: 1,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 3, 3, 2, 3], // AC + TT(ins) + GT
        };

        // Insertions are invisible in MD — just match counts
        let read_seq = vec![0, 1, 3, 3, 2, 3];
        let md = build_md_tag(&transcript, &read_seq, &genome, false);
        assert_eq!(md, "4");
    }

    #[test]
    fn test_build_md_tag_soft_clip() {
        // Genome: ACGTACGT
        let genome = make_test_genome();
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 2, // Starts at G
            genome_end: 6,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::SoftClip(2),
                CigarOp::Match(4),
                CigarOp::SoftClip(2),
            ],
            score: 4,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 0, 2, 3, 0, 1, 0, 0], // XX + GTAC + XX
        };

        // Soft clips don't appear in MD
        let read_seq = vec![0, 0, 2, 3, 0, 1, 0, 0];
        let md = build_md_tag(&transcript, &read_seq, &genome, false);
        assert_eq!(md, "4");
    }

    #[test]
    fn test_tags_jm_ji_md_in_record() {
        // Verify tags appear in a full transcript_to_record call
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(4)],
            score: 4,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        let read_seq = vec![0, 1, 2, 3];
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
            &all_attrs(),
        )
        .unwrap();

        let data = record.data();
        // No junctions → no jM/jI tags
        assert!(data.get(&Tag::new(b'j', b'M')).is_none());
        assert!(data.get(&Tag::new(b'j', b'I')).is_none());
        // MD should be present when in attrs
        assert_eq!(
            data.get(&Tag::new(b'M', b'D')),
            Some(&Value::String(BString::from("4")))
        );
    }

    #[test]
    fn test_build_paired_mate_record_cross_strand() {
        use crate::align::transcript::Exon;

        let genome = make_test_genome();

        // Mate1: forward, chr 0, pos 0
        let mate1_trans = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
            }],
            cigar: vec![CigarOp::Match(4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        // Mate2: reverse, chr 0, pos 4
        let mate2_trans = Transcript {
            chr_idx: 0,
            genome_start: 4,
            genome_end: 7,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 4,
                genome_end: 7,
                read_start: 0,
                read_end: 3,
            }],
            cigar: vec![CigarOp::Match(3)],
            score: 90,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2],
        };

        let seq = vec![0, 1, 2, 3];
        let qual = vec![30, 30, 30, 30];

        // Mate1 record: mate is reverse → 0x20 set, PNEXT=5
        let rec1 = build_paired_mate_record(
            "read1",
            &seq,
            &qual,
            &mate1_trans,
            &mate2_trans,
            &genome,
            255,
            true,
            true,
            7,
            1,
            1,
            &standard_attrs(),
        )
        .unwrap();

        assert!(rec1.flags().is_mate_reverse_complemented());
        assert!(!rec1.flags().is_reverse_complemented());
        assert_eq!(rec1.mate_alignment_start().map(|p| usize::from(p)), Some(5)); // genome_start=4, chr_start=0 → 5

        // Mate2 record: mate is forward → 0x20 NOT set, PNEXT=1
        let rec2 = build_paired_mate_record(
            "read1",
            &seq,
            &qual,
            &mate2_trans,
            &mate1_trans,
            &genome,
            255,
            false,
            true,
            -7,
            1,
            1,
            &standard_attrs(),
        )
        .unwrap();

        assert!(!rec2.flags().is_mate_reverse_complemented());
        assert!(rec2.flags().is_reverse_complemented());
        assert_eq!(rec2.mate_alignment_start().map(|p| usize::from(p)), Some(1)); // genome_start=0, chr_start=0 → 1
    }

    #[test]
    fn test_build_paired_mate_record_per_mate_tags() {
        use crate::align::transcript::Exon;

        let genome = make_test_genome();

        // Mate1: score=100, 0 mismatches, no junctions
        let mate1_trans = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
            }],
            cigar: vec![CigarOp::Match(4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        // Mate2: score=80, 2 mismatches, 1 deletion
        let mate2_trans = Transcript {
            chr_idx: 0,
            genome_start: 4,
            genome_end: 7,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 4,
                genome_end: 7,
                read_start: 0,
                read_end: 3,
            }],
            cigar: vec![CigarOp::Match(2), CigarOp::Del(1), CigarOp::Match(1)],
            score: 80,
            n_mismatch: 2,
            n_gap: 1,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2],
        };

        let seq1 = vec![0, 1, 2, 3];
        let qual1 = vec![30, 30, 30, 30];
        let seq2 = vec![0, 1, 2];
        let qual2 = vec![30, 30, 30];

        // Mate1 record should have mate1's tags
        let rec1 = build_paired_mate_record(
            "read1",
            &seq1,
            &qual1,
            &mate1_trans,
            &mate2_trans,
            &genome,
            255,
            true,
            true,
            7,
            1,
            1,
            &standard_attrs(),
        )
        .unwrap();
        assert_eq!(
            rec1.data().get(&Tag::ALIGNMENT_SCORE),
            Some(&Value::from(100_i32)),
            "Mate1 AS should be 100"
        );
        assert_eq!(
            rec1.data().get(&Tag::EDIT_DISTANCE),
            Some(&Value::from(0_i32)),
            "Mate1 NM should be 0"
        );

        // Mate2 record should have mate2's tags
        let rec2 = build_paired_mate_record(
            "read1",
            &seq2,
            &qual2,
            &mate2_trans,
            &mate1_trans,
            &genome,
            255,
            false,
            true,
            -7,
            1,
            1,
            &standard_attrs(),
        )
        .unwrap();
        assert_eq!(
            rec2.data().get(&Tag::ALIGNMENT_SCORE),
            Some(&Value::from(80_i32)),
            "Mate2 AS should be 80"
        );
        // NM = 2 mismatches + 1 deleted base = 3
        assert_eq!(
            rec2.data().get(&Tag::EDIT_DISTANCE),
            Some(&Value::from(3_i32)),
            "Mate2 NM should be 3 (2 mismatches + 1 del)"
        );
    }

    #[test]
    fn test_build_paired_mate_record_both_forward() {
        use crate::align::transcript::Exon;

        let genome = make_test_genome();

        // Both mates forward
        let mate1_trans = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
            }],
            cigar: vec![CigarOp::Match(4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        let mate2_trans = Transcript {
            chr_idx: 0,
            genome_start: 4,
            genome_end: 7,
            is_reverse: false, // Also forward
            exons: vec![Exon {
                genome_start: 4,
                genome_end: 7,
                read_start: 0,
                read_end: 3,
            }],
            cigar: vec![CigarOp::Match(3)],
            score: 90,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2],
        };

        let seq = vec![0, 1, 2, 3];
        let qual = vec![30, 30, 30, 30];

        // Both forward → neither should have 0x20 set
        let rec1 = build_paired_mate_record(
            "read1",
            &seq,
            &qual,
            &mate1_trans,
            &mate2_trans,
            &genome,
            255,
            true,
            true,
            7,
            1,
            1,
            &standard_attrs(),
        )
        .unwrap();
        assert!(!rec1.flags().is_mate_reverse_complemented());
        assert!(!rec1.flags().is_reverse_complemented());

        let rec2 = build_paired_mate_record(
            "read1",
            &seq,
            &qual,
            &mate2_trans,
            &mate1_trans,
            &genome,
            255,
            false,
            true,
            -7,
            1,
            1,
            &standard_attrs(),
        )
        .unwrap();
        assert!(!rec2.flags().is_mate_reverse_complemented());
        assert!(!rec2.flags().is_reverse_complemented());
    }

    #[test]
    fn test_out_sam_attributes_standard() {
        // Default (Standard) → NH, HI, AS, NM present; XS, jM, jI, MD absent
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(4)],
            score: 100,
            n_mismatch: 1,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        let read_seq = vec![0, 1, 2, 3];
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
            &standard_attrs(),
        )
        .unwrap();

        let data = record.data();
        // Standard tags present
        assert!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT).is_some(),
            "NH should be present"
        );
        assert!(data.get(&Tag::HIT_INDEX).is_some(), "HI should be present");
        assert!(
            data.get(&Tag::ALIGNMENT_SCORE).is_some(),
            "AS should be present"
        );
        assert!(
            data.get(&Tag::EDIT_DISTANCE).is_some(),
            "NM should be present"
        );
        assert!(
            data.get(&Tag::new(b'n', b'M')).is_some(),
            "nM should be present"
        );
        // Non-standard tags absent
        assert!(
            data.get(&Tag::new(b'X', b'S')).is_none(),
            "XS should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'M')).is_none(),
            "jM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'I')).is_none(),
            "jI should be absent"
        );
        assert!(
            data.get(&Tag::new(b'M', b'D')).is_none(),
            "MD should be absent"
        );
    }

    #[test]
    fn test_out_sam_attributes_none() {
        // None → no optional tags at all
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(4)],
            score: 100,
            n_mismatch: 1,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let empty_attrs: HashSet<String> = HashSet::new();
        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            &empty_attrs,
        )
        .unwrap();

        let data = record.data();
        assert!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT).is_none(),
            "NH should be absent"
        );
        assert!(data.get(&Tag::HIT_INDEX).is_none(), "HI should be absent");
        assert!(
            data.get(&Tag::ALIGNMENT_SCORE).is_none(),
            "AS should be absent"
        );
        assert!(
            data.get(&Tag::EDIT_DISTANCE).is_none(),
            "NM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'n', b'M')).is_none(),
            "nM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'X', b'S')).is_none(),
            "XS should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'M')).is_none(),
            "jM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'I')).is_none(),
            "jI should be absent"
        );
        assert!(
            data.get(&Tag::new(b'M', b'D')).is_none(),
            "MD should be absent"
        );
    }

    #[test]
    fn test_out_sam_attributes_explicit() {
        // Explicit ["NH", "MD"] → only NH and MD present
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let attrs: HashSet<String> = ["NH", "MD"].iter().map(|s| s.to_string()).collect();
        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            &attrs,
        )
        .unwrap();

        let data = record.data();
        // Only NH and MD present
        assert!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT).is_some(),
            "NH should be present"
        );
        assert!(
            data.get(&Tag::new(b'M', b'D')).is_some(),
            "MD should be present"
        );
        // Others absent
        assert!(data.get(&Tag::HIT_INDEX).is_none(), "HI should be absent");
        assert!(
            data.get(&Tag::ALIGNMENT_SCORE).is_none(),
            "AS should be absent"
        );
        assert!(
            data.get(&Tag::EDIT_DISTANCE).is_none(),
            "NM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'n', b'M')).is_none(),
            "nM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'X', b'S')).is_none(),
            "XS should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'M')).is_none(),
            "jM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'I')).is_none(),
            "jI should be absent"
        );
    }

    #[test]
    fn test_sam_attribute_set_expansion() {
        use crate::params::Parameters;

        // Standard
        let p = Parameters::parse_from(vec!["ruSTAR", "--readFilesIn", "r.fq"]);
        let attrs = p.sam_attribute_set();
        assert_eq!(attrs.len(), 5);
        assert!(attrs.contains("NH"));
        assert!(attrs.contains("HI"));
        assert!(attrs.contains("AS"));
        assert!(attrs.contains("NM"));
        assert!(attrs.contains("nM"));

        // All
        let p = Parameters::parse_from(vec![
            "ruSTAR",
            "--readFilesIn",
            "r.fq",
            "--outSAMattributes",
            "All",
        ]);
        let attrs = p.sam_attribute_set();
        assert_eq!(attrs.len(), 9);
        assert!(attrs.contains("nM"));
        assert!(attrs.contains("XS"));
        assert!(attrs.contains("MD"));
        assert!(attrs.contains("jM"));
        assert!(attrs.contains("jI"));

        // None
        let p = Parameters::parse_from(vec![
            "ruSTAR",
            "--readFilesIn",
            "r.fq",
            "--outSAMattributes",
            "None",
        ]);
        let attrs = p.sam_attribute_set();
        assert!(attrs.is_empty());

        // Explicit subset
        let p = Parameters::parse_from(vec![
            "ruSTAR",
            "--readFilesIn",
            "r.fq",
            "--outSAMattributes",
            "NH",
            "AS",
            "MD",
        ]);
        let attrs = p.sam_attribute_set();
        assert_eq!(attrs.len(), 3);
        assert!(attrs.contains("NH"));
        assert!(attrs.contains("AS"));
        assert!(attrs.contains("MD"));
        assert!(!attrs.contains("HI"));
    }

    #[test]
    fn test_nm_vs_nm_mismatch_difference() {
        // Verify NM (edit distance) ≠ nM (mismatches only) when indels present
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 60,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(20),
                CigarOp::Ins(5),
                CigarOp::Match(10),
                CigarOp::Del(3),
                CigarOp::Match(20),
            ],
            score: 80,
            n_mismatch: 2,
            n_gap: 2,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0, 1, 2, 3],
        };

        let read_seq = vec![0, 1, 2, 3];
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
            &all_attrs(),
        )
        .unwrap();

        let data = record.data();
        // NM = edit distance: 2 mismatches + 5 ins + 3 del = 10
        assert_eq!(
            data.get(&Tag::EDIT_DISTANCE),
            Some(&Value::from(10_i32)),
            "NM should be 10 (2 mismatches + 5 ins + 3 del)"
        );
        // nM = mismatches only: 2
        assert_eq!(
            data.get(&Tag::new(b'n', b'M')),
            Some(&Value::from(2_i32)),
            "nM should be 2 (mismatches only, excluding indels)"
        );
    }
}
