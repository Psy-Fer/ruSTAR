/// BAM output writer with noodles (streaming, unsorted)
use crate::error::Error;
use crate::genome::Genome;
use crate::params::Parameters;
use crate::quant::transcriptome::TranscriptomeIndex;
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as SamWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

/// Buffer for BAM records built by parallel threads
#[derive(Default)]
pub struct BufferedBamRecords {
    pub records: Vec<RecordBuf>,
}

impl BufferedBamRecords {
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

/// BAM file writer (streaming, unsorted)
///
/// This writer streams BAM records directly to disk as they're generated,
/// without buffering or sorting. The output is BGZF-compressed but unsorted.
/// Users can sort the output with `samtools sort` if needed.
pub struct BamWriter {
    writer: bam::io::Writer<noodles::bgzf::Writer<BufWriter<File>>>,
    header: sam::Header,
}

impl BamWriter {
    /// Create a BAM writer at `output_path` with the given prepared header.
    fn with_header(output_path: &Path, header: sam::Header) -> Result<Self, Error> {
        let buf_writer = BufWriter::new(File::create(output_path)?);
        let mut writer = bam::io::Writer::new(buf_writer);
        writer.write_header(&header)?;
        Ok(Self { writer, header })
    }

    /// Create a new BAM writer with header from genome index.
    pub fn create(output_path: &Path, genome: &Genome, params: &Parameters) -> Result<Self, Error> {
        Self::with_header(
            output_path,
            crate::io::sam::build_sam_header(genome, params)?,
        )
    }

    /// Create a BAM writer whose @SQ header lists the transcripts (not
    /// chromosomes) from `tr_idx`.  Used for
    /// `Aligned.toTranscriptome.out.bam`.
    pub fn create_transcriptome(
        output_path: &Path,
        tr_idx: &TranscriptomeIndex,
        params: &Parameters,
    ) -> Result<Self, Error> {
        let refs = tr_idx
            .tr_ids
            .iter()
            .zip(tr_idx.tr_length.iter())
            .map(|(id, len)| (id.as_str(), *len as usize));
        Self::with_header(
            output_path,
            crate::io::sam::build_sam_header_from_refs(refs, params)?,
        )
    }

    /// Access the underlying SAM header (needed by serial-write code paths
    /// that build records outside of `write_batch`).
    pub fn header(&self) -> &sam::Header {
        &self.header
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

    /// Flush and close BAM file
    pub fn finish(mut self) -> Result<(), Error> {
        self.writer.finish(&self.header)?;
        log::info!("BAM file written successfully");
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::{CigarOp, Exon, Transcript};
    use clap::Parser;
    use tempfile::NamedTempFile;

    fn create_test_genome() -> Genome {
        Genome {
            sequence: vec![0, 1, 2, 3, 0, 1, 2, 3], // ACGTACGT
            n_genome: 8,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![8],
            chr_start: vec![0, 8],
        }
    }

    fn create_test_params() -> Parameters {
        Parameters::parse_from(vec!["ruSTAR", "--readFilesIn", "test.fq"])
    }

    #[test]
    fn test_bam_writer_creation() {
        let genome = create_test_genome();
        let params = create_test_params();
        let temp_file = NamedTempFile::new().unwrap();

        let writer = BamWriter::create(temp_file.path(), &genome, &params);
        assert!(writer.is_ok(), "BAM writer creation should succeed");
    }

    #[test]
    fn test_bam_unmapped_write() {
        let genome = create_test_genome();
        let params = create_test_params();
        let temp_file = NamedTempFile::new().unwrap();

        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();

        let read_name = "read1";
        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let read_qual = vec![30, 30, 30, 30];

        // Build record using SAM builder
        let record = crate::io::sam::SamWriter::build_unmapped_record(
            read_name, &read_seq, &read_qual, None,
        )
        .unwrap();

        let result = writer.write_batch(&[record]);
        assert!(result.is_ok(), "Writing unmapped read should succeed");

        // Finish the file
        let result = writer.finish();
        assert!(result.is_ok(), "Finishing BAM file should succeed");
    }

    #[test]
    fn test_bam_alignment_write() {
        let genome = create_test_genome();
        let params = create_test_params();
        let temp_file = NamedTempFile::new().unwrap();

        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();

        // Create a simple transcript
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 104,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 100,
                genome_end: 104,
                read_start: 0,
                read_end: 4,
            }],
            cigar: vec![CigarOp::Match(4)],
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        let read_name = "read1";
        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let read_qual = vec![30, 30, 30, 30];

        // Build records using SAM builder
        let records = crate::io::sam::SamWriter::build_alignment_records(
            read_name,
            &read_seq,
            &read_qual,
            &[transcript],
            &genome,
            &params,
            1,
        )
        .unwrap();

        let result = writer.write_batch(&records);
        assert!(result.is_ok(), "Writing alignment should succeed");

        let result = writer.finish();
        assert!(result.is_ok(), "Finishing BAM file should succeed");
    }

    #[test]
    fn test_bam_transcriptome_writer_creation() {
        use crate::junction::gtf::GtfRecord;
        use crate::quant::transcriptome::TranscriptomeIndex;
        use std::collections::HashMap;

        let genome = create_test_genome();
        // Stretch genome / exon to fit the tiny test chromosome.
        let mut attrs = HashMap::new();
        attrs.insert("gene_id".to_string(), "G1".to_string());
        attrs.insert("transcript_id".to_string(), "T1".to_string());
        let exons = vec![GtfRecord {
            seqname: "chr1".to_string(),
            feature: "exon".to_string(),
            start: 1,
            end: 8,
            strand: '+',
            attributes: attrs,
        }];
        let tr_idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(tr_idx.n_transcripts(), 1);

        let params = create_test_params();
        let temp_file = NamedTempFile::new().unwrap();
        let writer = BamWriter::create_transcriptome(temp_file.path(), &tr_idx, &params);
        assert!(
            writer.is_ok(),
            "transcriptome BAM writer creation should succeed"
        );

        // Header should contain exactly 1 @SQ entry (matching n_transcripts).
        let writer = writer.unwrap();
        assert_eq!(writer.header().reference_sequences().len(), 1);
    }

    #[test]
    fn test_bam_batch_write() {
        let genome = create_test_genome();
        let params = create_test_params();
        let temp_file = NamedTempFile::new().unwrap();

        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();

        // Create a batch of unmapped records
        let records = vec![
            crate::io::sam::SamWriter::build_unmapped_record(
                "read1",
                &[0, 1, 2, 3],
                &[30, 30, 30, 30],
                None,
            )
            .unwrap(),
            crate::io::sam::SamWriter::build_unmapped_record(
                "read2",
                &[0, 1, 2, 3],
                &[30, 30, 30, 30],
                None,
            )
            .unwrap(),
        ];

        let result = writer.write_batch(&records);
        assert!(result.is_ok(), "Writing batch should succeed");

        let result = writer.finish();
        assert!(result.is_ok(), "Finishing BAM file should succeed");
    }
}
