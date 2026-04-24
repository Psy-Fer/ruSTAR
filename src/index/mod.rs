pub mod io;
pub mod packed_array;
pub mod sa_index;
pub mod suffix_array;

use std::fs;
use std::path::Path;

use crate::error::Error;
use crate::genome::Genome;
use crate::junction::SpliceJunctionDb;
use crate::params::Parameters;
use crate::quant::transcriptome::TranscriptomeIndex;
use sa_index::SaIndex;
use suffix_array::SuffixArray;

/// Complete genome index (genome + suffix array + SA index + junction database).
#[derive(Clone)]
pub struct GenomeIndex {
    pub genome: Genome,
    pub suffix_array: SuffixArray,
    pub sa_index: SaIndex,
    pub junction_db: SpliceJunctionDb,
    /// Populated when the index was built with a GTF (`--sjdbGTFfile`).
    /// Mirrors STAR's `Transcriptome` object and is written to disk as
    /// `transcriptInfo.tab` + friends at `genomeGenerate`, reloaded at
    /// `alignReads` from the same files.
    pub transcriptome: Option<TranscriptomeIndex>,
}

impl GenomeIndex {
    /// Build a complete genome index from FASTA files.
    pub fn build(params: &Parameters) -> Result<Self, Error> {
        log::info!("Loading FASTA files...");
        let genome = Genome::from_fasta(params)?;

        log::info!(
            "Loaded {} chromosomes, total padded genome size: {} bytes",
            genome.n_chr_real,
            genome.n_genome
        );

        log::info!("Building suffix array...");
        let suffix_array = SuffixArray::build(&genome)?;

        log::info!("Suffix array built: {} entries", suffix_array.len());

        log::info!("Building SA index...");
        let sa_index = SaIndex::build(&genome, &suffix_array, params.genome_sa_index_nbases)?;

        log::info!(
            "SA index built: nbases={}, {} indices",
            sa_index.nbases,
            sa_index.data.len()
        );

        // Load GTF annotations if provided.  This parses the GTF once and
        // shares the result between the junction database and the
        // transcriptome index so we don't pay the cost twice.
        let (junction_db, transcriptome) = if let Some(ref gtf_path) = params.sjdb_gtf_file {
            let jdb = SpliceJunctionDb::from_gtf(gtf_path, &genome)?;
            let exons = crate::junction::gtf::parse_gtf(gtf_path)?;
            let tr = TranscriptomeIndex::from_gtf_exons(&exons, &genome)?;
            log::info!(
                "Transcriptome index built from GTF: {} transcripts, {} genes",
                tr.n_transcripts(),
                tr.gene_ids.len()
            );
            (jdb, Some(tr))
        } else {
            log::info!("No GTF file provided, all junctions will be novel");
            (SpliceJunctionDb::empty(), None)
        };

        log::info!(
            "Junction database initialized: {} annotated junctions",
            junction_db.len()
        );

        Ok(GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db,
            transcriptome,
        })
    }

    /// Convert a raw SA position for a reverse-strand match to forward genome coordinates.
    ///
    /// The SA stores reverse-strand positions as offsets within the RC genome region.
    /// For chromosome identification and SAM output, we need the corresponding position
    /// in the forward genome (leftmost aligned base in forward coordinates).
    ///
    /// For forward strand: returns the position unchanged.
    /// For reverse strand: `forward_pos = n_genome - 1 - sa_pos - (match_length - 1)`
    ///
    /// The raw SA position is still needed for genome base access (add n_genome offset).
    pub fn sa_pos_to_forward(&self, sa_pos: u64, is_reverse: bool, match_length: usize) -> u64 {
        if is_reverse {
            self.genome.n_genome - sa_pos - match_length as u64
        } else {
            sa_pos
        }
    }

    /// Write index files to directory.
    pub fn write(&self, dir: &Path, params: &Parameters) -> Result<(), Error> {
        use std::io::Write;

        // Write genome files
        self.genome.write_index_files(dir, params)?;

        // Write SA file
        let sa_path = dir.join("SA");
        fs::write(&sa_path, self.suffix_array.data.data()).map_err(|e| Error::io(e, &sa_path))?;

        // Write SAindex file
        let sai_path = dir.join("SAindex");
        let mut sai_file = fs::File::create(&sai_path).map_err(|e| Error::io(e, &sai_path))?;

        // Write header: gSAindexNbases as u64
        sai_file
            .write_all(&(self.sa_index.nbases as u64).to_le_bytes())
            .map_err(|e| Error::io(e, &sai_path))?;

        // Write genomeSAindexStart array
        for &val in &self.sa_index.genome_sa_index_start {
            sai_file
                .write_all(&val.to_le_bytes())
                .map_err(|e| Error::io(e, &sai_path))?;
        }

        // Write packed SAindex data
        sai_file
            .write_all(self.sa_index.data.data())
            .map_err(|e| Error::io(e, &sai_path))?;

        // Update genomeParameters.txt with SA file size. Matches STAR's
        // `genomeFileSizes\t<n_genome> <sa_size>\n` pattern (tab before first
        // value, space between subsequent values) — written out in
        // Genome::write_genome_parameters_txt with `0` as the SA placeholder.
        let genome_params_path = dir.join("genomeParameters.txt");
        let sa_size = self.suffix_array.data.data().len();
        let content = fs::read_to_string(&genome_params_path)
            .map_err(|e| Error::io(e, &genome_params_path))?;
        let updated_content = content.replace(
            &format!("genomeFileSizes\t{} 0", self.genome.n_genome),
            &format!("genomeFileSizes\t{} {}", self.genome.n_genome, sa_size),
        );
        fs::write(&genome_params_path, updated_content)
            .map_err(|e| Error::io(e, &genome_params_path))?;

        // Write transcriptome index files (STAR-compatible) when the GTF
        // was supplied. Matches STAR's `GTF_transcriptGeneSJ.cpp` outputs.
        if let Some(tr) = &self.transcriptome {
            tr.write_transcript_info(dir)?;
            tr.write_exon_info(dir)?;
            tr.write_gene_info(dir)?;
            tr.write_exon_ge_tr_info(dir)?;
            tr.write_sjdb_list_from_gtf(dir, &self.genome)?;
            log::info!(
                "Wrote transcriptome index files: {} transcripts, {} genes",
                tr.n_transcripts(),
                tr.gene_ids.len()
            );
        }

        Ok(())
    }
}
