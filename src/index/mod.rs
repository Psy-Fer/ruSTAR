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
use sa_index::SaIndex;
use suffix_array::SuffixArray;

/// Complete genome index (genome + suffix array + SA index + junction database).
#[derive(Clone)]
pub struct GenomeIndex {
    pub genome: Genome,
    pub suffix_array: SuffixArray,
    pub sa_index: SaIndex,
    pub junction_db: SpliceJunctionDb,
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

        // Load GTF annotations if provided
        let junction_db = if let Some(ref gtf_path) = params.sjdb_gtf_file {
            SpliceJunctionDb::from_gtf(gtf_path, &genome)?
        } else {
            log::info!("No GTF file provided, all junctions will be novel");
            SpliceJunctionDb::empty()
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
        })
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

        // Update genomeParameters.txt with SA file size
        let genome_params_path = dir.join("genomeParameters.txt");
        let sa_size = self.suffix_array.data.data().len();

        // Read existing file, update genomeFileSizes line
        let content = fs::read_to_string(&genome_params_path)
            .map_err(|e| Error::io(e, &genome_params_path))?;

        let updated_content = content.replace(
            &format!("genomeFileSizes\t{}\t0", self.genome.n_genome),
            &format!("genomeFileSizes\t{}\t{}", self.genome.n_genome, sa_size),
        );

        fs::write(&genome_params_path, updated_content)
            .map_err(|e| Error::io(e, &genome_params_path))?;

        Ok(())
    }
}
