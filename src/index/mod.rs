pub mod packed_array;
pub mod suffix_array;

use std::fs;
use std::path::Path;

use crate::error::Error;
use crate::genome::Genome;
use crate::params::Parameters;
use suffix_array::SuffixArray;

// Future submodules for Phase 3:
// pub mod sa_index;

/// Complete genome index (genome + suffix array + metadata).
pub struct GenomeIndex {
    pub genome: Genome,
    pub suffix_array: SuffixArray,
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

        log::info!(
            "Suffix array built: {} entries",
            suffix_array.len()
        );

        Ok(GenomeIndex {
            genome,
            suffix_array,
        })
    }

    /// Write index files to directory.
    pub fn write(&self, dir: &Path, params: &Parameters) -> Result<(), Error> {
        // Write genome files
        self.genome.write_index_files(dir, params)?;

        // Write SA file
        let sa_path = dir.join("SA");
        fs::write(&sa_path, self.suffix_array.data.data())
            .map_err(|e| Error::io(e, &sa_path))?;

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
