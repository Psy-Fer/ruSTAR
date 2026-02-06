use std::fs::File;
use std::io::Read;
use std::path::Path;

use byteorder::{LittleEndian, ReadBytesExt};

use crate::error::Error;
use crate::genome::Genome;
use crate::index::GenomeIndex;
use crate::index::packed_array::PackedArray;
use crate::index::sa_index::SaIndex;
use crate::index::suffix_array::SuffixArray;
use crate::params::Parameters;

impl GenomeIndex {
    /// Load a genome index from disk.
    ///
    /// Reads Genome, SA, and SAindex files from the specified directory.
    pub fn load(genome_dir: &Path, params: &Parameters) -> Result<Self, Error> {
        log::info!("Loading genome from {}...", genome_dir.display());

        // Load Genome file
        let genome = load_genome(genome_dir, params)?;
        log::info!(
            "Loaded genome: {} chromosomes, {} bytes",
            genome.n_chr_real,
            genome.n_genome
        );

        // Load SA file
        let suffix_array = load_suffix_array(genome_dir, &genome)?;
        log::info!("Loaded suffix array: {} entries", suffix_array.len());

        // Load SAindex file
        let sa_index = load_sa_index(genome_dir, suffix_array.gstrand_bit)?;
        log::info!(
            "Loaded SA index: nbases={}, {} indices",
            sa_index.nbases,
            sa_index.data.len()
        );

        Ok(GenomeIndex {
            genome,
            suffix_array,
            sa_index,
        })
    }
}

/// Load genome from disk.
fn load_genome(genome_dir: &Path, _params: &Parameters) -> Result<Genome, Error> {
    // Read chromosome metadata
    let chr_name_path = genome_dir.join("chrName.txt");
    let chr_name_contents =
        std::fs::read_to_string(&chr_name_path).map_err(|e| Error::io(e, &chr_name_path))?;
    let chr_name: Vec<String> = chr_name_contents.lines().map(|s| s.to_string()).collect();

    let chr_length_path = genome_dir.join("chrLength.txt");
    let chr_length_contents =
        std::fs::read_to_string(&chr_length_path).map_err(|e| Error::io(e, &chr_length_path))?;
    let chr_length: Vec<u64> = chr_length_contents
        .lines()
        .map(|s| s.parse().unwrap())
        .collect();

    let chr_start_path = genome_dir.join("chrStart.txt");
    let chr_start_contents =
        std::fs::read_to_string(&chr_start_path).map_err(|e| Error::io(e, &chr_start_path))?;
    let chr_start: Vec<u64> = chr_start_contents
        .lines()
        .map(|s| s.parse().unwrap())
        .collect();

    let n_chr_real = chr_name.len();
    let n_genome = chr_start[n_chr_real]; // Last entry is total size

    // Load Genome sequence file
    let genome_path = genome_dir.join("Genome");
    let genome_data = std::fs::read(&genome_path).map_err(|e| Error::io(e, &genome_path))?;

    if genome_data.len() != n_genome as usize {
        return Err(Error::Index(format!(
            "Genome file size mismatch: expected {} bytes, got {}",
            n_genome,
            genome_data.len()
        )));
    }

    // Build full sequence buffer (forward + reverse complement)
    let mut sequence = vec![5u8; (n_genome * 2) as usize];
    sequence[..n_genome as usize].copy_from_slice(&genome_data);

    // Build reverse complement
    for i in 0..n_genome as usize {
        let base = sequence[i];
        let complement = if base < 4 { 3 - base } else { base };
        sequence[2 * n_genome as usize - 1 - i] = complement;
    }

    Ok(Genome {
        sequence,
        n_genome,
        n_chr_real,
        chr_name,
        chr_length,
        chr_start,
    })
}

/// Load suffix array from disk.
fn load_suffix_array(genome_dir: &Path, genome: &Genome) -> Result<SuffixArray, Error> {
    let sa_path = genome_dir.join("SA");
    let sa_data = std::fs::read(&sa_path).map_err(|e| Error::io(e, &sa_path))?;

    let gstrand_bit = SuffixArray::calculate_gstrand_bit(genome.n_genome);
    let word_length = gstrand_bit + 1;
    let gstrand_mask = (1u64 << gstrand_bit) - 1;

    // Calculate expected length from file size
    // Formula from STAR: lengthByte = (length-1)*wordLength/8 + 8
    // We need to solve for length, accounting for integer division:
    // total_bits = (lengthByte - 8) * 8
    // length = (total_bits / wordLength) + 1
    // BUT we need ceiling division to account for partial entries
    let length_byte = sa_data.len();
    let length = if length_byte < 8 {
        0
    } else {
        let total_bits = (length_byte - 8) * 8;
        let entries = total_bits.div_ceil(word_length as usize);
        entries + 1
    };

    let data = PackedArray::from_bytes(word_length, length, sa_data);

    Ok(SuffixArray {
        data,
        gstrand_bit,
        gstrand_mask,
    })
}

/// Load SA index from disk.
fn load_sa_index(genome_dir: &Path, gstrand_bit: u32) -> Result<SaIndex, Error> {
    let sai_path = genome_dir.join("SAindex");
    let mut file = File::open(&sai_path).map_err(|e| Error::io(e, &sai_path))?;

    // Read nbases (u64)
    let nbases = file
        .read_u64::<LittleEndian>()
        .map_err(|e| Error::io(e, &sai_path))? as u32;

    // Read genomeSAindexStart array (nbases + 1 entries)
    let mut genome_sa_index_start = Vec::with_capacity((nbases + 1) as usize);
    for _ in 0..=nbases {
        let val = file
            .read_u64::<LittleEndian>()
            .map_err(|e| Error::io(e, &sai_path))?;
        genome_sa_index_start.push(val);
    }

    // Read packed data
    let mut packed_data = Vec::new();
    file.read_to_end(&mut packed_data)
        .map_err(|e| Error::io(e, &sai_path))?;

    let word_length = gstrand_bit + 3;
    let num_indices = SaIndex::calculate_num_indices(nbases);

    let data = PackedArray::from_bytes(word_length, num_indices as usize, packed_data);

    Ok(SaIndex {
        nbases,
        genome_sa_index_start,
        data,
        word_length,
        gstrand_bit,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn load_generated_index() {
        // Create a simple genome
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap();

        let dir = tempfile::tempdir().unwrap();

        let args = vec![
            "ruSTAR",
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            file.path().to_str().unwrap(),
            "--genomeDir",
            dir.path().to_str().unwrap(),
            "--genomeChrBinNbits",
            "2",
            "--genomeSAindexNbases",
            "1",
        ];

        let params = Parameters::parse_from(args.clone());

        // Build index
        let index = GenomeIndex::build(&params).unwrap();
        index.write(dir.path(), &params).unwrap();

        // Load index back
        let loaded_index = GenomeIndex::load(dir.path(), &params).unwrap();

        // Verify
        assert_eq!(loaded_index.genome.n_genome, index.genome.n_genome);
        assert_eq!(loaded_index.genome.n_chr_real, index.genome.n_chr_real);
        assert_eq!(loaded_index.suffix_array.len(), index.suffix_array.len());
        assert_eq!(loaded_index.sa_index.nbases, index.sa_index.nbases);
        assert_eq!(loaded_index.sa_index.data.len(), index.sa_index.data.len());

        // Verify first few SA entries match
        for i in 0..loaded_index.suffix_array.len().min(5) {
            assert_eq!(loaded_index.suffix_array.get(i), index.suffix_array.get(i));
        }
    }
}
