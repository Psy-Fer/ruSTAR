pub mod fasta;

use std::path::Path;

use crate::error::Error;
use crate::params::Parameters;

use fasta::parse_fasta_files;

/// STAR's genome spacing character (used for inter-chromosome padding).
const GENOME_SPACING_CHAR: u8 = 5;

/// Packed genome with chromosome metadata.
///
/// The genome sequence is stored as one byte per base:
/// - A=0, C=1, G=2, T=3, N=4, padding=5
/// - Chromosomes are concatenated with padding to bin boundaries
/// - The reverse complement occupies the second half of the `sequence` buffer
pub struct Genome {
    /// Forward genome (0..n_genome) + reverse complement (n_genome..2*n_genome).
    /// Initialized to GENOME_SPACING_CHAR (5), then overwritten with actual bases.
    pub sequence: Vec<u8>,

    /// Total length of the forward (padded) genome.
    pub n_genome: u64,

    /// Number of real chromosomes (not including scaffold/contigs if excluded).
    pub n_chr_real: usize,

    /// Chromosome names.
    pub chr_name: Vec<String>,

    /// True (unpadded) chromosome lengths.
    pub chr_length: Vec<u64>,

    /// Padded start positions of each chromosome in the genome.
    /// Length = n_chr_real + 1; the last entry is n_genome (total size).
    pub chr_start: Vec<u64>,
}

impl Genome {
    /// Build a genome from FASTA files, matching STAR's layout.
    ///
    /// # Arguments
    /// - `params`: CLI parameters (genomeFastaFiles, genomeChrBinNbits)
    ///
    /// # Returns
    /// A `Genome` with forward + reverse complement sequences and metadata.
    pub fn from_fasta(params: &Parameters) -> Result<Self, Error> {
        let chromosomes = parse_fasta_files(&params.genome_fasta_files)?;

        // Compute padding bin size
        let bin_nbits = params.genome_chr_bin_nbits;
        let bin_size = 1u64 << bin_nbits;

        // First pass: compute padded positions and total genome size
        let mut chr_name = Vec::new();
        let mut chr_length = Vec::new();
        let mut chr_start = Vec::new();

        let mut n: u64 = 0; // current position in the padded genome

        for chrom in &chromosomes {
            let len = chrom.sequence.len() as u64;

            if len == 0 {
                return Err(Error::Fasta(format!(
                    "chromosome '{}' has zero length",
                    chrom.name
                )));
            }

            // Apply STAR's padding formula before this chromosome (except for the first)
            if n > 0 {
                n = ((n + 1) / bin_size + 1) * bin_size;
            }

            chr_name.push(chrom.name.clone());
            chr_length.push(len);
            chr_start.push(n);

            n += len;
        }

        // Final padding after the last chromosome
        n = ((n + 1) / bin_size + 1) * bin_size;
        let n_genome = n;
        chr_start.push(n_genome); // STAR adds this extra entry

        let n_chr_real = chromosomes.len();

        // Allocate buffer: forward (0..n_genome) + reverse (n_genome..2*n_genome)
        let total_len = (n_genome * 2) as usize;
        let mut sequence = vec![GENOME_SPACING_CHAR; total_len];

        // Second pass: copy actual chromosome sequences into the buffer
        for (i, chrom) in chromosomes.iter().enumerate() {
            let start = chr_start[i] as usize;
            let len = chrom.sequence.len();
            sequence[start..start + len].copy_from_slice(&chrom.sequence);
        }

        // Build reverse complement in the second half
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

    /// Access a base from the genome (forward or reverse strand).
    ///
    /// # Arguments
    /// - `pos`: Position in the genome (0..2*n_genome)
    ///
    /// # Returns
    /// The base value (0-3 for ACGT, 4 for N, 5 for padding), or None if out of bounds.
    pub fn get_base(&self, pos: u64) -> Option<u8> {
        if pos < self.sequence.len() as u64 {
            Some(self.sequence[pos as usize])
        } else {
            None
        }
    }

    /// Get the chromosome containing a given genomic position.
    ///
    /// # Returns
    /// `(chr_index, offset_within_chr)` or None if position is in padding.
    pub fn position_to_chr(&self, pos: u64) -> Option<(usize, u64)> {
        // Binary search to find the chromosome
        for i in 0..self.n_chr_real {
            let start = self.chr_start[i];
            let end = start + self.chr_length[i];
            if pos >= start && pos < end {
                return Some((i, pos - start));
            }
        }
        None
    }

    /// Write genome index files to the specified directory.
    ///
    /// Creates:
    /// - `Genome` — raw binary file (n_genome bytes, forward strand only)
    /// - `chrName.txt` — chromosome names, one per line
    /// - `chrLength.txt` — chromosome lengths, one per line
    /// - `chrStart.txt` — chromosome start positions + final n_genome entry
    /// - `chrNameLength.txt` — tab-separated name + length
    /// - `genomeParameters.txt` — key-value pairs of genome generation parameters
    pub fn write_index_files(&self, dir: &Path, params: &Parameters) -> Result<(), Error> {
        use std::fs;
        use std::io::Write;

        // Create output directory if needed
        fs::create_dir_all(dir).map_err(|e| Error::io(e, dir))?;

        // Write Genome file (forward strand only, n_genome bytes)
        let genome_path = dir.join("Genome");
        fs::write(&genome_path, &self.sequence[..self.n_genome as usize])
            .map_err(|e| Error::io(e, &genome_path))?;

        // Write chrName.txt
        let chr_name_path = dir.join("chrName.txt");
        let mut f = fs::File::create(&chr_name_path).map_err(|e| Error::io(e, &chr_name_path))?;
        for name in &self.chr_name {
            writeln!(f, "{}", name).map_err(|e| Error::io(e, &chr_name_path))?;
        }

        // Write chrLength.txt
        let chr_length_path = dir.join("chrLength.txt");
        let mut f =
            fs::File::create(&chr_length_path).map_err(|e| Error::io(e, &chr_length_path))?;
        for &len in &self.chr_length {
            writeln!(f, "{}", len).map_err(|e| Error::io(e, &chr_length_path))?;
        }

        // Write chrStart.txt (includes the extra n_genome entry)
        let chr_start_path = dir.join("chrStart.txt");
        let mut f = fs::File::create(&chr_start_path).map_err(|e| Error::io(e, &chr_start_path))?;
        for &start in &self.chr_start {
            writeln!(f, "{}", start).map_err(|e| Error::io(e, &chr_start_path))?;
        }

        // Write chrNameLength.txt (tab-separated)
        let chr_name_length_path = dir.join("chrNameLength.txt");
        let mut f = fs::File::create(&chr_name_length_path)
            .map_err(|e| Error::io(e, &chr_name_length_path))?;
        for (name, &len) in self.chr_name.iter().zip(&self.chr_length) {
            writeln!(f, "{}\t{}", name, len).map_err(|e| Error::io(e, &chr_name_length_path))?;
        }

        // Write genomeParameters.txt
        let genome_params_path = dir.join("genomeParameters.txt");
        let mut f =
            fs::File::create(&genome_params_path).map_err(|e| Error::io(e, &genome_params_path))?;

        writeln!(f, "### STAR genome generation parameters").unwrap();
        writeln!(f, "versionGenome\t2.7.11b").unwrap(); // ruSTAR version (we can bump this)

        // genomeFastaFiles (space-separated)
        write!(f, "genomeFastaFiles").unwrap();
        for path in &params.genome_fasta_files {
            write!(f, "\t{}", path.display()).unwrap();
        }
        writeln!(f).unwrap();

        writeln!(f, "genomeSAindexNbases\t{}", params.genome_sa_index_nbases).unwrap();
        writeln!(f, "genomeChrBinNbits\t{}", params.genome_chr_bin_nbits).unwrap();
        writeln!(f, "genomeSAsparseD\t{}", params.genome_sa_sparse_d).unwrap();

        // sjdb parameters
        if let Some(ref gtf) = params.sjdb_gtf_file {
            writeln!(f, "sjdbGTFfile\t{}", gtf.display()).unwrap();
        } else {
            writeln!(f, "sjdbGTFfile\t-").unwrap();
        }
        writeln!(f, "sjdbOverhang\t{}", params.sjdb_overhang).unwrap();

        // genomeFileSizes: <genome_size> <SA_size> (SA is 0 for now since we haven't built it)
        writeln!(f, "genomeFileSizes\t{}\t0", self.n_genome).unwrap();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_params(fasta_paths: Vec<std::path::PathBuf>, bin_nbits: u32) -> Parameters {
        use clap::Parser;
        let mut args = vec!["ruSTAR", "--runMode", "genomeGenerate"];

        for path in &fasta_paths {
            args.push("--genomeFastaFiles");
            args.push(path.to_str().unwrap());
        }

        let bin_str = bin_nbits.to_string();
        args.extend(["--genomeChrBinNbits", &bin_str]);

        Parameters::parse_from(args)
    }

    #[test]
    fn single_chromosome_padding() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap(); // 4 bases

        let params = make_params(vec![file.path().to_path_buf()], 3); // bin_size = 8
        let genome = Genome::from_fasta(&params).unwrap();

        // Padding formula: n=4, then ((4+1)/8 + 1)*8 = (0+1)*8 = 8
        assert_eq!(genome.n_genome, 8);
        assert_eq!(genome.n_chr_real, 1);
        assert_eq!(genome.chr_start, vec![0, 8]);
        assert_eq!(genome.chr_length, vec![4]);

        // Check bases
        assert_eq!(genome.get_base(0), Some(0)); // A
        assert_eq!(genome.get_base(1), Some(1)); // C
        assert_eq!(genome.get_base(2), Some(2)); // G
        assert_eq!(genome.get_base(3), Some(3)); // T
        assert_eq!(genome.get_base(4), Some(5)); // padding
        assert_eq!(genome.get_base(7), Some(5)); // padding

        // Check reverse complement
        // Formula: pos(2*n - 1 - i) = complement(pos i)
        // For n=8: pos 15 = complement(pos 0), pos 12 = complement(pos 3), etc.
        assert_eq!(genome.get_base(15), Some(3)); // T (complement of A at pos 0)
        assert_eq!(genome.get_base(14), Some(2)); // G (complement of C at pos 1)
        assert_eq!(genome.get_base(13), Some(1)); // C (complement of G at pos 2)
        assert_eq!(genome.get_base(12), Some(0)); // A (complement of T at pos 3)
        assert_eq!(genome.get_base(11), Some(5)); // padding (complement of padding at pos 4)
        assert_eq!(genome.get_base(8), Some(5)); // padding (complement of padding at pos 7)
    }

    #[test]
    fn two_chromosomes_padding() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "AA").unwrap(); // 2 bases
        writeln!(file, ">chr2").unwrap();
        writeln!(file, "TT").unwrap(); // 2 bases

        let params = make_params(vec![file.path().to_path_buf()], 2); // bin_size = 4
        let genome = Genome::from_fasta(&params).unwrap();

        // chr1 starts at 0, length 2
        // After chr1: n=2, padding ((2+1)/4 + 1)*4 = (0+1)*4 = 4
        // chr2 starts at 4, length 2
        // After chr2: n=6, padding ((6+1)/4 + 1)*4 = (1+1)*4 = 8
        assert_eq!(genome.n_genome, 8);
        assert_eq!(genome.chr_start, vec![0, 4, 8]);
        assert_eq!(genome.chr_length, vec![2, 2]);

        // chr1 bases
        assert_eq!(genome.get_base(0), Some(0)); // A
        assert_eq!(genome.get_base(1), Some(0)); // A
        assert_eq!(genome.get_base(2), Some(5)); // padding
        assert_eq!(genome.get_base(3), Some(5)); // padding

        // chr2 bases
        assert_eq!(genome.get_base(4), Some(3)); // T
        assert_eq!(genome.get_base(5), Some(3)); // T
        assert_eq!(genome.get_base(6), Some(5)); // padding
        assert_eq!(genome.get_base(7), Some(5)); // padding
    }

    #[test]
    fn reverse_complement_correctness() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">test").unwrap();
        writeln!(file, "ACGTN").unwrap(); // A=0, C=1, G=2, T=3, N=4

        let params = make_params(vec![file.path().to_path_buf()], 3); // bin_size = 8
        let genome = Genome::from_fasta(&params).unwrap();

        let n = genome.n_genome as usize;

        // Forward: A C G T N (then padding)
        assert_eq!(genome.sequence[0], 0); // A
        assert_eq!(genome.sequence[1], 1); // C
        assert_eq!(genome.sequence[2], 2); // G
        assert_eq!(genome.sequence[3], 3); // T
        assert_eq!(genome.sequence[4], 4); // N

        // Reverse complement should be at positions [2n-1, 2n-2, 2n-3, 2n-4, 2n-5]
        // which maps to the reverse of [0,1,2,3,4]
        assert_eq!(genome.sequence[2 * n - 1 - 0], 3); // T (complement of A at pos 0)
        assert_eq!(genome.sequence[2 * n - 1 - 1], 2); // G (complement of C at pos 1)
        assert_eq!(genome.sequence[2 * n - 1 - 2], 1); // C (complement of G at pos 2)
        assert_eq!(genome.sequence[2 * n - 1 - 3], 0); // A (complement of T at pos 3)
        assert_eq!(genome.sequence[2 * n - 1 - 4], 4); // N (complement of N at pos 4)
    }

    #[test]
    fn position_to_chr_mapping() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "AAA").unwrap();
        writeln!(file, ">chr2").unwrap();
        writeln!(file, "TTT").unwrap();

        let params = make_params(vec![file.path().to_path_buf()], 2); // bin_size = 4
        let genome = Genome::from_fasta(&params).unwrap();

        // chr1: positions 0-2 (starts at 0, length 3)
        // After chr1: n=3, padding ((3+1)/4 + 1)*4 = (1+1)*4 = 8
        // chr2: positions 8-10 (starts at 8, length 3)
        // After chr2: n=11, padding ((11+1)/4 + 1)*4 = (3+1)*4 = 16

        assert_eq!(genome.n_genome, 16);
        assert_eq!(genome.chr_start, vec![0, 8, 16]);

        assert_eq!(genome.position_to_chr(0), Some((0, 0)));
        assert_eq!(genome.position_to_chr(1), Some((0, 1)));
        assert_eq!(genome.position_to_chr(2), Some((0, 2)));
        assert_eq!(genome.position_to_chr(3), None); // padding
        assert_eq!(genome.position_to_chr(7), None); // padding
        assert_eq!(genome.position_to_chr(8), Some((1, 0)));
        assert_eq!(genome.position_to_chr(9), Some((1, 1)));
        assert_eq!(genome.position_to_chr(10), Some((1, 2)));
        assert_eq!(genome.position_to_chr(11), None); // padding
    }
}
