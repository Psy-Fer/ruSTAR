use std::path::PathBuf;

use clap::Parser;

// ---------------------------------------------------------------------------
// Run mode enum
// ---------------------------------------------------------------------------

/// STAR's `--runMode` values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RunMode {
    AlignReads,
    GenomeGenerate,
}

impl std::str::FromStr for RunMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "alignReads" => Ok(Self::AlignReads),
            "genomeGenerate" => Ok(Self::GenomeGenerate),
            _ => Err(format!(
                "unknown runMode '{s}'; expected 'alignReads' or 'genomeGenerate'"
            )),
        }
    }
}

impl std::fmt::Display for RunMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AlignReads => write!(f, "alignReads"),
            Self::GenomeGenerate => write!(f, "genomeGenerate"),
        }
    }
}

// ---------------------------------------------------------------------------
// SAM output type enums
// ---------------------------------------------------------------------------

/// STAR's `--outSAMtype` format component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutSamFormat {
    Sam,
    Bam,
    None,
}

/// STAR's `--outSAMtype` sort order component (only applies to BAM).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutSamSortOrder {
    Unsorted,
    SortedByCoordinate,
}

/// Combined `--outSAMtype` value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutSamType {
    pub format: OutSamFormat,
    pub sort_order: Option<OutSamSortOrder>,
}

impl Default for OutSamType {
    fn default() -> Self {
        Self {
            format: OutSamFormat::Sam,
            sort_order: None,
        }
    }
}

impl std::fmt::Display for OutSamType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match (&self.format, &self.sort_order) {
            (OutSamFormat::Sam, _) => write!(f, "SAM"),
            (OutSamFormat::None, _) => write!(f, "None"),
            (OutSamFormat::Bam, Some(OutSamSortOrder::SortedByCoordinate)) => {
                write!(f, "BAM SortedByCoordinate")
            }
            (OutSamFormat::Bam, _) => write!(f, "BAM Unsorted"),
        }
    }
}

// ---------------------------------------------------------------------------
// SAM unmapped output
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutSamUnmapped {
    None,
    Within,
    WithinKeepPairs,
}

impl Default for OutSamUnmapped {
    fn default() -> Self {
        Self::None
    }
}

impl std::str::FromStr for OutSamUnmapped {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(Self::None),
            "Within" => Ok(Self::Within),
            "Within KeepPairs" => Ok(Self::WithinKeepPairs),
            _ => Err(format!("unknown outSAMunmapped value: '{s}'")),
        }
    }
}

// ---------------------------------------------------------------------------
// Output filter type
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutFilterType {
    Normal,
    BySJout,
}

impl Default for OutFilterType {
    fn default() -> Self {
        Self::Normal
    }
}

impl std::str::FromStr for OutFilterType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Normal" => Ok(Self::Normal),
            "BySJout" => Ok(Self::BySJout),
            _ => Err(format!("unknown outFilterType value: '{s}'")),
        }
    }
}

// ---------------------------------------------------------------------------
// Two-pass mode
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TwopassMode {
    None,
    Basic,
}

impl Default for TwopassMode {
    fn default() -> Self {
        Self::None
    }
}

impl std::str::FromStr for TwopassMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(Self::None),
            "Basic" => Ok(Self::Basic),
            _ => Err(format!("unknown twopassMode value: '{s}'")),
        }
    }
}

// ---------------------------------------------------------------------------
// Parameters struct
// ---------------------------------------------------------------------------

/// ruSTAR command-line parameters, matching STAR's `--camelCase` argument names.
///
/// Only the ~40 most important parameters are included; more will be added
/// incrementally as later phases need them.
#[derive(Debug, Clone, Parser)]
#[command(
    name = "ruSTAR",
    about = "RNA-seq aligner (Rust reimplementation of STAR)",
    version
)]
pub struct Parameters {
    // ── Run ─────────────────────────────────────────────────────────────
    /// Run mode: alignReads or genomeGenerate
    #[arg(long = "runMode", default_value = "alignReads")]
    pub run_mode: RunMode,

    /// Number of threads
    #[arg(long = "runThreadN", default_value_t = 1)]
    pub run_thread_n: usize,

    // ── Genome ──────────────────────────────────────────────────────────
    /// Path to genome index directory
    #[arg(long = "genomeDir", default_value = "./GenomeDir")]
    pub genome_dir: PathBuf,

    /// FASTA file(s) with genome reference sequences (for genomeGenerate)
    #[arg(long = "genomeFastaFiles", num_args = 1..)]
    pub genome_fasta_files: Vec<PathBuf>,

    /// Length of SA pre-indexing string (log2-based)
    #[arg(long = "genomeSAindexNbases", default_value_t = 14)]
    pub genome_sa_index_nbases: u32,

    /// Log2(chromosome bin size) for genome storage
    #[arg(long = "genomeChrBinNbits", default_value_t = 18)]
    pub genome_chr_bin_nbits: u32,

    /// Suffix array sparsity (larger = less RAM, slower mapping)
    #[arg(long = "genomeSAsparseD", default_value_t = 1)]
    pub genome_sa_sparse_d: u32,

    // ── Read files ──────────────────────────────────────────────────────
    /// Input read file(s); second file is mate 2 for paired-end
    #[arg(long = "readFilesIn", num_args = 1..=2)]
    pub read_files_in: Vec<PathBuf>,

    /// Command to decompress input files (e.g. "zcat" for .gz)
    #[arg(long = "readFilesCommand")]
    pub read_files_command: Option<String>,

    /// Number of reads to map; -1 = all
    #[arg(long = "readMapNumber", default_value_t = -1, allow_hyphen_values = true)]
    pub read_map_number: i64,

    /// Bases to clip from 5' end of each mate
    #[arg(long = "clip5pNbases", default_value_t = 0)]
    pub clip5p_nbases: u32,

    /// Bases to clip from 3' end of each mate
    #[arg(long = "clip3pNbases", default_value_t = 0)]
    pub clip3p_nbases: u32,

    // ── Output ──────────────────────────────────────────────────────────
    /// Output file name prefix (including path)
    #[arg(long = "outFileNamePrefix", default_value = "./")]
    pub out_file_name_prefix: PathBuf,

    /// Output type: SAM, BAM Unsorted, BAM SortedByCoordinate, None.
    /// Provide as space-separated tokens, e.g. "BAM SortedByCoordinate".
    #[arg(long = "outSAMtype", num_args = 1..=2, default_values_t = vec!["SAM".to_string()])]
    pub out_sam_type_raw: Vec<String>,

    /// Strand field: None or intronMotif
    #[arg(long = "outSAMstrandField", default_value = "None")]
    pub out_sam_strand_field: String,

    /// SAM attributes to include (Standard, All, None, or explicit list)
    #[arg(long = "outSAMattributes", num_args = 1.., default_values_t = vec!["Standard".to_string()])]
    pub out_sam_attributes: Vec<String>,

    /// Unmapped reads in SAM output: None or Within
    #[arg(long = "outSAMunmapped", default_value = "None")]
    pub out_sam_unmapped: OutSamUnmapped,

    /// MAPQ value for unique mappers
    #[arg(long = "outSAMmapqUnique", default_value_t = 255)]
    pub out_sam_mapq_unique: u8,

    /// Output filter type: Normal or BySJout
    #[arg(long = "outFilterType", default_value = "Normal")]
    pub out_filter_type: OutFilterType,

    /// Max multimap loci (reads mapping to more loci are unmapped)
    #[arg(long = "outFilterMultimapNmax", default_value_t = 10)]
    pub out_filter_multimap_nmax: u32,

    /// Score range for multi-mapping (keep alignments within this range of best score)
    #[arg(long = "outFilterMultimapScoreRange", default_value_t = 1)]
    pub out_filter_multimap_score_range: i32,

    /// Max mismatches per pair
    #[arg(long = "outFilterMismatchNmax", default_value_t = 10)]
    pub out_filter_mismatch_nmax: u32,

    /// Max ratio of mismatches to mapped length
    #[arg(long = "outFilterMismatchNoverLmax", default_value_t = 0.3)]
    pub out_filter_mismatch_nover_lmax: f64,

    /// Min alignment score (absolute)
    #[arg(long = "outFilterScoreMin", default_value_t = 0)]
    pub out_filter_score_min: i32,

    /// Min alignment score normalized to read length
    #[arg(long = "outFilterScoreMinOverLread", default_value_t = 0.66)]
    pub out_filter_score_min_over_lread: f64,

    /// Min matched bases (absolute)
    #[arg(long = "outFilterMatchNmin", default_value_t = 0)]
    pub out_filter_match_nmin: u32,

    /// Min matched bases normalized to read length
    #[arg(long = "outFilterMatchNminOverLread", default_value_t = 0.66)]
    pub out_filter_match_nmin_over_lread: f64,

    // ── Alignment scoring ───────────────────────────────────────────────
    /// Min intron size (smaller gaps are deletions)
    #[arg(long = "alignIntronMin", default_value_t = 21)]
    pub align_intron_min: u32,

    /// Max intron size; 0 = auto
    #[arg(long = "alignIntronMax", default_value_t = 0)]
    pub align_intron_max: u32,

    /// Max genomic distance between mates; 0 = auto
    #[arg(long = "alignMatesGapMax", default_value_t = 0)]
    pub align_mates_gap_max: u32,

    /// Min overhang for novel spliced alignments
    #[arg(long = "alignSJoverhangMin", default_value_t = 5)]
    pub align_sj_overhang_min: u32,

    /// Min overhang for annotated splice junctions
    #[arg(long = "alignSJDBoverhangMin", default_value_t = 3)]
    pub align_sjdb_overhang_min: u32,

    /// Max mismatches for stitching SJs (4 ints: noncanonical, GC/AG, AT/AC, noncanonical)
    #[arg(long = "alignSJstitchMismatchNmax", num_args = 4,
          default_values_t = vec![0, -1, 0, 0], allow_hyphen_values = true)]
    pub align_sj_stitch_mismatch_nmax: Vec<i32>,

    /// Splice junction penalty (canonical)
    #[arg(long = "scoreGap", default_value_t = 0)]
    pub score_gap: i32,

    /// Non-canonical junction penalty
    #[arg(long = "scoreGapNoncan", default_value_t = -8, allow_hyphen_values = true)]
    pub score_gap_noncan: i32,

    /// GC/AG junction penalty
    #[arg(long = "scoreGapGCAG", default_value_t = -4, allow_hyphen_values = true)]
    pub score_gap_gcag: i32,

    /// AT/AC junction penalty
    #[arg(long = "scoreGapATAC", default_value_t = -8, allow_hyphen_values = true)]
    pub score_gap_atac: i32,

    /// Deletion open penalty
    #[arg(long = "scoreDelOpen", default_value_t = -2, allow_hyphen_values = true)]
    pub score_del_open: i32,

    /// Deletion extension penalty per base
    #[arg(long = "scoreDelBase", default_value_t = -2, allow_hyphen_values = true)]
    pub score_del_base: i32,

    /// Insertion open penalty
    #[arg(long = "scoreInsOpen", default_value_t = -2, allow_hyphen_values = true)]
    pub score_ins_open: i32,

    /// Insertion extension penalty per base
    #[arg(long = "scoreInsBase", default_value_t = -2, allow_hyphen_values = true)]
    pub score_ins_base: i32,

    /// Max score reduction for SJ stitching shift
    #[arg(long = "scoreStitchSJshift", default_value_t = 1)]
    pub score_stitch_sj_shift: i32,

    // ── Seed and anchor parameters ──────────────────────────────────────
    /// Max number of loci a seed can map to (seeds with more loci are discarded)
    #[arg(long = "seedMultimapNmax", default_value_t = 10000)]
    pub seed_multimap_nmax: usize,

    /// Max number of loci anchors are allowed to map to
    #[arg(long = "winAnchorMultimapNmax", default_value_t = 50)]
    pub win_anchor_multimap_nmax: usize,

    /// Max number of seed loci per window
    #[arg(long = "seedNoneLociPerWindow", default_value_t = 10)]
    pub seed_none_loci_per_window: usize,

    // ── Splice junction database ────────────────────────────────────────
    /// GTF file for splice junction annotations
    #[arg(long = "sjdbGTFfile")]
    pub sjdb_gtf_file: Option<PathBuf>,

    /// Overhang length for splice junction database
    #[arg(long = "sjdbOverhang", default_value_t = 100)]
    pub sjdb_overhang: u32,

    /// Extra score for alignments crossing annotated junctions
    #[arg(long = "sjdbScore", default_value_t = 2)]
    pub sjdb_score: i32,

    // ── Two-pass ────────────────────────────────────────────────────────
    /// Two-pass mode: None or Basic
    #[arg(long = "twopassMode", default_value = "None")]
    pub twopass_mode: TwopassMode,

    /// Reads to process in first pass; -1 = all
    #[arg(long = "twopass1readsN", default_value_t = -1, allow_hyphen_values = true)]
    pub twopass1_reads_n: i64,

    // ── Chimeric ────────────────────────────────────────────────────────
    /// Min chimeric segment length; 0 = disable chimeric detection
    #[arg(long = "chimSegmentMin", default_value_t = 0)]
    pub chim_segment_min: u32,

    /// Min total chimeric score
    #[arg(long = "chimScoreMin", default_value_t = 0)]
    pub chim_score_min: i32,

    /// Chimeric output type
    #[arg(long = "chimOutType", num_args = 1..=2, default_values_t = vec!["Junctions".to_string()])]
    pub chim_out_type: Vec<String>,
}

impl Parameters {
    /// Parse the raw `--outSAMtype` tokens into a structured `OutSamType`.
    pub fn out_sam_type(&self) -> Result<OutSamType, String> {
        match self
            .out_sam_type_raw
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>()
            .as_slice()
        {
            ["SAM"] => Ok(OutSamType {
                format: OutSamFormat::Sam,
                sort_order: None,
            }),
            ["None"] => Ok(OutSamType {
                format: OutSamFormat::None,
                sort_order: None,
            }),
            ["BAM", "Unsorted"] => Ok(OutSamType {
                format: OutSamFormat::Bam,
                sort_order: Some(OutSamSortOrder::Unsorted),
            }),
            ["BAM", "SortedByCoordinate"] => Ok(OutSamType {
                format: OutSamFormat::Bam,
                sort_order: Some(OutSamSortOrder::SortedByCoordinate),
            }),
            other => Err(format!("unknown outSAMtype: {:?}", other)),
        }
    }

    /// Validate parameter combinations that clap alone cannot enforce.
    pub fn validate(&self) -> Result<(), crate::error::Error> {
        // genomeGenerate requires FASTA files
        if self.run_mode == RunMode::GenomeGenerate && self.genome_fasta_files.is_empty() {
            return Err(crate::error::Error::Parameter(
                "--genomeFastaFiles is required when --runMode genomeGenerate".into(),
            ));
        }

        // alignReads requires read files
        if self.run_mode == RunMode::AlignReads && self.read_files_in.is_empty() {
            return Err(crate::error::Error::Parameter(
                "--readFilesIn is required when --runMode alignReads".into(),
            ));
        }

        // Validate outSAMtype
        self.out_sam_type()
            .map_err(crate::error::Error::Parameter)?;

        // Thread count must be at least 1
        if self.run_thread_n == 0 {
            return Err(crate::error::Error::Parameter(
                "--runThreadN must be >= 1".into(),
            ));
        }

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: parse a STAR-style command line (without program name).
    fn parse(args: &[&str]) -> Parameters {
        let mut full = vec!["ruSTAR"];
        full.extend_from_slice(args);
        Parameters::parse_from(full)
    }

    #[test]
    fn defaults() {
        let p = parse(&["--readFilesIn", "reads.fq"]);
        assert_eq!(p.run_mode, RunMode::AlignReads);
        assert_eq!(p.run_thread_n, 1);
        assert_eq!(p.genome_dir, PathBuf::from("./GenomeDir"));
        assert_eq!(p.genome_sa_index_nbases, 14);
        assert_eq!(p.genome_chr_bin_nbits, 18);
        assert_eq!(p.genome_sa_sparse_d, 1);
        assert_eq!(p.read_map_number, -1);
        assert_eq!(p.clip5p_nbases, 0);
        assert_eq!(p.clip3p_nbases, 0);
        assert_eq!(p.out_file_name_prefix, PathBuf::from("./"));
        assert_eq!(p.out_sam_type_raw, vec!["SAM".to_string()]);
        assert_eq!(p.out_sam_strand_field, "None");
        assert_eq!(p.out_sam_attributes, vec!["Standard".to_string()]);
        assert_eq!(p.out_sam_unmapped, OutSamUnmapped::None);
        assert_eq!(p.out_sam_mapq_unique, 255);
        assert_eq!(p.out_filter_type, OutFilterType::Normal);
        assert_eq!(p.out_filter_multimap_nmax, 10);
        assert_eq!(p.out_filter_mismatch_nmax, 10);
        assert!((p.out_filter_mismatch_nover_lmax - 0.3).abs() < f64::EPSILON);
        assert_eq!(p.out_filter_score_min, 0);
        assert!((p.out_filter_score_min_over_lread - 0.66).abs() < f64::EPSILON);
        assert_eq!(p.out_filter_match_nmin, 0);
        assert!((p.out_filter_match_nmin_over_lread - 0.66).abs() < f64::EPSILON);
        assert_eq!(p.align_intron_min, 21);
        assert_eq!(p.align_intron_max, 0);
        assert_eq!(p.align_mates_gap_max, 0);
        assert_eq!(p.align_sj_overhang_min, 5);
        assert_eq!(p.align_sjdb_overhang_min, 3);
        assert_eq!(p.align_sj_stitch_mismatch_nmax, vec![0, -1, 0, 0]);
        assert_eq!(p.score_gap, 0);
        assert_eq!(p.score_gap_noncan, -8);
        assert_eq!(p.score_gap_gcag, -4);
        assert_eq!(p.score_gap_atac, -8);
        assert_eq!(p.score_del_open, -2);
        assert_eq!(p.score_del_base, -2);
        assert_eq!(p.score_ins_open, -2);
        assert_eq!(p.score_ins_base, -2);
        assert_eq!(p.score_stitch_sj_shift, 1);
        assert_eq!(p.seed_multimap_nmax, 10000);
        assert_eq!(p.win_anchor_multimap_nmax, 50);
        assert_eq!(p.seed_none_loci_per_window, 10);
        assert!(p.sjdb_gtf_file.is_none());
        assert_eq!(p.sjdb_overhang, 100);
        assert_eq!(p.sjdb_score, 2);
        assert_eq!(p.twopass_mode, TwopassMode::None);
        assert_eq!(p.twopass1_reads_n, -1);
        assert_eq!(p.chim_segment_min, 0);
        assert_eq!(p.chim_score_min, 0);
        assert_eq!(p.chim_out_type, vec!["Junctions".to_string()]);
    }

    #[test]
    fn genome_generate_mode() {
        let p = parse(&[
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            "/data/genome",
            "--genomeFastaFiles",
            "chr1.fa",
            "chr2.fa",
            "--runThreadN",
            "8",
            "--genomeSAindexNbases",
            "11",
        ]);
        assert_eq!(p.run_mode, RunMode::GenomeGenerate);
        assert_eq!(p.genome_dir, PathBuf::from("/data/genome"));
        assert_eq!(
            p.genome_fasta_files,
            vec![PathBuf::from("chr1.fa"), PathBuf::from("chr2.fa")]
        );
        assert_eq!(p.run_thread_n, 8);
        assert_eq!(p.genome_sa_index_nbases, 11);
    }

    #[test]
    fn typical_align_command() {
        let p = parse(&[
            "--runMode",
            "alignReads",
            "--genomeDir",
            "/idx/hg38",
            "--readFilesIn",
            "R1.fq.gz",
            "R2.fq.gz",
            "--readFilesCommand",
            "zcat",
            "--runThreadN",
            "16",
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--outFileNamePrefix",
            "/out/sample1_",
            "--outFilterMultimapNmax",
            "20",
            "--alignIntronMax",
            "1000000",
            "--sjdbGTFfile",
            "gencode.gtf",
            "--twopassMode",
            "Basic",
        ]);
        assert_eq!(p.run_mode, RunMode::AlignReads);
        assert_eq!(p.genome_dir, PathBuf::from("/idx/hg38"));
        assert_eq!(
            p.read_files_in,
            vec![PathBuf::from("R1.fq.gz"), PathBuf::from("R2.fq.gz")]
        );
        assert_eq!(p.read_files_command, Some("zcat".to_string()));
        assert_eq!(p.run_thread_n, 16);
        assert_eq!(
            p.out_sam_type_raw,
            vec!["BAM".to_string(), "SortedByCoordinate".to_string()]
        );
        let sam_type = p.out_sam_type().unwrap();
        assert_eq!(sam_type.format, OutSamFormat::Bam);
        assert_eq!(
            sam_type.sort_order,
            Some(OutSamSortOrder::SortedByCoordinate)
        );
        assert_eq!(p.out_file_name_prefix, PathBuf::from("/out/sample1_"));
        assert_eq!(p.out_filter_multimap_nmax, 20);
        assert_eq!(p.align_intron_max, 1_000_000);
        assert_eq!(p.sjdb_gtf_file, Some(PathBuf::from("gencode.gtf")));
        assert_eq!(p.twopass_mode, TwopassMode::Basic);
    }

    #[test]
    fn scoring_overrides() {
        let p = parse(&[
            "--readFilesIn",
            "reads.fq",
            "--scoreGap",
            "0",
            "--scoreGapNoncan",
            "-12",
            "--scoreGapGCAG",
            "-6",
            "--scoreGapATAC",
            "-10",
            "--scoreDelOpen",
            "-3",
            "--scoreDelBase",
            "-1",
            "--scoreInsOpen",
            "-3",
            "--scoreInsBase",
            "-1",
        ]);
        assert_eq!(p.score_gap, 0);
        assert_eq!(p.score_gap_noncan, -12);
        assert_eq!(p.score_gap_gcag, -6);
        assert_eq!(p.score_gap_atac, -10);
        assert_eq!(p.score_del_open, -3);
        assert_eq!(p.score_del_base, -1);
        assert_eq!(p.score_ins_open, -3);
        assert_eq!(p.score_ins_base, -1);
    }

    #[test]
    fn validate_genome_generate_needs_fasta() {
        let p = parse(&["--runMode", "genomeGenerate"]);
        let err = p.validate().unwrap_err();
        assert!(err.to_string().contains("genomeFastaFiles"));
    }

    #[test]
    fn validate_align_needs_reads() {
        let p = parse(&["--runMode", "alignReads"]);
        let err = p.validate().unwrap_err();
        assert!(err.to_string().contains("readFilesIn"));
    }

    #[test]
    fn out_sam_type_parsing() {
        let p = parse(&["--readFilesIn", "r.fq", "--outSAMtype", "SAM"]);
        let t = p.out_sam_type().unwrap();
        assert_eq!(t.format, OutSamFormat::Sam);
        assert_eq!(t.sort_order, None);

        let p = parse(&["--readFilesIn", "r.fq", "--outSAMtype", "BAM", "Unsorted"]);
        let t = p.out_sam_type().unwrap();
        assert_eq!(t.format, OutSamFormat::Bam);
        assert_eq!(t.sort_order, Some(OutSamSortOrder::Unsorted));

        let p = parse(&["--readFilesIn", "r.fq", "--outSAMtype", "None"]);
        let t = p.out_sam_type().unwrap();
        assert_eq!(t.format, OutSamFormat::None);
    }

    #[test]
    fn chimeric_params() {
        let p = parse(&[
            "--readFilesIn",
            "r.fq",
            "--chimSegmentMin",
            "20",
            "--chimScoreMin",
            "10",
            "--chimOutType",
            "WithinBAM",
            "SoftClip",
        ]);
        assert_eq!(p.chim_segment_min, 20);
        assert_eq!(p.chim_score_min, 10);
        assert_eq!(
            p.chim_out_type,
            vec!["WithinBAM".to_string(), "SoftClip".to_string()]
        );
    }

    #[test]
    fn sj_stitch_mismatch() {
        let p = parse(&[
            "--readFilesIn",
            "r.fq",
            "--alignSJstitchMismatchNmax",
            "1",
            "-1",
            "2",
            "3",
        ]);
        assert_eq!(p.align_sj_stitch_mismatch_nmax, vec![1, -1, 2, 3]);
    }
}
