# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ruSTAR is a Rust reimplementation of [STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference), an RNA-seq aligner originally written in C++ by Alexander Dobin. Licensed under MIT to match the original STAR license.

The primary goal is a faithful port — matching the original STAR behavior as closely as possible. Extra features and divergences from the original will come in later releases/forks. When implementing, refer to the [STAR source code](https://github.com/alexdobin/STAR) to ensure correctness and behavioral parity.

## Build Commands

Rust 2024 edition. Standard Cargo commands:

```bash
cargo build            # Debug build
cargo build --release  # Release build
cargo test             # Run all tests
cargo test <name>      # Run a single test by name
cargo clippy           # Lint
cargo fmt              # Format code
```

Always run `cargo clippy`, `cargo fmt --check`, and `cargo test` before considering a phase complete.

## Current Implementation Status

See [ROADMAP.md](ROADMAP.md) for detailed phase tracking. **Currently on Phase 7** (GTF/Splice Junction Annotation).

**Completed**:
- Phase 1 (CLI/params)
- Phase 2 (FASTA loading + genome generation)
- Phase 3 (Suffix array + SAindex)
- Phase 4 (Index loading + seed finding)
- Phase 5 (Seed stitching + alignment scoring)
- Phase 6 (FASTQ reading + SAM output) ← **FIRST END-TO-END ALIGNMENT**

## Source Layout

```
src/
  main.rs          -- Thin entry: parse CLI (clap), init logging, call lib::run()
  lib.rs           -- run() dispatches on RunMode (AlignReads | GenomeGenerate)
  params.rs        -- ~40 STAR CLI params via clap derive, --camelCase long names
  error.rs         -- Error enum with thiserror (Parameter, Io, Fasta, Index, Alignment, Gtf)
  mapq.rs          -- ✅ MAPQ calculation (unique=255, multi=-10*log10(1-1/n))
  stats.rs         -- ✅ Alignment statistics tracking and reporting
  genome/
    mod.rs         -- ✅ Genome struct, padding logic, reverse complement, file writing
    fasta.rs       -- ✅ FASTA parser, base encoding (A=0,C=1,G=2,T=3,N=4)
  index/
    mod.rs         -- ✅ GenomeIndex (build + load + write)
    packed_array.rs-- ✅ Variable-width bit packing (1-64 bits per element)
    suffix_array.rs-- ✅ SA construction, custom comparator, strand encoding
    sa_index.rs    -- ✅ K-mer lookup table (35-bit entries with flags)
    io.rs          -- ✅ Load index from disk (Genome, SA, SAindex)
  align/
    mod.rs         -- ✅ Module definition
    seed.rs        -- ✅ Seed finding via MMP search + binary search on SA
    stitch.rs      -- ✅ Seed clustering + DP stitching
    score.rs       -- ✅ Scoring functions (gaps, mismatches, splice junctions)
    transcript.rs  -- ✅ Transcript struct (exon coords, CIGAR, scores)
    read_align.rs  -- ✅ Per-read alignment driver
  io/
    mod.rs         -- ✅ Module exports
    fastq.rs       -- ✅ FASTQ reader (plain + gzip, noodles wrapper)
    sam.rs         -- ✅ SAM writer (header + records, noodles wrapper)
  junction/mod.rs  -- (stub) Phase 7: Splice junctions, GTF parsing, motif detection
```

## Key Conventions

- **Crate name is `ruSTAR`** — `#![allow(non_snake_case)]` in lib.rs suppresses the crate name warning
- **STAR params use `--camelCase` naming** — clap `#[arg(long = "camelCase")]` maps to snake_case Rust fields
- **Multi-value params** (genomeFastaFiles, readFilesIn, outSAMtype, outSAMattributes, chimOutType, alignSJstitchMismatchNmax) need explicit `num_args`
- **Negative defaults** (scoreGapNoncan=-8, readMapNumber=-1, etc.) need `allow_hyphen_values = true`
- **`outSAMtype`** is parsed as raw `Vec<String>` then structured via `Parameters::out_sam_type()` method
- **Validation** beyond clap's type checking is in `Parameters::validate()` (e.g. genomeGenerate requires FASTA files)
- **No async** — CPU-bound work; async adds complexity with zero benefit
- **Error handling** — `thiserror` for `Error` enum, `anyhow` for top-level result propagation

## Dependencies

```toml
[dependencies]
clap = { version = "4", features = ["derive"] }
anyhow = "1"
thiserror = "2"
log = "0.4"
env_logger = "0.11"
memmap2 = "0.9"
byteorder = "1"
noodles = { version = "0.80", features = ["fastq", "sam"] }
flate2 = "1"

[dev-dependencies]
tempfile = "3"
assert_cmd = "2"
predicates = "3"
```

Future phases will add: `rayon` (threading), `noodles-bam` (BAM output).

## Testing Pattern

- Unit tests: `#[cfg(test)]` in each module
- Integration tests: `tests/` directory (future phases)
- Every phase uses differential testing against STAR where applicable
- Test data tiers: synthetic micro-genome → chr22 → full human genome

**Current test status**: 84/84 tests passing, zero critical clippy warnings

## Current Capabilities

ruSTAR can now perform **end-to-end single-end RNA-seq alignment**:
- Generate genome indices from FASTA files
- Read plain or gzipped FASTQ files
- Align single-end RNA-seq reads with splice junction detection
- Output valid SAM files with:
  - Proper headers (@HD, @SQ, @PG)
  - Correct CIGAR strings (M, I, D, N for junctions)
  - Appropriate FLAGS (unmapped, reverse complement)
  - MAPQ scores
  - 1-based genomic positions
- Print alignment statistics (unique/multi/unmapped percentages)

## Limitations (to be addressed in future phases)

- Single-end reads only (paired-end in Phase 8)
- SAM output only (BAM in Phase 10)
- No multithreading (Phase 9)
- No SAM optional tags (AS, NM, NH, HI) - noodles lifetime complexity
- No GTF-based junction scoring (Phase 7)
