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

See [ROADMAP.md](ROADMAP.md) for detailed phase tracking. Currently on **Phase 5** (Seed Stitching + Alignment Scoring).

**Completed**:
- Phase 1 (CLI/params)
- Phase 2 (FASTA loading + genome generation)
- Phase 3 (Suffix array + SAindex)
- Phase 4 (Index loading + seed finding)

## Source Layout

```
src/
  main.rs          -- Thin entry: parse CLI (clap), init logging, call lib::run()
  lib.rs           -- run() dispatches on RunMode (AlignReads | GenomeGenerate)
  params.rs        -- ~40 STAR CLI params via clap derive, --camelCase long names
  error.rs         -- Error enum with thiserror (Parameter, Io, Fasta, Index, Alignment, Gtf)
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
    mod.rs         -- Module definition
    seed.rs        -- ✅ Seed finding via MMP search + binary search on SA
  junction/mod.rs  -- (stub) Phase 7: Splice junctions, motif detection
  io/mod.rs        -- (stub) Phase 6+: FASTQ reader, SAM/BAM output, SJ.out.tab
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

[dev-dependencies]
tempfile = "3"
assert_cmd = "2"
predicates = "3"
```

Future phases will add: `rayon`, `flate2`, `memmap2`, `byteorder`, `noodles-bam`/`rust-htslib`.

## Testing Pattern

- Unit tests: `#[cfg(test)]` in each module
- Integration tests: `tests/` directory (from Phase 6)
- Every phase uses differential testing against STAR where applicable
- Test data tiers: synthetic micro-genome → chr22 → full human genome
