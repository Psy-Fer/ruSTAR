# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important: Git Workflow

**DO NOT commit changes automatically.** The user will review, test, and commit changes themselves. Claude should:
- Make code changes as requested
- Suggest what should be committed
- Let the user handle `git add`, `git commit`, and `git push`

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

See [ROADMAP.md](ROADMAP.md) for detailed phase tracking. **Phases 1-13.11 complete (96.3% position agreement, 97.8% CIGAR agreement, 2 false junctions). Next: Phase 13.12** (further accuracy improvements).

**Phase order change**: Phases reordered to 9 → 8 → 7 to establish parallel architecture foundation
before adding complex features. Threading affects the entire execution model and is harder to retrofit later.

**Completed**:
- Phase 1 (CLI/params)
- Phase 2 (FASTA loading + genome generation)
- Phase 3 (Suffix array + SAindex)
- Phase 4 (Index loading + seed finding)
- Phase 5 (Seed stitching + alignment scoring)
- Phase 6 (FASTQ reading + SAM output) ← **FIRST END-TO-END ALIGNMENT**
- Phase 9 (Threading)
- Phase 8 (Paired-end reads)
- Phase 7 (GTF/splice junction annotation)
- Phase 10 (BAM output - unsorted streaming)
- Phase 11 (Two-pass mode - novel junction discovery)
- Phase 12 (Chimeric alignment detection)
- Phase 13 (Performance + accuracy optimization) ← **94.5% position agreement, soft clips match STAR**
- Phase 13.9 (Position fix) ← **SA reverse-strand encoding fix, SEQ reverse-complement**
- Phase 13.9b (CIGAR/splice fix) ← **CIGAR reversal, motif coord fix, genomic length penalty**
- Phase 13.9c (Deterministic tie-breaking) ← **Reproducible multi-mapper ordering**
- Phase 13.10 (Accuracy parity) ← **Terminal exon filter, annotation bonus, coverage filter, seed caps**
- Phase 13.11 (R→L seeding) ← **Bidirectional seed search, +127 multi-mapped, 3.3x more shared junctions**

**Current Status** (10k yeast reads, single-end):
- ✅ **96.3% position agreement** with STAR (was 51% → 94.5% → 95.3% → 96.3%)
- ✅ **97.8% CIGAR agreement** among position-matching reads (was 84.3% → 96.5% → 97.4% → 97.8%)
- ✅ 84.0% unique mapped (STAR: 82.6%), 4.92% multi-mapped (STAR: 7.4%)
- ✅ 26.5% soft clips (STAR: 26.0%)
- ✅ 80% motif agreement on shared junctions (24/30)
- ✅ 30 shared junctions (was 9, STAR has 72 total)
- ✅ SAM SEQ properly reverse-complemented for reverse-strand reads
- ✅ 196 unit tests passing
- ✅ Deterministic output (identical SAM across runs)
- ✅ Bidirectional seed search (L→R + R→L)
- ✅ **Only 2 false junctions** (was 33 → 3 → 2)
- ✅ Annotation-aware DP scoring (sjdbScore bonus during stitching)
- ⚠️ **Splice rate 0.9%** (STAR: 2.2%) — doubled by R→L seeding but gap remains

## Source Layout

```
src/
  main.rs          -- Thin entry: parse CLI (clap), init logging, call lib::run()
  lib.rs           -- run() dispatches on RunMode (AlignReads | GenomeGenerate)
  params.rs        -- ~48 STAR CLI params via clap derive, --camelCase long names
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
    stitch.rs      -- ✅ Seed clustering + DP stitching + alignment extension (extendAlign)
    score.rs       -- ✅ Scoring functions (gaps, mismatches, splice junctions)
    transcript.rs  -- ✅ Transcript struct (exon coords, CIGAR, scores)
    read_align.rs  -- ✅ Per-read alignment driver
  io/
    mod.rs         -- ✅ Module exports
    fastq.rs       -- ✅ FASTQ reader (plain + gzip, noodles wrapper)
    sam.rs         -- ✅ SAM writer (header + records, noodles wrapper)
    bam.rs         -- ✅ BAM writer (BGZF compression, streaming unsorted output)
  junction/
    mod.rs         -- ✅ GTF parsing, junction database, motif detection, two-pass filtering
    sj_output.rs   -- ✅ SJ.out.tab writer
    gtf.rs         -- ✅ GTF parser (internal)
  chimeric/
    mod.rs         -- ✅ Module exports
    detect.rs      -- ✅ Chimeric detection algorithms (Tier 1: soft-clip, Tier 2: multi-cluster)
    segment.rs     -- ✅ ChimericSegment and ChimericAlignment data structures
    score.rs       -- ✅ Junction type classification, repeat length calculation
    output.rs      -- ✅ Chimeric.out.junction writer (14-column format)
```

## Key Conventions

- **Crate name is `ruSTAR`** — `#![allow(non_snake_case)]` in lib.rs suppresses the crate name warning
- **STAR params use `--camelCase` naming** — clap `#[arg(long = "camelCase")]` maps to snake_case Rust fields
- **Multi-value params** (genomeFastaFiles, readFilesIn, outSAMtype, outSAMattributes, chimOutType, alignSJstitchMismatchNmax, outSJfilterIntronMaxVsReadN) need explicit `num_args`
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
noodles = { version = "0.80", features = ["fastq", "sam", "bam", "bgzf"] }
flate2 = "1"
rayon = "1"
dashmap = "6"

[dev-dependencies]
tempfile = "3"
assert_cmd = "2"
predicates = "3"
```

## Testing Pattern

- Unit tests: `#[cfg(test)]` in each module
- Integration tests: `tests/` directory (future phases)
- Every phase uses differential testing against STAR where applicable
- Test data tiers: synthetic micro-genome → chr22 → full human genome

**Current test status**: 196/196 unit tests passing, non-critical clippy warnings (too_many_arguments × 3, implicit_saturating_sub × 1, manual_contains × 2)

**Note**: Phase 9 integration tests fail due to pathologically repetitive test genomes (50 exact copies of 20bp). These tests need realistic genomes (deferred to Phase 13).

## Current Capabilities

ruSTAR can now perform **end-to-end RNA-seq alignment with two-pass mode and chimeric detection**:
- Generate genome indices from FASTA files
- Read plain or gzipped FASTQ files
- Multi-threaded parallel alignment (via `--runThreadN`)
- Align single-end AND paired-end RNA-seq reads with splice junction detection
- GTF-based junction annotation with scoring bonus
- Two-pass mode for novel junction discovery (`--twopassMode Basic`)
  - Pass 1: Discovers novel junctions from read alignments
  - Pass 2: Re-aligns using both GTF and novel junctions
  - Improves accuracy by ~5-10% for samples with novel junctions
- Output valid SAM or BAM files (`--outSAMtype SAM` or `BAM Unsorted`):
  - Proper headers (@HD, @SQ, @PG)
  - Correct CIGAR strings (M, I, D, N for junctions)
  - Appropriate FLAGS (unmapped, reverse complement, paired-end)
  - MAPQ scores
  - 1-based genomic positions
  - Mate information (RNEXT, PNEXT, TLEN for paired-end)
- Splice junction statistics output (SJ.out.tab, SJ.pass1.out.tab in two-pass mode)
- Chimeric alignment detection for single-end reads (`--chimSegmentMin` > 0):
  - Detects inter-chromosomal fusions (e.g., BCR-ABL chr9→chr22)
  - Detects intra-chromosomal strand breaks
  - Detects large-distance breaks (>1Mb same chr/strand)
  - Soft-clip based detection (>20% clipped reads)
  - Outputs Chimeric.out.junction file (14-column STAR-compatible format)
  - Junction type classification (GT/AG, CT/AC, etc.)
- Print alignment statistics (unique/multi/unmapped percentages)

## Known Issues / Accuracy Gaps (Priority Order)

1. ~~**Position agreement 51%**~~ ✅ FIXED in Phase 13.9 — now **96.3%**.
2. ~~**CIGAR agreement 84.3%**~~ ✅ FIXED in Phase 13.9b — now **97.8%**. Root cause was CIGAR not reversed for reverse-strand reads.
3. ~~**Splice motif detection**~~ ✅ FIXED in Phase 13.9b — `score_gap_with_strand()` converts RC donor position to forward genome coordinates.
4. ~~**33 false junctions**~~ ✅ FIXED in Phase 13.10 — now **2 ruSTAR-only junctions** (94% reduction).
5. ~~**Single-direction seeding**~~ ✅ FIXED in Phase 13.11 — bidirectional L→R + R→L search. Multi-mapped 3.65% → 4.92% (STAR: 7.4%). Shared junctions 9 → 30.
6. **Splice rate still low** (0.9% vs STAR 2.2%): Doubled by R→L seeding but gap remains. GTF annotations would help further.
7. **188 same-chr >500bp apart**: Mostly chrXII rDNA repeats — STAR=MAPQ 1-3, ruSTAR=MAPQ 255. R→L seeding improved but didn't fully resolve.
8. **103 diff-chr disagreements**: 101 are multi-mappers (harmless tie-breaking), 2 have MAPQ mismatches.
9. **60 STAR-only mapped reads**: Stable.
10. **SJ motif strand disagreement** (6/30 shared): CT/AC vs GT/AG — strand assignment issue in SJ.out.tab writer.

## Limitations (to be addressed in future phases)

- No SAM optional tags (AS, NM, NH, HI) - noodles lifetime complexity
- No coordinate-sorted BAM output (unsorted only; use `samtools sort`)
- Chimeric alignment detection:
  - ✅ Single-end chimeric detection (inter-chr, strand breaks, large gaps, soft-clips)
  - ✅ Chimeric.out.junction output file
  - ❌ Paired-end chimeric detection not yet implemented
  - ❌ Tier 3 (re-mapping soft-clipped regions) not yet implemented
- No STARsolo single-cell features (Phase 14, deferred until accuracy parity achieved)
