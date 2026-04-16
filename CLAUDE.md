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

## Current Status

**268 tests passing, 0 clippy warnings.** SE: 8796/8926 compare_sam.py (98.5%), 2.2% splice rate (STAR: 2.2%), 66 shared junctions, **100.0% MAPQ agreement, MAPQ inflation: 0, deflation: 0**. 127 position disagreements (ALL verified as genuine ties). 1 CIGAR-only disagree (ERR12389696.13573895, insertion placement, seed-level tie). **0 STAR-only / 0 ruSTAR-only SE reads**. PE: **8390/8390 both-mapped (0 gap, exact STAR match)**, 0 half-mapped, **4 MAPQ inflations** (rDNA/repeat secondary loci), **98.903% PE faithfulness** (Phase 16.50). See [ROADMAP.md](ROADMAP.md) for detailed phase tracking and [docs/](docs/) for per-phase notes.

## Source Layout

```
src/
  main.rs          -- Thin entry: parse CLI (clap), init logging, call lib::run()
  lib.rs           -- run() dispatches on RunMode (AlignReads | GenomeGenerate)
  params.rs        -- ~52 STAR CLI params via clap derive, --camelCase long names
  error.rs         -- Error enum with thiserror (Parameter, Io, Fasta, Index, Alignment, Gtf)
  mapq.rs          -- MAPQ calculation (STAR lookup table: n=1→255, n=2→3, n≥5→0)
  stats.rs         -- Alignment statistics, Log.final.out writer, UnmappedReason enum
  genome/
    mod.rs         -- Genome struct, padding logic, reverse complement, file writing
    fasta.rs       -- FASTA parser, base encoding (A=0,C=1,G=2,T=3,N=4)
  index/
    mod.rs         -- GenomeIndex (build + load + write)
    packed_array.rs-- Variable-width bit packing (1-64 bits per element)
    suffix_array.rs-- SA construction, custom comparator, strand encoding
    sa_index.rs    -- K-mer lookup table (35-bit entries with flags)
    io.rs          -- Load index from disk (Genome, SA, SAindex)
  align/
    mod.rs         -- Module definition
    seed.rs        -- Seed finding via hierarchical SAindex lookup + MMP search
    stitch.rs      -- Seed clustering + DP stitching + alignment extension (extendAlign)
    score.rs       -- Scoring functions (gaps, mismatches, splice junctions)
    transcript.rs  -- Transcript struct (exon coords, CIGAR, scores)
    read_align.rs  -- Per-read alignment driver
  io/
    mod.rs         -- Module exports
    fastq.rs       -- FASTQ reader (plain + gzip, noodles wrapper)
    sam.rs         -- SAM writer (header + records, noodles wrapper)
    bam.rs         -- BAM writer (BGZF compression, streaming unsorted output)
  junction/
    mod.rs         -- GTF parsing, junction database, motif detection, two-pass filtering
    sj_output.rs   -- SJ.out.tab writer
    gtf.rs         -- GTF parser (internal)
  chimeric/
    mod.rs         -- Module exports
    detect.rs      -- Chimeric detection (Tier 1: soft-clip, Tier 2: multi-cluster)
    segment.rs     -- ChimericSegment and ChimericAlignment data structures
    score.rs       -- Junction type classification, repeat length calculation
    output.rs      -- Chimeric.out.junction writer (14-column format)
```

## Development Philosophy — Match STAR Exactly

**Never backtrack on changes that make the method closer to STAR's logic, even if metrics temporarily regress.** A regression in alignment statistics (position agreement, splice rate, etc.) after making the code more STAR-faithful does NOT mean the change was wrong — it means the implementation is incomplete. More STAR-matching changes are needed to close the gap.

- It is better to go backwards in metrics temporarily to eventually reach 100% matching results, than to keep approximate methods that can never get there.
- Every algorithm, threshold, and code path should match STAR's C++ source exactly. Any deviation — even a small one — should be corrected.
- Do not revert STAR-faithful changes. Instead, investigate what STAR does differently in the surrounding logic and implement that too.
- When a STAR-matching change causes regression, the response is: "What else does STAR do that we don't?" — not "revert and try something different."

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
bstr = "1"
flate2 = "1"
rayon = "1"
dashmap = "6"
chrono = "0.4"

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

**Current test status**: 268/268 tests passing (264 unit + 4 integration), 0 clippy warnings

## Known Issues — Disagreement Root Causes (10k SE yeast)

**127 total position disagreements — ALL verified as genuine ties** (confirmed via STAR debug tracing):

Both tools find identical alignment sets for all 127 disagreements. The primary difference is tie-breaking order (SA iteration order). Neither alignment is more correct than the other.

- **100 diff-chr ties** — same set of alignments, different repeat copy chosen as primary.
- **27 same-chr ties** — same alignment set, different primary due to tie-breaking (includes multi-intron reads where both tools find same 2 alignments but select different primaries).

**1 CIGAR-only disagreement (same position, different CIGAR):**
- `ERR12389696.13573895`: both tools align to XV:218357 MAPQ=255, but ruSTAR gives `100M1I45M4S` (insertion at read pos 100) while STAR gives `108M1I37M4S` (insertion at 108). Root cause: both alignments score AS=133. The 71-base seed is found at RC pos 29 (ruSTAR) vs RC pos 37 (STAR) due to different Lmapped chain paths through a long homopolymer region. Same diagonal, different starting position → different insertion placement. Seed-level tie.

**1 truly actionable SE issue:**

1. **CIGAR-only insertion placement (1 read)** — `ERR12389696.13573895` (see above). Requires matching STAR's exact Lmapped chain to fix.

Previously listed issues now resolved:
- `ERR12389696.18296181` (ruSTAR-only false splice): filtered out by score threshold (score=95/92 < 98 threshold for 150bp read). Both the 1-exon soft-clip AND the 2-exon false splice fail the minimum score gate — read correctly unmapped.
- `ERR12389696.13766843` (STAR-only NM=10 read): ruSTAR already maps this correctly (VII:24449, 28S121M1S, MAPQ=255) — it was incorrectly listed as "STAR-only".

**Phase 16.38 (STAR-faithful filter ordering):**
- Moved dedup + score-range filter (multMapSelect) to run BEFORE quality filters (mappedFilter). STAR's ordering: `multMapSelect → mappedFilter`. Old code ran quality filters first, which could remove the high-scoring primary leaving a lower-scoring secondary as apparent "best". No benchmark change (the specific false splice read was already filtered by score threshold), but structurally correct for edge cases.

**Phase 16.37 (alignIntronMax=0 fix) resolved:**
- `ERR12389696.5825571`: now aligns as `XV:80779 121M607028N13M16S` (exact match with STAR). Root cause: ruSTAR computed a finite intron limit of 589824 when `alignIntronMax=0`, blocking the 607kb intron. STAR uses `>0` guard in `stitchAlignToTranscript.cpp` line 100 — alignIntronMax=0 means no limit. Fix: use `u32::MAX` sentinel in score.rs.
- `ERR12389696.16030539`: MAPQ inflation resolved. Both tools now find two alignments (XV:121224 128M925400N10M12S and XV:598336 128M448288N10M12S), both MAPQ=3. Primary selection differs (tie). MAPQ inflate: 1→0.

See [ROADMAP.md](ROADMAP.md) and [docs/](docs/) for full issue tracking.

## PE Status (Updated 2026-04-09 — Phase 16.45)

**Phase 16.45** (split_working_transcript junction split fix): **PE both-mapped = 8390/8390 (exact STAR match)**. PE MAPQ deflations: **16 → 0**. PE faithfulness: 98.784%.

**Root cause fixed**: `split_working_transcript` used `wt2_junc_start = boundary_idx.min(n_junctions)` to partition junction_shifts between the two mate portions. Since the inter-mate boundary does NOT produce a junction entry, the index offset was wrong: any junctions belonging to the second portion with index < boundary_idx were silently discarded. Fix: use `wt2_junc_start = wt1_junc_end` (junction split point = number of junctions in the first portion). For the specific case boundary_idx=1 (one mate2 exon in reverse cluster), wt1_junc_end=0 but wt2_junc_start was incorrectly set to 1, causing the 145907N junction's shift (jj_r=8) to be lost. finalize_transcript's overhang check then saw empty junction_shifts and passed the spurious 139M145907N11M alignment (11 < threshold 5+8=13 → should have been rejected).

**Current PE parity**: ruSTAR=8390, STAR=8390, **0 gap**. 2 ruSTAR-only FPs (.17779410, .6302610). 2 STAR-only mates (.18919121 half-mapped, .6302610 half-mapped). 24 MAPQ inflations (rDNA/repeat multi-mappers). 0 MAPQ deflations.

## Remaining Limitations (Top 5)

- No coordinate-sorted BAM output (use `samtools sort`) — Phase 17.2
- No PE chimeric detection — Phase 17.3
- No `--quantMode GeneCounts` — Phase 17.8
- No `--outStd SAM/BAM` (stdout output) — Phase 17.6
- No STARsolo single-cell features — Phase 14 (deferred)

See [docs/phase17_features.md](docs/phase17_features.md) for full feature status.
