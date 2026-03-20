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

**268 tests passing, 0 clippy warnings.** SE: 99.7% position agreement (adjusted), 99.9% CIGAR (pos-agreeing), 2.2% splice rate (STAR: 2.2%), 67 shared junctions, 99.9% MAPQ agreement, 26 actionable disagreements, 0 STAR-only / 0 ruSTAR-only mapped. MAPQ inflation: 4, deflation: 4. PE: 8295/8390 both-mapped (−95 under STAR), 0 half-mapped, 32 ruSTAR-only false positives (down from 144 via Phase 16.31 overlap check fix), 127 STAR-only missed. See [ROADMAP.md](ROADMAP.md) for detailed phase tracking and [docs/](docs/) for per-phase notes.

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

126 total position disagreements, 26 actionable:

1. **Diff-chr multi-mapper ties (100 reads)** — Same CIGAR/MAPQ, different positions in repeat copies (rDNA). **Unavoidable** without matching STAR's internal SA iteration order.
2. **Same-chr repeat ties (~19 reads)** — rDNA/segmental dup copies on XII and IV. Same CIGAR/MAPQ, different position. **Unavoidable.**
3. **Wrong intron choice (4 reads)** — ruSTAR picks a different large intron than STAR (e.g., 170kb vs 84kb) at multi-mapping loci on VII/XV.
4. **MAPQ inflation (4 reads)** — ruSTAR=255 (unique) when STAR finds an additional spliced or indel alternative that ruSTAR misses.
5. **MAPQ deflation (4 reads)** — ruSTAR generates extra unspliced secondary that makes the primary multi-mapped. 2 pre-existing + 2 unmasked by Phase 16.29 scoring fix.
6. **ruSTAR-only (1 read)** — `ERR12389696.18296181` false splice with adapter contamination, 279kb intron.
7. **STAR-only (1 read)** — `ERR12389696.13766843` high-mismatch read (NM=10), not mapped by ruSTAR.

See [ROADMAP.md](ROADMAP.md) and [docs/](docs/) for full issue tracking.

## PE False Positives — Root Cause (Confirmed 2026-03-19)

**144 ruSTAR-only false positives** are non-palindromic overlapping PE pairs rejected by STAR via `outFilterMatchNminOverLread` applied to the **combined read length**, not per-mate.

- STAR `Lread` for a PE combined read = `len1 + SPACER(11) + len2` ≈ 311 for 150bp reads
- `outFilterMatchNminOverLread` default = 0.66 → threshold = `0.66 × 310 = 204.6 → 205` matched bases required
- Overlapping false-positive pairs have `nMatch ≈ 130-150` (only one mate's worth) → fail `mappedFilter` → `unmapped=short`
- ruSTAR's `combined_score_threshold = 0.66 × (len1 + SPACER + len2 - 1) = 198` is correct for the score, but `nMatch` filtering against the combined Lread is missing

**Fix needed**: After `split_working_transcript`, check that `t1.n_match + t2.n_match >= 0.66 × (combined_Lread - 1)`. This is the `outFilterMatchNminOverLread` check on the combined read.

**Confirmed by**: STAR debug tracing (`test/debug_star.sh` + instrumented STAR at `/home/jamfer/Dropbox/Bioinformatics/tools/repos/STAR/source/STAR`). Example: `ERR12389696.10060311` has `nMatch=137`, threshold=198 → `MAPPEDFILTER-RESULT: unmapped=short`.

## Remaining Limitations (Top 5)

- No coordinate-sorted BAM output (use `samtools sort`) — Phase 17.2
- No PE chimeric detection — Phase 17.3
- No `--quantMode GeneCounts` — Phase 17.8
- No `--outStd SAM/BAM` (stdout output) — Phase 17.6
- No STARsolo single-cell features — Phase 14 (deferred)

See [docs/phase17_features.md](docs/phase17_features.md) for full feature status.
