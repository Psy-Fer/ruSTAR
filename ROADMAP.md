# ruSTAR Implementation Roadmap

Tracks implementation progress across sessions. Each phase lists its deliverables, files touched, and completion status.

## Phase Dependency Graph

```
Phase 1 (CLI) ✅
  └→ Phase 2 (FASTA/genome)
       └→ Phase 3 (suffix array)
       └→ Phase 4 (seed finding) ← can load STAR index, no need to wait for Phase 3
            └→ Phase 5 (stitching/scoring)
                 └→ Phase 6 (SAM output) ← FIRST END-TO-END ALIGNMENT
                      ├→ Phase 7 (splice junctions)
                      ├→ Phase 8 (paired-end)
                      └→ Phase 9 (threading)
                           └→ Phase 10 (BAM output)
                                └→ Phase 11 (two-pass)
                                     └→ Phase 12 (chimeric)
                                          └→ Phase 13 (optimization)
                                               └→ Phase 14 (STARsolo)
```

---

## Phase 1: CLI + Parameters + Infrastructure ✅

**Status**: Complete

**Deliverables**:
- `src/params.rs` — `Parameters` struct with ~40 STAR CLI params via clap derive
  - Enums: `RunMode`, `OutSamFormat`, `OutSamSortOrder`, `OutSamType`, `OutSamUnmapped`, `OutFilterType`, `TwopassMode`
  - `Parameters::out_sam_type()` parses raw tokens into structured type
  - `Parameters::validate()` checks cross-field constraints
  - 9 unit tests covering defaults, genome generate, typical align, scoring overrides, validation, outSAMtype parsing, chimeric params, SJ stitch mismatch
- `src/error.rs` — `Error` enum (Parameter, Io, Fasta, Index, Alignment, Gtf) with thiserror
- `src/lib.rs` — `run()` dispatcher (AlignReads | GenomeGenerate), `#![allow(non_snake_case)]`
- `src/main.rs` — thin entry: parse CLI, init env_logger, call `lib::run()`
- Module stubs: `genome/`, `index/`, `align/`, `junction/`, `io/` (each with `mod.rs`)
- `Cargo.toml` — clap, anyhow, thiserror, log, env_logger; dev: tempfile, assert_cmd, predicates

**Verified**: `cargo build` clean, `cargo clippy` zero warnings, `cargo fmt --check` pass, `cargo test` 9/9 pass

**Key decisions**:
- Multi-value params need explicit `num_args` (genomeFastaFiles, readFilesIn, outSAMtype, etc.)
- Negative defaults need `allow_hyphen_values = true` (scoreGapNoncan, readMapNumber, etc.)
- `outSAMtype` stored as raw `Vec<String>`, parsed on demand via method

---

## Phase 2: FASTA Loading + Packed Genome ✅

**Status**: Complete

**Goal**: Parse multi-FASTA, concatenate chromosomes with padding (matching STAR's layout), encode bases, append reverse complement, write binary `Genome` file matching STAR's format.

**Files created**:
- `src/genome/mod.rs` — `Genome` struct with padding logic, reverse complement, file writing, chromosome lookup
- `src/genome/fasta.rs` — FASTA parser with base encoding (A=0, C=1, G=2, T=3, N=4)

**Key implementation details**:
- **NOT 2-bit packed** — STAR's Genome file uses 1 byte per base, not 2 bits
- Base encoding: A=0, C=1, G=2, T=3, N/other=4, padding=5 (GENOME_SPACING_CHAR)
- Padding formula: `n = ((n+1)/binSize + 1) * binSize` where `binSize = 1 << genomeChrBinNbits`
- This formula always advances to at least the next bin boundary (not a simple ceiling)
- Reverse complement stored in second half of buffer (positions `n_genome..2*n_genome-1`), NOT written to Genome file
- `chrStart.txt` has `n_chr_real + 1` entries (last entry = n_genome)
- All padding bytes initialized to value 5, not 0
- Metadata files: `chrName.txt`, `chrLength.txt`, `chrStart.txt`, `chrNameLength.txt`, `genomeParameters.txt`

**Tests**: 19 unit tests covering:
- FASTA parsing (single/multi-chromosome, case insensitivity, error handling)
- Padding calculations (verified with multiple bin sizes)
- Reverse complement correctness
- Chromosome position lookups
- End-to-end genome generation validated with synthetic 2-chromosome FASTA

**Verified**: `cargo test` (19/19 pass), `cargo clippy` (clean), manual inspection of output files matches STAR format

**New dependencies**: none

---

## Phase 3: Suffix Array Generation ✅

**Status**: Complete (including SAindex refinement)

**Goal**: Build suffix array from packed genome, matching STAR's SA format. Write `SA` and `SAindex` binary files.

**Files to create/modify**:
- `src/index/mod.rs` — `GenomeIndex` struct
- `src/index/suffix_array.rs` — SA construction (port STAR's comparator, Rayon parallel sort)
- `src/index/sa_index.rs` — Pre-computed prefix lookup table
- `src/index/io.rs` — Read/write index files in STAR's binary format

**STAR reference files**: `suffixArrayFuns.cpp`, `PackedArray.cpp`, `genomeGenerate.cpp`

**Key details**:
- STAR uses `PackedArray` for compressed SA storage (variable bits per element)
- SA is sorted using custom comparator that compares 8-byte words from packed genome
- Prefix-based bucketing (4-nt prefixes) for parallel chunk sorting
- `SAindex` is a lookup table: for a given k-mer prefix, gives SA range to search

**Tests**: Binary-diff SA/SAindex vs STAR-generated; load ruSTAR index in STAR

**New dependencies**: `rayon`, `byteorder`

---

## Phase 4: Index Loading + Seed Finding

**Status**: Not started

**Goal**: Memory-map STAR index files, implement MMP (Maximal Mappable Prefix) search via SA binary search.

**Files to create/modify**:
- `src/index/io.rs` — memory-mapped loading of SA, SAindex, Genome
- `src/align/seed.rs` — `Seed` struct, MMP search on both strands

**Key details**:
- Can start before Phase 3 by loading a STAR-generated index
- For each read position, binary search SA for longest exact match using SAindex to narrow range
- Search both forward and reverse complement

**Tests**: Compare seeds found vs STAR on synthetic + chr22 reads

**New dependencies**: `memmap2`

---

## Phase 5: Seed Stitching + Alignment Scoring

**Status**: Not started

**Goal**: Cluster seeds into windows, DP stitching, STAR's exact scoring model, build Transcript with CIGAR.

**Files to create/modify**:
- `src/align/stitch.rs` — Seed clustering, DP stitching
- `src/align/score.rs` — Scoring functions (gaps, mismatches, splice junctions)
- `src/align/transcript.rs` — `Transcript` struct (exon coords, CIGAR, scores)
- `src/align/read_align.rs` — Per-read alignment driver

**Tests**: Compare POS, CIGAR, AS tag for 1000 reads vs STAR

---

## Phase 6: SAM Output (First End-to-End Alignment)

**Status**: Not started

**Goal**: FASTQ reader, SAM writer, MAPQ calculation — first complete single-end alignment pipeline.

**Files to create/modify**:
- `src/io/fastq.rs` — FASTQ reader (plain + gzipped)
- `src/io/sam.rs` — SAM header + records + tags
- `src/mapq.rs` — MAPQ: 255 unique, -10*log10(1-1/n) multi
- `src/stats.rs` — Basic alignment statistics

**Tests**: Align 10K single-end reads to chr22, field-by-field SAM comparison vs STAR

**New dependencies**: `flate2`

---

## Phase 7: Splice Junction Handling

**Status**: Not started

**Goal**: GTF parsing, canonical motif detection, JunctionDb, SJ.out.tab output.

**Files**: `src/junction/mod.rs`, `src/junction/motif.rs`, `src/junction/gtf.rs`, `src/io/sj_out.rs`

---

## Phase 8: Paired-End Reads

**Status**: Not started

**Goal**: Paired FASTQ, concordant/discordant pairing, proper SAM FLAG/TLEN/mate fields.

---

## Phase 9: Threading

**Status**: Not started

**Goal**: Rayon chunk-based parallelism matching `--runThreadN`. Identical output regardless of thread count.

---

## Phase 10: BAM Output + Coordinate Sorting

**Status**: Not started

**New dependencies**: `noodles-bam` or `rust-htslib`

---

## Phase 11: Two-Pass Mode

**Status**: Not started

---

## Phase 12: Chimeric Detection

**Status**: Not started

---

## Phase 13: Performance Optimization

**Status**: Not started

---

## Phase 14: STARsolo (Single-Cell)

**Status**: Not started
