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

## Phase 2: FASTA Loading + Packed Genome

**Status**: Not started

**Goal**: Parse multi-FASTA, concatenate chromosomes with padding (matching STAR's layout), 2-bit pack, append reverse complement, write binary `Genome` file matching STAR's format.

**Files to create/modify**:
- `src/genome/mod.rs` — `Genome` struct: packed sequence, chromosome metadata (names, lengths, start positions)
- `src/genome/fasta.rs` — FASTA parser, chromosome concatenation with padding
- `src/genome/pack.rs` — 2-bit base encoding/decoding (A=0, C=1, G=2, T=3), 4 bases per byte

**STAR reference files**: `genomeFasta.cpp`, `Genome.cpp`, `PackedArray.cpp`

**Key details (from STAR source)**:
- Chromosomes are concatenated with padding of `2^genomeChrBinNbits` boundary alignment
- After forward genome, reverse complement is appended (total length = 2 * padded genome length)
- STAR stores packed genome in `Genome` binary file; `genomeParameters.txt` has metadata
- Also writes `chrName.txt`, `chrLength.txt`, `chrStart.txt`, `chrNameLength.txt`

**Tests**:
- Round-trip: pack then unpack every base combination
- Byte-for-byte comparison of `Genome` file vs STAR output on a small synthetic FASTA
- Verify chromosome start positions match STAR's padding logic

**New dependencies**: none expected (2-bit packing is trivial bit math)

---

## Phase 3: Suffix Array Generation

**Status**: Not started

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
