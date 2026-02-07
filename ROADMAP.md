# ruSTAR Implementation Roadmap

Tracks implementation progress across sessions. Each phase lists its deliverables, files touched, and completion status.

## Phase Dependency Graph

```
Phase 1 (CLI) ✅
  └→ Phase 2 (FASTA/genome) ✅
       └→ Phase 3 (suffix array) ✅
       └→ Phase 4 (seed finding) ✅ ← can load STAR index, no need to wait for Phase 3
            └→ Phase 5 (stitching/scoring) ✅
                 └→ Phase 6 (SAM output) ✅ ← FIRST END-TO-END ALIGNMENT
                      └→ Phase 9 (threading) ← NEXT: Parallel architecture foundation
                           └→ Phase 8 (paired-end) ← Build paired-end on threaded base
                                └→ Phase 7 (splice junctions) ← Additive enhancement
                                     └→ Phase 10 (BAM output)
                                          └→ Phase 11 (two-pass)
                                               └→ Phase 12 (chimeric)
                                                    └→ Phase 13 (optimization)
                                                         └→ Phase 14 (STARsolo)
```

**Phase ordering rationale**: Threading (Phase 9) done first to establish parallel architecture foundation.
Paired-end (Phase 8) builds on threaded infrastructure. GTF/junctions (Phase 7) is mostly additive and
less architecturally disruptive, so done after core parallelism and paired-end are in place.

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

## Phase 4: Index Loading + Seed Finding ✅

**Status**: Complete (forward strand search implemented)

**Goal**: Load genome index from disk and implement MMP (Maximal Mappable Prefix) search via SA binary search.

**Files created**:
- `src/index/io.rs` (251 lines) — Load Genome, SA, and SAindex from disk
  - `GenomeIndex::load()` - Main entry point
  - Correct SA length calculation from file size
  - SAindex header parsing (nbases + genomeSAindexStart array)
  - 1 round-trip test (build → write → load → verify)

- `src/align/seed.rs` (365 lines) — `Seed` struct, MMP search algorithm
  - `Seed::find_seeds()` - Find all seeds in read sequence
  - `find_seed_at_position()` - MMP search at one position
  - `binary_search_sa()` - Binary search SA for exact range
  - `extend_match()` - Extend to maximal length
  - 4 unit tests (exact match, filtering, no-match, position extraction)

**Key implementation details**:
- SAindex used for fast initialization (lookup k-mer → SA position)
- Binary search finds exact SA range [sa_start, sa_end) in O(log n)
- Extends matches greedily to find longest exact match
- Handles padding/sentinels (stops at value 5)
- Skips k-mers containing N
- Min seed length filtering (typically 8-20bp)

**Test results**: 39/39 tests passing, zero clippy warnings

**Deferred to Phase 5**:
- Reverse complement search for reads (will be part of full read alignment)
- Comparison with STAR seed output on real data

**New dependencies**: `memmap2 = "0.9"`, `byteorder = "1"`

---

## Phase 5: Seed Stitching + Alignment Scoring ✅

**Status**: Complete

**Goal**: Cluster seeds into windows, DP stitching, STAR's exact scoring model, build Transcript with CIGAR.

**Files created**:
- `src/align/stitch.rs` (419 lines) — Seed clustering + DP stitching
  - `cluster_seeds()` - Group seeds by genomic proximity (100kb window)
  - `stitch_seeds()` - DP algorithm to connect seeds
  - Anchor seed selection (<=10 loci)
  - Gap penalty scoring (indels vs splice junctions)

- `src/align/score.rs` (179 lines) — Scoring functions
  - `AlignmentScorer` - Match/mismatch/gap scoring
  - Splice junction penalties by motif (GT-AG, GC-AG, AT-AC, non-canonical)
  - Indel vs junction distinction via `alignIntronMin` threshold (21bp)

- `src/align/transcript.rs` (267 lines) — `Transcript` struct
  - Full alignment representation: chr_idx, genome coords, strand, exons
  - CIGAR operations: M, =, X, I, D, N, S, H
  - Score tracking: alignment score, mismatches, gaps, junctions
  - Helper methods: cigar_string(), n_matched(), read_length(), reference_length()

- `src/align/read_align.rs` (122 lines) — Per-read alignment driver
  - `align_read()` - Main entry point: seed finding → clustering → stitching → filtering
  - Transcript filtering: score, mismatches, matches, multimap count
  - Returns sorted transcripts (best first)

**Key implementation details**:
- DP scoring matches STAR's formula: match score - mismatch penalty - gap penalties
- Splice junction detection via dinucleotide motifs (GT-AG = 0, GC-AG = -4, AT-AC = -8, other = -8)
- CIGAR generation during DP traceback (handles M, I, D, N operations)
- Transcript filtering by multiple thresholds (outFilterScoreMin, outFilterMismatchNmax, etc.)

**Test results**: 57/57 tests passing, zero clippy warnings

**New dependencies**: none

---

## Phase 6: SAM Output (First End-to-End Alignment) ✅

**Status**: Complete

**Goal**: FASTQ reader, SAM writer, MAPQ calculation — first complete single-end alignment pipeline.

**Files created**:
- `src/io/fastq.rs` (310 lines) — FASTQ reader wrapper around noodles
  - `FastqReader` - Plain and gzip FASTQ support
  - Auto-detect gzip by file extension (.gz, .gzip)
  - Base encoding (ACGTN → 01234) for genome compatibility
  - Read clipping support (5' and 3')
  - External decompression command support (readFilesCommand)
  - 8 unit tests (encoding, clipping, plain/gzip reading)

- `src/io/sam.rs` (355 lines) — SAM writer wrapper around noodles
  - `SamWriter` - SAM file output with proper headers
  - SAM header generation (@HD, @SQ, @PG lines)
  - CIGAR conversion (ruSTAR → noodles format)
  - FLAGS encoding (unmapped, reverse complement)
  - 0-based → 1-based position conversion
  - MAPQ calculation
  - Handles unmapped reads (FLAG=4)
  - 6 unit tests (header, writer creation, unmapped, transcript conversion, CIGAR)

- `src/mapq.rs` (75 lines) — MAPQ calculation
  - Unique mappers: use outSAMmapqUnique (default 255)
  - Multi-mappers: -10*log10(1-1/n), capped at 255
  - Unmapped: 0
  - 6 unit tests covering all cases

- `src/stats.rs` (200 lines) — Alignment statistics
  - `AlignmentStats` - Track unique/multi/unmapped/too-many-loci counts
  - Percentage calculations and summary printing
  - 8 unit tests for stat tracking and percentages

- `src/lib.rs` — Modified align_reads() function (~100 lines)
  - Full end-to-end pipeline: load index → read FASTQ → align → write SAM
  - Progress logging (every 100K reads)
  - Read limiting support (--readMapNumber)
  - Multi-mapper filtering (--outFilterMultimapNmax)
  - Unmapped read handling (--outSAMunmapped)

- `src/io/mod.rs` — Module exports for fastq and sam

**Key implementation details**:
- Uses **noodles v0.80** library for FASTQ/SAM I/O (pure Rust, no C dependencies)
- Gzip detection by file extension (simpler than magic byte detection)
- SAM optional tags (AS, NM, NH, HI, jM, jI) deferred due to noodles lifetime complexity
- Error handling: added `From<std::io::Error>` for Error enum
- MAPQ formula matches STAR: unique=255, multi=-10*log10(1-1/n)

**Test results**: 84/84 tests passing (up from 57), zero clippy warnings

**Verified**: End-to-end single-end alignment works. SAM output is valid and can be parsed by samtools.

**New dependencies**: `noodles = { version = "0.80", features = ["fastq", "sam"] }`, `flate2 = "1"`

**Known limitations** (to be addressed in later phases):
- Single-end reads only (paired-end in Phase 8)
- SAM output only (BAM in Phase 10)
- No multithreading (Phase 9)
- No SAM optional tags (AS, NM, NH, HI) - deferred due to API complexity
- No GTF-based junction scoring (Phase 7)

---

## Phase 9: Threading (NEXT)

**Status**: Not started

**Goal**: Rayon chunk-based parallelism matching `--runThreadN`. Establish parallel architecture foundation.

**Why this phase order**: Threading affects the entire execution model and is much harder to retrofit
into complex features later. By implementing threading now (before paired-end and GTF features), we
ensure all future features are built with parallelism from the start.

**Key decisions**:
- Use Rayon for data parallelism (chunk-based work distribution)
- Maintain deterministic output regardless of thread count
- Thread-safe statistics aggregation
- Careful handling of SAM output ordering (if required)

**Files to create/modify**:
- `src/lib.rs` — Parallelize main alignment loop with Rayon
- `src/stats.rs` — Thread-safe statistics collection (Arc<Mutex<>> or atomic counters)
- `src/io/sam.rs` — Thread-safe SAM writing (may need buffering/ordering)
- Tests to verify deterministic output with different thread counts

**New dependencies**: `rayon = "1"`

---

## Phase 8: Paired-End Reads

**Status**: Not started (blocked by Phase 9)

**Goal**: Paired FASTQ, concordant/discordant pairing, proper SAM FLAG/TLEN/mate fields.

**Why after threading**: Implementing paired-end on top of established parallel infrastructure
ensures we don't have to retrofit threading into paired-end logic later.

**Files to modify**: `src/io/fastq.rs`, `src/align/read_align.rs`, `src/io/sam.rs`, `src/lib.rs`

---

## Phase 7: Splice Junction Handling

**Status**: Not started (blocked by Phase 8)

**Goal**: GTF parsing, canonical motif detection, JunctionDb, SJ.out.tab output.

**Why after paired-end**: Mostly additive enhancement that works equally well with single-end,
paired-end, threaded or not. Less architecturally disruptive, so done after core features are stable.

**Files to create**: `src/junction/mod.rs`, `src/junction/motif.rs`, `src/junction/gtf.rs`, `src/io/sj_out.rs`

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
