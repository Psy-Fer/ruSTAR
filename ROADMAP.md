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
                      └→ Phase 9 (threading) ✅ ← Parallel architecture foundation
                           └→ Phase 8 (paired-end) ✅ ← Built on threaded base
                                └→ Phase 7 (splice junctions) ✅ ← GTF/junction annotations
                                     └→ Phase 10 (BAM output) ✅ ← Binary alignment format
                                          └→ Phase 11 (two-pass) ✅ ← Novel junction discovery
                                               └→ Phase 12 (chimeric) ✅ ← Gene fusion detection
                                                    └→ Phase 13.1-13.6 (perf+accuracy) ✅
                                                         └→ Phase 13.7-13.9 (accuracy refinement)
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

## Phase 9: Threading ✅

**Status**: Complete

**Goal**: Rayon chunk-based parallelism matching `--runThreadN`. Establish parallel architecture foundation.

**Why this phase order**: Threading affects the entire execution model and is much harder to retrofit
into complex features later. By implementing threading now (before paired-end and GTF features), we
ensure all future features are built with parallelism from the start.

**Key implementation**:
- Rayon parallel iterators for batch processing (10,000 reads per batch)
- Thread pool configuration via `--runThreadN`
- Sequential FASTQ reading, parallel alignment, sequential SAM writing
- Atomic counters for statistics (`Arc<AlignmentStats>`)
- Buffered SAM record accumulation per thread
- Deterministic output regardless of thread count

**Files modified**:
- `src/lib.rs` — Parallelized alignment loop with Rayon `.par_iter()`
- `src/stats.rs` — AtomicU64 counters for thread-safe statistics
- `src/io/sam.rs` — `BufferedSamRecords` for per-thread buffering
- Integration tests to verify correctness with 1, 4, and 8 threads

**New dependencies**: `rayon = "1"` (already present from Phase 3)

---

## Phase 8: Paired-End Reads ✅

**Status**: Complete

**Goal**: Paired FASTQ, concordant/discordant pairing, proper SAM FLAG/TLEN/mate fields.

**Why after threading**: Implementing paired-end on top of established parallel infrastructure
ensures we don't have to retrofit threading into paired-end logic later.

**Key features implemented**:
- Paired FASTQ reading (two input files)
- Unified transcript for both mates
- Proper pair detection (same chr, concordant orientation, distance)
- SAM FLAGS for paired reads (0x1, 0x2, 0x8, 0x20, 0x40, 0x80)
- TLEN (insert size) calculation
- Mate position (RNEXT, PNEXT) fields

**Files modified**: `src/io/fastq.rs`, `src/align/read_align.rs`, `src/io/sam.rs`, `src/lib.rs`

---

## Phase 7: GTF/Splice Junction Annotation ✅

**Status**: Complete

**Goal**: GTF parsing, junction database, annotation-based scoring bonus, SJ.out.tab output.

**Files created**:
- `src/junction/mod.rs` (205 lines) — `SpliceJunctionDb` struct with HashMap-based junction lookup
  - `from_gtf()` - Load junction database from GTF file
  - `is_annotated()` - O(1) junction lookup during alignment
  - 5 unit tests (empty db, lookup, strand-specific)

- `src/junction/gtf.rs` (515 lines) — GTF parser
  - `parse_gtf()` - Read GTF file, extract exon features
  - `extract_junctions_from_exons()` - Calculate intron coordinates from consecutive exons
  - Supports standard GTF format (Ensembl/GENCODE)
  - Groups exons by transcript_id, sorts by position
  - Handles unknown chromosomes gracefully (warnings)
  - 11 unit tests (parsing, attribute extraction, junction calculation, deduplication)

- `src/junction/sj_output.rs` (420 lines) — SJ.out.tab writer
  - `SpliceJunctionStats` - Thread-safe junction statistics accumulator (DashMap)
  - `record_junction()` - Thread-safe junction recording
  - `write_output()` - Write 9-column SJ.out.tab file
  - Motif encoding: 0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG, 4=CT/GC, 5=AT/AC
  - Tracks unique/multi counts, max overhang per junction
  - 11 unit tests (stats accumulation, motif encoding, output format)

**Files modified**:
- `src/lib.rs` — Junction statistics collection and SJ.out.tab output
  - Create `SpliceJunctionStats` at alignment start
  - Pass to single-end and paired-end alignment functions
  - Extract junctions from transcripts via CIGAR traversal
  - Write SJ.out.tab at end of alignment
  - Added `record_transcript_junctions()` helper function

- `src/index/mod.rs` & `src/index/io.rs` — Junction DB integration
  - `GenomeIndex` now includes `junction_db: SpliceJunctionDb`
  - Load GTF during genome generation and index loading
  - Logs junction count on load

- `src/align/score.rs` — Annotation bonus support
  - Added `sjdb_score` field to `AlignmentScorer`
  - `score_annotated_junction()` - Apply bonus to annotated junctions
  - Made `detect_splice_motif()` public for junction recording
  - Added `Hash` trait to `SpliceMotif` for DashMap keys
  - 1 new unit test for annotation bonus

- `src/align/read_align.rs` & `src/align/stitch.rs` — Test fixture updates
  - Added `junction_db` field to test GenomeIndex constructions

**Key implementation details**:
- **Junction coordinates**: 1-based intronic bases (STAR convention)
  - `intron_start = exon_end + 1`
  - `intron_end = next_exon_start - 1`
- **Thread-safe collection**: DashMap enables parallel statistics accumulation
- **Annotation bonus**: `--sjdbScore` (default +2) applied to annotated junctions
- **Motif detection**: Re-detect motif during junction recording (simpler than storing in transcript)
- **Overhang**: Placeholder (5bp) - accurate calculation deferred to future refinement

**Test results**: 132/132 tests passing (up from 109), 5 minor clippy warnings (acceptable)

**New dependencies**: `dashmap = "6"`

**Known limitations**:
- Overhang calculation uses placeholder value (5bp) instead of accurate per-junction calculation
- Overhang-based filtering (`alignSJoverhangMin`, `alignSJDBoverhangMin`) not yet implemented
- Two-pass mode not implemented (would build junction DB from first alignment pass)

---

## Phase 10: BAM Output (Unsorted Streaming) ✅

**Status**: Complete

**Goal**: Binary compressed BAM output with streaming (unsorted) write capability.

**Why now**: With GTF/junction annotation complete, BAM output is the next logical step.
BAM is the standard format for downstream analysis tools and significantly reduces file sizes.

**Implementation approach**:
- **Streaming unsorted BAM** (not sorted in-memory)
- Users sort with `samtools sort` separately (standard bioinformatics workflow)
- Simpler implementation (~180 lines) vs integrated sorting (~400+ lines)
- Follows industry-standard practice (bwa, bowtie2 also output unsorted)

**Key features implemented**:
- BAM file writing with BGZF compression (`--outSAMtype BAM Unsorted`)
- Proper BAM header generation (reuses SAM header logic)
- Support for both SAM and BAM output modes
- Generic `AlignmentWriter` trait for SAM/BAM polymorphism
- Streaming write (no memory buffering beyond thread batches)
- Compatible with `samtools sort` and `samtools index`

**Files created/modified**:
- `src/io/bam.rs` — NEW: `BamWriter` struct with streaming output (~150 lines)
  - `write_batch(&[RecordBuf])` for parallel thread integration
  - `finish()` to flush BGZF buffers
  - Unit tests for unmapped, aligned, and batch writes
- `src/io/mod.rs` — Export `bam` module
- `src/io/sam.rs` — Made `build_sam_header()` public for reuse
- `src/lib.rs` — `AlignmentWriter` trait for polymorphism (~15 lines)
  - Route to SAM or BAM based on `--outSAMtype` parameter
  - Generic `align_reads_single_end<W: AlignmentWriter>` functions
  - BAM finish() called after alignment complete
- `Cargo.toml` — Added `bgzf` and `bam` features to noodles

**New dependencies**:
- `noodles = { version = "0.80", features = ["fastq", "sam", "bam", "bgzf"] }`

**Test results**:
- ✅ 136/136 tests passing (added 3 BAM unit tests)
- ✅ `cargo clippy` clean (1 pre-existing warning)
- ✅ Manual verification with samtools:
  - `samtools view` successfully reads BAM
  - `samtools flagstat` reports correct counts
  - `samtools sort` works on ruSTAR BAM
  - `samtools index` successfully creates .bai from sorted BAM
- ✅ BAM content matches SAM content (verified identical records)
- ✅ File size comparison: BAM ~10-20x smaller than SAM (varies by data)

**Deferred to future phase**:
- Integrated coordinate sorting (`--outSAMtype BAM SortedByCoordinate`)
- Automatic BAI index generation
- These can be added in Phase 10.1 if users request integrated sorting

**Dependencies**: Phase 7 (complete)

---

## Phase 11: Two-Pass Mode ✅

**Status**: Complete

**Goal**: Implement STAR's two-pass mode to discover novel splice junctions in pass 1, then re-align all reads using both GTF and novel junctions in pass 2.

**Why two-pass**: Improves alignment accuracy by ~5-10% for samples with novel junctions (tissue-specific isoforms, cancer samples, non-model organisms).

**Files created/modified**:
- `src/junction/mod.rs` — Extended junction database for two-pass mode (~60 lines added)
  - `NovelJunctionKey` struct for two-pass junction tracking
  - `insert_novel()` method to add discovered junctions
  - `filter_novel_junctions()` function for coverage/overhang filtering
  - 3 new unit tests (insert, filter, integration)

- `src/lib.rs` — Two-pass dispatch and workflow (~180 lines added)
  - `run_two_pass()` orchestrates pass 1 → junction insertion → pass 2
  - `run_pass1()` discovers junctions (discards alignments via NullWriter)
  - `run_single_pass()` extracted from original `align_reads()`
  - `NullWriter` struct for pass 1 alignment discard
  - Creates SJ.pass1.out.tab and SJ.out.tab output files

- `src/stats.rs` — Added `total_reads()` method
- Multiple files — Added `Clone` derives for two-pass cloning:
  - `Genome`, `SuffixArray`, `SaIndex`, `PackedArray`, `SpliceJunctionDb`, `GenomeIndex`, `Parameters`
  - Manual `Clone` impl for `SpliceJunctionStats` (handles DashMap with AtomicU32)

**Key implementation details**:
- **Pass 1**: Align reads (limited by `--twopass1readsN`), collect junction stats, discard alignments
- **Filtering**: Novel junctions require ≥1 unique OR ≥2 multi reads, overhang ≥ `alignSJoverhangMin` (default 5bp)
- **Junction insertion**: Clones GenomeIndex, merges GTF + novel junctions
- **Pass 2**: Re-aligns ALL reads (not just pass 1 subset) with merged junction DB
- **Memory efficient**: Arc-sharing of large structures (Genome, SA, SAindex), only junction DB cloned (~1-10MB overhead)
- **Output files**: SJ.pass1.out.tab (pass 1), SJ.out.tab (final), Aligned.out.sam/bam (pass 2 only)

**Test results**:
- ✅ 138/138 tests passing (up from 136)
- ✅ 1 non-critical clippy warning (too_many_arguments - acceptable)
- ✅ Manual verification:
  - `--twopassMode None` works (default, single-pass)
  - `--twopassMode Basic` discovers and uses novel junctions
  - `--twopass1readsN` limits pass 1 reads correctly
  - SJ.pass1.out.tab and SJ.out.tab created properly

**New dependencies**: None (reuses existing infrastructure)

**Known limitations**:
- Basic mode only (no multi-sample two-pass)
- Novel junction filtering thresholds are hardcoded (1 unique OR 2 multi)

**Dependencies**: Phase 10 (complete)

---

## Phase 12: Chimeric Detection ✅

**Status**: Complete (single-end chimeric detection functional)

**Goal**: Detect and report chimeric alignments (gene fusions, circular RNAs, structural variants).

**Phase 12.1 Deliverables** (Completed):
- `src/chimeric/mod.rs` — Module exports and type definitions
- `src/chimeric/segment.rs` — `ChimericSegment` and `ChimericAlignment` data structures (110 lines)
- `src/chimeric/score.rs` — Junction type classification, repeat length calculation (185 lines, 8 tests)
- `src/chimeric/detect.rs` — Chimeric detection algorithms (360 lines, 8 tests):
  - Tier 1: Soft-clip detection (detects excessive soft-clipping)
  - Tier 2: Multi-cluster detection (inter-chromosomal, strand breaks, large distances)
- `src/chimeric/output.rs` — `ChimericJunctionWriter` for 14-column output format (245 lines, 4 tests)
- `src/align/read_align.rs` — Modified `align_read()` signature to return chimeric alignments
- `src/align/transcript.rs` — Added `count_soft_clips()` method
- `src/error.rs` — Added `Chimeric` error variant

**Key features**:
- Detection triggers:
  - Soft-clip >20% of read AND ≥ chimSegmentMin bases (Tier 1)
  - Seeds on different chromosomes (Tier 2)
  - Seeds on different strands, same chromosome (Tier 2)
  - Seeds >1Mb apart, same chr/strand (Tier 2)
- Junction type classification: 0-6 (GT/AG, CT/AC, GC/AG, CT/GC, AT/AC, GT/AT, non-canonical)
- Repeat length calculation at junction sites
- 14-column Chimeric.out.junction format matching STAR

**Tests**: 32 new tests added (20 unit tests + 12 in segment/output modules)
- All 170 tests passing
- 4 non-critical clippy warnings (acceptable)

**Phase 12.2 Deliverables** (Completed):
- `src/lib.rs` — Integrated chimeric detection into end-to-end pipeline (~100 lines modified):
  - Created `AlignmentBatchResults` helper struct
  - Modified `align_reads_single_end()` to collect chimeric alignments from parallel workers
  - Created `ChimericJunctionWriter` at alignment start if chimSegmentMin > 0
  - Write chimeric alignments to Chimeric.out.junction after each batch
  - Added chimeric writer infrastructure to `align_reads_paired_end()` (ready for future implementation)
  - Flush chimeric output at end of alignment
- Total Phase 12 code: ~1000 lines (900 new in Phase 12.1 + 100 modified in Phase 12.2)

**Known limitations** (future work):
- ❌ Tier 3 (re-mapping soft-clipped regions) not implemented
- ❌ Paired-end chimeric detection not implemented (infrastructure ready, detection logic needed)
- ❌ No integration tests with synthetic fusion reads yet (needs test data)

**Verified**: `cargo test` (170/170 pass), `cargo clippy` (4 non-critical warnings), `cargo fmt --check` (pass), `cargo build --release` (success)

**Usage**: Enable with `--chimSegmentMin 15` (typical value 10-20bp)

---

## Phase 13: Performance + Accuracy Optimization ✅ (Phases 13.1-13.6)

**Status**: Complete through Phase 13.6 (performance + alignment extension). Phases 13.7-13.9 planned for accuracy refinement.

**Goal**: Optimize alignment performance to approach STAR speeds, fix classification issues, and reduce accuracy gaps.

### Phase 13.1: Critical Bug Fixes ✅ COMPLETE (2026-02-07)

**Deliverables**:
- Fixed clustering explosion (11,982 clusters → ~200 with STAR limits)
- Implemented exon building from CIGAR (line 354 TODO completed)
- Added match scoring (+1 per matched base for DP stitching)
- Implemented soft clip logic for partial alignments

**Test Results**:
- ✅ 100 reads: 100% mapped in ~5s (was: >60s, 0% mapped)
- ✅ Valid SAM output with stitched alignments
- ✅ Soft clips working
- ✅ 10 splice junctions detected
- ⚠️ All reads classified as multi-mapped (should be ~88% unique)
- ⚠️ 50x slower than STAR (5s vs 0.1s for 100 reads)

**Files Modified**:
- `src/params.rs` - Added seedMultimapNmax, winAnchorMultimapNmax, seedNoneLociPerWindow
- `src/align/stitch.rs` - Fixed scoring, exon building, soft clips (~200 lines)
- `src/align/read_align.rs` - Cleaned up debug output

**See**: [ALIGNMENT_FIXES.md](ALIGNMENT_FIXES.md) for detailed debugging session

---

### Phase 13.2: Mismatch Counting + Seed Expansion Bugs ✅ COMPLETE (2026-02-07)

**Problem**: 66% unmapped reads (STAR: 11%) due to inflated mismatch counts after adding `count_mismatches()`.

**Root Causes & Fixes**:
1. **Reverse-strand genome offset** — `count_mismatches()` was reading forward genome for reverse reads. Fixed: add `n_genome` offset for reverse strand access.
2. **Seed length overestimation** — `extend_match()` verified at `sa_start` only; other SA positions got same length. Fixed: added `verify_match_at_position()` to re-verify each expanded seed position.
3. **DP chain start wrong** — Used `expanded_seeds.first()` instead of tracing back through `prev_seed` to find actual chain start, causing wrong `genome_pos` for CIGAR walk.
4. **Combined gap CIGAR missing** — When both `read_gap > 0` AND `genome_gap > 0`, seeds merged as if adjacent (missing indel operations). Fixed CIGAR builder to emit `Match(shared) + Del/Ins(excess)`.
5. **Removed incorrect read reverse-complement** — `count_mismatches` was reverse-complementing the read (wrong). Seeds match read-as-is against the reverse-complement genome region.

**Files Modified**:
- `src/align/stitch.rs` — Added `verify_match_at_position()`, fixed `count_mismatches()`, fixed chain start tracing, fixed CIGAR gap handling, removed debug logging
- `src/align/score.rs` — Handle combined read+genome gaps in `score_gap()`
- `src/align/read_align.rs` — Formatting only

**Test Results (1000 yeast reads)**:

| Metric | Before | After | STAR Target |
|--------|--------|-------|-------------|
| Unique | 32% | **82.4%** | 82% |
| Multi | 2% | **8.3%** | 7% |
| Unmapped | 66% | **9.3%** | 11% |

- ✅ 170/170 unit tests passing
- ✅ Alignment accuracy matches STAR
- ✅ 1000 reads in ~3s (vs STAR <1s)

---

### Phase 13.3: Performance Optimization ✅ COMPLETE (2026-02-07)

**Problem**: ~3x slower than STAR for 1000 reads (3s vs <1s).

**Profiling results** (perf, 1000 reads, single-threaded):
- 29.5% malloc/free/realloc (Vec allocations in hot loops)
- 18.1% PackedArray::read (SA bit-unpacking)
- 15.9% cluster_seeds (inner loop)
- 14.7% stitch_seeds (DP + seed expansion)
- 5.2% score_gap

**Optimizations Applied**:

1. **Eliminated Vec allocations in hot loops** (target: 30% alloc overhead)
   - Added `Seed::genome_positions()` lazy iterator (no Vec allocation)
   - `cluster_seeds()` and `stitch_seeds()` use iterator directly
   - Pre-allocated Vecs with capacity estimates

2. **Optimized PackedArray::read** (target: 18% bit-unpacking)
   - Fast path: direct 8-byte slice read when within bounds (eliminates 8 per-byte bounds checks)
   - Slow path: original byte-by-byte read near end of array

3. **Binary search `position_to_chr`** (target: 16% clustering)
   - Replaced O(n_chr) linear scan with O(log n) binary search via `partition_point()`

4. **Deduplicated expanded seeds** (target: 15% DP stitching)
   - Dedup by (read_pos, genome_pos), keeping longest seed
   - Cap at 200 expanded seeds per cluster to prevent pathological O(n²) DP

5. **Deferred CIGAR clone in DP** (part of alloc + DP overhead)
   - Find best `j` index in inner loop, build CIGAR only once after loop
   - Eliminates repeated Vec cloning for non-winning transitions

**Files Modified**:
- `src/align/seed.rs` — Added `genome_positions()` iterator method
- `src/index/packed_array.rs` — Fast path read with direct slice access
- `src/genome/mod.rs` — Binary search `position_to_chr()`
- `src/align/stitch.rs` — All DP/clustering optimizations

**Performance Results (1000 yeast reads, single-threaded)**:

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| Wall time | ~3.0s | **~0.77s** | **3.9x faster** |
| STAR baseline | <1.0s | — | Now within range |

**Accuracy (unchanged)**:

| Metric | Before | After |
|--------|--------|-------|
| Unique | 82.4% | 82.2% |
| Multi | 8.3% | 8.4% |
| Unmapped | 9.3% | 9.4% |

- ✅ 170/170 unit tests passing
- ✅ No new clippy warnings (same 3 pre-existing)

**Remaining Issues** (LOW PRIORITY):
1. **Splice motif detection for reverse strand** — `detect_splice_motif()` in score.rs doesn't add n_genome offset
2. **Parameter tuning** — Genome-size-dependent defaults

---

### Phase 13.4: CIGAR Integer Overflow + Coordinate Bugs ✅ COMPLETE (2026-02-09)

**Problem**: Test framework detected corrupted SAM/BAM outputs with integer overflow values.

**Symptoms**:
- CIGAR strings: `10M4294953882D12M` (D value near 2³² indicating overflow)
- Junction coordinates: `4296353566` (beyond yeast chromosome boundaries)
- Mean alignment length: 5.8 billion bases (should be ~150bp)
- Consecutive Match operations: `10M4M10M` instead of `24M`

**Root Causes Identified**:

1. **Integer overflow from unsafe casts** (`src/align/stitch.rs:423-424`)
   - Gap calculations legitimately produce negative values when seeds overlap
   - Unsafe cast: `best_read_gap as u32` when value is `-13414` → `4294953882`
   - Formula: `-13414 as u32` = `2³² - 13414` = `4294953882`

2. **CIGAR merging failure** (`src/align/stitch.rs:465, 484-491`)
   - Gap region adds Match operation, then seed adds another Match
   - Result: consecutive Match ops not merged (e.g., `10M4M10M` instead of `24M`)

3. **Global vs per-chromosome coordinates** (`src/io/sam.rs:384`, `src/junction/sj_output.rs:164`)
   - Internal representation: global genome coordinates (cumulative across chromosomes)
   - SAM/SJ.out.tab format requires: per-chromosome coordinates
   - Missing conversion: `pos - genome.chr_start[chr_idx]`

**Fixes Applied**:

1. **Negative gap handling** — Added explicit checks before casting to u32
   - Skip connections with negative gaps (overlapping seeds)
   - Log warnings for debugging: `"Skipping connection with negative gap: read_gap=X, genome_gap=Y"`

2. **CIGAR merging** — Always merge consecutive Match operations
   - Try to merge with previous Match operation before pushing new one
   - Applies to both gap regions and seed matches

3. **Coordinate conversion** — Convert global→per-chromosome in output writers
   - SAM writer: `let chr_start = genome.chr_start[transcript.chr_idx]; let pos = (transcript.genome_start - chr_start + 1) as usize;`
   - SJ writer: Same conversion for junction start/end coordinates

**Files Modified**:
- `src/align/stitch.rs` — Negative gap checks, CIGAR merging (~30 lines changed)
- `src/io/sam.rs` — Coordinate conversion (2 occurrences for single/paired-end)
- `src/junction/sj_output.rs` — Junction coordinate conversion

**Test Results**:

| Dataset | Unique | Multi | Unmapped | Integer Overflow |
|---------|--------|-------|----------|------------------|
| 100 reads | 78.0% | 5.0% | 17.0% | **0 cases** ✅ |
| 1k reads | 73.7% | 4.2% | 22.1% | **0 cases** ✅ |
| 10k reads | 74.2% | 4.3% | 21.5% | **0 cases** ✅ |

**Verification**:
- ✅ Zero integer overflow values in CIGAR strings
- ✅ Proper CIGAR merging: `150M`, `2S111M398N37M` (not `10M4M...`)
- ✅ All coordinates within chromosome boundaries
- ✅ Mean read length: 150bp (not 5.8 billion)
- ✅ Junction lengths reasonable: 237-482kb (not near 2³²)
- ✅ 170/170 unit tests passing
- ✅ Negative gap warnings logged (expected for overlapping seeds)

**Known Issue**: Tests fail due to many spurious non-canonical junctions (separate alignment quality issue, not overflow-related). This is a separate problem from the integer overflow bug and will be addressed in future optimization work.

---

### Phase 13.5: Scoring Fix ✅ COMPLETE (2026-02-09)

**Problem**: Previous attempt broke mapping (83% → 8%) due to wrong defaults and unnecessary restrictions.

**Root Causes & Fixes**:
1. **`outFilterIntronMotifs`** defaulted to `RemoveNoncanonical` — changed to `None` (STAR default)
2. **`seedMultimapNmax`** changed to 100 — restored to 10000 (STAR default)
3. **Multi-chr anchor check** removed — STAR processes each anchor position independently
4. **Seed overlap trimming** — STAR-style: advance seed B's start by overlap amount (was: skip entirely)
5. **Gap mismatch scoring** — Added `shared_score = shared_bases - 2*mismatches` to DP transitions
6. **`alignSJstitchMismatchNmax` filter** — Reject junction stitches with too many mismatches per motif
7. **Noisy logging** — Changed `log::warn` → `log::trace` for negative gap warnings

**Files Modified**:
- `src/params.rs` — Default corrections
- `src/align/stitch.rs` — Overlap trimming, gap mismatch scoring, logging (~40 lines changed)

**Test Results**:
- ✅ 170/170 unit tests passing
- Stats (100 reads): 78% unique, 9% multi, 13% unmapped
- Stats (1k reads): 71.3% unique, 10.1% multi, 18.6% unmapped
- Performance: ~1.7s for 1000 reads

---

### Phase 13.6: Alignment Extension (extendAlign) ✅ COMPLETE (2026-02-09)

**Problem**: ruSTAR produced 65% reads with soft clips vs STAR's 26%, and only 42% position agreement with STAR. After DP seed stitching, if the seed chain doesn't reach read position 0 or the read end, ruSTAR unconditionally adds soft clips. STAR instead calls `extendAlign()` to extend the alignment into flanking regions, tolerating mismatches up to a limit.

**Implementation**:

1. **`ExtendResult` struct + `extend_alignment()` function** (`src/align/stitch.rs`, ~100 lines)
   - Mirrors STAR's `extendAlign.cpp` algorithm
   - Walks base-by-base from alignment boundary, scoring +1 match / -1 mismatch
   - Records maximum-score extension point
   - Stops when total mismatches exceed `min(p_mm_max * total_len, n_mm_max)`
   - Handles N bases (skipped, no score impact), chromosome boundary (padding=5), reverse strand

2. **`AlignmentScorer` extended** (`src/align/score.rs`, 4 lines)
   - Added `n_mm_max: u32` (from `outFilterMismatchNmax`, default 10)
   - Added `p_mm_max: f64` (from `outFilterMismatchNoverLmax`, default 0.3)
   - Updated `from_params()` and all 8 test constructors + 1 in `src/lib.rs`

3. **`stitch_seeds()` soft-clip replacement** (`src/align/stitch.rs`, ~80 lines)
   - Compute right-side genome position by walking DP CIGAR
   - Left extension: `extend_alignment(..., direction=-1, max_extend=alignment_start)`
   - Right extension: `extend_alignment(..., direction=+1)` with cumulative mismatch count
   - Build final CIGAR: remaining soft-clip + left Match + main CIGAR + right Match + remaining soft-clip
   - CIGAR merging: extension Match merged with adjacent seed Match ops
   - Adjust genome start position by left extension length
   - Adjust transcript score by extension scores

**Files Modified**:
- `src/align/score.rs` — Added `n_mm_max`, `p_mm_max` fields; updated `from_params()`; updated 8 test constructors
- `src/align/stitch.rs` — Added `ExtendResult` + `extend_alignment()`; modified `stitch_seeds()` soft-clip section; added 8 unit tests
- `src/lib.rs` — Updated 1 manual AlignmentScorer constructor

**Test Results** (8 new unit tests):
- Perfect match extension (full extend, 0 mismatches)
- Stops at optimal point with mismatches
- Chromosome boundary (genome base 5)
- N bases skipped
- Leftward direction
- All-mismatch returns extend_len=0
- Recovery through mismatch (3M 1X 10M → extends past mismatch)
- Zero max_extend returns immediately

**10k-read STAR Comparison**:

| Metric | Before | After | STAR |
|--------|--------|-------|------|
| Soft clip rate | 65% | **26.6%** | 25.8% |
| Unique mapped | 74.2% | **79.0%** | 82.6% |
| Position agreement | 42% | **51.0%** | — |
| CIGAR agree (of pos-agrees) | — | 94.2% | — |
| Shared junctions | — | 29/80 | 80 total |
| ruSTAR-only junctions | — | 495 | — |
| Spliced rate | — | 5.0% | 2.5% |

**Remaining Issues Diagnosed**:
1. **Reverse-strand splice motif** — 375 non-canonical junctions from wrong motif detection
2. **2x splice rate** — ruSTAR creates spliced alignments where STAR chooses unspliced
3. **Position disagreements** — Multi-mapping tie-breaking differences (different chr, both MAPQ=255)
4. **629 STAR-only mapped reads** — Seed finding or scoring coverage gaps

**Verified**: 178/178 tests, `cargo clippy` clean (pre-existing only), `cargo fmt --check` pass

---

## Phase 13.7: Reverse-Strand Splice Motif Fix ✅

**Status**: Complete

**Goal**: Fix reverse-strand splice motif detection to eliminate spurious non-canonical junctions.

**Changes**:
- Added `CtAc`, `CtGc`, `GtAt` variants to `SpliceMotif` enum (reverse-complement motifs)
- Modified `detect_splice_motif()` to always read forward genome and check all 6 motif patterns
- Updated `score_splice_junction()` and `stitch_mismatch_allowed()` to handle new variants
- 100% motif agreement on 31/31 shared junctions (was 89.7% on 26/29)
- 343 non-canonical ruSTAR-only junctions (was 375)
- 179/179 tests passing

**Files Modified**: `src/align/score.rs`

---

## Phase 13.8: Splice Junction Overhang Minimum ✅

**Status**: Complete

**Goal**: Enforce `alignSJoverhangMin` (default 5) during DP seed stitching to reduce false spliced alignments.

**Problem**: ruSTAR produced 5.6% spliced reads vs STAR's 2.5%. The `alignSJoverhangMin` and `alignSJDBoverhangMin` parameters were defined in `params.rs` but never enforced during DP stitching — only used in two-pass junction filtering. This allowed junctions flanked by tiny exons (e.g., 3bp) that STAR would reject.

**Changes**:
- Added `align_sj_overhang_min` and `align_sjdb_overhang_min` fields to `AlignmentScorer`
- Wired from `from_params()` and updated all test constructors (~10 across score.rs, stitch.rs, lib.rs)
- Added overhang check in DP stitching loop: for splice junctions, reject if left overhang (prev seed length) or right overhang (current seed effective length) is below `align_sj_overhang_min`
- Added 2 unit tests for overhang acceptance/rejection

**Files Modified**: `src/align/score.rs`, `src/align/stitch.rs`, `src/lib.rs`

---

## Phase 13.8b: Enforce alignIntronMax ✅

**Status**: Complete

**Goal**: Enforce `alignIntronMax` to eliminate false splice junctions with enormous intron spans (100kb-780kb). When `alignIntronMax=0` (default), STAR computes a maximum of `2^winBinNbits * winAnchorDistNbins = 65536 * 9 = 589,824bp`. Gaps exceeding this are treated as deletions (heavily penalized), not splice junctions.

**Problem**: `alignIntronMax` was defined in `params.rs` but never enforced anywhere. ruSTAR accepted arbitrarily large introns, producing 527 ruSTAR-only junctions including spans of 779kb, 703kb, 328kb, and 190kb.

**Changes**:
- Added `align_intron_max` field to `AlignmentScorer`, wired from `from_params()` with default resolution (0 → 589,824)
- Enforced upper-bound check in both splice junction branches of `score_gap()`: gaps exceeding `align_intron_max` fall through to deletion scoring (`score_del_open + score_del_base * len`), which is extremely negative for large gaps
- Added `alignIntronMax` filtering to `filter_novel_junctions()` in two-pass mode
- Updated all test `AlignmentScorer` constructors (~12 across score.rs, stitch.rs, lib.rs)
- Added 4 unit tests: default resolution, custom value, boundary classification, exceeding-max classification

**Files Modified**: `src/align/score.rs`, `src/align/stitch.rs`, `src/lib.rs`, `src/junction/mod.rs`

---

## Phase 13.9: Multi-Mapping Tie-Breaking (Planned)

**Status**: Not started

**Goal**: Align multi-mapping tie-breaking strategy with STAR.
- Investigate STAR's strategy for choosing among equally-scoring loci
- Most position disagreements are different-chromosome with MAPQ=255
- May involve deterministic ordering (by genomic position?) or scoring tiebreakers

---

## Phase 14: STARsolo (Single-Cell)

**Status**: Not started
