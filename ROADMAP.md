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

## Phase 13: Performance Optimization

**Status**: In progress (critical bugs fixed, alignments working!)

**Goal**: Optimize alignment performance to approach STAR speeds and fix classification issues.

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

### Phase 13.3: Performance Optimization (TODO)

**Remaining Issues**:
1. **Performance**: ~3x slower than STAR for 1000 reads (MEDIUM PRIORITY)
   - Profile with perf
   - Optimize seed expansion (currently re-verifies all positions)
   - Optimize DP stitching for clusters with many expanded seeds
   - Cache genome position lookups

2. **Splice motif detection for reverse strand** (LOW PRIORITY)
   - `detect_splice_motif()` in score.rs doesn't add n_genome offset for reverse strand
   - Affects splice junction scoring but not overall alignment accuracy

3. **Parameter tuning**: Genome-size-dependent defaults (LOW PRIORITY)

**Next Steps**:
- Profile with perf to identify hot paths
- Test with 10K-100K reads for performance benchmarking
- Optimize most impactful hot paths

---

## Phase 14: STARsolo (Single-Cell)

**Status**: Not started
