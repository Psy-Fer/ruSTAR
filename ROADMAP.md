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
                                                         └→ Phase 13.7-13.9 (accuracy refinement) ✅
                                                              └→ Phase 13.9b (CIGAR/splice fix) ✅
                                                                   └→ Phase 13.9c (deterministic tie-breaking) ✅
                                                                        └→ Phase 13.10 (accuracy parity) ✅
                                                                             └→ Phase 13.11 (R→L seeding) ✅
                                                                                  └→ Phase 13.12 (SJ motif/strand fix) ✅
                                                                                       └→ Phase 13.13 (splice rate fix) ✅
                                                                                            └→ Phase 13.14 (outFilterBySJout) ✅
                                                                                                 └→ Phase 15.1 (NH/HI/AS/NM tags) ✅
                                                                                                      └→ Phase 15.2+ (XS, jM/jI/MD, PE fixes) ← Next
                                                                                                      └→ Phase 16 (accuracy parity)
                                                                                                           └→ Phase 17 (features + polish)
                                                                                                                └→ Phase 14 (STARsolo) [DEFERRED]
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

## Phase 13.8c: Reduce False Non-Canonical Splice Junctions ✅

**Status**: Complete

**Goal**: Reduce ruSTAR's 2x splice rate (5.2% vs STAR's 2.5%) and eliminate hundreds of false non-canonical junctions by implementing STAR's multi-layer splice junction filtering.

**Problem**: ruSTAR produced 516 ruSTAR-only junctions (349 non-canonical) because:
1. Overhang calculation was hardcoded to `5u32` — all SJ filtering based on overhang was useless
2. Missing `outFilterIntronStrands = RemoveInconsistentStrands` (STAR default) — alignments with conflicting junction strand motifs were accepted
3. Missing `outSJfilter*` parameters — no motif-specific filtering thresholds in SJ.out.tab or two-pass mode

**Changes**:
- **P0**: Fixed overhang calculation in `record_transcript_junctions()` — walks CIGAR to compute `min(left_exon_len, right_exon_len)` per junction instead of hardcoded `5u32`. Added `AlignmentScorer::from_params_minimal()` for lightweight motif detection.
- **P1**: Added `IntronStrandFilter` enum and `--outFilterIntronStrands` param (default `RemoveInconsistentStrands`). Added `SpliceMotif::implied_strand()` method. Rejects alignments where junctions imply both + and - transcript strands.
- **P2**: Added 4 `outSJfilter*` params (`OverhangMin`, `CountUniqueMin`, `CountTotalMin`, `DistToOtherSJmin`) — 4-element Vecs indexed by motif category [noncanon=0, GT/AG=1, GC/AG=2, AT/AC=3]. Added `SpliceMotif::filter_category()` and `filter_category_from_encoded()`. Modified `write_output()` to filter junctions by motif-specific thresholds (annotated junctions bypass all filters). Added distance-to-nearest-neighbor computation for `DistToOtherSJmin`.
- **P3**: Replaced hardcoded thresholds in `filter_novel_junctions()` with motif-specific `outSJfilter*` values (non-canonical now requires 30bp overhang + 3 unique reads vs canonical 12bp + 1 read).

**Files Modified**: `src/lib.rs`, `src/params.rs`, `src/align/score.rs`, `src/align/read_align.rs`, `src/junction/sj_output.rs`, `src/junction/mod.rs`

**Tests**: 192/192 passing (+11 new tests). New tests: `implied_strand`, `filter_category`, `filter_category_from_encoded`, `strand_consistency_filter`, `sj_filter_noncanonical_needs_high_overhang`, `sj_filter_annotated_bypasses_filters`, `filter_novel_junctions_noncanonical_strict`, param defaults for new params.

---

## Phase 13.9: Fix Position Agreement ✅

**Status**: Complete

**Goal**: Fix the 51% position agreement rate (4066/8298 both-mapped reads disagree).

**Root Cause**: SA reverse-strand position encoding bug. The suffix array stores reverse-strand positions as offsets within the RC genome region `[0, n_genome)`, but ruSTAR was using them directly as forward-genome coordinates for chromosome identification and SAM output. ~94% of all disagreements were different-chromosome, caused by this single bug.

**Fixes Applied**:
1. **`sa_pos_to_forward()` method** added to `GenomeIndex` (`src/index/mod.rs`): converts raw SA position to forward genome coordinate via `n_genome - sa_pos - match_length` for reverse strand
2. **`cluster_seeds()`** (`src/align/stitch.rs`): convert SA positions to forward coordinates before `position_to_chr()` for both anchor and non-anchor seeds
3. **`stitch_seeds()`** (`src/align/stitch.rs`): convert positions for chromosome bounds checks; keep raw SA positions for internal DP genome access; convert at final transcript creation
4. **SAM SEQ/QUAL** (`src/io/sam.rs`): reverse-complement SEQ and reverse QUAL for FLAG & 16 reads (both single-end and paired-end paths)
5. **`complement_base()`** (`src/io/fastq.rs`): new public function for encoded base complementation

**Files Modified**: `src/index/mod.rs`, `src/align/stitch.rs`, `src/io/sam.rs`, `src/io/fastq.rs`

**Files Created**: `test/diagnose_disagreements.py` (diagnostic script that revealed root cause)

**10k-read STAR Comparison**:

| Metric | Before (13.8c) | After (13.9) | STAR |
|--------|----------------|--------------|------|
| Position agreement | 51.0% | **94.5%** | — |
| Unique mapped | 79.1% | **83.8%** | 82.6% |
| Multi mapped | — | 6.1% | 7.4% |
| Diff-chr both MAPQ=255 | 3497 | **2** | — |
| STAR-only mapped | 629 | **42** | — |
| Soft clip rate | 26.5% | 26.5% | 25.8% |
| CIGAR agree (of pos-agree) | 94.2% | 84.3% | — |
| Spliced rate | 5.1% | 5.8% | 2.5% |

**Remaining Disagreements** (492 of 8444 both-mapped):
- 337 same-chr >500bp: mostly chrXII rDNA repeats, chrII duplicated regions
- 113 diff-chr: 97 multi-mappers (tie-breaking), only 2 both MAPQ=255
- 42 same-chr <500bp: CIGAR differences (soft-clip vs indel)

**Verified**: 192/192 tests passing, `cargo clippy` clean (pre-existing only), `cargo fmt --check` pass

---

## Phase 13.9b: CIGAR Reversal + Splice Motif Fix + Genomic Length Penalty ✅

**Status**: Complete

**Goal**: Fix CIGAR agreement (84.3% → 96.5%) and reduce splice rate (5.8% → 4.1%) through three fixes targeting reverse-strand alignment issues and missing scoring logic.

**Root Causes Identified**:
1. **CIGAR not reversed for reverse strand** (70.6% of CIGAR disagreements): DP builds CIGARs in RC genome order; SAM requires forward genome order
2. **Splice motif detection at wrong coordinates** (13.9% of disagreements): `score_gap()` received raw SA positions for motif detection; for reverse strand these point to unrelated forward-genome locations
3. **Missing `scoreGenomicLengthLog2scale` penalty**: STAR penalizes long-spanning alignments with `ceil(log2(span) * -0.25 - 0.5)`, preventing huge false introns from outscoring compact alignments

**Fixes Applied**:

1. **CIGAR reversal** (`src/align/stitch.rs`):
   - Added `final_cigar.reverse()` for reverse-strand clusters after `count_mismatches()` (which needs RC genome order) but before `ref_len` computation
   - Fixed 956 reads with exactly-reversed CIGARs (e.g., `124M26S` → `26S124M`)

2. **Strand-aware splice motif detection** (`src/align/score.rs`):
   - New `score_gap_with_strand()` method accepts `is_reverse` and `n_genome` parameters
   - Converts RC donor position to forward: `forward_donor = n_genome - rc_donor - intron_len`
   - `detect_splice_motif()` unchanged — always reads forward genome (handles both GT-AG and CT-AC patterns)
   - Original `score_gap()` preserved as wrapper for backward compatibility (tests)

3. **Genomic length penalty** (`src/params.rs`, `src/align/score.rs`, `src/align/stitch.rs`):
   - Added `--scoreGenomicLengthLog2scale` parameter (STAR default: -0.25)
   - `genomic_length_penalty()` method: `ceil(log2(genomic_span) * scale - 0.5)`
   - Applied to transcript score after DP stitching, floored at 0
   - For 121kb intron: penalty ≈ -4 points; for 150bp match: penalty ≈ -2 points

**Files Modified**:
- `src/align/stitch.rs` — CIGAR reversal, `score_gap_with_strand()` call site, genomic length penalty application
- `src/align/score.rs` — `score_gap_with_strand()`, `genomic_length_penalty()`, `score_genomic_length_log2_scale` field
- `src/params.rs` — `scoreGenomicLengthLog2scale` parameter

**10k-read STAR Comparison**:

| Metric | Before (13.9) | After (13.9b) | STAR |
|--------|---------------|---------------|------|
| Position agreement | 94.5% | **95.3%** | — |
| CIGAR agree (of pos-agree) | 84.3% | **96.5%** | — |
| Unique mapped | 83.8% | **84.2%** | 82.6% |
| Multi mapped | 6.1% | **5.3%** | 7.4% |
| Spliced rate | 5.8% | **4.1%** | 2.5% |
| ruSTAR-only junctions | 40 | **33** | — |
| Shared junctions | 40 | 40 | 80 total |
| Position disagreements | 492 | **416** | — |
| Soft clip rate | 26.5% | **26.2%** | 25.8% |

**Remaining Issues**:
- 33 ruSTAR-only junctions: huge introns (10k-500k), all canonical motifs found by coincidence; caused by missing seed positions at the correct unspliced locus
- 5/40 shared junctions: ruSTAR=CT/AC vs STAR=GT/AG (strand assignment in SJ.out.tab)
- 266 same-chr >500bp apart: rDNA repeats and duplicated regions
- 116 diff-chr: 99 multi-mappers (tie-breaking)

**Verified**: 192/192 tests passing, `cargo clippy` clean (pre-existing only), `cargo fmt --check` pass

---

## Phase 13.9c: Deterministic Multi-Mapper Tie-Breaking ✅

**Status**: Complete

**Goal**: Ensure deterministic ordering of equal-score multi-mapped alignments for reproducible output across runs.

**Problem**: When multiple alignments have identical scores, `transcripts.sort_by(|a, b| b.score.cmp(&a.score))` provides no deterministic ordering for equal-score transcripts. The "primary" alignment depends on internal cluster processing order, making output non-reproducible.

**Fixes Applied**:
1. **Single-end sort** (`src/align/read_align.rs:250`): Added secondary tie-breaking keys: smallest chr_idx → smallest genome_start → forward strand first
2. **Paired-end sort** (`src/align/read_align.rs:392`): Same deterministic tie-breaking for paired alignments

**Files Modified**: `src/align/read_align.rs` (~8 lines changed total)

**10k-read STAR Comparison**: No change in aggregate metrics (position agreement 95.3%, CIGAR agreement 96.5%). The 116 diff-chr disagreements remain — STAR uses a different internal ordering (SA-enumeration-based) rather than the simple chr/pos/strand sort. These are all multi-mapper ambiguity (99/116 are both MAPQ < 255).

**Key Result**: Output is now **deterministic** — running ruSTAR twice on the same input produces identical SAM files (verified with diff).

**Verified**: 192/192 tests passing, `cargo clippy` clean (pre-existing only), `cargo fmt --check` pass

---

## Phase 13.10: Accuracy Parity with STAR ✅

**Status**: Complete

**Goal**: Resolve tractable STAR output disagreements through terminal exon filtering, annotation-aware DP scoring, coverage filtering, and seed/window caps.

### Sub-phases Completed

**Phase 13.10a: Extension Mismatch Boundary Investigation**
- Investigated `>` vs `>=` in `extend_alignment()` mismatch check
- **Result**: STAR uses `>` (not `>=`) — confirmed by unit test failures when changed to `>=`
- No change needed; existing `>` is correct

**Phase 13.10b: Terminal Exon Overhang Enforcement** (highest impact)
- **Problem**: 8-11bp seeds accepted as first/last exons, creating 50-500kb false introns (~125 reads)
- **Fix**: In DP splice junction overhang check, first/last exons in chain require 12bp minimum overhang for novel junctions, 3bp for annotated junctions (`alignSJDBoverhangMin`)
- New `stitch_seeds_with_jdb()` function accepts optional junction DB reference
- Junction annotation lookup converts SA positions to forward genome coordinates for DB queries

**Phase 13.10c: `outSJfilterIntronMaxVsReadN` Filter**
- Added parameter `--outSJfilterIntronMaxVsReadN` (3 values, default `[50000, 100000, 200000]`)
- Junctions with intron length exceeding threshold for their supporting read count are filtered from SJ.out.tab
- 1 read → 50kb max, 2 reads → 100kb max, 3+ reads → 200kb max

**Phase 13.10d: `winReadCoverageRelativeMin` Filter**
- Added parameter `--winReadCoverageRelativeMin` (default 0.5)
- After clustering, computes total seed coverage per cluster (union of read ranges / read length)
- Discards clusters below threshold before stitching, preventing sparse clusters from producing bad alignments

**Phase 13.10e: Annotation Bonus During DP Stitching**
- `sjdbScore` (+2) now applied during DP transitions for annotated junctions
- Combined with 13.10b in same code path — annotation lookup determines both overhang threshold and score bonus
- Annotated junctions get preferential scoring during DP, not just during post-filtering

**Phase 13.10f: Seed/Window Cap Parameters**
- Added `--seedPerReadNmax` (1000): caps total seeds per read in `find_seeds()`
- Added `--seedPerWindowNmax` (50): caps seeds per cluster
- Added `--alignWindowsPerReadNmax` (10000): caps total clusters per read
- Added `--alignTranscriptsPerWindowNmax` (100): declared (enforcement in future refinement)

### Files Modified

| File | Changes |
|------|---------|
| `src/align/stitch.rs` | Terminal exon overhang, annotation bonus, `stitch_seeds_with_jdb()` |
| `src/align/read_align.rs` | Coverage filter, seed/window caps, junction DB passing |
| `src/align/seed.rs` | `seedPerReadNmax` cap |
| `src/params.rs` | 7 new parameters, updated defaults test |
| `src/junction/sj_output.rs` | `outSJfilterIntronMaxVsReadN` filter |

### 10k-read STAR Comparison

| Metric | Before (13.9c) | After (13.10) | STAR |
|--------|----------------|---------------|------|
| Position agreement | 95.3% | **96.3%** | — |
| CIGAR agree (of pos-agree) | 96.5% | **97.4%** | — |
| Unique mapped | 84.2% | **85.1%** | 82.6% |
| Multi mapped | 5.3% | **3.65%** | 7.4% |
| Spliced rate | 4.1% | **0.4%** | 2.5% |
| ruSTAR-only junctions | 33 | **3** | — |
| Shared junctions | 40 | 9 | 80 total |
| Same-chr >500bp disagree | 266 | **179** | — |
| Diff-chr disagree | 116 | **102** | — |
| STAR-only mapped | 60 | **75** | — |
| Soft clip rate | 26.2% | **26.7%** | 25.8% |

### Remaining Issues

1. **Splice rate too low** (0.4% vs STAR 2.5%): 12bp terminal overhang filter is aggressive without GTF. With GTF loaded, annotated junctions use 3bp threshold + sjdbScore bonus → should recover most real junctions.
2. **rDNA multi-mapping** (~179 same-chr >500bp): chrXII rDNA repeats, STAR=MAPQ 1-3, ruSTAR=MAPQ 255. Root cause: single-direction seed search misses tandem repeat copies.
3. **102 diff-chr disagreements**: 100 are multi-mappers (harmless tie-breaking).
4. **75 STAR-only mapped reads**: slight increase from coverage filter.

### Remaining Algorithm Gaps vs STAR

| Category | Gap | Impact | Status |
|----------|-----|--------|--------|
| Seed search | L→R only (STAR: both directions) | rDNA multi-mapping, missed junctions | ✅ Fixed in 13.11 |
| Seed search | No `seedSearchStartLmax` | Missed seeds at optimal positions | Open |
| DP stitching | No junction position optimization (jR scanning) | Suboptimal junction placement | Open |
| Multi-mapper | Single-direction seeding misses repeat copies | Inflated MAPQ in tandem repeats | Improved in 13.11 |
| Output | `outFilterBySJout` not implemented | Minor filtering difference | Open |

**Verified**: 192/192 tests passing, clippy clean (6 pre-existing warnings), `cargo fmt --check` pass

---

## Phase 13.11: Bidirectional R→L Seed Search ✅

**Status**: Complete

**Goal**: Implement STAR's reverse-direction seed search to find seeds that L→R search misses, particularly in tandem repeats (rDNA) and at splice junction boundaries.

**Implementation** (`src/align/seed.rs` only — all production changes in one file):

1. **`Seed.search_rc: bool` field** — Marks seeds found via R→L search
2. **`reverse_complement_read()` helper** — RC read for R→L SA search, uses existing `complement_base`
3. **`genome_positions()` updated** — When `search_rc=true`, converts `(pos, is_rev)` → `(n_genome - pos - len, !is_rev)` with `filter_map` to handle underflows in small genomes
4. **`find_seeds()` updated** — After L→R loop, RC read searched L→R, `read_pos` converted back (`read_len - rc_pos - seed.length`), shared `seedPerReadNmax` cap (L→R takes priority)
5. **`stitch.rs` test update** — 2 Seed literals updated with `search_rc: false`

**Bug Found**: `n_genome - pos - length` underflows for small genomes when SA positions are near genome boundary. Fixed by using `filter_map` with `pos + length <= n_genome` guard.

**Unit Tests Added** (4 new, 196 total):
- `test_reverse_complement_read` — RC correctness including N bases
- `test_rl_seeds_found` — R→L seeds appear in find_seeds() output
- `test_shared_seed_cap` — Combined L→R + R→L respects seedPerReadNmax
- `test_rc_seed_genome_positions` — Converted positions are valid

### 10k-read STAR Comparison

| Metric | Before (13.10) | After (13.11) | STAR |
|--------|----------------|---------------|------|
| Position agreement | 96.3% | **96.3%** | — |
| CIGAR agree (of pos-agree) | 97.4% | **97.8%** | — |
| Unique mapped | 85.1% | **84.0%** | 82.6% |
| Multi mapped | 3.65% | **4.92%** | 7.4% |
| Spliced rate | 0.4% | **0.9%** | 2.2% |
| Shared junctions | 9 | **30** | 72 total |
| ruSTAR-only junctions | 3 | **2** | — |
| STAR-only mapped | 75 | **60** | — |
| Soft clip rate | 26.7% | **26.5%** | 26.0% |

**Key Improvements**:
- Multi-mapped +35% (365→492) — R→L seeds expose tandem repeat copies
- Shared junctions 3.3x (9→30) — R→L seeds at junction boundaries enable spliced alignments
- Spliced rate doubled (0.4%→0.9%) — closer to STAR's 2.2%
- Only 2 false junctions (was 3)

**Remaining Algorithm Gaps vs STAR**:
- No `seedSearchStartLmax` — STAR starts R→L search from specific positions, not all
- DP stitching: No junction position optimization (jR scanning)
- rDNA: R→L improved but still misses some tandem repeat copies → inflated MAPQ

**Verified**: 196/196 tests passing, clippy clean (6 pre-existing warnings), `cargo fmt --check` pass

---

## Phase 13.12: SJ.out.tab Motif/Strand Fix ✅

**Status**: Complete

**Goal**: Fix SJ.out.tab motif/strand assignment to match STAR convention. Strand should derive from splice motif dinucleotides, not from read alignment strand.

**Root Cause**: `record_transcript_junctions()` in `lib.rs` derived junction strand from `transcript.is_reverse` (read alignment strand), but STAR derives it from the splice motif dinucleotides (GT/AG → strand=1, CT/AC → strand=2). Additionally, `encode_motif()` incorrectly transformed the motif code based on strand, when it should map directly from the detected motif.

**Fixes Applied**:

1. **Strand from motif** (`src/lib.rs`):
   - Changed `let strand = if transcript.is_reverse { 2 } else { 1 }` to use `motif.implied_strand()`:
     - `Some('+')` → strand=1, `Some('-')` → strand=2, `None` → strand=0

2. **Direct motif encoding** (`src/junction/sj_output.rs`):
   - Simplified `encode_motif()` to direct mapping without strand parameter
   - GtAg→1, CtAc→2, GcAg→3, CtGc→4, AtAc→5, GtAt→6, NonCanonical→0
   - Updated call site to remove strand argument

3. **Full test coverage** — 7 `encode_motif` tests covering all motif variants (was 4)

**Files Modified**:
- `src/lib.rs` — ~3 lines: motif-derived strand
- `src/junction/sj_output.rs` — ~25 lines: simplified `encode_motif()`, updated call site + tests

### 10k-read STAR Comparison

| Metric | Before (13.11) | After (13.12) | STAR |
|--------|----------------|---------------|------|
| Position agreement | 96.3% | **96.3%** | — |
| CIGAR agree (of pos-agree) | 97.8% | **97.8%** | — |
| Motif agreement | 80% (24/30) | **100% (30/30)** | — |
| Unique mapped | 84.0% | **84.0%** | 82.6% |
| Multi mapped | 4.92% | **4.92%** | 7.4% |
| Spliced rate | 0.9% | **0.9%** | 2.2% |
| Shared junctions | 30 | **30** | 72 total |
| ruSTAR-only junctions | 2 | **2** | — |

**Key Result**: Motif agreement on shared junctions improved from 80% to **100%**. Duplicate junction entries (same coordinates, different strand/motif) eliminated.

**Verified**: 199/199 tests passing, clippy clean (6 pre-existing warnings), `cargo fmt --check` pass

---

## Phase 13.13: Relax Terminal Exon Overhang Filter ✅

**Status**: Complete

**Goal**: Improve splice rate (0.9% → STAR's 2.2%) by removing the 12bp terminal exon overhang floor added in Phase 13.10b, which was too aggressive for novel junctions during DP stitching.

**Root Cause**: Phase 13.10b added `base_min_overhang.max(12)` for novel junctions in DP, but STAR only applies 12bp filter at SJ.out.tab write time (`outSJfilterOverhangMin`), not during DP. STAR uses `alignSJoverhangMin` (5bp) during DP.

**Fix** (`src/align/stitch.rs` — ~20 lines → 3 lines):
- Removed terminal exon overhang logic: first/last seed chain detection, 12bp floor, `is_first_seed_in_chain`, `terminal_min_overhang`, `left_min`, `right_min`
- Simplified to: `if left_overhang < base_min_overhang || right_overhang < base_min_overhang { continue; }`
- `base_min_overhang` correctly computed at 5bp (novel) or 3bp (annotated)

**Why safe**: Existing guards still prevent bad junctions:
1. `outSJfilterOverhangMin` [30,12,12,12] filters SJ.out.tab output (sj_output.rs)
2. Two-pass filtering requires 12bp overhang before entering pass 2 junction DB
3. `alignIntronMax` (589,824bp) — huge gaps scored as deletions
4. `scoreGenomicLengthLog2scale` (-0.25) — penalizes long-spanning alignments
5. `alignSJstitchMismatchNmax` — rejects junctions with mismatches at non-GT/AG
6. `winReadCoverageRelativeMin` (0.5) — discards sparse clusters

### 10k-read STAR Comparison

| Metric | Before (13.12) | After (13.13) | STAR |
|--------|----------------|---------------|------|
| Position agreement | 96.3% | **95.7%** | — |
| CIGAR agree (of pos-agree) | 97.8% | **97.3%** | — |
| Unique mapped | 84.0% | **82.9%** | 82.6% |
| Multi mapped | 4.92% | **6.12%** | 7.4% |
| Spliced rate | 0.9% | **3.4%** | 2.2% |
| Shared junctions | 30 | **50** | 72 total |
| ruSTAR-only junctions | 2 | **6** | — |
| STAR-only junctions | 42 | **22** | — |
| Motif agreement | 100% (30/30) | **100% (50/50)** | — |
| STAR-only mapped | 60 | **57** | — |

**Key Improvements**:
- Splice rate 3.8x (0.9% → 3.4%) — now exceeds STAR's 2.2%
- Shared junctions +67% (30 → 50)
- Unique/multi ratios now near-identical to STAR
- Small position/CIGAR regression (~0.6%) from reads now splicing vs soft-clipping

**Trade-off**: Over-splicing (3.4% vs 2.2%) and 6 false junctions (was 2). May need `outFilterBySJout` in future phase.

**Verified**: 199/199 tests passing, clippy clean (6 pre-existing warnings), `cargo fmt --check` pass

---

## Phase 13.14: Implement `outFilterBySJout` ✅

**Status**: Complete

**Goal**: Implement STAR's `outFilterType BySJout` mode — after all reads are aligned, compute which junctions survive `outSJfilter*` thresholds, then filter reads whose primary alignment has any junction not in the surviving set.

**Implementation** (4 files, +468/-78 lines):

1. **`src/junction/sj_output.rs`** — Added `compute_surviving_junctions()` method
   - Factors out filtering logic from `write_output()` into reusable method returning `HashSet<SjKey>`
   - Refactored `write_output()` to use `compute_surviving_junctions()` internally
   - Made `encode_motif()` `pub(crate)` for use by `lib.rs`
   - 3 new tests

2. **`src/junction/mod.rs`** — Re-exported `SjKey` and `encode_motif` as `pub(crate)`

3. **`src/lib.rs`** — BySJout buffering and filtering
   - `extract_junction_keys()` — walks CIGAR for RefSkip ops, builds SjKey per junction
   - `primary_junction_keys` field added to `AlignmentBatchResults`
   - `align_reads_single_end()` and `align_reads_paired_end()` modified:
     - Normal mode: unchanged behavior
     - BySJout mode: buffers all results → computes surviving junctions → filters reads → writes survivors
   - 0 new tests (covered by existing pipeline tests)

4. **`src/stats.rs`** — Added `undo_mapped_record_bysj()` method
   - Atomically moves a mapped read (unique or multi) to unmapped using CAS loops
   - 3 new tests

**Test Results**: 205/205 tests passing (was 199, +6 new)

### 10k-read STAR Comparison

Normal mode (**unchanged** from 13.13):
| Metric | ruSTAR | STAR |
|--------|--------|------|
| Position agreement | 95.7% | — |
| CIGAR agree | 97.3% | — |
| Splice rate | 3.4% | 2.2% |

BySJout mode (`--outFilterType BySJout`):
| Metric | BySJout | Normal | STAR |
|--------|---------|--------|------|
| Position agreement | **96.7%** | 95.7% | — |
| CIGAR agreement | **98.3%** | 97.3% | — |
| Splice rate | **1.1%** | 3.4% | 2.2% |
| Reads filtered | **207** | 0 | — |
| STAR-only mapped | 259 | 57 | — |

**Key Insight**: BySJout works correctly but is too aggressive without GTF annotations — all junctions are novel so subject to strict `outSJfilter*` thresholds. With GTF, annotated junctions bypass filters so more spliced reads survive. This is faithful STAR behavior.

**Memory Note**: BySJout buffers all SAM records in memory (~5MB for 10k reads). For 100M+ reads, disk-based buffering would be needed (future optimization).

**Verified**: 205/205 tests passing, clippy clean (6 pre-existing warnings), `cargo fmt --check` pass

---

## Phase 15: SAM Tags + Output Correctness

**Status**: In Progress

**Goal**: Add all SAM optional tags required by downstream tools (featureCounts, RSEM, StringTie, GATK, Picard, samtools markdup). Fix paired-end output bugs. Implement `--outSAMattributes` enforcement.

### Phase 15.1: NH, HI, AS, NM Tags (Foundation) ✅ COMPLETE (2026-02-10)

**Problem**: All downstream tools require one or more of NH/HI/AS/NM. Without them, ruSTAR output is unusable in standard pipelines.

**Implementation** (`src/io/sam.rs` only):
1. Added imports for `Tag` and `Value` from noodles `record::data::field` and `record_buf::data::field`
2. Added `compute_edit_distance()` helper — sums `n_mismatch + Ins(n) + Del(n)` from CIGAR (excludes RefSkip/SoftClip)
3. Updated `transcript_to_record()` — removed `_n_alignments`/`_hit_index` underscore prefixes, replaced TODO with 4 tag insertions via `record.data_mut().insert()`
4. Updated `build_paired_mate_record()` — added `n_alignments`/`hit_index` params + same 4 tag insertions
5. Updated `build_paired_records()` — passes `n_alignments`/`hit_index` (1-based from enumerate) to both mate record builders
6. Updated existing test call sites for new `build_paired_mate_record()` signature
7. Added 3 new tests: `test_tags_nh_hi_as_nm`, `test_edit_distance_computation`, `test_transcript_to_record_has_tags`

**Tag Comparison vs STAR** (8484 position-matching reads, 10k yeast):

| Tag | Agreement | Rate | Notes |
|-----|-----------|------|-------|
| NH | 8340/8484 | 98.3% | 144 differ (multi-mapper count from different seeding; 122 off-by-1) |
| HI | 8484/8484 | 100.0% | Perfect agreement |
| AS | 8373/8484 | 98.7% | 2 same-CIGAR diffs (minor), 109 from different CIGARs |
| NM vs nM | 8287/8484 | 97.7% | **0 unexplained**: 40 = indel bases (NM counts, nM doesn't), 157 from different CIGARs |

**Identified Problems**:
1. **STAR uses `nM` (mismatches only), ruSTAR uses `NM` (edit distance)** — different tags with different semantics. Both valid SAM tags. May need to also emit `nM` for STAR-compatible output.
2. **2 same-CIGAR AS disagreements** — both 150M with 0 mismatches but STAR gives slightly lower AS (147/144 vs 148). Likely base-quality or annotation scoring penalty in STAR.
3. **NH disagrees on 1.7%** — 52 STAR finds more hits, 92 ruSTAR finds more. From different seed search strategies.

**Test Results**: 208/208 tests passing (+3 new), clippy clean (pre-existing only), `cargo fmt --check` pass

**Files**: `src/io/sam.rs`

### Phase 15.2: XS Tag + Secondary Alignment Output

**Problem**: XS tag required by StringTie/Cufflinks. `--outSAMstrandField intronMotif` parsed but unused. Multi-mappers only output primary alignment, missing SECONDARY flag.

**Fix**:
- For spliced reads, derive strand from junction motifs via `implied_strand()`. Add XS:A:+/- tag.
- Add `--outSAMmultNmax` parameter. Set `FLAGS |= SECONDARY` for hit_index > 1. Respect limit.

**Files**: `src/io/sam.rs`, `src/params.rs`
**Depends on**: 15.1

### Phase 15.3: jM, jI, MD Tags

**Problem**: jM/jI are STAR-specific junction tags used by QC pipelines. MD is required by GATK and variant calling.

**Fix**:
- Walk CIGAR for RefSkip ops, encode motifs as integers → jM (array of i8), jI (array of i32 pairs)
- New `build_md_string()` function: walk CIGAR, count consecutive matches, emit mismatch bases

**Files**: `src/io/sam.rs`, possibly `src/align/stitch.rs`
**Depends on**: 15.1

### Phase 15.4: Paired-End FLAG/PNEXT Fixes

**Problem**: Two bugs in `src/io/sam.rs`:
- Line 484-486: Mate reverse flag (0x20) always assumes opposite strand (should use actual mate strand)
- Line 527: PNEXT set to same transcript's start (should be mate's position)

**Fix**: Pass both mate transcripts into record builder. Set 0x20 from actual mate strand, PNEXT from actual mate position.

**Files**: `src/io/sam.rs`

### Phase 15.5: --outSAMattributes Enforcement

**Problem**: Parameter parsed but no tags generated until 15.1-15.3, then need to control which tags are emitted.

**Fix**: Parse attribute list into a set. Check the set before emitting each tag.

**Files**: `src/io/sam.rs`, `src/params.rs`
**Depends on**: 15.3

### Phase 15.6: nM Tag (STAR-compatible Mismatch Count)

**Problem**: STAR outputs `nM:i:N` (mismatches only, excluding indels). ruSTAR outputs `NM:i:N` (edit distance = mismatches + indels, SAM spec). Both are valid tags but some downstream tools expect STAR's `nM`. On same-CIGAR reads with indels, the values differ by exactly the indel base count.

**Fix**: Add `nM` tag alongside `NM`. `nM` = `transcript.n_mismatch` directly (no indel sum). Custom tag: `Tag::new(b'n', b'M')`.

**Files**: `src/io/sam.rs`
**Depends on**: 15.1

---

## Phase 16: Accuracy + Algorithm Parity

**Status**: Not Started

**Goal**: Close remaining accuracy gaps vs STAR. Fix over-splicing, rDNA MAPQ, missing seed parameters, and DP junction optimization.

### Phase 16.1: max_cluster_dist from winBinNbits

**Problem**: `src/align/read_align.rs:59` hardcodes 100kb. STAR default is `2^winBinNbits * winAnchorDistNbins` = 589,824bp.

**Fix**: Add `--winBinNbits` (16) and `--winAnchorDistNbins` (9) params. Compute max_cluster_dist from them.

### Phase 16.2: RemoveNoncanonicalUnannotated Filter

**Problem**: `src/align/read_align.rs:217` TODO — falls through to `RemoveNoncanonical`, incorrectly rejecting annotated non-canonical junctions.

**Fix**: Pass junction DB into filter check. Only reject unannotated non-canonical junctions.

### Phase 16.3: Junction Position Optimization (jR Scanning)

**Problem**: Splice rate 3.4% vs STAR 2.2%. 6 false junctions. STAR shifts junction boundaries by ±scoreStitchSJshift bases to prefer canonical motifs.

**Fix**: After detecting a splice junction in DP, try shifting boundary ± `scoreStitchSJshift` bases. Keep the shift that produces the best motif score.

### Phase 16.4: seedSearchStartLmax + seedSearchLmax

**Problem**: Neither parsed. `seedSearchStartLmax` (50) controls R-to-L search start position.

**Fix**: Add parameters. In `find_seeds()`, start R-to-L loop from `read_len - seedSearchStartLmax`.

### Phase 16.5: SAindex Hint Usage (rDNA MAPQ Fix)

**Problem**: ~157 chrXII rDNA reads get MAPQ=255 (ruSTAR) vs MAPQ 1-3 (STAR). SA binary search ignores hint_pos.

**Fix**: Use hint_pos to set initial binary search bounds, finding more repeat copies.

### Phase 16.6: Paired-End Joint DP Stitching

**Problem**: Mates aligned independently then combined. Misses rescue opportunities.

**Fix**: Consider seeds from both mates in DP with fragment length gap penalty.

---

## Phase 17: Features + Polish

**Status**: Not Started

**Goal**: Production-ready features and quality-of-life improvements.

| Sub-phase | Description |
|-----------|-------------|
| 17.1 | Log.final.out statistics file (MultiQC/RNA-SeQC compatibility) |
| 17.2 | Coordinate-sorted BAM output (`--outSAMtype BAM SortedByCoordinate`) |
| 17.3 | Paired-end chimeric detection |
| 17.4 | `--outReadsUnmapped Fastx` |
| 17.5 | Fix clippy warnings (config structs for too_many_arguments) |
| 17.6 | `--outStd SAM/BAM` (stdout output for piping) |
| 17.7 | GTF tag parameters (`sjdbGTFchrPrefix`, `sjdbGTFtagExonParentTranscript/Gene`) |
| 17.8 | `--quantMode GeneCounts` |
| 17.9 | `--outBAMcompression` / `--limitBAMsortRAM` |
| 17.10 | Chimeric Tier 3 (re-map soft-clipped regions) |
| 17.11 | `--chimOutType WithinBAM` (supplementary FLAG 0x800) |
| 17.12 | BySJout memory optimization (disk buffering for 100M+ reads) |
| 17.13 | Phase 9 integration test fixes (realistic test genomes) |

---

## Phase 14: STARsolo (Single-Cell) [DEFERRED]

**Status**: Deferred until accuracy parity achieved

**Prerequisite**: All accuracy gaps resolved, all alignment-affecting parameters implemented, position agreement >99%
