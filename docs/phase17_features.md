[← Back to ROADMAP](../ROADMAP.md)

# Phase 17: Features + Polish

**Status**: In Progress (17.1 complete)

**Goal**: Production-ready features and quality-of-life improvements.

## Sub-phase Status

| Sub-phase | Description | Status |
|-----------|-------------|--------|
| 17.1 | Log.final.out statistics file (MultiQC/RNA-SeQC) | ✅ Complete |
| 17.2 | Coordinate-sorted BAM (`--outSAMtype BAM SortedByCoordinate`) | Planned |
| 17.3 | Paired-end chimeric detection | Planned |
| 17.4 | `--outReadsUnmapped Fastx` | Planned |
| 17.5 | Fix clippy warnings (config structs for too_many_arguments) | Planned |
| 17.6 | `--outStd SAM/BAM` (stdout output for piping) | Planned |
| 17.7 | GTF tag parameters (`sjdbGTFchrPrefix`, etc.) | Planned |
| 17.8 | `--quantMode GeneCounts` | Planned |
| 17.9 | `--outBAMcompression` / `--limitBAMsortRAM` | Planned |
| 17.10 | Chimeric Tier 3 (re-map soft-clipped regions) | Planned |
| 17.11 | `--chimOutType WithinBAM` (supplementary FLAG 0x800) | Planned |
| 17.12 | BySJout memory optimization (disk buffering for 100M+ reads) | Planned |
| 17.13 | Phase 9 integration test fixes (realistic test genomes) | Planned |

---

## Phase 17.1: Log.final.out ✅

**Problem**: No `Log.final.out` — MultiQC and RNA-SeQC can't parse results.

**Implementation**:

1. **`Cargo.toml`** — Added `chrono = "0.4"` for timestamps

2. **`src/stats.rs`** — Major expansion (+673 lines):
   - `UnmappedReason` enum (Other, TooShort, TooManyMismatches)
   - 15 new `AtomicU64` counters: read_bases, mapped_bases/mismatches, ins/del count/bases, splices_by_motif[7], splices_annotated, unmapped breakdown, chimeric_reads
   - `record_transcript_stats(&Transcript)` — walks CIGAR for all counters
   - `write_log_final(path, time_start, time_map_start, time_finish)`

3. **`src/align/read_align.rs`** — `align_read()` returns 4-tuple with `Option<UnmappedReason>`

4. **`src/lib.rs`** — Timing via `chrono::Local::now()`, stats collected before BySJout (matches STAR)

**Format**: 47-char right-justified field names + ` |\t` separator. All 37 STAR fields present.

**Differential Test** (10k SE):

| Field | ruSTAR | STAR |
|-------|--------|------|
| Input reads | 10000 | 10000 |
| Avg read length | 150 | 150 |
| Uniquely mapped % | 83.11% | 82.65% |
| Avg mapped length | 146.99 | 146.99 |
| Mismatch rate | 0.40% | 0.40% |

**Files**: `Cargo.toml`, `src/stats.rs`, `src/align/read_align.rs`, `src/lib.rs`

---

## Phase 14: STARsolo (Single-Cell) — DEFERRED

**Prerequisite**: All accuracy gaps resolved, position agreement >99%. (Current: 97.4% pos, 99.1% MAPQ agree after Phase 16.10 multi-transcript DP)
