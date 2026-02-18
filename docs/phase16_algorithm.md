[← Back to ROADMAP](../ROADMAP.md)

# Phase 16: Accuracy + Algorithm Parity

**Status**: In Progress (Phases 16.1-16.8 complete)

**Goal**: Close remaining accuracy gaps vs STAR. Fix over-splicing, rDNA MAPQ, seed parameters, DP junction optimization, and PE mate rescue.

---

## Phase 16.1: max_cluster_dist from winBinNbits ✅ (2026-02-13)

**Problem**: Hardcoded `max_cluster_dist = 100kb`. STAR computes `2^winBinNbits * winAnchorDistNbins` = 589,824bp (6x larger).

**Fix**:
- Added `--winBinNbits` (16) and `--winAnchorDistNbins` (9) params
- `win_bin_window_dist()` helper: `2^winBinNbits * winAnchorDistNbins`
- Replaced hardcoded 100kb and literal 589,824 with computed value

**Impact**: Splice rate **3.4% → 2.2%** (matches STAR exactly). CIGAR agreement 97.3% → 97.8%.

**Files**: `src/params.rs`, `src/align/read_align.rs`, `src/align/score.rs`, `src/junction/mod.rs`

---

## Phase 16.2: RemoveNoncanonicalUnannotated + GTF Testing ✅ (2026-02-13)

**Problem**: `RemoveNoncanonicalUnannotated` fell through to `RemoveNoncanonical`, rejecting ALL non-canonical junctions even when annotated.

**Fix**: `zip(junction_motifs, junction_annotated)` — only reject when `NonCanonical && !annotated`.

**GTF Testing**: Established differential baseline. STAR gains 8 more junctions with GTF because it inserts splice sites into genome index at alignment time (ruSTAR does not).

**Files**: `src/align/read_align.rs`

---

## Phase 16.3: Junction Position Optimization (jR Scanning) ✅ (2026-02-13)

STAR's 3-phase jR scanning algorithm as post-DP optimization:
1. `find_best_junction_position()` — scan left, scan right (match quality + motif), repeat detection
2. `optimize_junction_positions()` — walks winning chain's CIGAR, applies to each RefSkip
3. `jr_shift` clamped to `[-prev_match_len, next_match_len]` to prevent CIGAR corruption

**Architecture**: Post-DP on winning chain's ~1-3 junctions. Zero performance impact.

**Results**: Neutral on yeast (+1 shared junction, 42 total). Expected benefit on mammalian genomes.

**Files**: `src/align/score.rs`, `src/align/stitch.rs`

---

## Phase 16.4: Seed Search Params + Sparse Infrastructure ✅ (2026-02-13)

- Added `seedSearchStartLmax` (50), `seedSearchStartLmaxOverLread` (1.0), `seedSearchLmax` (0), `seedMapMin` (5)
- Refactored `find_seed_at_position()` → `MmpResult` (always provides MMP advance length)
- `search_direction_sparse()` with STAR-matching Lmapped tracking — kept dormant (DP needs dense seeds)

**Files**: `src/params.rs`, `src/align/seed.rs`

---

## Phase 16.5: MAPQ Formula Fix ✅ (2026-02-16)

Replaced formula with STAR's lookup table: n=1→255, n=2→3, n=3-4→1, n≥5→0.

Threaded `n_for_mapq` through pipeline: `align_read()` → `align_paired_read()` → SAM builders. Currently `n_for_mapq = transcripts.len()` (no inflation).

**Files**: `src/mapq.rs`, `src/align/read_align.rs`, `src/io/sam.rs`, `src/io/bam.rs`, `src/lib.rs`

---

## Phase 16.5b: rDNA Window-Model Fix — DEFERRED

~157 chrXII rDNA reads get MAPQ=255 vs STAR 1-3. Root cause: 589kb clusters merge all tandem repeat copy seeds into ONE cluster → ONE transcript → MAPQ=255.

**Attempted**: Anchor-bin counting (3 iterations — raw, competitive, quality-filtered). All failed: 589kb max_cluster_dist means ALL clusters share the same seeds → identical transcripts.

**Correct fix**: Split clusters into ~65kb sub-windows before DP. Requires significant pipeline refactoring. `SeedCluster.anchor_bin` field kept for future use.

---

## Phase 16.6: Sparse Seed Bug Fixes ✅ (2026-02-17)

3 bugs fixed in `search_direction_sparse()`:
1. **While condition**: exit when `pos + seed_map_min >= read_len`
2. **Nstart**: `read_len.div_ceil(effective_start_lmax)` (not +1)
3. **RC read_pos**: convert `read_pos = original_read_len - read_pos - length` for R→L seeds

**Activation tested and reverted**: 91.1% pos (was 94.5%), 4.3% splice (was 2.1%). Root cause: sparse seeds at mismatch positions produce spurious locations without enough neighbors to vote them down. Bug-fixed function kept dormant.

**Files**: `src/align/seed.rs`

---

## Phase 16.7: Sparse Seed Activation — PLANNED

Adapt DP stitcher for sparse seeds (extension-based gap filling), then activate `search_direction_sparse()`. Study STAR's `stitchAlignToTranscript` for:
1. Extension-based gap filling between sparse seeds
2. Window/cluster formation with sparse anchors
3. DP gap scoring when seed gaps are 30-50bp

---

## Phase 16.8: PE Mate Rescue + Half-Mapped Output ✅ (2026-02-18)

**Problem**: 12.9% unmapped pairs because both mates must independently produce alignments.

**Fix**: 3-tier PE alignment with mate rescue:
1. Both map → pair concordantly
2. One fails → `rescue_unmapped_mate()` (genome-wide seeds → filter chr ± `alignMatesGapMax` → cluster → stitch) → pair or HalfMapped
3. Neither maps → unmapped

**`PairedAlignmentResult` enum**: `BothMapped(Box<PairedAlignment>)` | `HalfMapped { mapped_transcript, mate1_is_mapped }`

**Half-mapped output**: mapped mate FLAG 0x8, unmapped mate FLAG 0x4, co-located.

| Metric | With Rescue | Pre-Rescue | STAR |
|--------|-------------|------------|------|
| Both mapped | 8706 (96.6%) | 8714 | 8390 |
| Half-mapped | **311** | 0 (dropped) | 0 |
| Unmapped pairs | **0** | 286 | 0 |
| Per-mate pos agree | **95.6%** | 95.7% | — |
| STAR-only mates | **98** | 184 | — |
| Shared junctions | **76** | 72 | 90 total |

**Files**: `src/align/read_align.rs`, `src/io/sam.rs`, `src/stats.rs`, `src/lib.rs`

---

## Phase 16.9: PE Joint DP Stitching — PLANNED

311 pairs (3.4%) are half-mapped. STAR maps both via joint mate-aware DP with fragment length gap penalty. Seeds from both mates in same genomic window participate in DP together.
