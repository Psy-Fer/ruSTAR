[← Back to ROADMAP](../ROADMAP.md)

# Phase 16: Accuracy + Algorithm Parity

**Status**: In Progress (Phases 16.1-16.10 + 16.7b/16.7c complete)

**Goal**: Close remaining accuracy gaps vs STAR. Fix over-splicing, rDNA MAPQ, seed parameters, DP junction optimization, PE mate rescue, and multi-transcript DP.

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

## Phase 16.7: Bin-Based Windowing ✅ (2026-02-18)

Replaced proximity-based clustering with STAR's bin-based windowing architecture:
1. `cluster_seeds()` rewritten: identify anchors → create windows from anchor bins → merge nearby windows (±winAnchorDistNbins) → extend by ±winFlankNbins → assign seeds by bin HashMap lookup
2. `WindowAlignment` struct with verified (seed_idx, sa_pos, strand, read_pos, length) — no SA range re-expansion
3. `stitch_seeds_with_jdb()` converts WA entries directly to ExpandedSeeds
4. Added `--winFlankNbins` param (default 4)

**Performance**: 12× faster (10s vs 120s on 10k reads) — bin HashMap lookup vs O(n²) proximity checks.

**Files**: `src/align/stitch.rs`, `src/align/read_align.rs`, `src/params.rs`

---

## Phase 16.7b: Pre-DP Seed Extension ✅ (2026-02-20)

**Problem**: Bin-based windowing (16.7) caused splice rate regression: 2.1% → 3.6%. Wide ~589KB windows allow coincidental short seed matches to create false splice junctions. STAR prevents this via seed extension during DP initialization.

**STAR's approach** (ReadAlign_stitchWindowSeeds.cpp):
1. Pass 1: DP base score = `seed_length + left_extend_score` — true positions extend well (+40bp), coincidental extend ~0bp
2. Pass 2: Right extension for chain endpoint selection — `dp_score + right_ext_score`
3. Post-DP: authoritative extensions replace approximate pre-DP scores

**Implementation** (all in `stitch_seeds_with_jdb()`):
1. **Pre-DP left extension**: `left_ext_scores: Vec<i32>` computed via `extend_alignment()` for every expanded seed before DP loop
2. **DP base case**: `score = seed_length + left_ext_scores[i]` (was `seed_length` only)
3. **Right extension for chain selection**: Replaced `max_by_key(score)` with loop computing rightward extension per chain endpoint, selecting by `dp[i].score + right_ext_score`
4. **Post-DP score adjustment**: `adjusted_score = best_state.score - left_ext_scores[chain_start_idx] + left_extend + right_extend`

| Metric | Pre-16.7 | 16.7 (WA-DP) | 16.7b (pre-DP ext) | STAR |
|--------|----------|--------------|---------------------|------|
| Position agree | 94.5% | 94.2% | **96.2%** | — |
| CIGAR agree | 97.8% | 97.1% | **97.6%** | — |
| Splice rate | 2.1% | 3.6% | **2.7%** | 2.2% |
| Shared junctions | 42 | 48 | **49** | 72 |

**Files**: `src/align/stitch.rs`

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

## Phase 16.7c: Anchor Threshold Fix ✅ (2026-02-20)

**Problem**: `max_loci_for_anchor` was hardcoded to 10. STAR uses `winAnchorMultimapNmax` (default 50).

**Fix**: Changed `max_loci_for_anchor` from hardcoded 10 to `params.win_anchor_multimap_nmax` in `align_read()` and `rescue_unmapped_mate()`.

**Impact**: Position **96.2% → 97.4%**. Shared junctions **49 → 57**. Splice rate 2.7% → 3.0% (slight increase).

**Files**: `src/align/read_align.rs`

---

## Phase 16.9: MMP SA Range Narrowing ✅ (2026-02-20)

**Problem**: `extend_match()` only narrowed `sa_start`, overestimating match lengths for positions at the end of the SA range.

**Fix**: Replaced `extend_match()` with `max_mappable_length()` — binary searches within SA range to narrow both boundaries. Ports STAR's `maxMappableLength` + `compareSeqToGenome` + `findMultRange`.

New functions:
- `compare_seq_to_genome()` — starts from offset l_start, returns (match_len, is_read_greater)
- `max_mappable_length()` — binary search within SA range
- `find_mult_range()` — narrow to exact boundary
- `median_uint2()` — safe median for binary search

Removed anchor fallback in `cluster_seeds()` (no longer needed with accurate SA ranges).

**Impact**: Position **97.4% → 97.9%**. CIGAR **97.6% → 98.3%**. Splice rate **3.0% → 2.2%** (matches STAR!). Shared junctions **57 → 62**. False-positive junctions **3 → 0**.

**Files**: `src/align/seed.rs`, `src/align/stitch.rs`

---

## Phase 16.10: Multi-Transcript DP (MAPQ Fix) ✅ (2026-02-20)

**Problem**: ruSTAR produced exactly **one transcript per window** from DP stitching. For rDNA reads (chrXII tandem repeats), all ~6 copies landed in a single window. STAR found multiple transcripts → NH=6 → MAPQ=0. ruSTAR found one → NH=1 → MAPQ=255. 323 reads had inflated MAPQ.

**Fix**: Multi-endpoint DP — collect top-N DP endpoints instead of single best:

1. **`build_transcript_from_endpoint()`** — extracted transcript-building logic (chain traceback → junction optimization → extend → CIGAR → exons → Transcript) into helper function
2. **Top-N endpoint collection** — collect all `(total_score, endpoint_idx)` pairs, sort by score descending
3. **Chain-start dedup** — different endpoints in the same chain produce the same alignment, so skip duplicates via `seen_chain_starts` HashSet
4. **Score-range early termination** — `score < best_score - 1` breaks early (worse endpoints won't survive `outFilterMultimapScoreRange` filtering)
5. **`max_transcripts_per_window`** parameter — `alignTranscriptsPerWindowNmax` (default 100) controls the limit
6. **Chimeric detection** — `stitch_seeds()` convenience wrapper passes `max_transcripts_per_window=1` (chimeric detection doesn't need multi-transcript)

| Metric | Before (16.9) | After (16.10) | Improvement |
|--------|---------------|---------------|-------------|
| MAPQ inflation (ruSTAR=255, STAR<255) | **323** | **62** | -81% |
| MAPQ agreement | 96.1% | **99.1%** | +3.0pp |
| Position agreement | 97.9% | **97.4%** | see note |
| CIGAR agree (of pos-agree) | 98.3% | **98.5%** | +0.2pp |
| Multi-mapped reads | 355 | **627** | STAR: 661 |
| Splice rate | 2.2% | **1.9%** | 122 false splices fixed |
| Disagreements | 190 | **173** | -17 |
| Shared junctions | 62 | **62** | stable |

**Position agreement note**: The 0.5pp drop is expected — diff-chr disagreements are all MAPQ-tied multi-mappers where both tools pick different locations among equally-valid alternatives. All 100 diff-chr cases show matching MAPQ.

**Splice rate improvement**: 122 reads lost false splices (100kb+ spurious introns from wrong DP endpoint). Their CIGARs now match STAR exactly.

**rDNA MAPQ**: Was 304 reads all MAPQ=255, now 26 at MAPQ=1 + 220 at MAPQ=3 + 65 at MAPQ=255 — closely matches STAR's distribution (29 at 1, 225 at 3, 60 at 255).

**PE impact** (multi-transcript DP improves both SE and PE since each mate uses `align_read()`):

| PE Metric | Before (16.9) | After (16.10) | Improvement |
|-----------|---------------|---------------|-------------|
| Per-mate position agree | 95.6% | **97.8%** | +2.2pp |
| Per-mate CIGAR agree | 93.1% | **96.0%** | +2.9pp |
| Both mates mapped | 8706 (96.6%) | **8761 (97.1%)** | +55 pairs |
| Half-mapped pairs | 311 (3.4%) | **263 (2.9%)** | -48 pairs |
| Shared junctions | 76 | **85** | +9 |
| ruSTAR-only junctions | — | **3** | — |

**Files**: `src/align/stitch.rs`, `src/align/read_align.rs`

---

## Phase 16.11b: Fix extend_alignment() to Match STAR's extendAlign() ✅ (2026-02-20)

**Problem**: ruSTAR's `extend_alignment()` passed actual cumulative alignment length (e.g., seed length ~20) as `len_prev`, activating the proportional mismatch check (`pMMmax * total_length`). STAR always passes `Lprev=100000`, making only the absolute `nMMmax` limit apply.

**Example**: For a seed of length 20 extending leftward with `pMMmax=0.3`, `nMMmax=10`:
- STAR: allows up to 10 mismatches (proportional check disabled by Lprev=100000)
- ruSTAR: allows only `min(0.3*(20+i+1), 10)` ≈ 6 mismatches at start

**Fix**:
1. **All 4 call sites**: Changed `len_prev` to `100_000` (pre-DP left, endpoint right, post-DP left, post-DP right)
2. **Post-DP left extension**: Changed `n_mm_prev` from `best_state.n_mismatch` to `0` (matches STAR)
3. **Loop body restructured** to match STAR's `extendAlign()` exactly:
   - Record max score only on match (was: on every base)
   - Break only on mismatch, checking limit before incrementing nMM (was: after)
   - Break condition uses full `max_extend` length (was: current position `i+1`)

**Impact**: Neutral on 10k yeast dataset (all metrics unchanged). The fix ensures correctness for edge cases (very short seeds, high mismatch regions, stricter filter settings) and matches STAR's calling convention exactly.

**Files**: `src/align/stitch.rs`

---

## Phase 16.12: Diagnose and Fix Remaining SE Disagreements ✅ (2026-02-20)

**Problem**: 173 SE disagreements remain. Of these, **100 are diff-chr multi-mapper ties** — reads that map equally well to multiple chromosomes (all MAPQ 1 or 3), where the choice of primary alignment depends on implementation-specific processing order (SA iteration, window sequence). These are unfixable without matching STAR's low-level tie-breaking.

**Adjusted position agreement definition**: Excluding diff-chr multi-mapper ties (both mapped, different chromosome, both MAPQ < 255) from the denominator. These are not meaningful disagreements — both tools report valid alignments, just with different tie-breaking among equally-good locations.

**Infrastructure added**:
1. **`--readNameFilter` parameter** — when set, produces detailed alignment trace on stderr for a specific read: seed counts, cluster details, DP expanded seeds + endpoints + scores, transcript filtering decisions, final output
2. **`compare_sam_thorough.py` updated** — new "Adjusted Agreement" section counting diff-chr ties and reporting adjusted metrics
3. **`extract_disagreement_reads.py`** — parses SAM comparison, categorizes disagreements (false_splice, missed_splice, same_chr_close, same_chr_far, star_only, rustar_only, diff_chr_tie), extracts representative reads for targeted debugging

**Adjusted metrics**:
- Raw position agreement: 97.4% (8620/8793 both-mapped reads)
- Diff-chr multi-mapper ties: ~100
- **Adjusted position agreement: 99.2%** (excluding ties)
- Actionable disagreements: 73 same-chr + 33 STAR-only + 26 ruSTAR-only = 132

**Files**: `src/params.rs`, `src/align/read_align.rs`, `src/align/stitch.rs`, `test/compare_sam_thorough.py`, `test/extract_disagreement_reads.py`

---

## Phase 16.11: PE Joint DP Stitching — PLANNED

263 pairs (2.9%) are half-mapped. STAR maps both via joint mate-aware DP with fragment length gap penalty. Seeds from both mates in same genomic window participate in DP together.
