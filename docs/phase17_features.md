[← Back to ROADMAP](../ROADMAP.md)

# Phase 17: Features + Polish

**Status**: In Progress (17.1, 17.5, 17.8, 17.A, 17.B, 17.C, 17.D, 17.2, 17.4, 17.6 complete)

**Goal**: Production-ready features and quality-of-life improvements.

## Sub-phase Status

| Sub-phase | Description | Status |
|-----------|-------------|--------|
| 17.1 | Log.final.out statistics file (MultiQC/RNA-SeQC) | ✅ Complete |
| 17.A | `scoreSeedBest` pre-extension on WA entries (STAR faithful) | ✅ Complete |
| 17.B | Per-mate seeding (fix `.18919121`, `.6302610` arch failures) | ✅ Complete — `.18919121` fixed; regressions under investigation |
| 17.C | STAR-faithful SCORE-GATE + mappedFilter for PE (fix 4 MAPQ inflations) | ✅ Complete |
| 17.D | PE combined-span penalty + dedup-before-score-range ordering (248→236 half-mapped) | ✅ Complete |
| 17.2 | Coordinate-sorted BAM (`--outSAMtype BAM SortedByCoordinate`) | ✅ Complete |
| 17.3 | Paired-end chimeric detection | Planned |
| 17.4 | `--outReadsUnmapped Fastx` | ✅ Complete |
| 17.5 | Fix clippy warnings (0 warnings) | ✅ Complete |
| 17.6 | `--outStd SAM/BAM` (stdout output for piping) | ✅ Complete |
| 17.7 | GTF tag parameters (`sjdbGTFchrPrefix`, etc.) | Planned |
| 17.8 | `--quantMode GeneCounts` | ✅ Complete |
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

## Phase 17.5: Clippy Cleanup ✅

**Problem**: 13 clippy warnings creating noise during debugging.

**Changes**:
- Removed dead code: `verify_match_at_position()`, unused `read_seq` param from `cluster_seeds()`
- **`cluster_seeds()`**: 9 args → 3 args — now takes `&Parameters` instead of 6 windowing params
- **`search_direction_sparse()`**: 8 args → 7 args — folded `effective_start_lmax` into body
- **`ChimericSegment::new()`**: Removed constructor, all 6 call sites use struct literal syntax
- Added `AlignReadResult` type alias for `align_read()` return type
- Idiomatic fixes: `.contains()`, `.div_ceil()`, `.saturating_sub()`
- `#[allow(clippy::too_many_arguments)]` on 4 functions with genuinely many distinct args

**Result**: 0 clippy warnings, 264/264 tests passing.

---

## Phase 17.8: `--quantMode GeneCounts` ✅ (2026-04-17)

**Goal**: Output `ReadsPerGene.out.tab` matching STAR's HTSeq-union gene-level counting.

**Implementation**: New `src/quant/mod.rs` with:
- `GeneAnnotation`: per-chromosome sorted interval list (absolute genome coords) built from GTF exons
- `GeneCounts`: atomic per-gene counters + 3 independent N_noFeature/N_ambiguous arrays
- `QuantContext`: `Arc`-shared bundle for rayon parallel threads
- `--quantMode GeneCounts` + `--sjdbGTFfile` validation in `params.rs`
- SE and PE counting paths in `lib.rs`

**Three bugs fixed vs initial implementation**:
1. **Coordinate mismatch**: GTF exon positions were stored chr-relative; `Transcript.exon.genome_start` uses absolute concatenated-genome coords. Fix: add `genome.chr_start[chr_idx]` offset when converting GTF positions.
2. **Single counting pass**: All 3 columns were identical. STAR runs 3 INDEPENDENT passes — col1 (any strand), col2 (same strand as read), col3 (opposite strand) — each with separate N_noFeature and N_ambiguous.
3. **Too-many-loci bucket**: These were going to N_multimapping. STAR puts them in N_unmapped.

**Results vs STAR (10k SE yeast)**:

| Metric | STAR | ruSTAR |
|--------|------|--------|
| N_unmapped | 1073 | 1074 (+1) |
| N_multimapping | 661 | 661 |
| N_noFeature col1/col2/col3 | 131/3653/4240 | 131/3653/4240 |
| N_ambiguous col1 | 567 | 566 (-1) |
| Gene total col1 | 7568 | 7568 |
| Col1 gene disagreements | — | **0** |
| Col2/col3 gene disagreements | — | 1 each (boundary edge case) |

The ±1 discrepancies (N_unmapped + N_ambiguous) are a single read at a gene overlap boundary — likely a minor coordinate boundary difference.

**Files**: `src/quant/mod.rs` (new), `src/params.rs`, `src/junction/mod.rs` (pub(crate) gtf), `src/lib.rs`

**Tests**: 274/274 (added 6 new quant unit tests), 0 clippy warnings.

---

## Phase 17.A: scoreSeedBest Pre-Extension ✅ (2026-04-16)

**Goal**: Match STAR's `ReadAlign_stitchWindowSeeds.cpp` — pre-extend each seed left+right before the recursive DP and store the result as `pre_ext_score` on each `WindowAlignment` entry.

**What STAR does**: Before `stitchWindowAligns`, STAR computes `scoreSeedBest[iS]` for every seed in the window via a two-level DP: (1) base case: `length + left_ext`, (2) chain case: `stitchAlignToTranscript(iS2→iS1) + scoreSeedBest[iS2]`. Then adds `right_ext` universally. Used for seed ordering in the recursive aligner (start from highest-scoring seed).

**Implementation**:

1. **`src/align/stitch.rs`** — `WindowAlignment` struct: added `pub pre_ext_score: i32` field. All construction sites updated (`pre_ext_score: length as i32` default).

2. **`src/align/score.rs`** — `AlignmentScorer`: added `pub out_filter_score_min_over_lread: f64`. All constructor paths updated.

3. **`src/chimeric/detect.rs`** — `WindowAlignment` construction updated.

4. **`src/align/stitch.rs`** — `stitch_seeds_core`: inserted pre-extension block after seed dedup/sort, before `stitch_recurse`:
   - EXTEND_ORDER respected: left-first for forward clusters (`!stitch_is_reverse`), right-first for reverse clusters (matching `stitch_recurse` base case)
   - `right_len_prev = wa.length + first_ext.extend_len` (mirrors base case's `len_after_first`)
   - Chain DP: `dp[i] = max(dp[i], dp[j] + wa_entries[i].pre_ext_score)` with colinearity check
   - No hard pre-filter gate: STAR uses `scoreSeedBest` for ordering only, not window rejection

**Key finding during implementation**: A pre-filter gate at full `outFilterScoreMinOverLread * (Lread-1)` threshold caused 42 false rejections — reads with only short seeds (9-16bp) in low-quality windows, where the full WT extension (starting from leftmost seed) can reach the threshold even though no individual seed's pre-extension does. STAR does NOT apply this gate; `scoreSeedBest` is used for seed ordering in `stitchWindowAligns` only.

**Result**: 268/268 tests, 0 warnings, 8796/8926 SE (baseline maintained), 8390/8390 PE (baseline maintained). `pre_ext_score` ready for Phase 17.B seed ordering.

---

## Phase 17.B: Per-Mate Seeding ✅ (2026-04-17)

**What this fixes**: `.18919121` (was STAR-only) — adapter-RC at start of rc_read1 caused a 15bp Nstart shift in the combined read's mate1 seed position, triggering reverse-cluster rejection. Per-mate seeding finds mate1 seeds from `mate1_seq` directly, avoiding the adapter-RC contamination.

**Root cause of original failures** (combined-read approach):
- `.18919121`: Nstart positions 21, 63, 106 in the 301bp combined-read fell within `rc_mate2` = RC(adapter-contaminated mate2). The adapter RC at stitch_read[155:171] caused a 15bp seed shift for mate1, firing the reverse-cluster reject condition.
- `.6302610`: In the forward cluster, rc_mate2 seeds at sa_pos=126596 (inside mate1's genome range) slipped through `fwd_reject` because the combined read blurred the mate boundary.

**Implementation** (`src/align/read_align.rs`):
1. **Per-mate seeding**: `Seed::find_seeds(mate1_seq, ...)` and `Seed::find_seeds(mate2_seq, ...)` separately. Each mate seeded with its own Nstart positions (0, 37, 74, 112 for 150bp reads).
2. **Independent clustering**: `cluster_seeds()` called separately for each mate.
3. **Independent stitching**: `stitch_seeds_with_jdb_debug()` per mate-cluster. Reverse clusters receive `mate2_seq` directly; stitch internally does RC and sets `is_reverse=true`.
4. **Pairwise matching**: `try_pair_transcripts()` — checks same chr, opposite strands, within `win_bin_window_dist()` span, combined score gate.
5. **Half-mapped fallback**: if no valid pair but one mate individually passes quality threshold, report as HalfMapped.

**Removed from stitch.rs**: `stitch_seeds_working`, `find_mate_boundary`, `split_working_transcript`, `adjust_mate2_coords`, `adjust_wt_read_coords` — no longer needed.

**Result**: `.18919121` now mapped as VIII:452300 15S134M1S + VIII:452301 133M17S (STAR: 16S133M1S + 133M17S). 1bp CIGAR difference is a seed-level tie.

**Regressions from per-mate approach (known, to fix later)**:
- **15 rDNA inter-copy junction reads missed**: Reads spanning the boundary between two adjacent rDNA repeat units (yeast chr XII, ~9.1kb inserts). STAR's combined-read boundary seed at position ~171 uniquely identifies the inter-copy junction. Per-mate seeding generates 55 mate1 × 9 mate2 = ~76 candidate pairs, hitting the TooManyLoci limit (>20). Root fix: apply position-dedup before TooManyLoci check (STAR's actual ordering), or implement targeted cross-boundary rescue.
- **~366 extra both-mapped pairs**: Cross-copy pairings created by combining mate1 and mate2 transcripts from different repeat copies. These inflate NH counts for some multi-mappers.
- **248 half-mapped pairs**: New behavior — reads where one mate individually maps but cannot pair. STAR doesn't output these by default (--outSAMunmapped None).

**Test status**: 274/274, 0 clippy warnings, SE 8796/8926 maintained.

---

## Phase 17.C: STAR-faithful SCORE-GATE + mappedFilter ✅ (2026-04-17)

**Problem**: 4 PE MAPQ inflations for rDNA/repeat multi-mappers. ruSTAR NH=2 vs STAR NH=3 for reads with cross-rDNA-copy pairs (M1@copy1 + M2@copy2, 9037bp gap), causing MAPQ=3 vs STAR's MAPQ=1.

**Root cause**: Two distinct bugs:

1. **Per-WT absolute threshold too strict** (`read_align.rs` forward/reverse cluster processing):
   - ruSTAR used `if adjusted_score < combined_score_threshold { continue; }` (hard cutoff at `outFilterScoreMinOverLread * (Lread-1)`)
   - STAR's `stitchWindowAligns.cpp:324` SCORE-GATE uses a RELATIVE criterion: `Score + outFilterMultimapScoreRange >= wTr[0]->maxScore` (within `scoreRange=1` of window best)
   - For cross-copy pairs: same-copy score=198 (g_span=100bp, penalty=-2), cross-copy score=197 (g_span=9237bp, penalty=-3). ruSTAR rejected cross-copy (197 < 198); STAR accepted it (197+1 ≥ 198)

2. **filter_paired_transcripts applied absolute threshold per-pair** (not just to best):
   - ruSTAR checked every pair's `combined_wt_score < absolute_threshold` → removed cross-copy (197 < 198)
   - STAR's `ReadAlign_mappedFilter.cpp` checks only `trBest->maxScore >= threshold` — if the best passes, ALL pairs in the score window are kept

**Fix**:

1. **`src/align/read_align.rs`** — both forward and reverse cluster processing (lines 750, 972):
   ```rust
   // Old:
   if adjusted_score < combined_score_threshold { continue; }
   // New:
   if adjusted_score + params.out_filter_multimap_score_range < combined_score_threshold { continue; }
   ```

2. **`src/align/read_align.rs`** — `filter_paired_transcripts` (line 1373):
   - Changed from per-pair retain to best-pair quality check
   - Find best pair (max `combined_wt_score`); if best fails any threshold → clear all (read unmapped)
   - If best passes → keep all pairs (they already passed multMapSelect relative criterion)

**Verification**: STAR debug trace on `.19790508` confirmed Score=197 cross-copy pair is INSERTED (`TR-INSERTED`) with `global_pass=1` because `scoreRange=1` (`outFilterMultimapScoreRange`). STAR's `mappedFilter` only checks `trBest->maxScore=198 >= 198` — passes.

**Result**: 268/268 tests, 0 warnings, 8796/8926 SE (maintained), 8390/8390 PE (maintained), **0 MAPQ inflations** (was 4), **0 MAPQ deflations**, faithfulness 98.915% (was 98.903%).

---

## Phase 17.D: PE Combined-Span Penalty + Dedup Ordering ✅ (2026-04-17)

**Problem 1**: `try_pair_transcripts` used `combined_wt_score = t1.score + t2.score`, double-applying the genomic-length penalty (each mate's finalized score already includes its own span penalty). This inflated the combined score relative to STAR's single-span formula, causing AS tag disagreements (was 99.6% of PE reads).

**Fix**: STAR computes ONE genomic-length penalty over the full PE span. Per-mate approach must undo per-mate penalties and apply the combined penalty:
```rust
combined_wt_score = t1.score + t2.score - p1 - p2 + combined_p
```
where `p1 = genomic_length_penalty(t1_span)`, `p2 = genomic_length_penalty(t2_span)`, `combined_p = genomic_length_penalty(right.genome_end - left.genome_start)`.
`try_pair_transcripts` now takes `scorer: &AlignmentScorer` to call `scorer.genomic_length_penalty()`.

**Problem 2**: Decision tree ordering was score-range → dedup → TooManyLoci (wrong). STAR's ordering is multMapSelect → dedup → TooManyLoci. Also, dedup was running after score-range filter, meaning some duplicate pairs at the same position could escape the score-range window.

**Fix** (`src/align/read_align.rs`): Reordered to: (1) position dedup, (2) score-range filter (multMapSelect), (3) TooManyLoci check, (4) sort by score, (5) quality filter (mappedFilter).

**Result**: 278/278 tests, 0 warnings. **248 → 236 half-mapped** (12 pairs fixed by correct ordering). PE diff-AS dropped from 99.6% → 3.1%. PE both-mapped: 8767 (STAR: 8390). SE 8796/8926 maintained.

**Investigation note**: Attempted quality-filter fallback (retry with pre-score-range pool when score-range winner fails quality) to recover additional half-mapped reads. The specific root cause of remaining 236 half-mapped: STAR's combined-read DP finds correct mate2 alignment with 0 mismatches; ruSTAR's per-mate DP finds a different alignment at the same position with 8 mismatches (combined_nm=14 > outFilterMismatchNmax=10). The fallback recovers 35 pairs but introduces ~100 position regressions (pairs at wrong positions passing individual quality checks). Reverted.

---

## Phase 17.2: Coordinate-Sorted BAM — Planned

High user value. STAR outputs `--outSAMtype BAM SortedByCoordinate` natively. Options:
1. In-memory sort during write (requires buffering all records)
2. Wrapper calling `samtools sort` (simpler, already documented as workaround)

---

## Phase 14: STARsolo (Single-Cell) — DEFERRED

**Prerequisite**: All accuracy gaps resolved, position agreement >99%. (Current: 99.92% parity excluding unavoidable ties)
