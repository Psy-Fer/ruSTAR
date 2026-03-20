# ruSTAR Implementation Roadmap

Tracks implementation progress across sessions. Each phase lists its deliverables, files touched, and completion status. Detailed notes for later phases are in `docs/`.

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
                                                    └→ Phase 13.1-13.14 (perf+accuracy) ✅
                                                         └→ Phase 15.1-15.6 (SAM tags) ✅
                                                              └→ Phase 16.1-16.10+16.11b (algorithm parity) ✅
                                                                   └→ Phase 16.PE1-PE3 (recursive stitcher, PE joint DP, PE arch refactor) ✅
                                                                        └→ Phase 16.14 (Nstart fix, 99.5% pos) ✅
                                                                             └→ Phase 16.26-16.29 (SA range fix, rev-strand fix, extendAlign fix, STITCH-SJ fix) ✅
                                                              └→ Phase 17.1 (Log.final.out) ✅
                                                                   └→ Phase 17.2+ (features + polish)
                                                              └→ Phase 14 (STARsolo) [DEFERRED]
```

**Phase ordering rationale**: Threading (Phase 9) done first to establish parallel architecture.
Paired-end (Phase 8) builds on threaded infrastructure. GTF/junctions (Phase 7) done after core parallelism.

---

## Phase Summary Table

| Phase | Description | Status | Tests | Key Result |
|-------|-------------|--------|-------|------------|
| 1 | CLI + Parameters | ✅ | 9 | clap derive, ~52 STAR params, validation |
| 2 | FASTA Loading + Genome | ✅ | 19 | 1 byte/base encoding, padding, RC |
| 3 | Suffix Array + SAindex | ✅ | — | PackedArray, SA construction, k-mer lookup |
| 4 | Index Loading + Seeds | ✅ | 39 | MMP search, binary search SA |
| 5 | Seed Stitching + Scoring | ✅ | 57 | DP clustering, CIGAR, transcript filtering |
| 6 | SAM Output (E2E) | ✅ | 84 | FASTQ reader, SAM writer, MAPQ, first pipeline |
| 9 | Threading | ✅ | — | Rayon parallel, `--runThreadN` |
| 8 | Paired-End | ✅ | — | Independent mate align + pairing |
| 7 | GTF/Splice Junctions | ✅ | 132 | GTF parsing, junction DB, SJ.out.tab |
| 10 | BAM Output | ✅ | 136 | BGZF streaming, `--outSAMtype BAM Unsorted` |
| 11 | Two-Pass Mode | ✅ | 138 | Novel junction discovery, pass1→pass2 |
| 12 | Chimeric Detection | ✅ | 170 | SE chimeric, Chimeric.out.junction |
| [13](docs/phase13_accuracy.md) | Performance + Accuracy | ✅ | 205 | 94.5% pos, 97.8% CIGAR, 2.1% splice |
| [15](docs/phase15_sam_tags.md) | SAM Tags + PE Fix | ✅ | 235 | NH/HI/AS/NM/nM/XS/jM/jI/MD, PE fix |
| [16](docs/phase16_algorithm.md) | Algorithm Parity | ✅* | 268 | SE: 99.7% pos, 2.2% splice, 26 actionable, 0 STAR-only; PE: 8295/8390 (−95 ruSTAR under), 32 FP (down from 144, Phase 16.31 overlap check fix), 127 STAR-only |
| [17](docs/phase17_features.md) | Features + Polish | ✅* | 268 | Log.final.out, clippy cleanup, sorted BAM planned |
| 14 | STARsolo | DEFERRED | — | Waiting for accuracy parity |

*Partially complete — see linked docs for sub-phase status.

---

## Phase 1: CLI + Parameters ✅

- `src/params.rs` — `Parameters` struct with ~52 STAR CLI params via clap derive
- `src/error.rs` — `Error` enum with thiserror
- `src/lib.rs` — `run()` dispatcher, `src/main.rs` — thin entry
- Multi-value params need explicit `num_args`; negative defaults need `allow_hyphen_values = true`

---

## Phase 2: FASTA Loading + Packed Genome ✅

- `src/genome/mod.rs` — Genome struct, padding, RC, file writing
- `src/genome/fasta.rs` — FASTA parser, base encoding (A=0, C=1, G=2, T=3, N=4)
- 1 byte per base (not 2-bit), padding=5, RC in second half of buffer

---

## Phase 3: Suffix Array Generation ✅

- `src/index/suffix_array.rs` — SA construction, Rayon parallel sort
- `src/index/sa_index.rs` — Pre-computed prefix lookup (35-bit entries)
- `src/index/packed_array.rs` — Variable-width bit packing

---

## Phase 4: Index Loading + Seed Finding ✅

- `src/index/io.rs` — Load Genome, SA, SAindex from disk
- `src/align/seed.rs` — MMP search, binary search SA, seed expansion

---

## Phase 5: Seed Stitching + Alignment Scoring ✅

- `src/align/stitch.rs` — Seed clustering (100kb window), DP stitching
- `src/align/score.rs` — Match/mismatch/gap scoring, splice motif penalties
- `src/align/transcript.rs` — Transcript struct, CIGAR ops, score tracking
- `src/align/read_align.rs` — Per-read alignment driver

---

## Phase 6: SAM Output (First End-to-End) ✅

- `src/io/fastq.rs` — FASTQ reader (plain + gzip, noodles)
- `src/io/sam.rs` — SAM writer (header + records)
- `src/mapq.rs` — MAPQ calculation
- `src/stats.rs` — Alignment statistics
- `src/lib.rs` — Full pipeline: load index → read FASTQ → align → write SAM

---

## Phase 9: Threading ✅

- Rayon parallel iterators, 10,000 reads per batch
- Sequential FASTQ reading, parallel alignment, sequential SAM writing
- `Arc<AlignmentStats>` with atomic counters

---

## Phase 8: Paired-End Reads ✅

- Independent SE alignment per mate then pairing by chr + distance
- SAM FLAGS (0x1, 0x2, 0x8, 0x20, 0x40, 0x80), TLEN, RNEXT/PNEXT
- Original seed-pooling approach (0% mapped) replaced by PE alignment fix

---

## Phase 7: GTF/Splice Junction Annotation ✅

- `src/junction/mod.rs` — `SpliceJunctionDb`, junction lookup
- `src/junction/gtf.rs` — GTF parser, exon→intron conversion
- `src/junction/sj_output.rs` — SJ.out.tab writer, motif encoding
- Junction coords: 1-based intronic bases (STAR convention)

---

## Phase 10: BAM Output ✅

- `src/io/bam.rs` — Streaming unsorted BAM with BGZF compression
- `AlignmentWriter` trait for SAM/BAM polymorphism
- Compatible with `samtools sort` and `samtools index`

---

## Phase 11: Two-Pass Mode ✅

- Pass 1: Discover junctions (NullWriter discards alignments)
- Filtering: novel junctions require ≥1 unique OR ≥2 multi reads, overhang ≥ 5bp
- Pass 2: Re-align ALL reads with merged GTF + novel junction DB
- Output: SJ.pass1.out.tab (pass 1) + SJ.out.tab + SAM/BAM (pass 2)

---

## Phase 12: Chimeric Detection ✅

- `src/chimeric/detect.rs` — Tier 1 (soft-clip) + Tier 2 (multi-cluster)
- `src/chimeric/score.rs` — Junction type classification, repeat length
- `src/chimeric/output.rs` — 14-column Chimeric.out.junction format
- Detects inter-chr fusions, strand breaks, large-distance breaks
- PE chimeric detection not yet implemented (Phase 17.3)

---

## Phase 13: Performance + Accuracy ✅

See [docs/phase13_accuracy.md](docs/phase13_accuracy.md) for detailed sub-phase notes (13.1-13.14).

**Summary**: From 42% to 94.5% position agreement through SA position encoding fix, CIGAR reversal, splice motif fix, extendAlign, bidirectional seeding, BySJout filtering, and scoring fixes.

---

## Phase 15: SAM Tags + Output Correctness ✅

See [docs/phase15_sam_tags.md](docs/phase15_sam_tags.md) for detailed sub-phase notes (15.1-15.6 + PE fix).

**Summary**: NH/HI/AS/NM/nM/XS/jM/jI/MD tags, SECONDARY flag, outSAMmultNmax, outSAMattributes enforcement, PE FLAG/PNEXT fixes, independent mate alignment.

---

## Phase 16: Algorithm Parity ✅ (partial)

See [docs/phase16_algorithm.md](docs/phase16_algorithm.md) for sub-phase notes (16.1-16.13), [docs/phase16_14_nstart_fix.md](docs/phase16_14_nstart_fix.md) for the Nstart fix.

**Summary**: Bin-based windowing, pre-DP seed extension, MMP SA range narrowing, multi-transcript DP, recursive combinatorial stitcher, STAR-faithful scoring (scoreStitchSJshift removed), sparse bidirectional seed search with Nstart +1 fix, WALrec persistent threshold, post-jR shared base scoring, hierarchical SAindex lookup, nWA reset + overlap detection, coverage filter removal, Lread-1 filter fix, too-many-loci filter, mate rescue, SA range narrowing fix (find_mult_range + max_mappable_length), reverse-strand stitcher coordinate fix (RC read + forward genome coords), PE joint DP stitching via combined-read path, STAR-faithful PE architecture (no cross-product), combined-read score threshold fix (pre-split check prevents double-counting), extendAlign EXTEND_ORDER fix (5' of read first; reverse-strand reads extend right before left) + float comparison fix.

**SE parity (10k yeast, post Phase 16.29):**

| Category | Count | % | Fixable? |
|----------|-------|---|----------|
| Exact match (chr + pos + CIGAR) | 8799 | 98.57% | — |
| Splice match (chr + pos + introns, CIGAR differs) | 1 | 0.01% | — |
| **Total match** | **8800** | **98.57%** | — |
| Unavoidable ties (repeat copy tiebreaking, same score) | 126 | 1.41% | No |
| Fixable algorithm differences | 26 | 0.29% | Yes |
| **Parity excl. unavoidable ties** | **8800/8826** | **99.70%** | — |

**Adjusted SE summary (post Phase 16.29)**: 99.7% position agreement, 99.9% CIGAR, 2.2% splice rate (= STAR), 99.9% MAPQ, 26 actionable disagreements, 1 STAR-only / 1 ruSTAR-only. MAPQ inflation: 4 reads, MAPQ deflation: 4 reads.

**PE parity (10k yeast pairs, 150 bp, post palindromic-early-exit removal + D5):**

| Metric | ruSTAR | STAR |
|--------|--------|------|
| Both-mapped pairs | 8412 | 8390 |
| Half-mapped pairs | 0 | 0 |
| Net gap | +22 (ruSTAR over) | — |
| Per-mate position agreement | ~97.8% | — |
| Per-mate CIGAR agreement | ~95.9% | — |
| ruSTAR-only false positives | 144 | — |
| STAR-only missed | 122 | — |

**PE implementation path:**
- Phase 16.PE1: Recursive combinatorial stitcher (`stitch_recurse`) replacing forward DP
- Phase 16.PE2: STAR-faithful combined-read approach `[mate1_fwd][SPACER][RC(mate2)]`; `mate_id` tagging; `split_working_transcript`; overlap consistency checks
- Phase 16.PE3: Removed non-STAR independent SE + cross-product path; decision tree: joint pairs only → TooShort/unmapped
- Score threshold fix: `split_working_transcript` copies `wt.score` to both halves; checking `wt1.score+wt2.score` doubled the threshold. Fixed to check `wt.score < combined_score_threshold` before split.
- Phase 16.28: extendAlign EXTEND_ORDER fix — for reverse-strand alignments, extend right (5' of read) before left.
- Phase 16.30: PE overlap check fix — forward cluster check 1 uses post-extension estimate for mate1 genome start.
- D5 (2026-03-12): Overlapping-mate junction consistency check (`pe_junctions_consistent`). Removed 2 false-positive pairs where splice junctions in the overlapping genome region disagreed between mates.
- Palindromic early-exit removed (2026-03-16): `if mate1_seq == rc_mate2 { return TooShort }` was wrong — STAR maps palindromic zero-insert RF pairs. Recovered ~28 pairs.

**PE false positives — root cause (confirmed 2026-03-19 via STAR debug tracing):**

The 144 ruSTAR-only false positives are rejected by STAR's `outFilterMatchNminOverLread` check in `mappedFilter`, applied against the **combined read length** (`len1 + SPACER + len2`). Overlapping non-palindromic pairs only have ~130-150 matched bases (one mate's worth); the combined threshold = `0.66 × 310 = 205` → `unmapped=short`. ruSTAR's score threshold is applied correctly but the `nMatch` check against combined Lread is missing.

**Fix**: After `split_working_transcript`, add: `if t1.n_match + t2.n_match < (0.66 * (combined_lread - 1) as f64) as i32 { reject }` in `read_align.rs`.

**Remaining fixable SE issues (deferred):**

| Issue | Count | Difficulty |
|-------|-------|------------|
| Wrong intron choice (same chr, different large intron) | 4 | High |
| MAPQ inflation (missed splice/indel secondary) | 4 | Medium |
| MAPQ deflation (extra unspliced secondary) | 4 | Medium |
| ruSTAR false splice (adapter contamination, 279 kb) | 1 | Medium |
| STAR-only mapped (high-mismatch read NM=10) | 1 | Unknown |

**Remaining PE gap:** 144 ruSTAR-only false positives (fix: combined nMatch threshold), 122 STAR-only missed (root cause TBD). See `docs/star_comparison/DIFFERENCES.md`.

---

## Debugging Tools

**STAR debug tracing** (added 2026-03-19): Instrumented STAR binary at `/home/jamfer/Dropbox/Bioinformatics/tools/repos/STAR/source/STAR` with read-name-filtered trace points.

Usage:
```bash
export STAR_DEBUG_READS="ERR12389696.12345,ERR12389696.67890"
STAR ... 2>star_debug.log
```

Helper script: `test/debug_star.sh`
```bash
./debug_star.sh pe <rustar.sam> <star.sam> [n_reads]  # extract & trace false positives
./debug_star.sh reads "read1,read2"                    # trace specific reads
```

Instrumented locations (all gated on read name match, no performance impact on non-target reads):
- `stitchWindowAligns.cpp`: FINALIZE, PE-MATE-CHECK, PE-OVERLAP, PE-CHECK1/2, PE-REJECT, STITCH-maxScoreMate, SCORE-GATE, TR-INSERTED
- `ReadAlign_multMapSelect.cpp`: MULTMAPSELECT, per-transcript pass/fail, MULTMAPSELECT-RESULT
- `ReadAlign_mappedFilter.cpp`: MAPPEDFILTER (all thresholds + result)
- `ReadAlign_stitchWindowSeeds.cpp`: SEEDSTITCH-maxScoreMate

---

## Phase 17: Features + Polish ✅ (partial)

See [docs/phase17_features.md](docs/phase17_features.md) for sub-phase table and 17.1 details.

**Summary**: Log.final.out complete. Clippy cleanup (0 warnings). Sorted BAM, PE chimeric, quantMode planned.

---

## Phase 14: STARsolo (Single-Cell) — DEFERRED

Waiting for accuracy parity (position agreement >99%).
