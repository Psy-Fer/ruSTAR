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
                                                              └→ Phase 16.1-16.8 (algorithm parity) ✅
                                                                   └→ Phase 16.9 (PE joint DP) ← Next
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
| [16](docs/phase16_algorithm.md) | Algorithm Parity | ✅* | 264 | 96.2% pos, 2.7% splice, mate rescue, MAPQ fix |
| [17](docs/phase17_features.md) | Features + Polish | ✅* | 264 | Log.final.out, sorted BAM planned |
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

See [docs/phase16_algorithm.md](docs/phase16_algorithm.md) for detailed sub-phase notes (16.1-16.9).

**Summary**: Bin-based windowing, pre-DP seed extension (96.2% pos, 2.7% splice), MAPQ lookup table, jR scanning, mate rescue (0 unmapped pairs), sparse seed infrastructure (dormant).

---

## Phase 17: Features + Polish ✅ (partial)

See [docs/phase17_features.md](docs/phase17_features.md) for sub-phase table and 17.1 details.

**Summary**: Log.final.out complete. Sorted BAM, PE chimeric, quantMode planned.

---

## Phase 14: STARsolo (Single-Cell) — DEFERRED

Waiting for accuracy parity (position agreement >99%).
