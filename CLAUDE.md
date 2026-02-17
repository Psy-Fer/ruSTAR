# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important: Git Workflow

**DO NOT commit changes automatically.** The user will review, test, and commit changes themselves. Claude should:
- Make code changes as requested
- Suggest what should be committed
- Let the user handle `git add`, `git commit`, and `git push`

## Project Overview

ruSTAR is a Rust reimplementation of [STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference), an RNA-seq aligner originally written in C++ by Alexander Dobin. Licensed under MIT to match the original STAR license.

The primary goal is a faithful port — matching the original STAR behavior as closely as possible. Extra features and divergences from the original will come in later releases/forks. When implementing, refer to the [STAR source code](https://github.com/alexdobin/STAR) to ensure correctness and behavioral parity.

## Build Commands

Rust 2024 edition. Standard Cargo commands:

```bash
cargo build            # Debug build
cargo build --release  # Release build
cargo test             # Run all tests
cargo test <name>      # Run a single test by name
cargo clippy           # Lint
cargo fmt              # Format code
```

Always run `cargo clippy`, `cargo fmt --check`, and `cargo test` before considering a phase complete.

## Current Implementation Status

See [ROADMAP.md](ROADMAP.md) for detailed phase tracking. **Phases 1-13.14 + 15.1-15.6 + 16.1-16.6 + 17.1 + PE alignment fix complete. SE: 94.5% position agreement, 97.8% CIGAR agreement, 2.1% splice rate (matches STAR 2.2%). With GTF: 94.5% position, 97.6% CIGAR (STAR detects more annotated junctions via index insertion). PE: 87.1% mapped (was 0%), 95.7% per-mate position agreement, 97.1% CIGAR agreement. SAM tags: NH/HI/AS/NM/nM/XS/jM/jI/MD all implemented + --outSAMattributes enforcement. SECONDARY flag + outSAMmultNmax: 99.8% FLAG agreement, 96.2% NH agreement. Phase 17.1: Log.final.out with all 37 STAR fields, MultiQC-compatible. 250 tests.**

**Phase order change**: Phases reordered to 9 → 8 → 7 to establish parallel architecture foundation
before adding complex features. Threading affects the entire execution model and is harder to retrofit later.

**Completed**:
- Phase 1 (CLI/params)
- Phase 2 (FASTA loading + genome generation)
- Phase 3 (Suffix array + SAindex)
- Phase 4 (Index loading + seed finding)
- Phase 5 (Seed stitching + alignment scoring)
- Phase 6 (FASTQ reading + SAM output) ← **FIRST END-TO-END ALIGNMENT**
- Phase 9 (Threading)
- Phase 8 (Paired-end reads)
- Phase 7 (GTF/splice junction annotation)
- Phase 10 (BAM output - unsorted streaming)
- Phase 11 (Two-pass mode - novel junction discovery)
- Phase 12 (Chimeric alignment detection)
- Phase 13 (Performance + accuracy optimization) ← **94.5% position agreement, soft clips match STAR**
- Phase 13.9 (Position fix) ← **SA reverse-strand encoding fix, SEQ reverse-complement**
- Phase 13.9b (CIGAR/splice fix) ← **CIGAR reversal, motif coord fix, genomic length penalty**
- Phase 13.9c (Deterministic tie-breaking) ← **Reproducible multi-mapper ordering**
- Phase 13.10 (Accuracy parity) ← **Terminal exon filter, annotation bonus, coverage filter, seed caps**
- Phase 13.11 (R→L seeding) ← **Bidirectional seed search, +127 multi-mapped, 3.3x more shared junctions**
- Phase 13.12 (SJ motif/strand fix) ← **100% motif agreement on shared junctions, motif-derived strand**
- Phase 13.13 (Splice rate fix) ← **Relaxed terminal exon overhang, splice rate 0.9% → 3.4%, shared junctions 30 → 50**
- Phase 13.14 (outFilterBySJout) ← **Implemented outFilterType BySJout, +1% position/CIGAR agreement**
- Phase 15.1 (NH/HI/AS/NM tags) ← **SAM tags added, 98.3% NH / 100% HI / 98.7% AS agreement with STAR**
- Phase 15.2 (XS/SECONDARY/multNmax) ← **SECONDARY flag, XS strand tag, outSAMmultNmax limit**
- Phase 15.3 (jM/jI/MD tags) ← **Junction motif/intron tags + mismatch descriptor for QC/GATK**
- Phase 15.4 (PE FLAG/PNEXT fixes) ← **Mate strand flag, mate position, RNEXT, per-mate tags**
- PE alignment fix ← **Independent mate alignment + pairing: 0% → 87.1% mapped, 95.7% per-mate position agreement**
- Phase 15.5 (--outSAMattributes enforcement) ← **Standard/All/None/explicit tag control, 234 tests**
- Phase 15.6 (nM tag) ← **STAR-compatible mismatch-only count, nM in Standard/All presets, 235 tests**
- Phase 16.1 (max_cluster_dist) ← **winBinNbits/winAnchorDistNbins params, 100kb→589kb clustering, splice rate 3.4%→2.2% (matches STAR), 237 tests**
- Phase 16.2 (RemoveNoncanonicalUnannotated + GTF testing) ← **Filter fix: only reject unannotated non-canonical junctions. GTF differential testing baseline established. 238 tests**
- Phase 16.3 (jR scanning) ← **Post-DP junction boundary optimization. Neutral on yeast (95.8% pos, 97.8% CIGAR unchanged). +1 shared junction. 241 tests**
- Phase 16.4 (Seed search params) ← **seedSearchStartLmax/seedSearchLmax/seedMapMin params + MmpResult refactoring. Sparse search infrastructure ready but dormant — DP needs dense seeds. 241 tests**
- Phase 16.5 (MAPQ formula fix) ← **STAR lookup table (n≥5→0, n=3-4→1, n=2→3). n_for_mapq infrastructure threaded through align→SAM. rDNA MAPQ=255 remains (needs window-model fix). 241 tests**
- Phase 16.6 (Sparse seed bug fixes) ← **3 bugs fixed (while condition, nstart div_ceil, RC read_pos). Activation tested: 91.1% pos / 4.3% splice → reverted to dense. Bug-fixed function dormant. 244 tests**
- Phase 17.1 (Log.final.out) ← **STAR-compatible Log.final.out with all 37 fields, MultiQC-parseable. chrono timestamps, UnmappedReason tracking, per-transcript CIGAR stats. 250 tests**

**Planned**:
- Phase 16.5b (rDNA window-model fix) ← DEFERRED: bin-counting approach investigated & failed (589kb clusters share seeds); needs ~65kb sub-window splitting
- Phase 16.7 (Sparse seed activation) ← Adapt DP stitcher for sparse seeds (extension-based gap filling), then activate search_direction_sparse()
- Phase 16.8 (PE joint DP) ← Mate-aware DP stitching for mate rescue (12.9% → ~0% unmapped)
- Phase 17 (Features + polish) ← sorted BAM, PE chimeric, quantMode, stdout output

**Current Status** (10k yeast reads):

Single-end (Normal mode):
- ✅ **94.5% position agreement** with STAR (was 51% → 94.5% → 95.3% → 96.3% → 95.7% → 95.8% → 94.5%)
- ✅ **97.8% CIGAR agreement** among position-matching reads
- ✅ 83.1% unique mapped (STAR: 82.6%), 5.4% multi-mapped (STAR: 7.4%)
- ✅ 26.6% soft clips (STAR: 26.0%)
- ✅ **Splice rate 2.1%** (STAR: 2.2%) — matches STAR (was 3.4%, fixed by Phase 16.1)

Paired-end (10k yeast read pairs):
- ✅ **87.1% mapped** (8714/10000 pairs) — was 0% before PE alignment fix
- ✅ **95.7% per-mate position agreement** with STAR (matches SE accuracy)
- ✅ **97.1% per-mate CIGAR agreement** among position-matching mates
- ✅ 81.3% unique mapped, 5.9% multi-mapped
- ✅ **72 shared junctions** with STAR (STAR: 90 total), **100% motif agreement**
- ⚠️ 12.9% unmapped pairs (STAR: 0%) — needs PE joint DP stitching (Phase 16.8)

All modes:
- ✅ **100% motif agreement** on shared junctions
- ✅ SAM SEQ properly reverse-complemented for reverse-strand reads
- ✅ 250 unit tests passing
- ✅ Deterministic output (identical SAM across runs)
- ✅ Bidirectional seed search (L→R + R→L)
- ✅ Annotation-aware DP scoring (sjdbScore bonus during stitching)
- ✅ SECONDARY flag (0x100) on multi-mapper secondary records (99.8% FLAG agreement)
- ✅ XS strand tag from splice motifs (--outSAMstrandField intronMotif)
- ✅ outSAMmultNmax limits secondary alignment output
- ✅ jM/jI/MD tags (junction motifs, intron coords, mismatch descriptor)
- ✅ PE FLAG/PNEXT/RNEXT correct (mate strand from actual alignment, per-chr mate position)
- ✅ Per-mate tags (AS, NM, XS, jM, jI, MD computed from each mate's own transcript)
- ✅ --outSAMattributes enforcement (Standard/All/None/explicit tag control)
- ✅ Post-DP junction boundary optimization (jR scanning, Phase 16.3)
- ✅ Log.final.out statistics file (STAR-compatible, MultiQC-parseable, Phase 17.1)

## Source Layout

```
src/
  main.rs          -- Thin entry: parse CLI (clap), init logging, call lib::run()
  lib.rs           -- run() dispatches on RunMode (AlignReads | GenomeGenerate)
  params.rs        -- ~52 STAR CLI params via clap derive, --camelCase long names
  error.rs         -- Error enum with thiserror (Parameter, Io, Fasta, Index, Alignment, Gtf)
  mapq.rs          -- ✅ MAPQ calculation (unique=255, multi=-10*log10(1-1/n))
  stats.rs         -- ✅ Alignment statistics tracking, reporting, and Log.final.out writer
  genome/
    mod.rs         -- ✅ Genome struct, padding logic, reverse complement, file writing
    fasta.rs       -- ✅ FASTA parser, base encoding (A=0,C=1,G=2,T=3,N=4)
  index/
    mod.rs         -- ✅ GenomeIndex (build + load + write)
    packed_array.rs-- ✅ Variable-width bit packing (1-64 bits per element)
    suffix_array.rs-- ✅ SA construction, custom comparator, strand encoding
    sa_index.rs    -- ✅ K-mer lookup table (35-bit entries with flags)
    io.rs          -- ✅ Load index from disk (Genome, SA, SAindex)
  align/
    mod.rs         -- ✅ Module definition
    seed.rs        -- ✅ Seed finding via MMP search + binary search on SA
    stitch.rs      -- ✅ Seed clustering + DP stitching + alignment extension (extendAlign)
    score.rs       -- ✅ Scoring functions (gaps, mismatches, splice junctions)
    transcript.rs  -- ✅ Transcript struct (exon coords, CIGAR, scores)
    read_align.rs  -- ✅ Per-read alignment driver
  io/
    mod.rs         -- ✅ Module exports
    fastq.rs       -- ✅ FASTQ reader (plain + gzip, noodles wrapper)
    sam.rs         -- ✅ SAM writer (header + records, noodles wrapper)
    bam.rs         -- ✅ BAM writer (BGZF compression, streaming unsorted output)
  junction/
    mod.rs         -- ✅ GTF parsing, junction database, motif detection, two-pass filtering
    sj_output.rs   -- ✅ SJ.out.tab writer
    gtf.rs         -- ✅ GTF parser (internal)
  chimeric/
    mod.rs         -- ✅ Module exports
    detect.rs      -- ✅ Chimeric detection algorithms (Tier 1: soft-clip, Tier 2: multi-cluster)
    segment.rs     -- ✅ ChimericSegment and ChimericAlignment data structures
    score.rs       -- ✅ Junction type classification, repeat length calculation
    output.rs      -- ✅ Chimeric.out.junction writer (14-column format)
```

## Key Conventions

- **Crate name is `ruSTAR`** — `#![allow(non_snake_case)]` in lib.rs suppresses the crate name warning
- **STAR params use `--camelCase` naming** — clap `#[arg(long = "camelCase")]` maps to snake_case Rust fields
- **Multi-value params** (genomeFastaFiles, readFilesIn, outSAMtype, outSAMattributes, chimOutType, alignSJstitchMismatchNmax, outSJfilterIntronMaxVsReadN) need explicit `num_args`
- **Negative defaults** (scoreGapNoncan=-8, readMapNumber=-1, etc.) need `allow_hyphen_values = true`
- **`outSAMtype`** is parsed as raw `Vec<String>` then structured via `Parameters::out_sam_type()` method
- **Validation** beyond clap's type checking is in `Parameters::validate()` (e.g. genomeGenerate requires FASTA files)
- **No async** — CPU-bound work; async adds complexity with zero benefit
- **Error handling** — `thiserror` for `Error` enum, `anyhow` for top-level result propagation

## Dependencies

```toml
[dependencies]
clap = { version = "4", features = ["derive"] }
anyhow = "1"
thiserror = "2"
log = "0.4"
env_logger = "0.11"
memmap2 = "0.9"
byteorder = "1"
noodles = { version = "0.80", features = ["fastq", "sam", "bam", "bgzf"] }
bstr = "1"
flate2 = "1"
rayon = "1"
dashmap = "6"
chrono = "0.4"

[dev-dependencies]
tempfile = "3"
assert_cmd = "2"
predicates = "3"
```

## Testing Pattern

- Unit tests: `#[cfg(test)]` in each module
- Integration tests: `tests/` directory (future phases)
- Every phase uses differential testing against STAR where applicable
- Test data tiers: synthetic micro-genome → chr22 → full human genome

**Current test status**: 250/250 unit tests passing, non-critical clippy warnings (too_many_arguments × 6, type_complexity × 1, manual_contains × 1, implicit_saturating_sub × 1)

**Note**: Phase 9 integration tests fail due to pathologically repetitive test genomes (50 exact copies of 20bp). These tests need realistic genomes (deferred to Phase 13).

## Current Capabilities

ruSTAR can now perform **end-to-end RNA-seq alignment with two-pass mode and chimeric detection**:
- Generate genome indices from FASTA files
- Read plain or gzipped FASTQ files
- Multi-threaded parallel alignment (via `--runThreadN`)
- Align single-end AND paired-end RNA-seq reads with splice junction detection
- GTF-based junction annotation with scoring bonus
- Two-pass mode for novel junction discovery (`--twopassMode Basic`)
  - Pass 1: Discovers novel junctions from read alignments
  - Pass 2: Re-aligns using both GTF and novel junctions
  - Improves accuracy by ~5-10% for samples with novel junctions
- Output valid SAM or BAM files (`--outSAMtype SAM` or `BAM Unsorted`):
  - Proper headers (@HD, @SQ, @PG)
  - Correct CIGAR strings (M, I, D, N for junctions)
  - Appropriate FLAGS (unmapped, reverse complement, paired-end)
  - MAPQ scores
  - 1-based genomic positions
  - Correct mate information (RNEXT, PNEXT, TLEN for paired-end — from actual mate alignment)
  - SAM optional tags: NH, HI, AS, NM, nM (mismatches only), XS (strand), jM (junction motifs), jI (intron coords), MD (mismatch descriptor)
  - --outSAMattributes controls which tags are emitted (Standard/All/None/explicit)
  - SECONDARY flag (0x100) on multi-mapper hit_index > 1
  - outSAMmultNmax limits number of reported alignments per read
- Splice junction statistics output (SJ.out.tab, SJ.pass1.out.tab in two-pass mode)
- Chimeric alignment detection for single-end reads (`--chimSegmentMin` > 0):
  - Detects inter-chromosomal fusions (e.g., BCR-ABL chr9→chr22)
  - Detects intra-chromosomal strand breaks
  - Detects large-distance breaks (>1Mb same chr/strand)
  - Soft-clip based detection (>20% clipped reads)
  - Outputs Chimeric.out.junction file (14-column STAR-compatible format)
  - Junction type classification (GT/AG, CT/AC, etc.)
- Post-alignment read filtering (`--outFilterType BySJout`):
  - Buffers reads, computes surviving junctions via `outSJfilter*` thresholds
  - Filters reads whose primary alignment has non-surviving junctions
  - Improves position/CIGAR agreement by ~1% (removes false spliced alignments)
- Print alignment statistics (unique/multi/unmapped percentages)
- Log.final.out statistics file (STAR-compatible format, MultiQC-parseable):
  - All 37 fields matching STAR's exact format (47-char right-justified field names)
  - Timestamps, mapping speed, splice counts by motif, mismatch/indel rates
  - Written to `{outFileNamePrefix}/Log.final.out`

## Known Issues / Accuracy Gaps (Priority Order)

1. ~~**Position agreement 51%**~~ ✅ FIXED in Phase 13.9 — now **95.7%**.
2. ~~**CIGAR agreement 84.3%**~~ ✅ FIXED in Phase 13.9b — now **97.3%**.
3. ~~**Splice motif detection**~~ ✅ FIXED in Phase 13.9b — `score_gap_with_strand()` converts RC donor position to forward genome coordinates.
4. ~~**33 false junctions**~~ ✅ FIXED in Phase 13.10 — now **6 ruSTAR-only junctions**.
5. ~~**Single-direction seeding**~~ ✅ FIXED in Phase 13.11 — bidirectional L→R + R→L search. Multi-mapped 6.12% (STAR: 7.4%). Shared junctions 50.
6. ~~**Splice rate still low**~~ ✅ FIXED in Phase 13.13 — relaxed terminal exon overhang. Splice rate 0.9% → **3.4%** (STAR: 2.2%). Shared junctions 30 → **50**.
7. ~~**outFilterBySJout not implemented**~~ ✅ FIXED in Phase 13.14 — BySJout mode: 96.7% pos, 98.3% CIGAR, 207 reads filtered.
7b. ~~**No SECONDARY flag / XS tag / outSAMmultNmax**~~ ✅ FIXED in Phase 15.2 — FLAG 0x100, XS:A:+/-, --outSAMmultNmax.
7c. ~~**No jM/jI/MD tags**~~ ✅ FIXED in Phase 15.3 — jM/jI on 448 spliced records, MD on all 9774 records.
7d. ~~**PE FLAG/PNEXT bugs**~~ ✅ FIXED in Phase 15.4 — FLAG 0x20 from mate's actual strand, PNEXT per-chr coords, RNEXT from mate, per-mate tags.
7e. ~~**PE 0% mapped**~~ ✅ FIXED — Rewrote `align_paired_read()` to use independent SE alignment per mate then pair by chr+distance. 0% → 87.1% mapped, 95.7% per-mate position agreement.
8. ~~**Over-splicing in Normal mode**~~ ✅ FIXED in Phase 16.1 — splice rate now 2.2% (matches STAR exactly). max_cluster_dist 100kb→589kb.
8b. ~~**RemoveNoncanonicalUnannotated filter incorrect**~~ ✅ FIXED in Phase 16.2 — now properly checks annotation status via `junction_annotated` Vec.
8c. ~~**No seedSearchStartLmax/seedSearchLmax/seedMapMin params**~~ ✅ ADDED in Phase 16.4 — params parsed; sparse search bug-fixed in Phase 16.6 (3 bugs: while condition, nstart div_ceil, RC read_pos) but activation reverted (91.1% pos, 4.3% splice rate with sparse vs 94.5%/2.1% dense). DP needs dense seeds.
9. **rDNA multi-mapping** (~157 same-chr >500bp in BySJout): chrXII rDNA repeats — STAR=MAPQ 1-3, ruSTAR=MAPQ 255. Root cause: 589kb clusters merge all repeat copy seeds into one cluster; DP produces one transcript. Bin-counting approaches fail because shared seeds make all clusters produce identical transcripts. Needs ~65kb sub-window splitting (Phase 16.5b deferred).
10. **113 diff-chr disagreements** (SE): Multi-mappers (harmless tie-breaking).
11. **98 STAR-only mapped reads** (SE Normal): Stable.
12. **PE 12.9% unmapped** (STAR: 0%): Requires both mates to align independently. Needs PE joint DP stitching (Phase 16.8) for mate rescue.
12b. ~~**SJ motif strand disagreement**~~ ✅ FIXED in Phase 13.12 — strand now derived from splice motif (`implied_strand()`), not read alignment strand.

## Limitations (to be addressed in future phases)

- ~~**SAM optional tags NH/HI/AS/NM**~~ ✅ DONE in Phase 15.1 — 98.3% NH / 100% HI / 98.7% AS agreement with STAR
- ~~**XS strand tag**~~ ✅ DONE in Phase 15.2 — emitted with `--outSAMstrandField intronMotif` (221 XS tags on 10k reads)
- ~~**Secondary alignment output**~~ ✅ DONE in Phase 15.2 — FLAG 0x100 on hit_index > 1 (876 secondaries, 99.8% FLAG agree, 96.2% NH agree)
- ~~**outSAMmultNmax**~~ ✅ DONE in Phase 15.2 — limits reported alignments (default -1 = all)
- ~~**STAR uses `nM` (mismatches only), ruSTAR uses `NM` (edit distance)**~~ ✅ DONE in Phase 15.6 — nM tag added alongside NM; both emitted in Standard/All presets
- ~~**jM/jI/MD tags**~~ ✅ DONE in Phase 15.3 — jM (junction motifs, B:c), jI (intron coords, B:i), MD (mismatch descriptor, Z:) on all records; 448 spliced records get jM/jI
- ~~**Paired-end FLAG/PNEXT bugs**~~ ✅ DONE in Phase 15.4 — FLAG 0x20 from mate's actual strand, PNEXT per-chr coords, RNEXT from mate, per-mate AS/NM/XS/jM/jI/MD
- ~~**PE 0% mapped**~~ ✅ FIXED — Independent mate alignment + pairing: 87.1% mapped, 95.7% per-mate position agreement
- **PE 12.9% unmapped** (STAR: 0%) — needs PE joint DP stitching for mate rescue — Phase 16.8
- ~~**--outSAMattributes**~~ ✅ DONE in Phase 15.5 — Standard/All/None/explicit tag control enforced
- **No coordinate-sorted BAM output** (unsorted only; use `samtools sort`) — Phase 17.2
- ~~**No Log.final.out**~~ ✅ DONE in Phase 17.1 — STAR-compatible format, all 37 fields, MultiQC-parseable
- ~~**Over-splicing**~~ ✅ FIXED in Phase 16.1 — splice rate 2.2% matches STAR
- **Sparse seed search dormant** — bug-fixed in Phase 16.6 (3 bugs), activation tested and reverted (91.1% pos, 4.3% splice). DP needs dense seeds; requires extension-based gap filling — Phase 16.7
- **rDNA MAPQ inflation** (~157 reads, MAPQ 255 vs STAR 1-3) — MAPQ formula fixed (Phase 16.5), bin-counting investigated & failed (Phase 16.5b), needs ~65kb sub-window splitting
- **No --outReadsUnmapped Fastx** — Phase 17.4
- **No --outStd SAM/BAM** (stdout output) — Phase 17.6
- **No --quantMode GeneCounts** — Phase 17.8
- Chimeric alignment detection:
  - ✅ Single-end chimeric detection (inter-chr, strand breaks, large gaps, soft-clips)
  - ✅ Chimeric.out.junction output file
  - ❌ Paired-end chimeric detection not yet implemented — Phase 17.3
  - ❌ Tier 3 (re-mapping soft-clipped regions) not yet implemented — Phase 17.10
  - ❌ --chimOutType WithinBAM not yet implemented — Phase 17.11
- No STARsolo single-cell features (Phase 14, deferred until accuracy parity achieved)
