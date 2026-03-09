# ruSTAR

A Rust reimplementation of [STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference), the widely-used RNA-seq aligner originally written in C++ by Alexander Dobin.

## Overview

ruSTAR aims to be a faithful port of STAR, matching the original behavior as closely as possible. It uses the same genome index format, accepts the same `--camelCase` command-line parameters, and produces compatible SAM/BAM output.

**Current status**: End-to-end single-end and paired-end RNA-seq alignment with splice junction detection, two-pass mode, chimeric alignment detection, and multi-threaded parallel processing. 268 tests passing.

## Quick Start

### Build

```bash
cargo build --release
```

### Generate genome index

```bash
target/release/ruSTAR --runMode genomeGenerate \
  --genomeDir /path/to/genome_index \
  --genomeFastaFiles /path/to/genome.fa
```

### Align reads

```bash
target/release/ruSTAR \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --outSAMtype SAM \
  --outSAMstrandField intronMotif \
  --outFileNamePrefix /path/to/output_
```

### Paired-end alignment

```bash
target/release/ruSTAR \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq reads_2.fq \
  --outSAMtype SAM \
  --outFileNamePrefix /path/to/output_
```

### BAM output

```bash
target/release/ruSTAR \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix /path/to/output_
```

### Two-pass mode

```bash
target/release/ruSTAR \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --twopassMode Basic \
  --outFileNamePrefix /path/to/output_
```

## Accuracy Comparison vs STAR

Benchmarked on 10,000 yeast RNA-seq reads (150 bp SE, ERR12389696), compared to STAR 2.7.x with identical parameters and genome index.

### Single-End Alignment Rates

| Metric | ruSTAR | STAR |
|--------|--------|------|
| Unique mapped | 92.6% | 92.6% |
| Multi-mapped | 7.4% | 7.4% |
| Soft-clipped reads | 26.0% | 26.0% |
| Splice rate | 2.2% | 2.2% |
| Shared splice junctions | 67 / 72 STAR junctions | — |
| Motif agreement (shared junctions) | 100% | — |

### Strict Per-Read Comparison (SE)

A read is counted as a **match** only if it aligns to the exact same chromosome, exact same start position, and has identical splice junctions (intron coordinates). Any difference in any of these is a mismatch.

| Result | Count | % |
|--------|-------|---|
| Exact match (chr + pos + CIGAR identical) | 8792 | 98.49% |
| Splice match (chr + pos + introns match, CIGAR differs) | 1 | 0.01% |
| **Total match** | **8793** | **98.50%** |
| Mismatch — unavoidable tie-breaking | 127 | 1.42% |
| Mismatch — fixable algorithm differences | 7 | 0.08% |
| **Parity (excluding unavoidable ties)** | **8920 / 8927** | **99.92%** |

#### Mismatch Classification

| Category | Count | Fixable? |
|----------|-------|----------|
| Diff chromosome, both multi-mapper (repeat copy tie-breaking) | 100 | No — same score, different copy chosen |
| Same chr, identical CIGAR, different position (repeat copy tie-breaking) | 21 | No — same score, different copy chosen |
| Same chr, same score, different intron at repeat locus | 4 | No — same AS/NH, tiebreaking only |
| Same chr, same splice, position differs 3–4 bp (soft-clip boundary) | 2 | No — same AS/NH, different soft-clip split |
| Same chr + pos, different splice junctions | 4 | Partial (2 missed splices + 1 jR scan off-by-4bp + 1 STAR false splice) |
| Same chr, STAR spliced / ruSTAR not (missed splice) | 1 | Yes |
| STAR mapped, ruSTAR unmapped | 1 | Maybe (high-mismatch read, NM=10) |
| ruSTAR mapped, STAR unmapped (false splice) | 1 | Yes (adapter contamination → 279 kb intron) |
| **Total mismatches** | **134** | |

> **Unavoidable ties (127 reads):** Both tools find the same set of equally-scored alignments but choose different ones as primary due to internal processing order. Neither alignment is more correct than the other.

> **STAR false splice (1 read):** ERR12389696.5825571 — STAR creates a 607 kb intron from adapter-contaminated bases, scoring 2 points higher than the correct soft-clipped alignment. ruSTAR correctly soft-clips this read.

### MAPQ Agreement (SE)

| Metric | Value |
|--------|-------|
| MAPQ agreement (position-matched reads) | 99.9% |
| MAPQ inflation (ruSTAR=255, STAR<255) | 5 reads |
| MAPQ deflation (ruSTAR<255, STAR=255) | 2 reads |

### Paired-End (10k yeast read pairs, 150 bp)

| Metric | ruSTAR | STAR |
|--------|--------|------|
| Both mates mapped | 8382 (99.9%) | 8390 (100%) |
| Half-mapped pairs | 0 | 0 |
| Unmapped pairs | 0 | 0 |
| Per-mate position agreement | 98.3% | — |
| Per-mate CIGAR agreement | 97.4% | — |

> **8-pair gap vs STAR**: ruSTAR uses STAR's combined-read PE path (`[mate1_fwd][SPACER][RC(mate2)]`), producing near-identical output. The remaining 8-pair difference stems from scoring edge cases in short-insert overlapping reads.

## Supported Features

- Single-end and paired-end alignment with mate rescue
- SAM and unsorted BAM output (`--outSAMtype SAM` or `BAM Unsorted`)
- Multi-threaded parallel alignment (`--runThreadN`)
- GTF-based junction annotation with scoring bonus (`--sjdbGTFfile`)
- Two-pass mode for novel junction discovery (`--twopassMode Basic`)
- Chimeric alignment detection for single-end reads (`--chimSegmentMin`)
- Post-alignment read filtering (`--outFilterType BySJout`)
- Splice junction output (SJ.out.tab)
- Gzip-compressed FASTQ input (`--readFilesCommand zcat`)
- SAM optional tags: NH, HI, AS, NM, nM, XS, jM, jI, MD
- `--outSAMattributes` control (Standard/All/None/explicit)
- SECONDARY flag (0x100) on multi-mapper alignments
- Configurable output limits (`--outSAMmultNmax`)
- Bidirectional seed search (L-to-R and R-to-L)
- Junction boundary optimization (jR scanning)
- Deterministic output (identical SAM across runs)
- Log.final.out statistics file (STAR-compatible, MultiQC-parseable)

## Known Limitations

- No coordinate-sorted BAM output (use `samtools sort` post-alignment)
- No paired-end chimeric detection
- No `--quantMode GeneCounts`
- No `--outReadsUnmapped Fastx`
- No `--outStd SAM/BAM` (stdout output)
- Residual MAPQ inflation (5 reads in 10k SE benchmark) — missed splice/indel secondary alignments
- No STARsolo single-cell features

See [ROADMAP.md](ROADMAP.md) for detailed implementation tracking.

## Building from Source

Requires Rust 2024 edition (rustc 1.85+).

```bash
cargo build --release    # Release build
cargo test               # Run tests
cargo clippy             # Lint
cargo fmt                # Format
```

## Development

The majority of ruSTAR's code was written by [Claude Code](https://claude.ai/code) (Anthropic's AI coding assistant), with technical direction, architecture decisions, and validation by the project maintainer.

## License

MIT (matching the original STAR license)
