# ruSTAR

A Rust reimplementation of [STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference), the widely-used RNA-seq aligner originally written in C++ by Alexander Dobin.

## Overview

ruSTAR aims to be a faithful port of STAR, matching the original behavior as closely as possible. It uses the same genome index format, accepts the same `--camelCase` command-line parameters, and produces compatible SAM/BAM output.

**Current status**: End-to-end single-end and paired-end RNA-seq alignment with splice junction detection, two-pass mode, chimeric alignment detection, and multi-threaded parallel processing. 264 tests passing.

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

### Single-End (10k yeast reads, 150bp)

#### Alignment Rates

| Metric | ruSTAR | STAR |
|--------|--------|------|
| Unique mapped | 93.0% | 92.6% |
| Multi-mapped | 7.0% | 7.4% |
| Soft-clipped reads | 26.0% | 26.0% |
| Splice rate | 1.9% | 2.2% |

#### Position and CIGAR Agreement

| Mode | Position agree | CIGAR agree | Splice rate |
|------|---------------|-------------|-------------|
| Normal (default) | 97.4% | 98.5% | 1.9% |

#### MAPQ Agreement

| Metric | Value |
|--------|-------|
| MAPQ agreement | 99.1% |
| MAPQ inflation (ruSTAR=255, STAR<255) | 62 reads (was 323 pre-16.10) |

#### Junction Statistics (SE)

| Metric | ruSTAR | STAR |
|--------|--------|------|
| Shared junctions | 62 | 72 total |
| ruSTAR-only junctions | 0 | -- |
| Motif agreement (shared) | 100% | -- |

### Paired-End (10k yeast read pairs, 150bp)

#### Alignment Rates

| Metric | ruSTAR | STAR |
|--------|--------|------|
| Both mates mapped | 8761 (97.1%) | 8390 (100%) |
| Half-mapped pairs | 263 (2.9%) | 0 |
| Unmapped pairs | 0 | 0 |

#### Per-Mate Agreement

| Metric | Value |
|--------|-------|
| Per-mate position agree | 97.8% |
| Per-mate CIGAR agree | 96.0% |

#### Junction Statistics (PE)

| Metric | ruSTAR | STAR |
|--------|--------|------|
| Shared junctions | 85 | 90 total |
| ruSTAR-only junctions | 3 | -- |
| Motif agreement (shared) | 100% | -- |

> **Note**: 263 half-mapped pairs (2.9%) are cases where one mate maps but the other fails even with mate rescue. STAR uses joint DP stitching which recovers these.

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
- Residual MAPQ inflation (~62 reads MAPQ 255 vs STAR <255) — mostly multi-mapper tie-breaking
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
