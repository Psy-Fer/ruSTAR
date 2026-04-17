//! Transcript-level quantification (`--quantMode TranscriptomeSAM`).
//!
//! Builds a per-transcript exon map from GTF records and projects genome-space
//! alignments onto the transcriptome coordinate system so downstream tools
//! (Salmon, RSEM) can read `Aligned.toTranscriptome.out.bam`.
//!
//! This module mirrors STAR's `Transcriptome` / `alignToTranscript` /
//! `quantAlign` logic (see `source/Transcriptome.cpp` and
//! `source/Transcriptome_quantAlign.cpp` in the upstream repo), with the
//! only substantive divergence that ruSTAR builds the transcript tables on
//! the fly from the input GTF instead of loading persisted
//! `transcriptInfo.tab`/`exonInfo.tab` files.
use std::collections::HashMap;

use crate::align::transcript::{CigarOp, Exon, Transcript};
use crate::error::Error;
use crate::genome::Genome;
use crate::junction::gtf::GtfRecord;

/// Per-transcript exon in absolute genome coordinates (0-based half-open),
/// paired with the cumulative transcript-space length of all preceding exons
/// (STAR's `exLenCum`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrExon {
    /// Absolute genome start (0-based, inclusive).
    pub genome_start: u64,
    /// Absolute genome end (0-based, exclusive).
    pub genome_end: u64,
    /// Cumulative transcript-space length of all prior exons (0 for the first).
    pub ex_len_cum: u32,
}

/// Transcriptome metadata built from GTF exon records.
///
/// Each transcript is a contiguous set of exons on a single chromosome and
/// strand, ordered by ascending genome start. `tr_order` / `tr_starts_sorted`
/// / `tr_end_max_sorted` support STAR's `binarySearch1a` + running-max
/// early-exit in `quantAlign`.
#[derive(Debug, Clone)]
pub struct TranscriptomeIndex {
    /// Transcript IDs in insertion order (index = transcript_idx).
    pub tr_ids: Vec<String>,
    /// Genome chromosome index per transcript.
    pub tr_chr_idx: Vec<usize>,
    /// Strand per transcript: 1 = forward/`+`, 2 = reverse/`-`.
    pub tr_strand: Vec<u8>,
    /// gene_id string per transcript.
    pub tr_gene_id: Vec<String>,
    /// Absolute genome start of the first exon (0-based).
    pub tr_start: Vec<u64>,
    /// Absolute genome end of the last exon (0-based half-open).
    pub tr_end: Vec<u64>,
    /// Exons per transcript (sorted by genome_start, with `ex_len_cum`).
    pub tr_exons: Vec<Vec<TrExon>>,
    /// Total transcript-space length (sum of exon spans).
    pub tr_length: Vec<u32>,
    /// Permutation of `0..n_transcripts` sorted by `(tr_start, tr_end)`.
    pub tr_order: Vec<usize>,
    /// `tr_start[tr_order[i]]` — for `binary_search`.
    pub tr_starts_sorted: Vec<u64>,
    /// Running max of `tr_end` along `tr_order` (STAR's `trEmax`).
    pub tr_end_max_sorted: Vec<u64>,
}

impl TranscriptomeIndex {
    /// Build from already-parsed GTF exon records.
    ///
    /// Records are grouped by `transcript_id`.  Transcripts with inconsistent
    /// strands / chromosomes are skipped with a warning.  Transcripts on
    /// unknown chromosomes are skipped with a warning.
    pub fn from_gtf_exons(exons: &[GtfRecord], genome: &Genome) -> Result<Self, Error> {
        // Group exons by transcript_id, preserving first-seen insertion order.
        let mut order: Vec<String> = Vec::new();
        let mut groups: HashMap<String, Vec<&GtfRecord>> = HashMap::new();
        for rec in exons {
            let tid = match rec.attributes.get("transcript_id") {
                Some(id) => id.clone(),
                None => continue,
            };
            if !groups.contains_key(&tid) {
                order.push(tid.clone());
            }
            groups.entry(tid).or_default().push(rec);
        }

        let mut tr_ids: Vec<String> = Vec::new();
        let mut tr_chr_idx: Vec<usize> = Vec::new();
        let mut tr_strand: Vec<u8> = Vec::new();
        let mut tr_gene_id: Vec<String> = Vec::new();
        let mut tr_start: Vec<u64> = Vec::new();
        let mut tr_end: Vec<u64> = Vec::new();
        let mut tr_exons_vec: Vec<Vec<TrExon>> = Vec::new();
        let mut tr_length: Vec<u32> = Vec::new();

        for tid in &order {
            let mut recs = groups.remove(tid).unwrap();
            if recs.is_empty() {
                continue;
            }

            // Validate consistent chromosome + strand.
            let first = recs[0];
            let chr_name = &first.seqname;
            let strand_char = first.strand;
            let mut inconsistent = false;
            for r in &recs {
                if &r.seqname != chr_name || r.strand != strand_char {
                    inconsistent = true;
                    break;
                }
            }
            if inconsistent {
                log::warn!(
                    "quantMode TranscriptomeSAM: transcript {} has inconsistent chromosome/strand across exons — skipping",
                    tid
                );
                continue;
            }

            // Map chromosome name to genome index.
            let chr_idx = match genome.chr_name.iter().position(|n| n == chr_name) {
                Some(i) if i < genome.n_chr_real => i,
                _ => {
                    log::warn!(
                        "quantMode TranscriptomeSAM: transcript {} on unknown chromosome {} — skipping",
                        tid,
                        chr_name
                    );
                    continue;
                }
            };
            let chr_offset = genome.chr_start[chr_idx];

            // Sort exons by GTF start ASC.
            recs.sort_by_key(|r| r.start);

            // Build TrExon list in absolute genome coords (0-based half-open).
            let mut tr_exons: Vec<TrExon> = Vec::with_capacity(recs.len());
            let mut ex_len_cum: u32 = 0;
            for r in &recs {
                // GTF is 1-based inclusive; ruSTAR uses 0-based half-open.
                let abs_start = chr_offset + r.start.saturating_sub(1);
                let abs_end = chr_offset + r.end;
                if abs_end <= abs_start {
                    log::warn!(
                        "quantMode TranscriptomeSAM: transcript {} has invalid exon {}-{} — skipping",
                        tid,
                        r.start,
                        r.end
                    );
                    tr_exons.clear();
                    break;
                }
                let len = (abs_end - abs_start) as u32;
                tr_exons.push(TrExon {
                    genome_start: abs_start,
                    genome_end: abs_end,
                    ex_len_cum,
                });
                ex_len_cum = ex_len_cum.saturating_add(len);
            }
            if tr_exons.is_empty() {
                continue;
            }

            let first_ex = tr_exons.first().unwrap();
            let last_ex = tr_exons.last().unwrap();
            let start_abs = first_ex.genome_start;
            let end_abs = last_ex.genome_end;
            let total_len = ex_len_cum;

            let strand_u8 = match strand_char {
                '+' => 1u8,
                '-' => 2u8,
                _ => 1u8, // unknown → treat as forward (STAR default)
            };

            let gene_id = first
                .attributes
                .get("gene_id")
                .cloned()
                .unwrap_or_default();

            tr_ids.push(tid.clone());
            tr_chr_idx.push(chr_idx);
            tr_strand.push(strand_u8);
            tr_gene_id.push(gene_id);
            tr_start.push(start_abs);
            tr_end.push(end_abs);
            tr_exons_vec.push(tr_exons);
            tr_length.push(total_len);
        }

        let n_tr = tr_ids.len();

        // Sorted view by (tr_start, tr_end) for binary-search +
        // running-max early-exit.
        let mut tr_order: Vec<usize> = (0..n_tr).collect();
        tr_order.sort_by(|&a, &b| {
            tr_start[a]
                .cmp(&tr_start[b])
                .then_with(|| tr_end[a].cmp(&tr_end[b]))
        });

        let tr_starts_sorted: Vec<u64> = tr_order.iter().map(|&i| tr_start[i]).collect();
        let mut tr_end_max_sorted: Vec<u64> = Vec::with_capacity(n_tr);
        let mut running_max: u64 = 0;
        for &i in &tr_order {
            running_max = running_max.max(tr_end[i]);
            tr_end_max_sorted.push(running_max);
        }

        Ok(TranscriptomeIndex {
            tr_ids,
            tr_chr_idx,
            tr_strand,
            tr_gene_id,
            tr_start,
            tr_end,
            tr_exons: tr_exons_vec,
            tr_length,
            tr_order,
            tr_starts_sorted,
            tr_end_max_sorted,
        })
    }

    /// Number of transcripts indexed.
    pub fn n_transcripts(&self) -> usize {
        self.tr_ids.len()
    }
}

// ---------------------------------------------------------------------------
// Projection — STAR's `Transcriptome::quantAlign` + `alignToTranscript`.
// ---------------------------------------------------------------------------

/// Project a genome-space `Transcript` onto every transcript whose coordinates
/// fully contain the alignment.  Returns transcript-space alignments ready for
/// BAM emission.
///
/// Mirrors STAR `Transcriptome_quantAlign.cpp`:
///   * binary-search `tr_starts_sorted` for the greatest `tr_start <=
///     align.genome_start`,
///   * walk back while `tr_end_max_sorted[i] >= align.genome_end`,
///   * for each candidate whose `[tr_start, tr_end]` fully contains the align,
///     call `align_to_one_transcript`.
///
/// `lread` is the total read length and is needed to flip read-space
/// coordinates on reverse-strand transcripts (`read_pos' = Lread - (read_pos +
/// len)`).  For paired-end inputs pass the sum of per-mate lengths.
pub fn align_to_transcripts(
    align: &Transcript,
    idx: &TranscriptomeIndex,
    lread: u32,
) -> Vec<Transcript> {
    let mut out: Vec<Transcript> = Vec::new();
    if idx.n_transcripts() == 0 || align.exons.is_empty() {
        return out;
    }

    let a_start = align.genome_start;
    let a_end = align.genome_end;

    // Binary-search `tr_starts_sorted` for greatest tr_start <= a_start.
    // partition_point returns the first index with tr_start > a_start; subtract 1.
    let upper = idx
        .tr_starts_sorted
        .partition_point(|&s| s <= a_start);
    if upper == 0 {
        return out; // a_start is to the left of all transcripts
    }
    let mut sorted_i = upper - 1;

    // Walk backwards while the running-max end still covers this alignment.
    loop {
        if idx.tr_end_max_sorted[sorted_i] < a_end {
            break;
        }
        let tr_idx = idx.tr_order[sorted_i];
        if idx.tr_chr_idx[tr_idx] == align.chr_idx
            && idx.tr_start[tr_idx] <= a_start
            && idx.tr_end[tr_idx] >= a_end
            && let Some(projected) = align_to_one_transcript(align, tr_idx, idx, lread)
        {
            out.push(projected);
        }
        if sorted_i == 0 {
            break;
        }
        sorted_i -= 1;
    }
    out
}

/// Project a single alignment onto a single transcript.  Returns `None` if the
/// alignment is inconsistent with the transcript's exon structure (block
/// extends past an exon boundary, splice boundary doesn't match a transcript
/// junction, etc.).
///
/// Direct port of STAR's `alignToTranscript` (see
/// `source/Transcriptome_quantAlign.cpp:5-89`).
fn align_to_one_transcript(
    align: &Transcript,
    tr_idx: usize,
    idx: &TranscriptomeIndex,
    lread: u32,
) -> Option<Transcript> {
    let tr_exons = &idx.tr_exons[tr_idx];
    let tr_strand = idx.tr_strand[tr_idx];
    let tr_length = idx.tr_length[tr_idx];

    // Find exon containing first block's start.
    let first_block = &align.exons[0];
    let g1 = first_block.genome_start;

    let mut ex = match find_containing_exon(tr_exons, g1) {
        Some(idx) => idx,
        None => return None,
    };

    // Build projected exons (in t-space) as we walk the alignment blocks.
    let mut proj_exons: Vec<Exon> = Vec::new();

    for iab in 0..align.exons.len() {
        let block = &align.exons[iab];
        let block_start = block.genome_start;
        let block_end = block.genome_end; // half-open exclusive

        // Block must not extend past the current transcript exon's end.
        if block_end > tr_exons[ex].genome_end {
            return None;
        }
        if block_start < tr_exons[ex].genome_start {
            return None;
        }

        // Determine if this block starts a new projected exon or extends the
        // previous one.  STAR starts a new projected exon on: (a) the first
        // block, or (b) when the preceding canonSJ was a junction (>= 0).  We
        // handle those boundaries via `mate_or_junction_before`.
        let start_new = iab == 0 || is_splice_boundary_before(align, iab);

        if start_new {
            // t-space position = ex_len_cum + (block_start - exon_start)
            let tr_offset =
                tr_exons[ex].ex_len_cum as u64 + (block_start - tr_exons[ex].genome_start);
            let len = (block_end - block_start) as usize;
            proj_exons.push(Exon {
                genome_start: tr_offset,
                genome_end: tr_offset + len as u64,
                read_start: block.read_start,
                read_end: block.read_start + len,
            });
        } else {
            // Coalesce: extend the last projected exon by this block's length
            // (STAR adds block EX_L directly — does NOT update start position).
            if let Some(last) = proj_exons.last_mut() {
                let len = (block_end - block_start) as u64;
                last.genome_end += len;
                last.read_end = block.read_end;
            } else {
                return None;
            }
        }

        // Advance `ex` across any splice boundary BEFORE the next block.
        if iab + 1 < align.exons.len() && is_splice_boundary_before(align, iab + 1) {
            // Require the junction to match a transcript junction.
            let next_block = &align.exons[iab + 1];
            if ex + 1 >= tr_exons.len() {
                return None;
            }
            if block_end != tr_exons[ex].genome_end {
                return None;
            }
            if next_block.genome_start != tr_exons[ex + 1].genome_start {
                return None;
            }
            ex += 1;
        }
    }

    // Apply reverse-strand coordinate flip (STAR canonSJ == -999 branch).
    let projected_is_reverse = if tr_strand == 2 {
        !align.is_reverse
    } else {
        align.is_reverse
    };

    if tr_strand == 2 {
        let tr_len = tr_length as u64;
        let lread_u = lread as u64;
        for e in proj_exons.iter_mut() {
            let len = e.genome_end - e.genome_start;
            let new_g = tr_len - (e.genome_start + len);
            e.genome_start = new_g;
            e.genome_end = new_g + len;

            let read_len = e.read_end - e.read_start;
            let new_r = (lread_u as usize).saturating_sub(e.read_start + read_len);
            e.read_start = new_r;
            e.read_end = new_r + read_len;
        }
        proj_exons.reverse();
    }

    // Build projected CIGAR: drop N operations (splices collapse in t-space);
    // for reverse-strand transcripts, reverse the resulting op sequence.
    let mut proj_cigar: Vec<CigarOp> = align
        .cigar
        .iter()
        .filter(|op| !matches!(op, CigarOp::RefSkip(_)))
        .copied()
        .collect();
    if tr_strand == 2 {
        proj_cigar.reverse();
    }

    // Projected genome bounds = outermost t-space exon positions.
    let proj_start = proj_exons.first().map(|e| e.genome_start).unwrap_or(0);
    let proj_end = proj_exons.last().map(|e| e.genome_end).unwrap_or(0);

    Some(Transcript {
        chr_idx: tr_idx,
        genome_start: proj_start,
        genome_end: proj_end,
        is_reverse: projected_is_reverse,
        exons: proj_exons,
        cigar: proj_cigar,
        score: align.score,
        n_mismatch: align.n_mismatch,
        n_gap: align.n_gap,
        // Splices collapse in t-space; leave junction metadata empty.
        n_junction: 0,
        junction_motifs: Vec::new(),
        junction_annotated: Vec::new(),
        read_seq: align.read_seq.clone(),
    })
}

/// Find the transcript exon (by index in `tr_exons`) that contains position
/// `pos` (0-based genome coord).  Returns `None` if `pos` is in an intron or
/// outside the transcript.
fn find_containing_exon(tr_exons: &[TrExon], pos: u64) -> Option<usize> {
    let mut lo = 0usize;
    let mut hi = tr_exons.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if tr_exons[mid].genome_end <= pos {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo < tr_exons.len() && pos >= tr_exons[lo].genome_start && pos < tr_exons[lo].genome_end {
        Some(lo)
    } else {
        None
    }
}

/// Return true if the boundary between `align.exons[iab-1]` and
/// `align.exons[iab]` is a splice (`RefSkip`) rather than an indel.
///
/// ruSTAR's `Transcript.exons` is a list of read-contiguous match blocks;
/// splices / insertions / deletions all create block boundaries.  We
/// discriminate based on the read-side gap:
///   * `read_end_prev == read_start_curr` AND `genome gap` → deletion OR splice.
///     We call it a splice — the caller verifies the boundary matches a
///     transcript junction and rejects the alignment if it does not.
///   * read gap (insertion) → not a splice (coalesce).
///   * Pure deletion (read contiguous, small genome gap that does NOT match a
///     transcript junction) → handled by the caller rejecting the alignment
///     when the junction check fails.  In practice ruSTAR produces `Del` ops
///     inside a single exon (no block split for pure deletions because the
///     stitch merge coalesces across Del), so this branch rarely fires.
fn is_splice_boundary_before(align: &Transcript, iab: usize) -> bool {
    if iab == 0 {
        return false;
    }
    let prev = &align.exons[iab - 1];
    let cur = &align.exons[iab];
    // Insertion: read gap between blocks with no genome gap → coalesce.
    if prev.read_end < cur.read_start && prev.genome_end == cur.genome_start {
        return false;
    }
    // Splice / large gap on the genome side, with read-contiguous: treat as
    // potential splice junction (caller validates).
    prev.genome_end < cur.genome_start
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap as StdHashMap;

    fn make_genome() -> Genome {
        Genome {
            sequence: vec![0u8; 3000],
            n_genome: 3000,
            n_chr_real: 2,
            chr_start: vec![0, 1000, 3000],
            chr_length: vec![1000, 2000],
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
        }
    }

    fn make_exon(
        seqname: &str,
        start: u64,
        end: u64,
        strand: char,
        gene_id: &str,
        transcript_id: &str,
    ) -> GtfRecord {
        let mut attrs: StdHashMap<String, String> = StdHashMap::new();
        attrs.insert("gene_id".to_string(), gene_id.to_string());
        attrs.insert("transcript_id".to_string(), transcript_id.to_string());
        GtfRecord {
            seqname: seqname.to_string(),
            feature: "exon".to_string(),
            start,
            end,
            strand,
            attributes: attrs,
        }
    }

    #[test]
    fn single_exon_transcript_metadata() {
        let genome = make_genome();
        let exons = vec![make_exon("chr1", 101, 200, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 1);
        assert_eq!(idx.tr_ids[0], "T1");
        assert_eq!(idx.tr_gene_id[0], "G1");
        assert_eq!(idx.tr_strand[0], 1);
        assert_eq!(idx.tr_chr_idx[0], 0);
        // GTF 1-based 101 → absolute 0-based 100; end 200 → absolute 200 (exclusive)
        assert_eq!(idx.tr_start[0], 100);
        assert_eq!(idx.tr_end[0], 200);
        assert_eq!(idx.tr_length[0], 100);
        assert_eq!(idx.tr_exons[0].len(), 1);
        assert_eq!(idx.tr_exons[0][0].ex_len_cum, 0);
    }

    #[test]
    fn multi_exon_transcript_ex_len_cum() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),   // len 100
            make_exon("chr1", 301, 400, '+', "G1", "T1"),   // len 100
            make_exon("chr1", 501, 650, '+', "G1", "T1"),   // len 150
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 1);
        assert_eq!(idx.tr_length[0], 350);

        let ex = &idx.tr_exons[0];
        assert_eq!(ex[0].ex_len_cum, 0);
        assert_eq!(ex[1].ex_len_cum, 100);
        assert_eq!(ex[2].ex_len_cum, 200);

        // Ex_len_cum must be monotonically non-decreasing.
        for w in ex.windows(2) {
            assert!(w[0].ex_len_cum <= w[1].ex_len_cum);
        }

        // tr_start / tr_end from first/last absolute exon.
        assert_eq!(idx.tr_start[0], 100);
        assert_eq!(idx.tr_end[0], 650);
    }

    #[test]
    fn reverse_strand_transcript() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '-', "G2", "T2"),
            make_exon("chr1", 301, 400, '-', "G2", "T2"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 1);
        assert_eq!(idx.tr_strand[0], 2);
        assert_eq!(idx.tr_length[0], 200);
    }

    #[test]
    fn two_transcripts_same_gene() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '+', "G1", "T1"),
            make_exon("chr1", 101, 200, '+', "G1", "T2"),
            make_exon("chr1", 301, 500, '+', "G1", "T2"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 2);
        assert_eq!(idx.tr_gene_id[0], "G1");
        assert_eq!(idx.tr_gene_id[1], "G1");
        // Different lengths: T1 = 100+100 = 200, T2 = 100+200 = 300
        assert_eq!(idx.tr_length[0], 200);
        assert_eq!(idx.tr_length[1], 300);
    }

    #[test]
    fn unknown_chromosome_skipped() {
        let genome = make_genome();
        let exons = vec![make_exon("chrX", 101, 200, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 0);
    }

    #[test]
    fn inconsistent_strand_skipped() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '-', "G1", "T1"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 0);
    }

    fn make_align(
        chr_idx: usize,
        is_reverse: bool,
        exons: Vec<(u64, u64, usize, usize)>,
        cigar: Vec<CigarOp>,
    ) -> Transcript {
        let proj_exons: Vec<Exon> = exons
            .into_iter()
            .map(|(gs, ge, rs, re)| Exon {
                genome_start: gs,
                genome_end: ge,
                read_start: rs,
                read_end: re,
            })
            .collect();
        let gs = proj_exons.first().map(|e| e.genome_start).unwrap_or(0);
        let ge = proj_exons.last().map(|e| e.genome_end).unwrap_or(0);
        Transcript {
            chr_idx,
            genome_start: gs,
            genome_end: ge,
            is_reverse,
            exons: proj_exons,
            cigar,
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        }
    }

    #[test]
    fn project_single_exon_align_into_single_exon_transcript() {
        let genome = make_genome();
        // Transcript: chr1 [100, 200) forward.
        let gtf = vec![make_exon("chr1", 101, 200, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align fully inside exon: genome [110, 150), read [0, 40).
        let align = make_align(
            0,
            false,
            vec![(110, 150, 0, 40)],
            vec![CigarOp::Match(40)],
        );
        let results = align_to_transcripts(&align, &idx, 40);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.chr_idx, 0); // transcript index
        assert_eq!(r.is_reverse, false);
        assert_eq!(r.exons.len(), 1);
        // t-space offset = 0 (ex_len_cum) + (110 - 100) = 10
        assert_eq!(r.exons[0].genome_start, 10);
        assert_eq!(r.exons[0].genome_end, 50);
        assert_eq!(r.exons[0].read_start, 0);
        assert_eq!(r.exons[0].read_end, 40);
        // CIGAR must have no N
        assert!(r.cigar.iter().all(|op| !matches!(op, CigarOp::RefSkip(_))));
        assert_eq!(r.cigar.len(), 1);
    }

    #[test]
    fn project_two_exon_align_matching_junction() {
        let genome = make_genome();
        // Transcript T1: chr1 [100, 200) + [300, 400) forward. tr_length = 200.
        let gtf = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '+', "G1", "T1"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align: two blocks (50 M each), junction matches transcript junction.
        // Genome: [150, 200) + [300, 350); read [0, 50) + [50, 100).
        let align = make_align(
            0,
            false,
            vec![(150, 200, 0, 50), (300, 350, 50, 100)],
            vec![CigarOp::Match(50), CigarOp::RefSkip(100), CigarOp::Match(50)],
        );
        let results = align_to_transcripts(&align, &idx, 100);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        // Two t-space exons, no N in CIGAR.
        assert_eq!(r.exons.len(), 2);
        // First exon: t-space [50, 100)
        assert_eq!(r.exons[0].genome_start, 50);
        assert_eq!(r.exons[0].genome_end, 100);
        // Second exon: t-space [100, 150) — starts right after prev (splice collapsed)
        assert_eq!(r.exons[1].genome_start, 100);
        assert_eq!(r.exons[1].genome_end, 150);
        // CIGAR: no N
        assert!(r.cigar.iter().all(|op| !matches!(op, CigarOp::RefSkip(_))));
        assert_eq!(r.cigar.len(), 2);
    }

    #[test]
    fn project_mismatched_junction_fails() {
        let genome = make_genome();
        let gtf = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '+', "G1", "T1"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align with junction NOT matching transcript — splice ends at 195 instead of 200.
        let align = make_align(
            0,
            false,
            vec![(150, 195, 0, 45), (305, 350, 45, 90)],
            vec![CigarOp::Match(45), CigarOp::RefSkip(110), CigarOp::Match(45)],
        );
        let results = align_to_transcripts(&align, &idx, 90);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn project_onto_reverse_strand_transcript() {
        let genome = make_genome();
        // Reverse transcript: chr1 [100, 200) - strand. tr_length = 100.
        let gtf = vec![make_exon("chr1", 101, 200, '-', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align forward on genome [120, 160), read [0, 40).
        let align = make_align(
            0,
            false,
            vec![(120, 160, 0, 40)],
            vec![CigarOp::Match(40)],
        );
        let results = align_to_transcripts(&align, &idx, 40);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        // is_reverse flipped (transcript strand == 2, align was false → result true)
        assert_eq!(r.is_reverse, true);
        // t-space position after flip:
        //   pre-flip: genome_start=20, length=40 → new_g = 100 - (20 + 40) = 40
        //   read_start=0, read_len=40 → new_r = 40 - (0 + 40) = 0
        assert_eq!(r.exons.len(), 1);
        assert_eq!(r.exons[0].genome_start, 40);
        assert_eq!(r.exons[0].genome_end, 80);
        assert_eq!(r.exons[0].read_start, 0);
        assert_eq!(r.exons[0].read_end, 40);
    }

    #[test]
    fn project_multi_exon_align_onto_longer_transcript() {
        let genome = make_genome();
        // 3-exon transcript: [100,200) + [300,400) + [500,600). tr_length = 300.
        let gtf = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '+', "G1", "T1"),
            make_exon("chr1", 501, 600, '+', "G1", "T1"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align spans just exon 2 and 3 (skipping first exon): genome [350, 400) + [500, 550).
        let align = make_align(
            0,
            false,
            vec![(350, 400, 0, 50), (500, 550, 50, 100)],
            vec![CigarOp::Match(50), CigarOp::RefSkip(100), CigarOp::Match(50)],
        );
        let results = align_to_transcripts(&align, &idx, 100);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.exons.len(), 2);
        // First t-space exon: ex_len_cum[1]=100, offset within exon = 350-300=50 → t-space start 150
        assert_eq!(r.exons[0].genome_start, 150);
        assert_eq!(r.exons[0].genome_end, 200);
        // Second t-space exon: ex_len_cum[2]=200, offset = 500-500=0 → t-space start 200
        assert_eq!(r.exons[1].genome_start, 200);
        assert_eq!(r.exons[1].genome_end, 250);
    }

    #[test]
    fn project_past_transcript_end_fails() {
        let genome = make_genome();
        let gtf = vec![make_exon("chr1", 101, 200, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align extends past transcript end (genome [150, 250), transcript ends at 200).
        let align = make_align(
            0,
            false,
            vec![(150, 250, 0, 100)],
            vec![CigarOp::Match(100)],
        );
        let results = align_to_transcripts(&align, &idx, 100);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn project_before_all_transcripts_returns_empty() {
        let genome = make_genome();
        let gtf = vec![make_exon("chr1", 501, 600, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        let align = make_align(
            0,
            false,
            vec![(100, 150, 0, 50)],
            vec![CigarOp::Match(50)],
        );
        let results = align_to_transcripts(&align, &idx, 50);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn project_onto_multiple_overlapping_transcripts() {
        let genome = make_genome();
        // Two transcripts both containing the align:
        //  T1: [100, 400) single exon
        //  T2: [100, 400) single exon, different ID
        let gtf = vec![
            make_exon("chr1", 101, 400, '+', "G1", "T1"),
            make_exon("chr1", 101, 400, '+', "G2", "T2"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        let align = make_align(
            0,
            false,
            vec![(200, 250, 0, 50)],
            vec![CigarOp::Match(50)],
        );
        let results = align_to_transcripts(&align, &idx, 50);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn tr_end_max_sorted_is_running_max() {
        let genome = make_genome();
        let exons = vec![
            // T1 chr1 starts at 100, ends at 200
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            // T2 chr1 starts at 150, ends at 500 (encloses T1)
            make_exon("chr1", 151, 500, '+', "G1", "T2"),
            // T3 chr1 starts at 200, ends at 300 (nested)
            make_exon("chr1", 201, 300, '+', "G1", "T3"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 3);

        // tr_order must be sorted by (start, end)
        let sorted_starts: Vec<u64> = idx.tr_starts_sorted.clone();
        let mut check = sorted_starts.clone();
        check.sort();
        assert_eq!(sorted_starts, check);

        // tr_end_max_sorted must be monotonically non-decreasing
        for w in idx.tr_end_max_sorted.windows(2) {
            assert!(w[0] <= w[1]);
        }
        // Last entry equals overall maximum tr_end
        let overall_max = *idx.tr_end.iter().max().unwrap();
        assert_eq!(*idx.tr_end_max_sorted.last().unwrap(), overall_max);
    }
}
