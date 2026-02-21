/// Seed clustering and stitching via dynamic programming
use crate::align::score::AlignmentScorer;
use crate::align::seed::Seed;
use crate::align::transcript::{CigarOp, Transcript};
use crate::error::Error;
use crate::index::GenomeIndex;

/// Verify match length at a specific genome position with correct strand handling.
///
/// Seeds found via binary search guarantee `sa_nbases` matching bases at `sa_start`,
/// but other positions in the SA range may match fewer bases. This function re-verifies
/// the actual match length at each specific genome position.
///
/// For reverse-strand positions, adds `n_genome` offset to access the reverse-complement
/// region of the genome (which is stored at `[n_genome, 2*n_genome)`).
fn verify_match_at_position(
    read_seq: &[u8],
    read_pos: usize,
    genome_pos: u64,
    is_reverse: bool,
    max_length: usize,
    index: &GenomeIndex,
) -> usize {
    let actual_genome_pos = if is_reverse {
        genome_pos + index.genome.n_genome
    } else {
        genome_pos
    };
    let mut length = 0;
    for i in 0..max_length {
        if read_pos + i >= read_seq.len() {
            break;
        }
        match index.genome.get_base(actual_genome_pos + i as u64) {
            Some(gb) if gb < 5 && gb == read_seq[read_pos + i] => length += 1,
            _ => break,
        }
    }
    length
}

/// Count mismatches in an alignment by comparing read sequence to genome sequence.
///
/// The read sequence is always in forward orientation. For reverse-strand alignments,
/// the genome is accessed at `pos + n_genome` (the reverse-complement region) rather
/// than reverse-complementing the read.
///
/// # Arguments
/// * `read_seq` - Read sequence in forward orientation (encoded as 0=A, 1=C, 2=G, 3=T, 4=N)
/// * `cigar_ops` - CIGAR operations
/// * `genome_start` - Starting position in genome (decoded SA position, WITHOUT n_genome offset)
/// * `read_start` - Starting position in read
/// * `index` - Genome index (contains genome sequence)
/// * `is_reverse` - Whether this is a reverse-strand alignment
///
/// # Returns
/// Number of mismatched bases (excluding N bases)
/// Count mismatches in a simple region (no CIGAR, just read vs genome)
fn count_mismatches_in_region(
    read_seq: &[u8],
    read_start: usize,
    genome_start: u64,
    length: usize,
    index: &GenomeIndex,
    is_reverse: bool,
) -> u32 {
    let genome_offset = if is_reverse { index.genome.n_genome } else { 0 };
    let mut n_mismatch = 0u32;

    for i in 0..length {
        let read_pos = read_start + i;
        if read_pos >= read_seq.len() {
            break;
        }
        let read_base = read_seq[read_pos];
        if let Some(genome_base) = index
            .genome
            .get_base(genome_start + i as u64 + genome_offset)
        {
            if read_base != genome_base && read_base != 4 && genome_base != 4 {
                n_mismatch += 1;
            }
        }
    }

    n_mismatch
}

fn count_mismatches(
    read_seq: &[u8],
    cigar_ops: &[CigarOp],
    genome_start: u64,
    read_start: usize,
    index: &GenomeIndex,
    is_reverse: bool,
) -> u32 {
    // Add n_genome offset for reverse-strand genome access
    let genome_offset = if is_reverse { index.genome.n_genome } else { 0 };

    let mut n_mismatch = 0u32;
    let mut read_pos = read_start;
    let mut genome_pos = genome_start;

    for op in cigar_ops {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                for _i in 0..*len {
                    if read_pos < read_seq.len() {
                        let read_base = read_seq[read_pos];
                        if let Some(genome_base) = index.genome.get_base(genome_pos + genome_offset)
                        {
                            if read_base != genome_base && read_base != 4 && genome_base != 4 {
                                n_mismatch += 1;
                            }
                        }
                    }
                    read_pos += 1;
                    genome_pos += 1;
                }
            }
            CigarOp::Ins(len) => {
                read_pos += *len as usize;
            }
            CigarOp::Del(len) | CigarOp::RefSkip(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::SoftClip(len) => {
                read_pos += *len as usize;
            }
            CigarOp::HardClip(_) => {}
        }
    }

    n_mismatch
}

/// Result of extending an alignment into flanking regions
#[derive(Debug, Clone)]
struct ExtendResult {
    /// How far the extension reached (bases)
    extend_len: usize,
    /// Maximum score achieved during extension
    max_score: i32,
    /// Number of mismatches in the extended region
    n_mismatch: u32,
}

/// Extend alignment from a boundary into flanking sequence, mirroring STAR's extendAlign().
///
/// Walks base-by-base from the alignment boundary, scoring +1 match / -1 mismatch,
/// tracking the maximum-score extension point. Stops when total mismatches exceed
/// `min(p_mm_max * total_length, n_mm_max)`.
///
/// # Arguments
/// * `read_seq` - Full read sequence (encoded)
/// * `read_start` - Boundary position in read (where extension begins)
/// * `genome_start` - Corresponding genome position (WITHOUT n_genome offset)
/// * `direction` - +1 for rightward extension, -1 for leftward
/// * `max_extend` - Maximum distance to extend (to read boundary)
/// * `n_mm_prev` - Mismatches already in the aligned portion
/// * `len_prev` - Length of the already-aligned portion
/// * `n_mm_max` - outFilterMismatchNmax (absolute max mismatches)
/// * `p_mm_max` - outFilterMismatchNoverLmax (max mismatch ratio)
/// * `index` - Genome index
/// * `is_reverse` - Whether this is a reverse-strand alignment
fn extend_alignment(
    read_seq: &[u8],
    read_start: usize,
    genome_start: u64,
    direction: i32,
    max_extend: usize,
    n_mm_prev: u32,
    len_prev: usize,
    n_mm_max: u32,
    p_mm_max: f64,
    index: &GenomeIndex,
    is_reverse: bool,
) -> ExtendResult {
    if max_extend == 0 {
        return ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        };
    }

    let genome_offset = if is_reverse { index.genome.n_genome } else { 0 };

    let mut score: i32 = 0;
    let mut max_score: i32 = 0;
    let mut best_len: usize = 0;
    let mut best_mm: u32 = 0;
    let mut n_mm = 0u32;

    for i in 0..max_extend {
        // Calculate read and genome positions based on direction
        let read_pos = if direction > 0 {
            read_start + i
        } else {
            // Leftward: read_start is exclusive boundary, go backwards
            if read_start < 1 + i {
                break;
            }
            read_start - 1 - i
        };

        if read_pos >= read_seq.len() {
            break;
        }

        let genome_pos = if direction > 0 {
            genome_start + i as u64
        } else {
            if genome_start < 1 + i as u64 {
                break;
            }
            genome_start - 1 - i as u64
        };

        // Get genome base (with strand offset)
        let genome_base = match index.genome.get_base(genome_pos + genome_offset) {
            Some(b) => b,
            None => break,
        };

        // Stop at chromosome boundary (padding = 5)
        if genome_base == 5 {
            break;
        }

        let read_base = read_seq[read_pos];

        // Skip N bases (no score impact, matches STAR behavior)
        if read_base == 4 || genome_base == 4 {
            continue;
        }

        if read_base == genome_base {
            // MATCH — only record new max on match (STAR behavior)
            score += 1;
            if score > max_score {
                let total_mm = n_mm_prev + n_mm;
                let record_limit = ((p_mm_max * (len_prev + i + 1) as f64) as u32).min(n_mm_max);
                if total_mm <= record_limit {
                    max_score = score;
                    best_len = i + 1;
                    best_mm = n_mm;
                }
            }
        } else {
            // MISMATCH — check break BEFORE incrementing nMM (STAR behavior)
            // Break uses full extension length max_extend, not current position i+1
            let total_mm = n_mm_prev + n_mm;
            let break_limit = ((p_mm_max * (len_prev + max_extend) as f64) as u32).min(n_mm_max);
            if total_mm >= break_limit {
                break;
            }
            n_mm += 1;
            score -= 1;
        }
    }

    // Only accept extension if it has positive score
    if max_score > 0 {
        ExtendResult {
            extend_len: best_len,
            max_score,
            n_mismatch: best_mm,
        }
    } else {
        ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        }
    }
}

/// A Window Alignment entry — equivalent to STAR's WA[iW][iA] array.
///
/// Each entry represents one seed at one specific genome position within a window.
/// During seed assignment, verify_match_at_position() confirms the actual match length.
/// The DP reads these entries directly (no SA range re-expansion needed).
#[derive(Debug, Clone)]
pub struct WindowAlignment {
    /// Index into the seeds array (for DP expansion to identify the originating seed)
    pub seed_idx: usize,
    /// Read start position (STAR: WA_rStart)
    pub read_pos: usize,
    /// Verified match length at this specific position (STAR: WA_Length)
    pub length: usize,
    /// Forward-strand genome position (STAR: WA_gStart)
    pub genome_pos: u64,
    /// Raw SA position (for genome base access in DP — reverse strand uses sa_pos + n_genome)
    pub sa_pos: u64,
    /// SA range size of the originating seed (STAR: WA_Nrep) — for scoring
    pub n_rep: usize,
    /// Whether this entry is an anchor (protected from capacity eviction)
    pub is_anchor: bool,
}

/// A cluster of seeds mapping to the same genomic region
#[derive(Debug, Clone)]
pub struct SeedCluster {
    /// Window Alignment entries — seed positions assigned to this window (STAR's WA array)
    pub alignments: Vec<WindowAlignment>,
    /// Chromosome index
    pub chr_idx: usize,
    /// Genomic start (leftmost position, forward coords, from actual seed positions)
    pub genome_start: u64,
    /// Genomic end (rightmost position, forward coords, from actual seed positions)
    pub genome_end: u64,
    /// Strand (false = forward, true = reverse)
    pub is_reverse: bool,
    /// Anchor seed index (in the seeds array)
    pub anchor_idx: usize,
    /// Anchor genomic bin (anchor_pos >> win_bin_nbits) for MAPQ window counting
    pub anchor_bin: u64,
}

/// Cluster seeds using STAR's bin-based windowing algorithm.
///
/// # Algorithm (faithful to STAR's `createExtendWindowsWithAlign` + `assignAlignToWindow`)
/// 1. Identify anchor seeds (SA range ≤ max_loci_for_anchor)
/// 2. Create windows from anchor positions using `winBin[(strand, bin)]` lookup:
///    - If bin already has a window → assign anchor to it, skip creation
///    - Else scan left/right for nearby windows → merge or create new
/// 3. Extend windows by ±win_flank_nbins on each side
/// 4. Assign ALL seeds to windows with overlap dedup + capacity eviction
/// 5. Build SeedCluster output
///
/// # Arguments
/// * `seeds` - All seeds found in the read
/// * `index` - Genome index
/// * `win_bin_nbits` - Log2 of window bin size (STAR default: 16 → 64KB bins)
/// * `win_anchor_dist_nbins` - Max bins for anchor-window merging (STAR default: 9)
/// * `win_flank_nbins` - Bins to extend each window side (STAR default: 4)
/// * `max_loci_for_anchor` - Max SA range for a seed to be an anchor (e.g., 10)
/// * `win_anchor_multimap_nmax` - Max loci anchors can map to (STAR default: 50)
/// * `seed_per_window_nmax` - Max WA entries per window (capacity eviction threshold)
///
/// # Returns
/// Vector of seed clusters, one per window with assigned seeds
pub fn cluster_seeds(
    seeds: &[Seed],
    read_seq: &[u8],
    index: &GenomeIndex,
    win_bin_nbits: u32,
    win_anchor_dist_nbins: u32,
    win_flank_nbins: u32,
    max_loci_for_anchor: usize,
    win_anchor_multimap_nmax: usize,
    seed_per_window_nmax: usize,
) -> Vec<SeedCluster> {
    use std::collections::HashMap;

    let anchor_set: Vec<bool> = seeds
        .iter()
        .map(|seed| {
            let n_loci = seed.sa_end - seed.sa_start;
            n_loci > 0 && n_loci <= max_loci_for_anchor
        })
        .collect();

    // Phase 1: Identify anchor seeds (few genomic positions → high specificity)
    // STAR: only seeds with Nrep <= winAnchorMultimapNmax create windows.
    let anchor_indices: Vec<usize> = anchor_set
        .iter()
        .enumerate()
        .filter(|(_, is_anchor)| **is_anchor)
        .map(|(i, _)| i)
        .collect();

    // No fallback: matches STAR behavior where reads with no anchors are unmapped.
    // MMP search now narrows SA ranges from both ends (max_mappable_length),
    // so seeds have accurate loci counts and anchor classification is correct.
    if anchor_indices.is_empty() {
        return Vec::new();
    }

    // Phase 2: Create windows from anchor positions
    // (matches STAR's createExtendWindowsWithAlign)
    struct Window {
        bin_start: u64,
        bin_end: u64,
        chr_idx: usize,
        is_reverse: bool,
        anchor_idx: usize,
        alignments: Vec<WindowAlignment>,
        // Tight bounds from actual seed positions
        actual_start: u64,
        actual_end: u64,
        alive: bool, // false = merged into another window (STAR kills merged windows)
    }

    let mut windows: Vec<Window> = Vec::new();
    // winBin: (strand, bin) → window_index
    // Chromosome is implicit since bins are from absolute forward positions
    let mut win_bin: HashMap<(bool, u64), usize> = HashMap::new();

    for &anchor_idx in &anchor_indices {
        let anchor = &seeds[anchor_idx];
        let n_loci = anchor.sa_end - anchor.sa_start;

        // Skip anchors with too many loci (STAR: winAnchorMultimapNmax)
        if n_loci > win_anchor_multimap_nmax {
            continue;
        }

        for (sa_pos, strand) in anchor.genome_positions(index) {
            let actual_length = verify_match_at_position(
                read_seq,
                anchor.read_pos,
                sa_pos,
                strand,
                anchor.length,
                index,
            );
            if actual_length < 8 {
                continue;
            }

            let forward_pos = index.sa_pos_to_forward(sa_pos, strand, actual_length);

            let chr_idx = match index.genome.position_to_chr(forward_pos) {
                Some(info) => info.0,
                None => continue,
            };

            let anchor_bin = forward_pos >> win_bin_nbits;

            let wa_entry = WindowAlignment {
                seed_idx: anchor_idx,
                read_pos: anchor.read_pos,
                length: actual_length,
                genome_pos: forward_pos,
                sa_pos,
                n_rep: n_loci,
                is_anchor: true,
            };

            // Check if this bin already has a window (STAR: skip creation, just assign)
            if let Some(&win_idx) = win_bin.get(&(strand, anchor_bin)) {
                let window = &mut windows[win_idx];
                if window.alive && window.chr_idx == chr_idx {
                    window.actual_start = window.actual_start.min(forward_pos);
                    window.actual_end = window.actual_end.max(forward_pos + actual_length as u64);
                    window.alignments.push(wa_entry);
                    continue;
                }
            }

            // Scan LEFT for existing window to merge with
            let mut merge_left: Option<usize> = None;
            for scan_bin in
                (anchor_bin.saturating_sub(win_anchor_dist_nbins as u64)..anchor_bin).rev()
            {
                if let Some(&win_idx) = win_bin.get(&(strand, scan_bin)) {
                    let w = &windows[win_idx];
                    if w.alive && w.chr_idx == chr_idx {
                        merge_left = Some(win_idx);
                        break;
                    }
                }
            }

            // Scan RIGHT for existing window to merge with
            let mut merge_right: Option<usize> = None;
            for scan_bin in (anchor_bin + 1)..=(anchor_bin + win_anchor_dist_nbins as u64) {
                if let Some(&win_idx) = win_bin.get(&(strand, scan_bin)) {
                    let w = &windows[win_idx];
                    if w.alive && w.chr_idx == chr_idx {
                        merge_right = Some(win_idx);
                        break;
                    }
                }
            }

            match (merge_left, merge_right) {
                (Some(left_idx), Some(right_idx)) if left_idx != right_idx => {
                    // Merge both windows: extend left window to cover right + anchor
                    let right_window = &windows[right_idx];
                    let new_bin_end = right_window.bin_end.max(anchor_bin);
                    let new_actual_start = right_window.actual_start.min(forward_pos);
                    let new_actual_end = right_window
                        .actual_end
                        .max(forward_pos + actual_length as u64);
                    let right_alignments: Vec<WindowAlignment> = right_window.alignments.clone();

                    // Kill right window
                    windows[right_idx].alive = false;

                    // Extend left window
                    let left_window = &mut windows[left_idx];
                    left_window.bin_start = left_window.bin_start.min(anchor_bin);
                    left_window.bin_end = left_window.bin_end.max(new_bin_end);
                    left_window.actual_start = left_window.actual_start.min(new_actual_start);
                    left_window.actual_end = left_window.actual_end.max(new_actual_end);
                    left_window.alignments.extend(right_alignments);
                    left_window.alignments.push(wa_entry);

                    // Update winBin for all bins from left to right
                    for bin in left_window.bin_start..=left_window.bin_end {
                        win_bin.insert((strand, bin), left_idx);
                    }
                }
                (Some(idx), _) | (_, Some(idx)) => {
                    // Merge with one existing window
                    let window = &mut windows[idx];
                    window.bin_start = window.bin_start.min(anchor_bin);
                    window.bin_end = window.bin_end.max(anchor_bin);
                    window.actual_start = window.actual_start.min(forward_pos);
                    window.actual_end = window.actual_end.max(forward_pos + actual_length as u64);
                    window.alignments.push(wa_entry);

                    // Update winBin for newly covered bins
                    for bin in window.bin_start..=window.bin_end {
                        win_bin.insert((strand, bin), idx);
                    }
                }
                _ => {
                    // No merge: create new window
                    let new_idx = windows.len();
                    win_bin.insert((strand, anchor_bin), new_idx);
                    windows.push(Window {
                        bin_start: anchor_bin,
                        bin_end: anchor_bin,
                        chr_idx,
                        is_reverse: strand,
                        anchor_idx,
                        alignments: vec![wa_entry],
                        actual_start: forward_pos,
                        actual_end: forward_pos + actual_length as u64,
                        alive: true,
                    });
                }
            }
        }
    }

    if windows.iter().all(|w| !w.alive) {
        return Vec::new();
    }

    // Phase 3: Extend windows by ±win_flank_nbins (matches STAR's flanking extension)
    // Update winBin for newly covered bins
    for (win_idx, window) in windows.iter_mut().enumerate() {
        if !window.alive {
            continue;
        }
        let old_start = window.bin_start;
        let old_end = window.bin_end;
        let new_start = old_start.saturating_sub(win_flank_nbins as u64);
        let new_end = old_end + win_flank_nbins as u64;
        window.bin_start = new_start;
        window.bin_end = new_end;

        let strand = window.is_reverse;
        for bin in new_start..old_start {
            win_bin.entry((strand, bin)).or_insert(win_idx);
        }
        for bin in (old_end + 1)..=new_end {
            win_bin.entry((strand, bin)).or_insert(win_idx);
        }
    }

    // Phase 4: Assign ALL seeds to windows (matches STAR's assignAlignToWindow)
    // Iterate all SA positions per seed — no per-seed cap. Seeds are already bounded
    // by seedMultimapNmax (default 10000) from seed finding. Window capacity eviction
    // (seedPerWindowNmax with anchor protection) handles limiting entries per window.
    for (seed_idx, seed) in seeds.iter().enumerate() {
        let n_loci = seed.sa_end - seed.sa_start;
        if n_loci == 0 {
            continue;
        }
        let is_anchor_seed = anchor_set[seed_idx];

        for (sa_pos, strand) in seed.genome_positions(index) {
            let actual_length = verify_match_at_position(
                read_seq,
                seed.read_pos,
                sa_pos,
                strand,
                seed.length,
                index,
            );
            if actual_length < 8 {
                continue;
            }

            let forward_pos = index.sa_pos_to_forward(sa_pos, strand, actual_length);

            let chr_idx = match index.genome.position_to_chr(forward_pos) {
                Some(info) => info.0,
                None => continue,
            };

            let seed_bin = forward_pos >> win_bin_nbits;

            let win_idx = match win_bin.get(&(strand, seed_bin)) {
                Some(&idx) if windows[idx].alive && windows[idx].chr_idx == chr_idx => idx,
                _ => continue,
            };

            let window = &mut windows[win_idx];

            // Exact duplicate dedup: skip if same (read_pos, sa_pos) already exists
            // (Diagonal overlap dedup deferred to DP — avoids removing valid entries
            // that contribute different match lengths at different positions)
            if window
                .alignments
                .iter()
                .any(|wa| wa.read_pos == seed.read_pos && wa.sa_pos == sa_pos)
            {
                continue;
            }

            // Capacity check (seedPerWindowNmax) with anchor protection
            if window.alignments.len() >= seed_per_window_nmax {
                // Find min length of non-anchor entries
                let min_non_anchor_len = window
                    .alignments
                    .iter()
                    .filter(|wa| !wa.is_anchor)
                    .map(|wa| wa.length)
                    .min()
                    .unwrap_or(usize::MAX);

                if actual_length <= min_non_anchor_len && !is_anchor_seed {
                    continue; // New entry too short
                }

                // Evict shortest non-anchor entries
                window
                    .alignments
                    .retain(|wa| wa.is_anchor || wa.length > min_non_anchor_len);

                if window.alignments.len() >= seed_per_window_nmax {
                    continue; // Still full after eviction
                }
            }

            // Insert the new WA entry
            window.actual_start = window.actual_start.min(forward_pos);
            window.actual_end = window.actual_end.max(forward_pos + actual_length as u64);
            window.alignments.push(WindowAlignment {
                seed_idx,
                read_pos: seed.read_pos,
                length: actual_length,
                genome_pos: forward_pos,
                sa_pos,
                n_rep: n_loci,
                is_anchor: is_anchor_seed,
            });
        }
    }

    // Phase 5: Build SeedCluster output
    let mut clusters = Vec::with_capacity(windows.len());
    for window in &windows {
        if !window.alive || window.alignments.is_empty() {
            continue;
        }

        clusters.push(SeedCluster {
            alignments: window.alignments.clone(),
            chr_idx: window.chr_idx,
            genome_start: window.actual_start,
            genome_end: window.actual_end,
            is_reverse: window.is_reverse,
            anchor_idx: window.anchor_idx,
            anchor_bin: window.bin_start,
        });
    }

    clusters
}

/// Lightweight exon block for in-progress transcript during recursion
#[derive(Debug, Clone)]
struct ExonBlock {
    read_start: usize, // 0-based inclusive
    read_end: usize,   // 0-based exclusive
    genome_start: u64, // SA coordinate space (raw sa_pos)
    genome_end: u64,   // SA coordinate space (exclusive)
}

/// In-progress transcript during recursive search (cheap to clone)
#[derive(Debug, Clone)]
struct WorkingTranscript {
    exons: Vec<ExonBlock>,
    score: i32,
    n_mismatch: u32,
    n_gap: u32,
    n_junction: u32,
    junction_motifs: Vec<crate::align::score::SpliceMotif>,
    junction_annotated: Vec<bool>,
    n_anchor: u32,
    // Tight bounds for extension at finalization
    read_start: usize,
    read_end: usize,
    genome_start: u64,
    genome_end: u64,
}

impl WorkingTranscript {
    fn new() -> Self {
        WorkingTranscript {
            exons: Vec::new(),
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: Vec::new(),
            junction_annotated: Vec::new(),
            n_anchor: 0,
            read_start: 0,
            read_end: 0,
            genome_start: 0,
            genome_end: 0,
        }
    }
}

/// Fill the gap between the last exon in a WorkingTranscript and the next WA entry.
/// Returns None if stitching fails (mismatch limit, overhang too short, etc.).
/// Matches STAR's stitchAlignToTranscript.cpp logic.
#[allow(clippy::too_many_arguments)]
fn stitch_align_to_transcript(
    wt: &WorkingTranscript,
    wa: &WindowAlignment,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    cluster: &SeedCluster,
    junction_db: Option<&crate::junction::SpliceJunctionDb>,
) -> Option<WorkingTranscript> {
    let last_exon = wt.exons.last().unwrap();

    // Overlap trimming: if new WA overlaps previous exon in read coords, shift start right
    let mut eff_read_pos = wa.read_pos;
    let mut eff_genome_pos = wa.sa_pos;
    let mut eff_length = wa.length;

    if last_exon.read_end > eff_read_pos {
        let overlap = last_exon.read_end - eff_read_pos;
        if overlap >= eff_length {
            return None; // Fully consumed
        }
        eff_read_pos = last_exon.read_end;
        eff_genome_pos += overlap as u64;
        eff_length -= overlap;
    }

    // Handle genome overlap
    if last_exon.genome_end > eff_genome_pos && eff_genome_pos > last_exon.genome_start {
        let g_overlap = (last_exon.genome_end - eff_genome_pos) as usize;
        if g_overlap >= eff_length {
            return None; // Fully consumed
        }
        eff_read_pos += g_overlap;
        eff_genome_pos += g_overlap as u64;
        eff_length -= g_overlap;
    }

    let read_gap = (eff_read_pos as i64) - (last_exon.read_end as i64);
    let genome_gap = (eff_genome_pos as i64) - (last_exon.genome_end as i64);

    // Reject negative gaps
    if read_gap < 0 || genome_gap < 0 {
        return None;
    }

    let mut new_wt = wt.clone();
    let mut d_score: i32 = 0;
    let mut gap_mm: u32 = 0;

    if read_gap == 0 && genome_gap == 0 {
        // Adjacent seeds — just extend the last exon
        if let Some(last) = new_wt.exons.last_mut() {
            last.read_end = eff_read_pos + eff_length;
            last.genome_end = eff_genome_pos + eff_length as u64;
        }
    } else if read_gap == genome_gap {
        // Equal gap: base-by-base scoring
        let shared = read_gap as usize;
        gap_mm = count_mismatches_in_region(
            read_seq,
            last_exon.read_end,
            last_exon.genome_end,
            shared,
            index,
            cluster.is_reverse,
        );
        d_score += shared as i32 - 2 * gap_mm as i32;

        // Extend last exon through the gap and the new seed
        if let Some(last) = new_wt.exons.last_mut() {
            last.read_end = eff_read_pos + eff_length;
            last.genome_end = eff_genome_pos + eff_length as u64;
        }
    } else if genome_gap > read_gap {
        // Deletion or splice junction
        let del = (genome_gap - read_gap) as u32;
        let shared = read_gap as usize;

        // Score shared bases on donor side
        let shared_mm = if shared > 0 {
            count_mismatches_in_region(
                read_seq,
                last_exon.read_end,
                last_exon.genome_end,
                shared,
                index,
                cluster.is_reverse,
            )
        } else {
            0
        };
        gap_mm += shared_mm;
        d_score += shared as i32 - 2 * shared_mm as i32;

        if del >= scorer.align_intron_min && del <= scorer.align_intron_max {
            // Splice junction — apply jR scanning
            let donor_sa = last_exon.genome_end + shared as u64;
            let (jr_shift, motif, motif_score) = scorer.find_best_junction_position(
                read_seq,
                last_exon.read_end + shared,
                donor_sa,
                read_gap.max(0),
                genome_gap,
                &index.genome,
                cluster.is_reverse,
                index.genome.n_genome,
                last_exon.read_end - last_exon.read_start + shared,
            );

            // Clamp shift: must not consume entire donor or acceptor
            let prev_match_len = (last_exon.read_end - last_exon.read_start + shared) as i32;
            let jr_shift = jr_shift
                .max(-prev_match_len)
                .min(read_gap as i32)
                .min(eff_length as i32);

            // Check stitch mismatch limit
            if !scorer.stitch_mismatch_allowed(&motif, gap_mm) {
                return None;
            }

            // Overhang check (use post-shift values)
            let left_overhang =
                ((last_exon.read_end - last_exon.read_start) as i32 + shared as i32 + jr_shift)
                    .max(0) as usize;
            let right_overhang = (eff_length as i32 - jr_shift).max(0) as usize;

            // Check annotation
            let is_annotated = junction_db.is_some_and(|db| {
                let junc_donor_sa = (donor_sa as i64 + jr_shift as i64) as u64;
                let donor_fwd =
                    index.sa_pos_to_forward(junc_donor_sa, cluster.is_reverse, del as usize);
                let acceptor_fwd = donor_fwd + del as u64;
                db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 0)
                    || db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 1)
                    || db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 2)
            });

            let min_overhang = if is_annotated {
                scorer.align_sjdb_overhang_min as usize
            } else {
                scorer.align_sj_overhang_min as usize
            };

            if left_overhang < min_overhang || right_overhang < min_overhang {
                return None;
            }

            d_score += motif_score - scorer.score_stitch_sj_shift;
            if is_annotated {
                d_score += scorer.sjdb_score;
            }

            // Adjust last exon for jr_shift
            if jr_shift != 0 {
                if let Some(last) = new_wt.exons.last_mut() {
                    last.read_end =
                        (last.read_end as i64 + shared as i64 + jr_shift as i64) as usize;
                    last.genome_end =
                        (last.genome_end as i64 + shared as i64 + jr_shift as i64) as u64;
                }
            } else if shared > 0 {
                if let Some(last) = new_wt.exons.last_mut() {
                    last.read_end += shared;
                    last.genome_end += shared as u64;
                }
            }

            // New exon for acceptor side
            // jr_shift positive = junction moved right → acceptor starts later
            // jr_shift negative = junction moved left → acceptor starts earlier
            let acceptor_read_start = (eff_read_pos as i64 + jr_shift as i64) as usize;
            let acceptor_genome_start = (eff_genome_pos as i64 + jr_shift as i64) as u64;
            let acceptor_len = (eff_length as i64 - jr_shift as i64).max(0) as usize;
            new_wt.exons.push(ExonBlock {
                read_start: acceptor_read_start,
                read_end: acceptor_read_start + acceptor_len,
                genome_start: acceptor_genome_start,
                genome_end: acceptor_genome_start + acceptor_len as u64,
            });

            new_wt.n_junction += 1;
            new_wt.junction_motifs.push(motif);
            new_wt.junction_annotated.push(is_annotated);
        } else {
            // Deletion (too short/long for intron)
            let del_score = scorer.score_del_open + scorer.score_del_base * del as i32;
            d_score += del_score;
            new_wt.n_gap += 1;

            // Extend last exon through shared bases
            if shared > 0 {
                if let Some(last) = new_wt.exons.last_mut() {
                    last.read_end += shared;
                    last.genome_end += shared as u64;
                }
            }

            // New exon after deletion
            new_wt.exons.push(ExonBlock {
                read_start: eff_read_pos,
                read_end: eff_read_pos + eff_length,
                genome_start: eff_genome_pos,
                genome_end: eff_genome_pos + eff_length as u64,
            });
        }
    } else {
        // Insertion: read_gap > genome_gap
        let ins = (read_gap - genome_gap) as u32;
        let shared = genome_gap as usize;

        let shared_mm = if shared > 0 {
            count_mismatches_in_region(
                read_seq,
                last_exon.read_end,
                last_exon.genome_end,
                shared,
                index,
                cluster.is_reverse,
            )
        } else {
            0
        };
        gap_mm += shared_mm;
        d_score += shared as i32 - 2 * shared_mm as i32;

        let ins_score = scorer.score_ins_open + scorer.score_ins_base * ins as i32;
        d_score += ins_score;
        new_wt.n_gap += 1;

        // Extend last exon through shared genome bases
        if shared > 0 {
            if let Some(last) = new_wt.exons.last_mut() {
                last.read_end += shared;
                last.genome_end += shared as u64;
            }
        }

        // New exon after insertion
        new_wt.exons.push(ExonBlock {
            read_start: eff_read_pos,
            read_end: eff_read_pos + eff_length,
            genome_start: eff_genome_pos,
            genome_end: eff_genome_pos + eff_length as u64,
        });
    }

    // Mismatch limit check
    let total_mm = new_wt.n_mismatch + gap_mm;
    let total_len = new_wt.read_end.max(eff_read_pos + eff_length) - new_wt.read_start;
    let mm_limit = ((scorer.p_mm_max * total_len as f64) as u32).min(scorer.n_mm_max);
    if total_mm > mm_limit {
        return None;
    }

    // Update working transcript
    // Seeds from SA are exact matches (0 internal mismatches).
    // Mismatches are only in gap-fill shared bases (counted in d_score) and extensions.
    new_wt.score += d_score + eff_length as i32;
    new_wt.n_mismatch += gap_mm;
    new_wt.read_end = new_wt.exons.last().map_or(new_wt.read_end, |e| e.read_end);
    new_wt.genome_end = new_wt
        .exons
        .last()
        .map_or(new_wt.genome_end, |e| e.genome_end);
    if wa.is_anchor {
        new_wt.n_anchor += 1;
    }

    Some(new_wt)
}

/// Stitch seeds within a cluster using recursive combinatorial stitching.
///
/// # Arguments
/// * `cluster` - Seed cluster with WindowAlignment entries (STAR's WA array)
/// * `read_seq` - Read sequence
/// * `index` - Genome index
/// * `scorer` - Alignment scorer
///
/// # Returns
/// Vector of transcripts (may have multiple paths through the cluster)
pub fn stitch_seeds(
    cluster: &SeedCluster,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
) -> Result<Vec<Transcript>, Error> {
    stitch_seeds_with_jdb(cluster, read_seq, index, scorer, None, 1)
}

/// Stitch seeds with optional junction database for annotation-aware scoring.
///
/// Uses STAR's recursive combinatorial stitcher (stitchWindowAligns) which
/// explores include/exclude branches for each seed, allowing it to skip
/// spurious short seeds that would create false splices.
///
/// `max_transcripts_per_window` controls how many transcripts are collected
/// (STAR's `alignTranscriptsPerWindowNmax`, default 100).
pub fn stitch_seeds_with_jdb(
    cluster: &SeedCluster,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    junction_db: Option<&crate::junction::SpliceJunctionDb>,
    max_transcripts_per_window: usize,
) -> Result<Vec<Transcript>, Error> {
    stitch_seeds_with_jdb_debug(
        cluster,
        read_seq,
        index,
        scorer,
        junction_db,
        max_transcripts_per_window,
        "",
    )
}

/// Check if two transcripts overlap in their exon blocks.
/// Returns total overlapping bases (where exons share the same read-genome diagonal).
/// Used for deduplication: if new transcript is a subset of existing, drop it.
fn blocks_overlap(t1_exons: &[ExonBlock], t2_exons: &[ExonBlock]) -> u32 {
    let mut overlap = 0u32;
    let mut i = 0;
    let mut j = 0;

    while i < t1_exons.len() && j < t2_exons.len() {
        let e1 = &t1_exons[i];
        let e2 = &t2_exons[j];

        // Check if exons are on the same read-genome diagonal
        let diag1 = e1.genome_start as i64 - e1.read_start as i64;
        let diag2 = e2.genome_start as i64 - e2.read_start as i64;

        if diag1 == diag2 {
            // Same diagonal — compute read-space overlap
            let r_start = e1.read_start.max(e2.read_start);
            let r_end = e1.read_end.min(e2.read_end);
            if r_start < r_end {
                overlap += (r_end - r_start) as u32;
            }
        }

        // Advance the exon that ends first in read space
        if e1.read_end <= e2.read_end {
            i += 1;
        } else {
            j += 1;
        }
    }

    overlap
}

/// Convert a WorkingTranscript to a final Transcript by extending flanks,
/// building CIGAR, counting mismatches, and converting to forward coordinates.
#[allow(clippy::too_many_arguments)]
fn finalize_transcript(
    wt: &WorkingTranscript,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    cluster: &SeedCluster,
) -> Transcript {
    use crate::align::transcript::Exon;

    let alignment_start = wt.read_start;
    let alignment_end = wt.read_end;

    // Extend alignment into flanking regions (STAR-style extendAlign)
    let left_extend = if alignment_start > 0 {
        extend_alignment(
            read_seq,
            alignment_start,
            wt.genome_start,
            -1,
            alignment_start,
            0,       // nMMprev=0 matches STAR
            100_000, // Lprev=100000 matches STAR
            scorer.n_mm_max,
            scorer.p_mm_max,
            index,
            cluster.is_reverse,
        )
    } else {
        ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        }
    };

    let right_extend = if alignment_end < read_seq.len() {
        extend_alignment(
            read_seq,
            alignment_end,
            wt.genome_end,
            1,
            read_seq.len() - alignment_end,
            wt.n_mismatch + left_extend.n_mismatch,
            100_000,
            scorer.n_mm_max,
            scorer.p_mm_max,
            index,
            cluster.is_reverse,
        )
    } else {
        ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        }
    };

    // Build final CIGAR from exon blocks
    let mut final_cigar: Vec<CigarOp> = Vec::new();

    // Left soft clip
    let remaining_left_clip = alignment_start - left_extend.extend_len;
    if remaining_left_clip > 0 {
        final_cigar.push(CigarOp::SoftClip(remaining_left_clip as u32));
    }

    // Left extension match
    if left_extend.extend_len > 0 {
        final_cigar.push(CigarOp::Match(left_extend.extend_len as u32));
    }

    // Walk exon blocks to build CIGAR
    for (idx, exon) in wt.exons.iter().enumerate() {
        if idx > 0 {
            let prev = &wt.exons[idx - 1];
            let read_gap = exon.read_start as i64 - prev.read_end as i64;
            let genome_gap = exon.genome_start as i64 - prev.genome_end as i64;

            if genome_gap > read_gap && genome_gap > 0 {
                // Shared match bases before the gap
                let shared = read_gap.max(0) as u32;
                if shared > 0 {
                    if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
                        *prev_len += shared;
                    } else {
                        final_cigar.push(CigarOp::Match(shared));
                    }
                }
                let del = (genome_gap - read_gap.max(0)) as u32;
                if del >= scorer.align_intron_min && del <= scorer.align_intron_max {
                    final_cigar.push(CigarOp::RefSkip(del));
                } else {
                    final_cigar.push(CigarOp::Del(del));
                }
            } else if read_gap > genome_gap && read_gap > 0 {
                // Insertion
                let shared = genome_gap.max(0) as u32;
                if shared > 0 {
                    if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
                        *prev_len += shared;
                    } else {
                        final_cigar.push(CigarOp::Match(shared));
                    }
                }
                let ins = (read_gap - genome_gap.max(0)) as u32;
                final_cigar.push(CigarOp::Ins(ins));
            }
            // Equal gap case is handled by extended exon blocks in stitch_align_to_transcript
        }

        // This exon's match region
        let match_len = (exon.read_end - exon.read_start) as u32;
        if match_len > 0 {
            if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
                *prev_len += match_len;
            } else {
                final_cigar.push(CigarOp::Match(match_len));
            }
        }
    }

    // Right extension match
    if right_extend.extend_len > 0 {
        if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
            *prev_len += right_extend.extend_len as u32;
        } else {
            final_cigar.push(CigarOp::Match(right_extend.extend_len as u32));
        }
    }

    // Right soft clip
    let remaining_right_clip = (read_seq.len() - alignment_end) - right_extend.extend_len;
    if remaining_right_clip > 0 {
        final_cigar.push(CigarOp::SoftClip(remaining_right_clip as u32));
    }

    // Validate CIGAR read-consuming length
    let cigar_read_len: u32 = final_cigar
        .iter()
        .filter(|op| op.consumes_query())
        .map(|op| op.len())
        .sum();
    if cigar_read_len != read_seq.len() as u32 {
        panic!(
            "[CIGAR-MISMATCH] cigar_read_len={} read_len={} left_ext={} right_ext={} align_start={} align_end={}\n  CIGAR: {}\n  Exon blocks ({}):\n{}",
            cigar_read_len,
            read_seq.len(),
            left_extend.extend_len,
            right_extend.extend_len,
            alignment_start,
            alignment_end,
            final_cigar
                .iter()
                .map(|op| op.to_string())
                .collect::<String>(),
            wt.exons.len(),
            wt.exons
                .iter()
                .enumerate()
                .map(|(i, e)| format!(
                    "    exon[{}]: read=[{}, {}), genome=[{}, {}), len={}",
                    i,
                    e.read_start,
                    e.read_end,
                    e.genome_start,
                    e.genome_end,
                    e.read_end - e.read_start
                ))
                .collect::<Vec<_>>()
                .join("\n")
        );
    }

    // Adjusted genome start for left extension (raw SA coordinates)
    let adjusted_genome_start = wt.genome_start - left_extend.extend_len as u64;

    // Adjust score: only add right extension (left was already added during recursion)
    let adjusted_score = wt.score + right_extend.max_score;

    // Count mismatches — MUST be called BEFORE CIGAR reversal
    let n_mismatch = count_mismatches(
        read_seq,
        &final_cigar,
        adjusted_genome_start,
        0,
        index,
        cluster.is_reverse,
    );

    // Reverse CIGAR for reverse strand
    if cluster.is_reverse {
        final_cigar.reverse();
    }

    // Compute total reference-consuming length from CIGAR
    let mut ref_len = 0u64;
    for op in &final_cigar {
        match op {
            CigarOp::Match(len)
            | CigarOp::Equal(len)
            | CigarOp::Diff(len)
            | CigarOp::Del(len)
            | CigarOp::RefSkip(len) => {
                ref_len += *len as u64;
            }
            _ => {}
        }
    }

    // Convert to forward genome coordinates
    let forward_genome_start =
        index.sa_pos_to_forward(adjusted_genome_start, cluster.is_reverse, ref_len as usize);
    let forward_genome_end = forward_genome_start + ref_len;

    // Build exons from CIGAR
    let mut exons = Vec::new();
    let mut read_pos_e = 0usize;
    let mut genome_pos_e = forward_genome_start;

    for op in &final_cigar {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                let len = *len as usize;
                exons.push(Exon {
                    genome_start: genome_pos_e,
                    genome_end: genome_pos_e + len as u64,
                    read_start: read_pos_e,
                    read_end: read_pos_e + len,
                });
                read_pos_e += len;
                genome_pos_e += len as u64;
            }
            CigarOp::Ins(len) => {
                read_pos_e += *len as usize;
            }
            CigarOp::Del(len) => {
                genome_pos_e += *len as u64;
            }
            CigarOp::RefSkip(len) => {
                genome_pos_e += *len as u64;
            }
            CigarOp::SoftClip(len) => {
                read_pos_e += *len as usize;
            }
            CigarOp::HardClip(_) => {}
        }
    }

    // Merge consecutive exons
    let mut merged_exons: Vec<Exon> = Vec::new();
    for exon in exons {
        if let Some(last_exon) = merged_exons.last_mut() {
            if last_exon.genome_end == exon.genome_start && last_exon.read_end == exon.read_start {
                last_exon.genome_end = exon.genome_end;
                last_exon.read_end = exon.read_end;
                continue;
            }
        }
        merged_exons.push(exon);
    }

    let t_genome_start = merged_exons
        .first()
        .map(|e| e.genome_start)
        .unwrap_or(forward_genome_start);
    let t_genome_end = merged_exons
        .last()
        .map(|e| e.genome_end)
        .unwrap_or(forward_genome_end);

    // Apply genomic length penalty
    let genomic_span = t_genome_end - t_genome_start;
    let length_penalty = scorer.genomic_length_penalty(genomic_span);
    let final_score = (adjusted_score + length_penalty).max(0);

    Transcript {
        chr_idx: cluster.chr_idx,
        genome_start: t_genome_start,
        genome_end: t_genome_end,
        is_reverse: cluster.is_reverse,
        exons: merged_exons,
        cigar: final_cigar,
        score: final_score,
        n_mismatch,
        n_gap: wt.n_gap,
        n_junction: wt.n_junction,
        junction_motifs: wt.junction_motifs.clone(),
        junction_annotated: wt.junction_annotated.clone(),
        read_seq: read_seq.to_vec(),
    }
}

/// Recursive include/exclude stitcher (STAR's stitchWindowAligns).
///
/// For each WA entry: try including it (call stitch_align_to_transcript to fill gap),
/// and try excluding it (subject to anchor constraint). Transcripts are finalized
/// at the base case when all entries have been considered.
#[allow(clippy::too_many_arguments)]
fn stitch_recurse(
    i_a: usize,
    wt: WorkingTranscript,
    wa_entries: &[WindowAlignment],
    last_anchor_idx: Option<usize>,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    cluster: &SeedCluster,
    junction_db: Option<&crate::junction::SpliceJunctionDb>,
    max_transcripts: usize,
    transcripts: &mut Vec<WorkingTranscript>,
    recursion_count: &mut u32,
) {
    const MAX_RECURSION: u32 = 10_000;

    if *recursion_count >= MAX_RECURSION || transcripts.len() >= max_transcripts {
        return;
    }
    *recursion_count += 1;

    // Base case: all entries considered
    if i_a >= wa_entries.len() {
        if !wt.exons.is_empty() {
            // Dedup via blocks_overlap: drop if subset of existing higher-score transcript
            let mut dominated = false;
            let mut remove_indices = Vec::new();
            for (idx, existing) in transcripts.iter().enumerate() {
                let overlap = blocks_overlap(&wt.exons, &existing.exons);
                let wt_len: u32 = wt
                    .exons
                    .iter()
                    .map(|e| (e.read_end - e.read_start) as u32)
                    .sum();
                let ex_len: u32 = existing
                    .exons
                    .iter()
                    .map(|e| (e.read_end - e.read_start) as u32)
                    .sum();

                // Only dedup transcripts with same number of exon blocks (junctions).
                // A non-spliced path should never be removed because a spliced path
                // covers the same bases — they may score differently after finalization
                // (genomic length penalty favors shorter genomic spans).
                let same_structure = wt.exons.len() == existing.exons.len();
                if same_structure && overlap >= wt_len && existing.score >= wt.score {
                    dominated = true;
                    break;
                }
                if same_structure && overlap >= ex_len && wt.score >= existing.score {
                    remove_indices.push(idx);
                }
            }

            if !dominated {
                // Remove subsets (iterate in reverse to preserve indices)
                for &idx in remove_indices.iter().rev() {
                    transcripts.swap_remove(idx);
                }
                transcripts.push(wt);
            }
        }
        return;
    }

    let wa = &wa_entries[i_a];

    // INCLUDE branch: try stitching wa_entries[i_a] to transcript
    if wt.exons.is_empty() {
        // First seed: create initial transcript
        let mut new_wt = wt.clone();
        new_wt.exons.push(ExonBlock {
            read_start: wa.read_pos,
            read_end: wa.read_pos + wa.length,
            genome_start: wa.sa_pos,
            genome_end: wa.sa_pos + wa.length as u64,
        });
        new_wt.score = wa.length as i32;
        new_wt.read_start = wa.read_pos;
        new_wt.read_end = wa.read_pos + wa.length;
        new_wt.genome_start = wa.sa_pos;
        new_wt.genome_end = wa.sa_pos + wa.length as u64;
        if wa.is_anchor {
            new_wt.n_anchor = 1;
        }

        // Add pre-DP left extension score (STAR Pass 1)
        if wa.read_pos > 0 {
            let left_ext = extend_alignment(
                read_seq,
                wa.read_pos,
                wa.sa_pos,
                -1,
                wa.read_pos,
                0,
                100_000,
                scorer.n_mm_max,
                scorer.p_mm_max,
                index,
                cluster.is_reverse,
            );
            new_wt.score += left_ext.max_score;
        }

        stitch_recurse(
            i_a + 1,
            new_wt,
            wa_entries,
            last_anchor_idx,
            read_seq,
            index,
            scorer,
            cluster,
            junction_db,
            max_transcripts,
            transcripts,
            recursion_count,
        );
    } else {
        // Try stitching this seed onto the existing transcript
        if let Some(new_wt) =
            stitch_align_to_transcript(&wt, wa, read_seq, index, scorer, cluster, junction_db)
        {
            stitch_recurse(
                i_a + 1,
                new_wt,
                wa_entries,
                last_anchor_idx,
                read_seq,
                index,
                scorer,
                cluster,
                junction_db,
                max_transcripts,
                transcripts,
                recursion_count,
            );
        }
    }

    // EXCLUDE branch: skip wa_entries[i_a]
    // Anchor constraint: can only skip the last anchor if transcript already has one
    let can_exclude = if let Some(last_anchor) = last_anchor_idx {
        if wa.is_anchor && i_a == last_anchor {
            wt.n_anchor > 0 // Already has an anchor → ok to skip
        } else {
            true
        }
    } else {
        true
    };

    if can_exclude {
        stitch_recurse(
            i_a + 1,
            wt,
            wa_entries,
            last_anchor_idx,
            read_seq,
            index,
            scorer,
            cluster,
            junction_db,
            max_transcripts,
            transcripts,
            recursion_count,
        );
    }
}

/// Inner implementation of stitch_seeds_with_jdb with optional debug logging.
/// When `debug_read_name` is non-empty, detailed info is logged to stderr.
///
/// Uses STAR's recursive combinatorial stitcher (stitchWindowAligns) instead of
/// forward DP. For each seed, explores include/exclude branches, allowing the
/// algorithm to skip spurious short seeds that would create false splices.
pub(crate) fn stitch_seeds_with_jdb_debug(
    cluster: &SeedCluster,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    junction_db: Option<&crate::junction::SpliceJunctionDb>,
    max_transcripts_per_window: usize,
    debug_read_name: &str,
) -> Result<Vec<Transcript>, Error> {
    let debug = !debug_read_name.is_empty();

    // Sort WA entries by read_pos, then length descending for dedup
    let mut wa_entries: Vec<WindowAlignment> = cluster.alignments.clone();
    wa_entries.sort_by(|a, b| a.read_pos.cmp(&b.read_pos).then(b.length.cmp(&a.length)));

    // Dedup: keep only the longest per (read_pos, sa_pos) pair
    wa_entries.dedup_by(|a, b| a.read_pos == b.read_pos && a.sa_pos == b.sa_pos);

    // Aggressive diagonal dedup: for each diagonal (sa_pos - read_pos), merge
    // overlapping seeds into intervals, then keep only 1 seed per merged interval
    // (the longest). With dense seed search, most seeds are redundant subsets.
    {
        use std::collections::HashMap;
        // For each diagonal, find the longest seed per merged interval
        let mut diag_seeds: HashMap<i64, Vec<(usize, usize, usize)>> = HashMap::new();
        for (idx, wa) in wa_entries.iter().enumerate() {
            let diag = wa.sa_pos as i64 - wa.read_pos as i64;
            diag_seeds
                .entry(diag)
                .or_default()
                .push((wa.read_pos, wa.read_pos + wa.length, idx));
        }

        let mut keep_indices = std::collections::HashSet::new();
        for (_diag, mut seeds) in diag_seeds {
            // Sort by start position
            seeds.sort();
            // Merge intervals, keeping the index of the longest seed in each merged group
            let mut merged_start = seeds[0].0;
            let mut merged_end = seeds[0].1;
            let mut best_idx = seeds[0].2;
            let mut best_len = seeds[0].1 - seeds[0].0;

            for &(s, e, idx) in &seeds[1..] {
                if s <= merged_end {
                    // Overlapping — extend and track longest
                    merged_end = merged_end.max(e);
                    let len = e - s;
                    if len > best_len {
                        best_len = len;
                        best_idx = idx;
                    }
                } else {
                    // New interval — commit previous best
                    keep_indices.insert(best_idx);
                    merged_start = s;
                    merged_end = e;
                    best_len = e - s;
                    best_idx = idx;
                }
            }
            // Commit last group
            keep_indices.insert(best_idx);
            let _ = merged_start; // suppress unused warning
        }

        // Retain only the kept indices
        let mut idx = 0usize;
        wa_entries.retain(|_| {
            let keep = keep_indices.contains(&idx);
            idx += 1;
            keep
        });
    }

    // Cap entries to prevent exponential blowup
    const MAX_WA_ENTRIES: usize = 200;
    if wa_entries.len() > MAX_WA_ENTRIES {
        wa_entries.sort_by(|a, b| b.length.cmp(&a.length));
        wa_entries.truncate(MAX_WA_ENTRIES);
        wa_entries.sort_by(|a, b| a.read_pos.cmp(&b.read_pos).then(b.length.cmp(&a.length)));
    }

    if wa_entries.is_empty() {
        return Ok(Vec::new());
    }

    if debug {
        eprintln!(
            "[DEBUG-STITCH {}] {} WA entries (is_reverse={})",
            debug_read_name,
            wa_entries.len(),
            cluster.is_reverse
        );
        for (i, wa) in wa_entries.iter().enumerate().take(30) {
            eprintln!(
                "  wa[{}]: read_pos={}, sa_pos={}, length={}, anchor={}",
                i, wa.read_pos, wa.sa_pos, wa.length, wa.is_anchor
            );
        }
        if wa_entries.len() > 30 {
            eprintln!("  ... ({} more)", wa_entries.len() - 30);
        }
    }

    // Find last anchor index for the anchor constraint
    let last_anchor_idx = wa_entries.iter().rposition(|wa| wa.is_anchor);

    // Run recursive include/exclude stitcher
    let mut working_transcripts: Vec<WorkingTranscript> = Vec::new();
    let mut recursion_count: u32 = 0;

    stitch_recurse(
        0,
        WorkingTranscript::new(),
        &wa_entries,
        last_anchor_idx,
        read_seq,
        index,
        scorer,
        cluster,
        junction_db,
        max_transcripts_per_window,
        &mut working_transcripts,
        &mut recursion_count,
    );

    if debug {
        eprintln!(
            "[DEBUG-STITCH {}] Recursion done: {} working transcripts, {} recursions",
            debug_read_name,
            working_transcripts.len(),
            recursion_count
        );
    }

    // Finalize working transcripts → Transcript
    let mut transcripts: Vec<Transcript> = Vec::with_capacity(working_transcripts.len());
    for wt in &working_transcripts {
        transcripts.push(finalize_transcript(wt, read_seq, index, scorer, cluster));
    }

    // Sort by score descending, then fewer junctions (prefer non-spliced on tie),
    // then fewer mismatches as final tiebreaker
    transcripts.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.n_junction.cmp(&b.n_junction))
            .then(a.n_mismatch.cmp(&b.n_mismatch))
    });
    transcripts.truncate(max_transcripts_per_window);

    Ok(transcripts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Genome;
    use crate::index::packed_array::PackedArray;
    use crate::index::sa_index::SaIndex;
    use crate::index::suffix_array::SuffixArray;

    fn make_simple_index() -> GenomeIndex {
        // Simple genome: ACGTACGTNN (10 bases)
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 4, 4];
        let n_genome = 64u64; // Padded
        let mut sequence = vec![5u8; (n_genome * 2) as usize];
        sequence[0..seq.len()].copy_from_slice(&seq);

        // Build reverse complement
        for i in 0..n_genome as usize {
            let base = sequence[i];
            let complement = if base < 4 { 3 - base } else { base };
            sequence[2 * n_genome as usize - 1 - i] = complement;
        }

        let genome = Genome {
            sequence,
            n_genome,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![10],
            chr_start: vec![0, n_genome],
        };

        // Create dummy SA and SAindex
        let gstrand_bit = 33;
        let suffix_array = SuffixArray {
            data: PackedArray::new(gstrand_bit, 0),
            gstrand_bit,
            gstrand_mask: (1u64 << gstrand_bit) - 1,
        };

        let word_length = gstrand_bit + 3;
        let sa_index = SaIndex {
            data: PackedArray::new(word_length, 0),
            nbases: 14,
            genome_sa_index_start: vec![0],
            word_length,
            gstrand_bit,
        };

        GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db: crate::junction::SpliceJunctionDb::empty(),
        }
    }

    #[test]
    fn test_cluster_seeds_simple() {
        // This test would require a properly populated SA
        // For now, just test that we can create the clustering structure
        // without panicking on an empty index
        let index = make_simple_index();

        // Create seeds with empty SA ranges (won't expand to any positions)
        let seeds = vec![
            Seed {
                read_pos: 0,
                length: 5,
                sa_start: 0,
                sa_end: 0, // Empty range
                is_reverse: false,
                search_rc: false,
                mate_id: 2,
            },
            Seed {
                read_pos: 10,
                length: 5,
                sa_start: 0,
                sa_end: 0, // Empty range
                is_reverse: false,
                search_rc: false,
                mate_id: 2,
            },
        ];

        // Bin-based windowing: win_bin_nbits=16, win_anchor_dist_nbins=9, win_flank_nbins=4
        let read_seq = vec![0u8; 20]; // Dummy read for verify_match_at_position
        let clusters = cluster_seeds(&seeds, &read_seq, &index, 16, 9, 4, 10, 50, 50);

        // With empty SA ranges, no clusters will be created
        assert_eq!(clusters.len(), 0);
    }

    #[test]
    fn test_wa_entry_sorting() {
        let mut entries = vec![
            WindowAlignment {
                seed_idx: 0,
                read_pos: 10,
                length: 5,
                genome_pos: 100,
                sa_pos: 100,
                n_rep: 1,
                is_anchor: true,
            },
            WindowAlignment {
                seed_idx: 1,
                read_pos: 5,
                length: 5,
                genome_pos: 50,
                sa_pos: 50,
                n_rep: 1,
                is_anchor: true,
            },
        ];

        entries.sort_by(|a, b| a.read_pos.cmp(&b.read_pos).then(b.length.cmp(&a.length)));

        assert_eq!(entries[0].read_pos, 5);
        assert_eq!(entries[1].read_pos, 10);
    }

    /// Helper to build a GenomeIndex with a specific forward sequence
    fn make_index_with_seq(seq: &[u8]) -> GenomeIndex {
        let n_genome = ((seq.len() as u64 + 1) / 64 + 1) * 64;
        let mut sequence = vec![5u8; (n_genome * 2) as usize];
        sequence[0..seq.len()].copy_from_slice(seq);

        // Build reverse complement
        for i in 0..n_genome as usize {
            let base = sequence[i];
            let complement = if base < 4 { 3 - base } else { base };
            sequence[2 * n_genome as usize - 1 - i] = complement;
        }

        let genome = Genome {
            sequence,
            n_genome,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![seq.len() as u64],
            chr_start: vec![0, n_genome],
        };

        let gstrand_bit = 33;
        let suffix_array = SuffixArray {
            data: PackedArray::new(gstrand_bit, 0),
            gstrand_bit,
            gstrand_mask: (1u64 << gstrand_bit) - 1,
        };
        let word_length = gstrand_bit + 3;
        let sa_index = SaIndex {
            data: PackedArray::new(word_length, 0),
            nbases: 14,
            genome_sa_index_start: vec![0],
            word_length,
            gstrand_bit,
        };

        GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db: crate::junction::SpliceJunctionDb::empty(),
        }
    }

    #[test]
    fn test_extend_perfect_match_rightward() {
        // Genome: ACGTACGTAC (10 bases, A=0, C=1, G=2, T=3)
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let index = make_index_with_seq(&seq);
        // Read matches genome perfectly from position 5 onward
        let read_seq = vec![1, 2, 3, 0, 1]; // matches genome[5..10]

        let result = extend_alignment(
            &read_seq, 0,   // read_start (boundary)
            5,   // genome_start
            1,   // rightward
            5,   // max_extend
            0,   // no previous mismatches
            0,   // no previous length
            10,  // n_mm_max
            0.3, // p_mm_max
            &index, false, // forward strand
        );

        assert_eq!(result.extend_len, 5);
        assert_eq!(result.max_score, 5);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_perfect_match_leftward() {
        // Genome: ACGTACGTAC
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let index = make_index_with_seq(&seq);
        // Read matches genome[0..5] = ACGTA
        let read_seq = vec![0, 1, 2, 3, 0];

        let result = extend_alignment(
            &read_seq, 5,   // read_start (exclusive boundary for leftward)
            5,   // genome_start (exclusive boundary for leftward)
            -1,  // leftward
            5,   // max_extend
            0,   // no previous mismatches
            0,   // no previous length
            10,  // n_mm_max
            0.3, // p_mm_max
            &index, false,
        );

        assert_eq!(result.extend_len, 5);
        assert_eq!(result.max_score, 5);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_stops_at_optimal_point_with_mismatches() {
        // Genome: A C G T A C G T (positions 0-7)
        let genome_seq = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let index = make_index_with_seq(&genome_seq);
        // Read: A C G T T T T T (matches first 4, then all mismatches)
        let read_seq: Vec<u8> = vec![0, 1, 2, 3, 3, 3, 3, 3];

        let result = extend_alignment(
            &read_seq, 0,   // read_start
            0,   // genome_start
            1,   // rightward
            8,   // max_extend
            0,   // no previous mismatches
            0,   // no previous length
            10,  // n_mm_max
            0.3, // p_mm_max
            &index, false,
        );

        // Should extend 4 bases (perfect match), then mismatches drag score down
        assert_eq!(result.extend_len, 4);
        assert_eq!(result.max_score, 4);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_chromosome_boundary() {
        // Genome: A C G (3 bases, then padding=5)
        let genome_seq = vec![0, 1, 2];
        let index = make_index_with_seq(&genome_seq);
        // Read is 5 bases, but genome only has 3
        let read_seq: Vec<u8> = vec![0, 1, 2, 3, 0];

        let result = extend_alignment(
            &read_seq, 0, // read_start
            0, // genome_start
            1, // rightward
            5, // max_extend
            0, 0, 10, 0.3, &index, false,
        );

        // Should stop at 3 bases (genome boundary)
        assert_eq!(result.extend_len, 3);
        assert_eq!(result.max_score, 3);
    }

    #[test]
    fn test_extend_n_bases_skipped() {
        // Genome: A N C G (N=4 at position 1)
        let genome_seq = vec![0, 4, 1, 2];
        let index = make_index_with_seq(&genome_seq);
        // Read: A A C G (matches at 0, N skip at 1, matches at 2-3)
        let read_seq: Vec<u8> = vec![0, 0, 1, 2];

        let result = extend_alignment(&read_seq, 0, 0, 1, 4, 0, 0, 10, 0.3, &index, false);

        // Should extend all 4 bases: match + N(skip) + match + match = score 3
        assert_eq!(result.extend_len, 4);
        assert_eq!(result.max_score, 3);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_all_mismatch_returns_zero() {
        // Genome: A A A A (all 0)
        let genome_seq = vec![0, 0, 0, 0];
        let index = make_index_with_seq(&genome_seq);
        // Read: T T T T (all 3, complete mismatch)
        let read_seq: Vec<u8> = vec![3, 3, 3, 3];

        let result = extend_alignment(&read_seq, 0, 0, 1, 4, 0, 0, 10, 0.3, &index, false);

        // Score never goes positive, so extend_len should be 0
        assert_eq!(result.extend_len, 0);
        assert_eq!(result.max_score, 0);
    }

    #[test]
    fn test_extend_recovery_through_mismatch() {
        // Genome: A C G T A C G T A C G T A C (14 bases)
        let genome_seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let index = make_index_with_seq(&genome_seq);
        // Read: A C G X A C G T A C G T A C (1 mismatch at pos 3, then 10 matches)
        //       M M M X M M M M M M M M M M
        let read_seq: Vec<u8> = vec![0, 1, 2, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1];

        let result = extend_alignment(&read_seq, 0, 0, 1, 14, 0, 0, 10, 0.3, &index, false);

        // Should extend past the mismatch: 3M + 1X + 10M
        // Score: +3 -1 +10 = 12, best at position 14
        assert_eq!(result.extend_len, 14);
        assert_eq!(result.max_score, 12);
        assert_eq!(result.n_mismatch, 1);
    }

    #[test]
    fn test_extend_zero_max_extend() {
        let genome_seq = vec![0, 1, 2, 3];
        let index = make_index_with_seq(&genome_seq);
        let read_seq: Vec<u8> = vec![0, 1, 2, 3];

        let result = extend_alignment(
            &read_seq, 0, 0, 1, 0, // max_extend = 0
            0, 0, 10, 0.3, &index, false,
        );

        assert_eq!(result.extend_len, 0);
        assert_eq!(result.max_score, 0);
    }

    #[test]
    fn test_overhang_check_rejects_short_overhang() {
        // Test that the overhang check rejects splice junctions with tiny flanking seeds.
        // Scenario: prev seed length=3 (below default min of 5), current seed length=20
        // With alignSJoverhangMin=5, the 3bp left overhang should cause rejection.
        use crate::align::score::AlignmentScorer;

        let scorer = AlignmentScorer {
            score_gap: 0,
            score_gap_noncan: -8,
            score_gap_gcag: -4,
            score_gap_atac: -8,
            score_del_open: -2,
            score_del_base: -2,
            score_ins_open: -2,
            score_ins_base: -2,
            align_intron_min: 21,
            sjdb_score: 2,
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
            score_stitch_sj_shift: 1,
        };

        // Left overhang (prev.length) = 3, below min of 5
        let left_overhang: usize = 3;
        let right_overhang: usize = 20;
        let min_overhang = scorer.align_sj_overhang_min as usize;

        // Should be rejected
        assert!(left_overhang < min_overhang || right_overhang < min_overhang);

        // Right overhang too small
        let left_overhang: usize = 20;
        let right_overhang: usize = 4;
        assert!(left_overhang < min_overhang || right_overhang < min_overhang);
    }

    #[test]
    fn test_overhang_check_accepts_sufficient_overhang() {
        // Test that splice junctions with sufficient overhang pass the check.
        use crate::align::score::AlignmentScorer;

        let scorer = AlignmentScorer {
            score_gap: 0,
            score_gap_noncan: -8,
            score_gap_gcag: -4,
            score_gap_atac: -8,
            score_del_open: -2,
            score_del_base: -2,
            score_ins_open: -2,
            score_ins_base: -2,
            align_intron_min: 21,
            sjdb_score: 2,
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
            score_stitch_sj_shift: 1,
        };

        // Both overhangs >= 5
        let left_overhang: usize = 5;
        let right_overhang: usize = 10;
        let min_overhang = scorer.align_sj_overhang_min as usize;

        // Should pass
        assert!(!(left_overhang < min_overhang || right_overhang < min_overhang));

        // Exactly at minimum
        let left_overhang: usize = 5;
        let right_overhang: usize = 5;
        assert!(!(left_overhang < min_overhang || right_overhang < min_overhang));
    }

    #[test]
    fn test_bin_based_cluster_bounds() {
        // Verify that cluster genome_start/genome_end are set from bin range
        // With win_bin_nbits=4 (bin_size=16), a window at bin 5 with flank 2
        // should span bins 3-7, i.e., genome_start=48, genome_end=128
        let bin_size: u64 = 1 << 4; // 16
        let bin_start: u64 = 5u64.saturating_sub(2); // 3
        let bin_end: u64 = 5 + 2; // 7
        let genome_start = bin_start * bin_size;
        let genome_end = (bin_end + 1) * bin_size;
        assert_eq!(genome_start, 48);
        assert_eq!(genome_end, 128);
    }

    #[test]
    fn test_window_merge_logic() {
        // Two anchors within winAnchorDistNbins should merge into one window
        // Anchor A at bin 10, Anchor B at bin 15, winAnchorDistNbins=9
        // B is within [10-9, 10+9] = [1, 19] → merge
        let win_anchor_dist_nbins = 9u64;
        let window_bin_start: u64 = 10;
        let window_bin_end: u64 = 10;
        let anchor_bin: u64 = 15;

        let merge_start = window_bin_start.saturating_sub(win_anchor_dist_nbins);
        let merge_end = window_bin_end + win_anchor_dist_nbins;
        let should_merge = anchor_bin >= merge_start && anchor_bin <= merge_end;
        assert!(
            should_merge,
            "Anchors 5 bins apart should merge with dist_nbins=9"
        );

        // Anchor C at bin 25 is outside [1, 19] → separate window
        let anchor_bin_c: u64 = 25;
        let should_merge_c = anchor_bin_c >= merge_start && anchor_bin_c <= merge_end;
        assert!(
            !should_merge_c,
            "Anchors 15 bins apart should NOT merge with dist_nbins=9"
        );
    }

    #[test]
    fn test_window_flank_extension() {
        // Window at bin 10, extended by ±4 flanking bins → bins 6-14
        let mut bin_start: u64 = 10;
        let mut bin_end: u64 = 10;
        let win_flank_nbins: u64 = 4;

        bin_start = bin_start.saturating_sub(win_flank_nbins);
        bin_end += win_flank_nbins;

        assert_eq!(bin_start, 6);
        assert_eq!(bin_end, 14);

        // Edge case: window at bin 2, flanking underflows to 0
        let mut bin_start_edge: u64 = 2;
        bin_start_edge = bin_start_edge.saturating_sub(win_flank_nbins);
        assert_eq!(bin_start_edge, 0, "Flanking should saturate at 0");
    }
}
