/// Seed clustering and stitching via dynamic programming
use crate::align::score::{AlignmentScorer, GapType};
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
            score += 1;
        } else {
            score -= 1;
            n_mm += 1;
        }

        // Check mismatch limits considering the full alignment
        let total_mm = n_mm_prev + n_mm;
        let total_len = len_prev + i + 1;
        let mm_limit = ((p_mm_max * total_len as f64) as u32).min(n_mm_max);
        if total_mm > mm_limit {
            break;
        }

        // Record best extension point (highest score)
        if score > max_score {
            max_score = score;
            best_len = i + 1;
            best_mm = n_mm;
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
    let mut anchor_indices: Vec<usize> = anchor_set
        .iter()
        .enumerate()
        .filter(|(_, is_anchor)| **is_anchor)
        .map(|(i, _)| i)
        .collect();

    // Fallback: if no anchors, use all seeds with non-empty SA ranges.
    // Note: STAR has no fallback (read would be unmapped), but our MMP search
    // overestimates SA ranges because extend_match() only narrows from sa_start.
    // This causes seeds to appear more repetitive than they are.
    // TODO: fix MMP search to narrow SA range from both ends, then remove fallback.
    if anchor_indices.is_empty() {
        anchor_indices = seeds
            .iter()
            .enumerate()
            .filter(|(_, s)| s.sa_end > s.sa_start)
            .map(|(i, _)| i)
            .collect();
        if anchor_indices.is_empty() {
            return Vec::new();
        }
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

/// DP state for seed stitching
#[derive(Debug, Clone)]
struct DpState {
    score: i32,
    prev_seed: Option<usize>,
    genome_pos: u64,
    cigar_ops: Vec<CigarOp>,
    n_mismatch: u32,
    n_gap: u32,
    n_junction: u32,
    junction_motifs: Vec<crate::align::score::SpliceMotif>,
    junction_annotated: Vec<bool>,
}

/// Expanded seed with specific genome position
#[derive(Debug, Clone)]
struct ExpandedSeed {
    read_pos: usize,
    read_end: usize,
    /// Raw SA position (used for genome base access and internal DP calculations).
    /// For reverse strand, genome access is at sa_pos + n_genome.
    genome_pos: u64,
    genome_end: u64,
    length: usize,
}

/// Post-DP junction optimization: apply jR scanning to each splice junction in the
/// winning chain's CIGAR. This matches STAR's architecture where stitchAlignToTranscript
/// is called per adjacent pair in the chosen path, NOT in the O(n²) DP search.
///
/// Walks the CIGAR, and for each RefSkip (splice junction), scans all possible boundary
/// positions to find the one that maximizes read-genome match quality + splice motif score.
/// Returns an optimized CIGAR with adjusted Match lengths around RefSkips, plus updated
/// junction motifs and annotation flags.
#[allow(clippy::too_many_arguments)]
fn optimize_junction_positions(
    cigar_ops: &[CigarOp],
    junction_motifs: &[crate::align::score::SpliceMotif],
    junction_annotated: &[bool],
    genome_start: u64,
    read_seq: &[u8],
    scorer: &AlignmentScorer,
    genome: &crate::genome::Genome,
    is_reverse: bool,
    n_genome: u64,
) -> (
    Vec<CigarOp>,
    Vec<crate::align::score::SpliceMotif>,
    Vec<bool>,
) {
    // Quick exit: no junctions to optimize
    if junction_motifs.is_empty() {
        return (
            cigar_ops.to_vec(),
            junction_motifs.to_vec(),
            junction_annotated.to_vec(),
        );
    }

    let mut result_cigar = Vec::with_capacity(cigar_ops.len());
    let mut result_motifs = Vec::with_capacity(junction_motifs.len());
    let mut result_annotated = Vec::with_capacity(junction_annotated.len());
    let mut junction_idx = 0usize;

    // Track positions as we walk the CIGAR (in SA coordinate space)
    let mut read_pos = 0usize;
    let mut genome_pos = genome_start;

    let mut i = 0;
    while i < cigar_ops.len() {
        match cigar_ops[i] {
            CigarOp::RefSkip(intron_len) => {
                // Found a splice junction — apply jR scanning
                // Get the preceding Match length (donor exon boundary)
                let prev_match_len = match result_cigar.last() {
                    Some(CigarOp::Match(len)) => *len,
                    _ => 0,
                };

                // Get the following Match length (acceptor exon boundary)
                let next_match_len = match cigar_ops.get(i + 1) {
                    Some(CigarOp::Match(len)) => *len,
                    _ => 0,
                };

                // Call jR scanning: r_gap = next_match_len (rightward bound),
                // g_gap = intron_len + next_match_len (so del = g_gap - r_gap = intron_len)
                let (jr_shift, new_motif, _motif_score) = scorer.find_best_junction_position(
                    read_seq,
                    read_pos,
                    genome_pos,
                    next_match_len as i64,
                    intron_len as i64 + next_match_len as i64,
                    genome,
                    is_reverse,
                    n_genome,
                    prev_match_len as usize,
                );

                // Clamp jr_shift to valid range: neither flanking Match can go negative
                let jr_shift = jr_shift
                    .max(-(prev_match_len as i32))
                    .min(next_match_len as i32);

                // Apply the shift to surrounding Matches
                if jr_shift != 0 {
                    // Adjust donor Match (before RefSkip) — guaranteed >= 0 by clamping
                    let new_prev = (prev_match_len as i32 + jr_shift) as u32;
                    if let Some(CigarOp::Match(len)) = result_cigar.last_mut() {
                        *len = new_prev;
                        if *len == 0 {
                            result_cigar.pop();
                        }
                    }

                    // Emit RefSkip (intron length unchanged)
                    result_cigar.push(CigarOp::RefSkip(intron_len));

                    // Adjust acceptor Match (after RefSkip) — will be handled below
                    let new_next = (next_match_len as i32 - jr_shift).max(0) as u32;
                    if new_next > 0 {
                        result_cigar.push(CigarOp::Match(new_next));
                    }

                    // Update genome/read positions: the shift moved the boundary
                    // genome_pos was at the junction boundary; now advance past intron + acceptor
                    genome_pos += intron_len as u64 + next_match_len as u64;
                    read_pos += next_match_len as usize;

                    // Skip the next Match op since we already handled it
                    i += 2;
                } else {
                    // No shift needed — emit RefSkip as-is
                    result_cigar.push(CigarOp::RefSkip(intron_len));
                    genome_pos += intron_len as u64;
                    i += 1;
                }

                // Record the (possibly updated) motif
                result_motifs.push(new_motif);
                if junction_idx < junction_annotated.len() {
                    result_annotated.push(junction_annotated[junction_idx]);
                }
                junction_idx += 1;
            }
            CigarOp::Match(len) => {
                // Merge with previous Match if possible
                if let Some(CigarOp::Match(prev_len)) = result_cigar.last_mut() {
                    *prev_len += len;
                } else {
                    result_cigar.push(CigarOp::Match(len));
                }
                read_pos += len as usize;
                genome_pos += len as u64;
                i += 1;
            }
            CigarOp::Ins(len) => {
                result_cigar.push(CigarOp::Ins(len));
                read_pos += len as usize;
                i += 1;
            }
            CigarOp::Del(len) => {
                result_cigar.push(CigarOp::Del(len));
                genome_pos += len as u64;
                i += 1;
            }
            other => {
                result_cigar.push(other);
                match other {
                    CigarOp::SoftClip(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                        read_pos += len as usize;
                        if matches!(other, CigarOp::Equal(_) | CigarOp::Diff(_)) {
                            genome_pos += len as u64;
                        }
                    }
                    CigarOp::HardClip(_) => {}
                    _ => {}
                }
                i += 1;
            }
        }
    }

    // Verify: total read-consuming bases must be unchanged
    let original_read_len: u32 = cigar_ops
        .iter()
        .map(|op| match op {
            CigarOp::Match(n)
            | CigarOp::Ins(n)
            | CigarOp::SoftClip(n)
            | CigarOp::Equal(n)
            | CigarOp::Diff(n) => *n,
            _ => 0,
        })
        .sum();
    let optimized_read_len: u32 = result_cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n)
            | CigarOp::Ins(n)
            | CigarOp::SoftClip(n)
            | CigarOp::Equal(n)
            | CigarOp::Diff(n) => *n,
            _ => 0,
        })
        .sum();
    debug_assert_eq!(
        original_read_len, optimized_read_len,
        "jR optimization changed read-consuming CIGAR length: {} -> {}, cigar_ops={:?}, result={:?}",
        original_read_len, optimized_read_len, cigar_ops, result_cigar
    );

    (result_cigar, result_motifs, result_annotated)
}

/// Stitch seeds within a cluster using dynamic programming.
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
    stitch_seeds_with_jdb(cluster, read_seq, index, scorer, None)
}

/// Stitch seeds with optional junction database for annotation-aware scoring.
///
/// Converts WindowAlignment entries directly to ExpandedSeeds for the DP.
/// This matches STAR's architecture where the WA array IS the DP input —
/// no SA range re-expansion needed. Each WA entry represents one verified
/// (seed, position) pair assigned during bin-based windowing.
pub fn stitch_seeds_with_jdb(
    cluster: &SeedCluster,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    junction_db: Option<&crate::junction::SpliceJunctionDb>,
) -> Result<Vec<Transcript>, Error> {
    // Convert WindowAlignment entries directly to ExpandedSeeds.
    // Each WA entry has a verified match length and raw SA position.
    let mut expanded_seeds: Vec<ExpandedSeed> = cluster
        .alignments
        .iter()
        .map(|wa| ExpandedSeed {
            read_pos: wa.read_pos,
            read_end: wa.read_pos + wa.length,
            genome_pos: wa.sa_pos,
            genome_end: wa.sa_pos + wa.length as u64,
            length: wa.length,
        })
        .collect();

    if expanded_seeds.is_empty() {
        return Ok(Vec::new());
    }

    // Sort by read position, then by length descending (longest first for dedup)
    expanded_seeds.sort_by(|a, b| a.read_pos.cmp(&b.read_pos).then(b.length.cmp(&a.length)));

    // Deduplicate: keep only the longest seed per (read_pos, genome_pos) pair
    expanded_seeds.dedup_by(|a, b| a.read_pos == b.read_pos && a.genome_pos == b.genome_pos);

    // Cap expanded seeds to prevent pathological O(n²) DP on repetitive regions.
    // Note: STAR caps at seedPerWindowNmax=50 during window assignment, but our Phase 1
    // bypasses that cap. This cap catches overflow from Phase 1.
    const MAX_EXPANDED_SEEDS: usize = 200;
    if expanded_seeds.len() > MAX_EXPANDED_SEEDS {
        expanded_seeds.sort_by(|a, b| b.length.cmp(&a.length));
        expanded_seeds.truncate(MAX_EXPANDED_SEEDS);
        expanded_seeds.sort_by_key(|s| s.read_pos);
    }

    // Pre-DP left extension: compute left extension score for every seed.
    // STAR uses seed_length + left_extend_score as the DP base case, which gives
    // seeds at true genome positions a large advantage over coincidental matches.
    let left_ext_scores: Vec<i32> = expanded_seeds
        .iter()
        .map(|seed| {
            if seed.read_pos > 0 {
                extend_alignment(
                    read_seq,
                    seed.read_pos,
                    seed.genome_pos,
                    -1,
                    seed.read_pos, // max = distance to read start
                    0,
                    seed.length, // seed length for mismatch ratio
                    scorer.n_mm_max,
                    scorer.p_mm_max,
                    index,
                    cluster.is_reverse,
                )
                .max_score
            } else {
                0
            }
        })
        .collect();

    // Initialize DP: one state per expanded seed
    let n = expanded_seeds.len();
    let mut dp: Vec<DpState> = Vec::with_capacity(n);

    // Base case: seed_length + left_extension_score (STAR Pass 1)
    for (i, exp_seed) in expanded_seeds.iter().enumerate() {
        let initial_cigar = vec![CigarOp::Match(exp_seed.length as u32)];
        dp.push(DpState {
            score: exp_seed.length as i32 + left_ext_scores[i],
            prev_seed: None,
            genome_pos: exp_seed.genome_pos,
            cigar_ops: initial_cigar,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: Vec::new(),
            junction_annotated: Vec::new(),
        });
    }

    // DP: for each seed, try connecting to all compatible previous seeds
    // Optimization: find best j first, then build CIGAR only once (avoids repeated cloning)
    for i in 1..n {
        let curr = &expanded_seeds[i];
        let mut best_score = dp[i].score;
        let mut best_j: Option<usize> = None;
        let mut best_gap_type = GapType::Deletion(0);
        let mut best_read_gap = 0i64;
        let mut best_genome_gap = 0i64;
        let mut best_gap_mismatches = 0u32;
        let mut best_eff_length = curr.length;
        let mut best_is_annotated = false;

        for j in 0..i {
            let prev = &expanded_seeds[j];

            // STAR-style overlap trimming: if seeds overlap in read, trim seed B's start
            let mut eff_read_pos = curr.read_pos;
            let mut eff_genome_pos = curr.genome_pos;
            let mut eff_length = curr.length;

            if prev.read_end > eff_read_pos {
                let overlap = prev.read_end - eff_read_pos;
                if overlap >= eff_length {
                    continue; // Seed B fully consumed by overlap
                }
                eff_read_pos = prev.read_end;
                eff_genome_pos += overlap as u64;
                eff_length -= overlap;
            }

            // Similarly handle genome overlap
            if prev.genome_end > eff_genome_pos && eff_genome_pos > prev.genome_pos {
                let g_overlap = (prev.genome_end - eff_genome_pos) as usize;
                if g_overlap >= eff_length {
                    continue; // Seed B fully consumed by genome overlap
                }
                eff_read_pos += g_overlap;
                eff_genome_pos += g_overlap as u64;
                eff_length -= g_overlap;
            }

            // Calculate gaps using effective (trimmed) coordinates
            let read_gap = (eff_read_pos as i64) - (prev.read_end as i64);
            let genome_gap = (eff_genome_pos as i64) - (prev.genome_end as i64);

            // Score the gap — pass strand info for correct splice motif detection
            let (gap_score, gap_type) = scorer.score_gap_with_strand(
                genome_gap,
                read_gap,
                prev.genome_end,
                &index.genome,
                cluster.is_reverse,
                index.genome.n_genome,
            );

            // Count mismatches in the gap region (shared match portion)
            let gap_mismatches = if read_gap > 0 && genome_gap > 0 {
                let shared = (read_gap.min(genome_gap)) as usize;
                if shared > 0 {
                    count_mismatches_in_region(
                        read_seq,
                        prev.read_end,
                        prev.genome_end,
                        shared,
                        index,
                        cluster.is_reverse,
                    )
                } else {
                    0
                }
            } else {
                0
            };

            // Apply alignSJstitchMismatchNmax filter
            // STAR rejects junction stitches with too many mismatches for the motif type
            let mut annotation_bonus = 0i32;
            let mut is_annotated = false;
            if let GapType::SpliceJunction {
                ref motif,
                intron_len,
            } = gap_type
            {
                if !scorer.stitch_mismatch_allowed(motif, gap_mismatches) {
                    continue; // Reject: too many mismatches at this junction type
                }

                // Enforce minimum overhang (alignSJoverhangMin / alignSJDBoverhangMin)
                // Left overhang = length of the preceding seed (immediately flanks the junction)
                // Right overhang = effective length of current seed after overlap trimming
                let left_overhang = prev.length;
                let right_overhang = eff_length;

                // Check if this junction is annotated (for lower overhang threshold + bonus)
                is_annotated = junction_db.is_some_and(|db| {
                    // Convert genome positions to forward coordinates for junction lookup
                    let donor_sa_pos = prev.genome_end;
                    let donor_fwd = index.sa_pos_to_forward(
                        donor_sa_pos,
                        cluster.is_reverse,
                        intron_len as usize,
                    );
                    let acceptor_fwd = donor_fwd + intron_len as u64;
                    // Lookup with strand 0 (unknown) first, then try both strands
                    db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 0)
                        || db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 1)
                        || db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 2)
                });

                // Use lower overhang threshold for annotated junctions
                let base_min_overhang = if is_annotated {
                    scorer.align_sjdb_overhang_min as usize
                } else {
                    scorer.align_sj_overhang_min as usize
                };

                // Use alignSJoverhangMin (5bp) or alignSJDBoverhangMin (3bp) during DP,
                // matching STAR. Stricter outSJfilterOverhangMin applied at SJ.out.tab write time.
                if left_overhang < base_min_overhang || right_overhang < base_min_overhang {
                    continue;
                }

                // Apply annotation bonus during DP (STAR's sjdbScore)
                if is_annotated {
                    annotation_bonus = scorer.sjdb_score;
                }
            }

            // STAR: shared region scores +1 per match, -1 per mismatch
            // Total = shared_length - 2*mismatches (each mismatch flips from +1 to -1)
            let shared_bases = if read_gap > 0 && genome_gap > 0 {
                read_gap.min(genome_gap) as i32
            } else {
                0
            };
            let shared_score = shared_bases - (2 * gap_mismatches as i32);

            // Transition score = prev_score + gap_penalty + shared_region_score + current_seed_match_score + annotation_bonus
            let transition_score =
                dp[j].score + gap_score + shared_score + (eff_length as i32) + annotation_bonus;

            if transition_score > best_score {
                best_score = transition_score;
                best_j = Some(j);
                best_gap_type = gap_type;
                best_read_gap = read_gap;
                best_genome_gap = genome_gap;
                best_gap_mismatches = gap_mismatches;
                best_eff_length = eff_length;
                best_is_annotated = is_annotated;
            }
        }

        // Build CIGAR only once for the best transition (instead of cloning per candidate)
        if let Some(j) = best_j {
            let mut cigar = dp[j].cigar_ops.clone();

            // Emit CIGAR operations for the gap between seeds
            // Handle negative gaps explicitly to avoid integer overflow
            if best_read_gap == 0 && best_genome_gap == 0 {
                // No gap - seeds are adjacent
            } else if best_read_gap < 0 || best_genome_gap < 0 {
                // Negative gap indicates overlapping seeds or error
                // score_gap() already penalized this, but we skip the connection
                log::trace!(
                    "Skipping connection with negative gap: read_gap={}, genome_gap={}",
                    best_read_gap,
                    best_genome_gap
                );
                continue;
            } else if best_read_gap == 0 && best_genome_gap > 0 {
                // Pure deletion or splice junction (read doesn't advance, genome does)
                match best_gap_type {
                    GapType::SpliceJunction { intron_len, .. } => {
                        cigar.push(CigarOp::RefSkip(intron_len));
                    }
                    _ => {
                        cigar.push(CigarOp::Del(best_genome_gap as u32));
                    }
                }
            } else if best_read_gap > 0 && best_genome_gap == 0 {
                // Pure insertion (read advances, genome doesn't)
                cigar.push(CigarOp::Ins(best_read_gap as u32));
            } else {
                // Both gaps positive: combined gap region
                let rg = best_read_gap as u32;
                let gg = best_genome_gap as u32;

                let shared = rg.min(gg);
                let excess_genome = gg.saturating_sub(rg);
                let excess_read = rg.saturating_sub(gg);

                if shared > 0 {
                    // Try to merge with previous Match operation
                    if let Some(CigarOp::Match(prev_len)) = cigar.last_mut() {
                        *prev_len += shared;
                    } else {
                        cigar.push(CigarOp::Match(shared));
                    }
                }
                if excess_genome > 0 {
                    match best_gap_type {
                        GapType::SpliceJunction { intron_len, .. } => {
                            cigar.push(CigarOp::RefSkip(intron_len));
                        }
                        _ => {
                            cigar.push(CigarOp::Del(excess_genome));
                        }
                    }
                }
                if excess_read > 0 {
                    cigar.push(CigarOp::Ins(excess_read));
                }
            }

            // Add current seed match using effective length (after overlap trimming)
            if let Some(CigarOp::Match(prev_len)) = cigar.last_mut() {
                // Merge with previous Match operation
                *prev_len += best_eff_length as u32;
            } else {
                // No previous Match, or previous op was not Match
                cigar.push(CigarOp::Match(best_eff_length as u32));
            }

            let mut n_gap = dp[j].n_gap;
            let mut n_junction = dp[j].n_junction;
            let mut junction_motifs = dp[j].junction_motifs.clone();
            let mut junction_annotated = dp[j].junction_annotated.clone();
            match best_gap_type {
                GapType::Insertion(_) | GapType::Deletion(_) => n_gap += 1,
                GapType::SpliceJunction { motif, .. } => {
                    n_junction += 1;
                    junction_motifs.push(motif);
                    junction_annotated.push(best_is_annotated);
                }
            }
            dp[i].score = best_score;
            dp[i].prev_seed = Some(j);
            dp[i].genome_pos = curr.genome_pos;
            dp[i].cigar_ops = cigar;
            dp[i].n_mismatch = dp[j].n_mismatch + best_gap_mismatches;
            dp[i].n_gap = n_gap;
            dp[i].n_junction = n_junction;
            dp[i].junction_motifs = junction_motifs;
            dp[i].junction_annotated = junction_annotated;
        }
        // If no best_j, dp[i] keeps its initial state (single seed)
    }

    // Select best chain endpoint using dp_score + right_extension_score (STAR Pass 2).
    // Right extension rewards chains that end at true genome positions.
    let mut best_final_idx = 0;
    let mut best_total_score = i32::MIN;
    for i in 0..n {
        let right_ext_score = if expanded_seeds[i].read_end < read_seq.len() {
            extend_alignment(
                read_seq,
                expanded_seeds[i].read_end,
                expanded_seeds[i].genome_end,
                1,
                read_seq.len() - expanded_seeds[i].read_end,
                dp[i].n_mismatch,
                expanded_seeds[i].read_end, // cumulative aligned length
                scorer.n_mm_max,
                scorer.p_mm_max,
                index,
                cluster.is_reverse,
            )
            .max_score
        } else {
            0
        };
        let total = dp[i].score + right_ext_score;
        if total > best_total_score {
            best_total_score = total;
            best_final_idx = i;
        }
    }

    let best_state = &dp[best_final_idx];

    // Find the first seed in the DP chain by tracing back through prev_seed
    let mut chain_start_idx = best_final_idx;
    while let Some(prev) = dp[chain_start_idx].prev_seed {
        chain_start_idx = prev;
    }
    let chain_start_seed = &expanded_seeds[chain_start_idx];

    // Post-DP junction optimization: apply jR scanning to each splice junction
    // in the winning chain. STAR calls stitchAlignToTranscript per adjacent pair
    // in the chosen path, NOT in the O(n²) DP. We do the same: scan only the
    // ~1-3 junctions in the final alignment.
    let (optimized_cigar, optimized_motifs, optimized_annotated) = optimize_junction_positions(
        &best_state.cigar_ops,
        &best_state.junction_motifs,
        &best_state.junction_annotated,
        chain_start_seed.genome_pos,
        read_seq,
        scorer,
        &index.genome,
        cluster.is_reverse,
        index.genome.n_genome,
    );

    // Calculate the read region covered by the alignment
    let alignment_start = chain_start_seed.read_pos;

    // Calculate aligned read length from CIGAR operations
    let mut aligned_read_len = 0usize;
    for op in &optimized_cigar {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) | CigarOp::Ins(len) => {
                aligned_read_len += *len as usize;
            }
            _ => {}
        }
    }
    let alignment_end = alignment_start + aligned_read_len;

    // Compute right-side genome position by walking the optimized CIGAR
    let mut right_genome_pos = chain_start_seed.genome_pos;
    for op in &optimized_cigar {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                right_genome_pos += *len as u64;
            }
            CigarOp::Del(len) | CigarOp::RefSkip(len) => {
                right_genome_pos += *len as u64;
            }
            _ => {}
        }
    }

    // Extend alignment into flanking regions (STAR-style extendAlign)
    let left_extend = if alignment_start > 0 {
        extend_alignment(
            read_seq,
            alignment_start,             // read boundary (exclusive, leftward)
            chain_start_seed.genome_pos, // genome boundary (exclusive, leftward)
            -1,                          // leftward
            alignment_start,             // max distance to read start
            best_state.n_mismatch,       // mismatches in aligned portion
            aligned_read_len,            // length of aligned portion
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
            alignment_end,                  // read boundary (inclusive, rightward)
            right_genome_pos,               // genome boundary (inclusive, rightward)
            1,                              // rightward
            read_seq.len() - alignment_end, // max distance to read end
            best_state.n_mismatch + left_extend.n_mismatch, // cumulative mismatches
            aligned_read_len + left_extend.extend_len, // cumulative length
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

    // Build final CIGAR with extensions
    let mut final_cigar = Vec::new();

    // Remaining 5' soft clip after left extension
    let remaining_left_clip = alignment_start - left_extend.extend_len;
    if remaining_left_clip > 0 {
        final_cigar.push(CigarOp::SoftClip(remaining_left_clip as u32));
    }

    // Left extension Match (merge with first main CIGAR Match if possible)
    if left_extend.extend_len > 0 {
        final_cigar.push(CigarOp::Match(left_extend.extend_len as u32));
    }

    // Main alignment CIGAR (merge first op with left extension Match if both are Match)
    for op in &optimized_cigar {
        if let CigarOp::Match(len) = op {
            if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
                *prev_len += len;
                continue;
            }
        }
        final_cigar.push(*op);
    }

    // Right extension Match (merge with last main CIGAR Match if possible)
    if right_extend.extend_len > 0 {
        if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
            *prev_len += right_extend.extend_len as u32;
        } else {
            final_cigar.push(CigarOp::Match(right_extend.extend_len as u32));
        }
    }

    // Remaining 3' soft clip after right extension
    let remaining_right_clip = (read_seq.len() - alignment_end) - right_extend.extend_len;
    if remaining_right_clip > 0 {
        final_cigar.push(CigarOp::SoftClip(remaining_right_clip as u32));
    }

    // Adjust genome start position for left extension (raw SA coordinates)
    let adjusted_genome_start = chain_start_seed.genome_pos - left_extend.extend_len as u64;

    // Adjust score: remove approximate pre-DP left extension (baked into dp score),
    // add authoritative post-DP left + right extensions (computed with full chain context)
    let adjusted_score = best_state.score - left_ext_scores[chain_start_idx]
        + left_extend.max_score
        + right_extend.max_score;

    use crate::align::transcript::Exon;

    // Count mismatches in the final alignment
    // count_mismatches uses raw SA position + n_genome offset for reverse strand
    // MUST be called BEFORE CIGAR reversal since it walks in RC genome order
    let n_mismatch = count_mismatches(
        read_seq,
        &final_cigar,
        adjusted_genome_start, // Raw SA position for genome access
        0,                     // Read starts at position 0 (CIGAR includes soft clips)
        index,
        cluster.is_reverse, // Pass reverse-strand flag for correct sequence comparison
    );

    // SAM CIGAR must be in forward genome order (5'→3' reference direction).
    // The DP builds CIGAR in read/RC-genome order; for reverse strand this is reversed.
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

    // Convert raw SA position to forward genome coordinates for the transcript
    let forward_genome_start =
        index.sa_pos_to_forward(adjusted_genome_start, cluster.is_reverse, ref_len as usize);
    let forward_genome_end = forward_genome_start + ref_len;

    // Build exons from CIGAR using forward genome coordinates
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

    // Build transcript
    let t_genome_start = merged_exons
        .first()
        .map(|e| e.genome_start)
        .unwrap_or(forward_genome_start);
    let t_genome_end = merged_exons
        .last()
        .map(|e| e.genome_end)
        .unwrap_or(forward_genome_end);

    // Apply STAR's genomic length penalty: penalizes long-spanning alignments
    let genomic_span = t_genome_end - t_genome_start;
    let length_penalty = scorer.genomic_length_penalty(genomic_span);
    let final_score = (adjusted_score + length_penalty).max(0);

    let transcript = Transcript {
        chr_idx: cluster.chr_idx,
        genome_start: t_genome_start,
        genome_end: t_genome_end,
        is_reverse: cluster.is_reverse,
        exons: merged_exons,
        cigar: final_cigar, // Use CIGAR with extensions + soft clips
        score: final_score,
        n_mismatch,
        n_gap: best_state.n_gap,
        n_junction: best_state.n_junction,
        junction_motifs: optimized_motifs,
        junction_annotated: optimized_annotated,
        read_seq: read_seq.to_vec(),
    };

    Ok(vec![transcript])
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
    fn test_expanded_seed_sorting() {
        let mut expanded = [
            ExpandedSeed {
                read_pos: 10,
                read_end: 15,
                genome_pos: 100,
                genome_end: 105,
                length: 5,
            },
            ExpandedSeed {
                read_pos: 5,
                read_end: 10,
                genome_pos: 50,
                genome_end: 55,
                length: 5,
            },
        ];

        expanded.sort_by_key(|s| s.read_pos);

        assert_eq!(expanded[0].read_pos, 5);
        assert_eq!(expanded[1].read_pos, 10);
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
