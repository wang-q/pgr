//! Two-index merge alignment: seed hits -> anti-diagonal chains -> PSL blocks.

use super::dist::validate_compatible;
use super::PgiIndex;
use crate::libs::alignment::align_banded_local;
use crate::libs::alignment::coords::reverse_range_pair;
use crate::libs::fmt::psl::Psl;
use crate::libs::nt::rev_comp;
use crate::libs::poa::align::AlignmentParams;
use rayon::prelude::*;

/// Alignment parameters; defaults follow FastGA (`-f 10`, `-c 85`, `-s 1000`).
#[derive(Debug, Clone, Copy)]
pub struct AlignParams {
    /// Maximum k-mer frequency on either side to keep as a seed.
    pub freq: u32,
    /// Minimum per-axis seed span (bp) for a chain to be kept.
    pub min_span: u32,
    /// Maximum bp gap between consecutive seeds in a chain.
    pub max_gap: u32,
    /// Diagonal band half-width (bp) around the chain mean.
    pub band: u32,
    /// Maximum gap (bp) between adjacent colinear chains to merge.
    pub merge_gap: u32,
    /// Minimum shared seed length (bp); `None` means exact k-mers (default).
    /// Lower values enable adaptamer-style partial seeds, which are
    /// experimental: they degrade block structure (see pgi-align.md §5.9).
    pub min_shared: Option<usize>,
}

impl Default for AlignParams {
    fn default() -> Self {
        Self {
            freq: 10,
            min_span: 85,
            max_gap: 1000,
            band: 128,
            merge_gap: 5000,
            min_shared: None,
        }
    }
}

/// A seed hit between two indexes, orientation resolved into a shared
/// (contig_a, contig_b, strand) coordinate space.
#[derive(Debug, Clone, Copy)]
pub struct SeedHit {
    /// Reference (a) contig id.
    pub a_contig: u32,
    /// Reference position, forward coordinates.
    pub a_pos: u32,
    /// Query (b) contig id.
    pub b_contig: u32,
    /// Query position in orientation space (RC space when `strand` is 1).
    pub b_pos: u32,
    /// Shared prefix length (bp) of the two k-mers.
    pub shared: u32,
    /// 0 = forward, 1 = reverse (query window is the RC of the reference).
    pub strand: u8,
}

/// Merge two sorted k-mer indexes, emitting adaptamer-style seed hits for
/// every k-mer pair sharing at least `min_shared` bases (mirroring FastGA's
/// lcp-driven merge; exact k-mers are the `min_shared == k` special case).
///
/// Both-strand entries are emitted by `pgi build`, so shared keys with equal
/// strand flags mean a forward hit, and differing flags a reverse hit. A
/// k-mer whose matching prefix range in the other index exceeds `freq`
/// entries is skipped (FastGA's frequency filter).
pub fn merge_seed_hits(
    a: &PgiIndex,
    b: &PgiIndex,
    freq: u32,
    min_shared: usize,
) -> anyhow::Result<Vec<SeedHit>> {
    validate_compatible(a, b)?;
    let k = a.k;
    anyhow::ensure!(
        min_shared >= 1 && min_shared <= k,
        "min_shared must be in 1..={k}"
    );
    let k_bits = 2 * k;
    let drop_bits = k_bits - 2 * min_shared;
    let range = 1u128 << drop_bits;
    let mask = if k_bits >= 128 {
        u128::MAX
    } else {
        (1u128 << k_bits) - 1
    };
    let mut hits = Vec::new();
    let mut i = 0usize;
    while i < a.entries.len() {
        let ea = &a.entries[i];
        i += 1;
        if ea.freq > freq {
            continue;
        }
        // Prefix range in b: keys sharing the first `min_shared` bases.
        let lo = ea.kmer & !(range - 1) & mask;
        let hi = lo + range;
        let j0 = b.entries.partition_point(|e| e.kmer < lo);
        let mut j = j0;
        while j < b.entries.len() && b.entries[j].kmer < hi {
            j += 1;
        }
        let range_size = (j - j0) as u32;
        if range_size == 0 || range_size > freq {
            continue;
        }
        let ap = &a.positions[ea.pos_start as usize..(ea.pos_start + ea.freq) as usize];
        for eb in &b.entries[j0..j] {
            if eb.freq > freq {
                continue;
            }
            let shared = shared_prefix(ea.kmer, eb.kmer, k);
            if shared < min_shared as u32 {
                continue;
            }
            let bp = &b.positions[eb.pos_start as usize..(eb.pos_start + eb.freq) as usize];
            for &(ac, apos, astrand) in ap {
                for &(bc, bpos, bstrand) in bp {
                    let fwd = astrand == bstrand;
                    let b_len = b.contigs[bc as usize].1;
                    let oriented = if fwd {
                        bpos as u64
                    } else {
                        b_len - k as u64 - bpos as u64
                    };
                    hits.push(SeedHit {
                        a_contig: ac,
                        a_pos: apos,
                        b_contig: bc,
                        b_pos: oriented as u32,
                        shared,
                        strand: u8::from(!fwd),
                    });
                }
            }
        }
    }
    Ok(hits)
}

/// Longest common prefix (in bases) of two 2-bit k-mer keys.
fn shared_prefix(a: u128, b: u128, k: usize) -> u32 {
    let x = a ^ b;
    if x == 0 {
        return k as u32;
    }
    let bitlen = 128 - x.leading_zeros();
    ((2 * k).saturating_sub(bitlen as usize) / 2) as u32
}

/// Effective minimum shared seed length: `params.min_shared` or `k` (exact
/// k-mers), capped at `k`.
fn effective_min_shared(idx: &PgiIndex, params: &AlignParams) -> usize {
    params.min_shared.unwrap_or(idx.k).min(idx.k)
}

/// One chained diagonal segment on a single (contig_a, contig_b, strand) pair.
#[derive(Debug, Clone, Copy)]
pub struct Chain {
    /// Reference contig id.
    pub a_contig: u32,
    /// Query contig id.
    pub b_contig: u32,
    /// 0 = forward, 1 = reverse.
    pub strand: u8,
    /// Seed span start on the reference.
    pub a_start: u32,
    /// Seed span end (exclusive, includes `k`) on the reference.
    pub a_end: u32,
    /// Seed span start on the query, in orientation space.
    pub b_start: u32,
    /// Seed span end (exclusive, includes `k`) on the query, orientation space.
    pub b_end: u32,
    /// Number of seed hits in the chain.
    pub seeds: usize,
    /// Mean diagonal (a_pos - b_pos in orientation space) of the seeds.
    pub diag: i64,
    /// Seed diagonal spread (max - min).
    pub diag_span: i64,
}

/// Greedy anti-diagonal chaining of seed hits.
///
/// Hits are scanned per (contig_a, contig_b, strand) group sorted by diagonal
/// then reference position. A chain extends while the next hit stays within
/// `band` of the group mean diagonal and within `max_gap` on both axes; an
/// in-gap hit outside the band is ignored (it belongs to another tube and
/// must not split the chain), and only a group change or an over-`max_gap`
/// jump closes the chain. Chains covering less than `min_span` on either
/// axis are dropped.
pub fn chain_hits(
    hits: &[SeedHit],
    k: u32,
    min_span: u32,
    max_gap: u32,
    band: u32,
    merge_gap: u32,
) -> Vec<Chain> {
    let mut sorted = hits.to_vec();
    sorted.sort_unstable_by_key(|h| {
        let diag = h.a_pos as i64 - h.b_pos as i64;
        (h.a_contig, h.b_contig, h.strand, diag, h.a_pos)
    });

    let mut cur: Option<ChainCursor> = None;
    let mut chains = Vec::new();
    for h in &sorted {
        let diag = h.a_pos as i64 - h.b_pos as i64;
        if let Some(c) = cur {
            if h.a_contig == c.a_contig && h.b_contig == c.b_contig && h.strand == c.strand {
                let mean = c.diag_sum / c.count as i64;
                let gap_a = h.a_pos as i64 - c.last_a as i64;
                let gap_b = h.b_pos as i64 - c.last_b as i64;
                let in_gap =
                    (0..=max_gap as i64).contains(&gap_a) && (0..=max_gap as i64).contains(&gap_b);
                if in_gap && (diag - mean).abs() <= band as i64 {
                    cur = Some(ChainCursor {
                        last_a: h.a_pos,
                        last_b: h.b_pos,
                        diag_sum: c.diag_sum + diag,
                        count: c.count + 1,
                        diag_min: c.diag_min.min(diag),
                        diag_max: c.diag_max.max(diag),
                        ..c
                    });
                    continue;
                }
                if in_gap {
                    continue; // outside the diagonal band: ignore, keep chain
                }
            }
            push_chain(&mut chains, &cur, k, min_span);
        }
        cur = Some(ChainCursor {
            a_contig: h.a_contig,
            b_contig: h.b_contig,
            strand: h.strand,
            first_a: h.a_pos,
            last_a: h.a_pos,
            first_b: h.b_pos,
            last_b: h.b_pos,
            diag_sum: diag,
            diag_min: diag,
            diag_max: diag,
            count: 1,
        });
    }
    push_chain(&mut chains, &cur, k, min_span);
    merge_adjacent_chains(&mut chains, band, merge_gap);
    chains
}

/// Merge adjacent chains on the same contig pair and strand whose spans are
/// within `merge_gap` and whose diagonals are within `band` of each other.
///
/// Insertions (IS elements etc.) shift the local diagonal past the chaining
/// band and split one syntenic block into several greedy chains; this pass
/// stitches them back together.
fn merge_adjacent_chains(chains: &mut Vec<Chain>, band: u32, merge_gap: u32) {
    if chains.len() < 2 {
        return;
    }
    chains.sort_by_key(|c| (c.a_contig, c.b_contig, c.strand, c.a_start));
    let mut out: Vec<Chain> = Vec::with_capacity(chains.len());
    for c in chains.drain(..) {
        if let Some(last) = out.last_mut() {
            let same_group = last.a_contig == c.a_contig
                && last.b_contig == c.b_contig
                && last.strand == c.strand;
            let gap_a = c.a_start as i64 - last.a_end as i64;
            let gap_b = c.b_start as i64 - last.b_end as i64;
            if same_group
                && (0..=merge_gap as i64).contains(&gap_a)
                && (0..=merge_gap as i64).contains(&gap_b)
                && (c.diag - last.diag).abs() <= band as i64
            {
                let span = last.diag_span.max(c.diag_span) + (c.diag - last.diag).abs();
                last.a_end = last.a_end.max(c.a_end);
                last.b_end = last.b_end.max(c.b_end);
                last.seeds += c.seeds;
                last.diag = (last.diag + c.diag) / 2;
                last.diag_span = span;
                continue;
            }
        }
        out.push(c);
    }
    *chains = out;
}

/// Greedy chain cursor: the current chain under construction.
#[derive(Debug, Clone, Copy)]
struct ChainCursor {
    a_contig: u32,
    b_contig: u32,
    strand: u8,
    first_a: u32,
    last_a: u32,
    first_b: u32,
    last_b: u32,
    diag_sum: i64,
    diag_min: i64,
    diag_max: i64,
    count: usize,
}

fn push_chain(chains: &mut Vec<Chain>, cur: &Option<ChainCursor>, k: u32, min_span: u32) {
    if let Some(c) = cur {
        let span_a = (c.last_a - c.first_a).saturating_add(k);
        let span_b = (c.last_b - c.first_b).saturating_add(k);
        if span_a >= min_span && span_b >= min_span {
            let mean = c.diag_sum / c.count as i64;
            chains.push(Chain {
                a_contig: c.a_contig,
                b_contig: c.b_contig,
                strand: c.strand,
                a_start: c.first_a,
                a_end: c.last_a.saturating_add(k),
                b_start: c.first_b,
                b_end: c.last_b.saturating_add(k),
                seeds: c.count,
                diag: mean,
                diag_span: c.diag_max - c.diag_min,
            });
        }
    }
}

/// Convert one chain into a single-block PSL record (q = query, t = reference).
///
/// Reverse-strand blocks carry `q_start`/`q_end` in original query coordinates
/// (converted from orientation space), matching the pgr/UCSC PSL convention.
pub fn chain_to_psl(chain: &Chain, a: &PgiIndex, b: &PgiIndex) -> Psl {
    let (a_name, a_len) = {
        let (name, len) = &a.contigs[chain.a_contig as usize];
        (name, *len as u32)
    };
    let (b_name, b_len) = {
        let (name, len) = &b.contigs[chain.b_contig as usize];
        (name, *len as u32)
    };
    let (q_start, q_end, strand) = if chain.strand == 0 {
        (chain.b_start, chain.b_end, "+")
    } else {
        let (s, e) = reverse_range_pair(chain.b_start, chain.b_end, b_len);
        (s, e, "-")
    };
    let mut psl = Psl::new();
    psl.q_name = b_name.clone();
    psl.q_size = b_len;
    psl.q_start = q_start as i32;
    psl.q_end = q_end as i32;
    psl.t_name = a_name.clone();
    psl.t_size = a_len;
    psl.t_start = chain.a_start as i32;
    psl.t_end = chain.a_end as i32;
    psl.strand = strand.to_string();
    psl.block_count = 1;
    psl.block_sizes.push(q_end - q_start);
    psl.q_starts.push(q_start);
    psl.t_starts.push(chain.a_start);
    psl
}

/// Align two compatible indexes: merge seeds, chain, and emit PSL blocks.
pub fn align_to_psl(a: &PgiIndex, b: &PgiIndex, params: &AlignParams) -> anyhow::Result<Vec<Psl>> {
    let hits = merge_seed_hits(a, b, params.freq, effective_min_shared(a, params))?;
    let chains = chain_hits(
        &hits,
        a.k as u32,
        params.min_span,
        params.max_gap,
        params.band,
        params.merge_gap,
    );
    Ok(chains.iter().map(|c| chain_to_psl(c, a, b)).collect())
}

/// Default window size (bp) for chain extension.
pub const EXTEND_WINDOW: usize = 16_000;
/// Default step (bp) between extension windows (2 kb overlap).
pub const EXTEND_STEP: usize = 14_000;

/// Extend one chain into scored PSL records via windowed banded alignment.
///
/// Chains longer than `window` are split into overlapping windows placed on
/// the chain diagonal, so giant syntenic chains still get real identity
/// counts. Returns empty when no window scores; callers fall back to
/// [`chain_to_psl`].
pub fn extend_chain(
    chain: &Chain,
    a: &PgiIndex,
    b: &PgiIndex,
    a_seqs: &[(String, Vec<u8>)],
    b_seqs: &[(String, Vec<u8>)],
    window: usize,
    step: usize,
) -> Vec<Psl> {
    let jobs = chain_windows(0, chain, a_seqs, b_seqs, window, step);
    let dp_band = (chain.diag_span as usize + 32).min(128);
    jobs.par_iter()
        .filter_map(|job| extend_window(job, chain, a, b, a_seqs, b_seqs, dp_band))
        .collect()
}

/// One banded extension window (chain id + oriented query range).
#[derive(Debug, Clone, Copy)]
struct WindowJob {
    chain_id: usize,
    q_win: usize,
    win_end: usize,
}

/// Build the bounds-checked extension window list for one chain.
fn chain_windows(
    chain_id: usize,
    chain: &Chain,
    a_seqs: &[(String, Vec<u8>)],
    b_seqs: &[(String, Vec<u8>)],
    window: usize,
    step: usize,
) -> Vec<WindowJob> {
    let Some((_, a_int)) = a_seqs.get(chain.a_contig as usize) else {
        return Vec::new();
    };
    let Some((_, b_int)) = b_seqs.get(chain.b_contig as usize) else {
        return Vec::new();
    };
    let a_len = a_int.len();
    let b_len = b_int.len();
    let q0 = chain.b_start as usize;
    let q1 = chain.b_end as usize;
    if q1 > b_len || q0 >= q1 || step == 0 {
        return Vec::new();
    }
    let mut jobs = Vec::new();
    let mut q_win = q0;
    while q_win < q1 {
        let win_end = (q_win + window).min(q1);
        // Expected target window start from the chain diagonal (t = q + diag).
        let t_start = q_win as i64 + chain.diag;
        let t_end = t_start + (win_end - q_win) as i64;
        if t_start >= 0 && t_end <= a_len as i64 {
            jobs.push(WindowJob {
                chain_id,
                q_win,
                win_end,
            });
        }
        q_win += step;
    }
    jobs
}

/// Extend one window into a scored PSL record, if the banded alignment scores.
fn extend_window(
    job: &WindowJob,
    chain: &Chain,
    a: &PgiIndex,
    b: &PgiIndex,
    a_seqs: &[(String, Vec<u8>)],
    b_seqs: &[(String, Vec<u8>)],
    dp_band: usize,
) -> Option<Psl> {
    let (_, a_int) = a_seqs.get(chain.a_contig as usize)?;
    let (_, b_int) = b_seqs.get(chain.b_contig as usize)?;
    let b_len = b_int.len();
    let (q_win, win_end) = (job.q_win, job.win_end);
    let t_start = q_win as i64 + chain.diag;
    let (t_win_start, t_win_end) = (
        t_start as usize,
        (t_start + (win_end - q_win) as i64) as usize,
    );
    let t = &a_int[t_win_start..t_win_end];
    let q: Vec<u8> = if chain.strand == 0 {
        b_int[q_win..win_end].to_vec()
    } else {
        let orig = (b_len - win_end, b_len - q_win);
        rev_comp(&b_int[orig.0..orig.1]).collect()
    };
    // Window starts were placed on the chain diagonal, so the expected
    // within-window diagonal (q_i - t_j) is 0.
    let aln = align_banded_local(&q, t, dp_band, 0, &AlignmentParams::default())?;
    let q_aln = aln.q_aln;
    let t_aln = aln.t_aln;
    let q_start = aln.q_start;
    let t_start = aln.t_start;
    let q_covered = q_aln.iter().filter(|&&c| c != b'-').count();
    let t_covered = t_aln.iter().filter(|&&c| c != b'-').count();
    if q_covered == 0 || t_covered == 0 {
        return None;
    }
    let q_abs = q_win + q_start;
    // Reverse-strand PSLs use whole-contig RC-space coordinates here;
    // `Psl::from_align` converts them to original ascending coordinates.
    let strand = if chain.strand == 0 { "+" } else { "-" };
    Psl::from_align(
        &b.contigs[chain.b_contig as usize].0,
        b.contigs[chain.b_contig as usize].1 as u32,
        q_abs as i32,
        (q_abs + q_covered) as i32,
        &String::from_utf8_lossy(&q_aln),
        &a.contigs[chain.a_contig as usize].0,
        a.contigs[chain.a_contig as usize].1 as u32,
        (t_win_start + t_start) as i32,
        (t_win_start + t_start + t_covered) as i32,
        &String::from_utf8_lossy(&t_aln),
        strand,
    )
}

/// Align two indexes, extending chains when sequences are provided and
/// falling back to plain blocks otherwise.
pub fn align_to_psl_ext(
    a: &PgiIndex,
    b: &PgiIndex,
    params: &AlignParams,
    a_seqs: &[(String, Vec<u8>)],
    b_seqs: &[(String, Vec<u8>)],
) -> anyhow::Result<Vec<Psl>> {
    let hits = merge_seed_hits(a, b, params.freq, effective_min_shared(a, params))?;
    let chains = chain_hits(
        &hits,
        a.k as u32,
        params.min_span,
        params.max_gap,
        params.band,
        params.merge_gap,
    );
    // Flatten all extension windows across chains into one parallel stream so
    // a giant chain cannot serialize behind smaller chains on a single task.
    let jobs: Vec<WindowJob> = chains
        .iter()
        .enumerate()
        .flat_map(|(id, c)| chain_windows(id, c, a_seqs, b_seqs, EXTEND_WINDOW, EXTEND_STEP))
        .collect();
    let records: Vec<Option<Psl>> = jobs
        .par_iter()
        .map(|job| {
            let chain = &chains[job.chain_id];
            let dp_band = (chain.diag_span as usize + 32).min(128);
            extend_window(job, chain, a, b, a_seqs, b_seqs, dp_band)
        })
        .collect();
    let mut covered = vec![false; chains.len()];
    let mut out = Vec::with_capacity(records.len());
    for (job, rec) in jobs.iter().zip(records) {
        if let Some(psl) = rec {
            covered[job.chain_id] = true;
            out.push(psl);
        }
    }
    for (id, chain) in chains.iter().enumerate() {
        if !covered[id] {
            out.push(chain_to_psl(chain, a, b));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::pgi::build::build_from_seqs;

    fn pseudo_random_seq(len: usize, seed: u64) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut x = seed;
        (0..len)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                bases[(x >> 33) as usize & 3]
            })
            .collect()
    }

    fn build(seq: &[u8]) -> PgiIndex {
        build_from_seqs(vec![(String::from("c"), seq.to_vec())], 10, 4, 2, false).unwrap()
    }

    #[test]
    fn merge_forward_and_reverse_hits() {
        let seq_a = pseudo_random_seq(300, 42);
        let (ia, ib) = (build(&seq_a), build(&seq_a));
        let hits = merge_seed_hits(&ia, &ib, 10, 10).unwrap();
        assert!(!hits.is_empty());
        assert!(
            hits.iter().all(|h| h.strand == 0),
            "identical sequences must be forward"
        );
        assert!(
            hits.iter().all(|h| h.a_pos == h.b_pos),
            "identical coordinates"
        );

        // Reverse-complemented query: every hit is reverse-oriented and the
        // oriented query position equals the reference position.
        let rc: Vec<u8> = crate::libs::nt::rev_comp(&seq_a).collect();
        let ir = build(&rc);
        let rev = merge_seed_hits(&ia, &ir, 10, 10).unwrap();
        assert!(!rev.is_empty());
        assert!(
            rev.iter().all(|h| h.strand == 1),
            "RC query must be reverse"
        );
        assert!(
            rev.iter().all(|h| h.a_pos == h.b_pos),
            "RC-space position mismatch"
        );
    }

    #[test]
    fn freq_filter_drops_repetitive_keys() {
        // Periodic sequence: each 10-mer repeats ~30x, so freq=1 drops all.
        let seq: Vec<u8> = (0..128u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let idx = build(&seq);
        let hits = merge_seed_hits(&idx, &idx, 1, 10).unwrap();
        assert!(hits.is_empty(), "repetitive k-mers must be filtered");
        let hits = merge_seed_hits(&idx, &idx, 200, 10).unwrap();
        assert!(!hits.is_empty(), "high threshold keeps repeats");

        // Random sequence: k-mers are mostly unique, so freq=1 keeps hits.
        let rnd = pseudo_random_seq(300, 7);
        let rnd_idx = build(&rnd);
        let hits = merge_seed_hits(&rnd_idx, &rnd_idx, 1, 10).unwrap();
        assert!(!hits.is_empty(), "unique k-mers pass freq=1");
    }

    #[test]
    fn chain_groups_by_direction_and_contig() {
        let hits = vec![
            SeedHit {
                a_contig: 0,
                a_pos: 100,
                b_contig: 0,
                b_pos: 200,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 110,
                b_contig: 0,
                b_pos: 210,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 120,
                b_contig: 0,
                b_pos: 220,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 100,
                b_contig: 0,
                b_pos: 200,
                shared: 10,
                strand: 1,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 110,
                b_contig: 0,
                b_pos: 210,
                shared: 10,
                strand: 1,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 120,
                b_contig: 0,
                b_pos: 220,
                shared: 10,
                strand: 1,
            },
        ];
        let chains = chain_hits(&hits, 10, 20, 1000, 8, 0);
        assert_eq!(
            chains.len(),
            2,
            "forward and reverse strands must chain separately"
        );
        assert!(chains.iter().all(|c| c.seeds == 3));
    }

    #[test]
    fn chain_breaks_on_gap_and_span_filter() {
        let hits = vec![
            SeedHit {
                a_contig: 0,
                a_pos: 100,
                b_contig: 0,
                b_pos: 100,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 130,
                b_contig: 0,
                b_pos: 130,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 5000,
                b_contig: 0,
                b_pos: 5000,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 5030,
                b_contig: 0,
                b_pos: 5030,
                shared: 10,
                strand: 0,
            },
        ];
        // max_gap=1000 splits the two diagonal runs; each run has span 40 > 20.
        let chains = chain_hits(&hits, 10, 20, 1000, 8, 0);
        assert_eq!(chains.len(), 2);
        // min_span=50 drops both runs (span 40 each).
        let chains = chain_hits(&hits, 10, 50, 1000, 8, 0);
        assert!(chains.is_empty());
    }

    #[test]
    fn off_band_hit_does_not_split_main_chain() {
        let hits = vec![
            SeedHit {
                a_contig: 0,
                a_pos: 100,
                b_contig: 0,
                b_pos: 100,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 300,
                b_contig: 0,
                b_pos: 300,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 200,
                b_contig: 0,
                b_pos: 120,
                shared: 10,
                strand: 0,
            },
        ];
        // The diag-80 hit backtracks after the main diagonal ran: it must not
        // split the main chain (which covers 100..310), and its own tiny chain
        // (span 10) is dropped by min_span=200.
        let chains = chain_hits(&hits, 10, 200, 1000, 8, 0);
        assert_eq!(chains.len(), 1);
        assert_eq!((chains[0].a_start, chains[0].a_end), (100, 310));
        assert_eq!(chains[0].seeds, 2);
    }

    #[test]
    fn in_gap_hits_respect_band() {
        let hits = vec![
            SeedHit {
                a_contig: 0,
                a_pos: 100,
                b_contig: 0,
                b_pos: 100,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 150,
                b_contig: 0,
                b_pos: 150,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 200,
                b_contig: 0,
                b_pos: 190,
                shared: 10,
                strand: 0,
            },
        ];
        // Third hit (diag 10) is in-gap but outside band=8: ignored.
        let chains = chain_hits(&hits, 10, 5, 1000, 8, 0);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].seeds, 2);
        assert_eq!(chains[0].a_end, 160);
        // A wider band admits it.
        let chains = chain_hits(&hits, 10, 5, 1000, 20, 0);
        assert_eq!(chains[0].seeds, 3);
        assert_eq!(chains[0].a_end, 210);
    }

    #[test]
    fn merge_adjacent_chains_stitches_syntenic_blocks() {
        // Two colinear diagonal runs separated by a 1.2 kb insertion: greedy
        // chaining (max_gap 1000) splits them, the merge pass joins them back.
        let hits = vec![
            SeedHit {
                a_contig: 0,
                a_pos: 100,
                b_contig: 0,
                b_pos: 100,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 130,
                b_contig: 0,
                b_pos: 130,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 1300,
                b_contig: 0,
                b_pos: 1300,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 1330,
                b_contig: 0,
                b_pos: 1330,
                shared: 10,
                strand: 0,
            },
        ];
        let chains = chain_hits(&hits, 10, 20, 1000, 8, 0);
        assert_eq!(chains.len(), 2, "large gap must split chains");
        let chains = chain_hits(&hits, 10, 20, 1000, 8, 2000);
        assert_eq!(chains.len(), 1, "merge gap must stitch chains");
        assert_eq!((chains[0].a_start, chains[0].a_end), (100, 1340));
        assert_eq!(chains[0].seeds, 4);

        // A diagonal shift beyond the band must not merge.
        let shifted = vec![
            SeedHit {
                a_contig: 0,
                a_pos: 100,
                b_contig: 0,
                b_pos: 100,
                shared: 10,
                strand: 0,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 2000,
                b_contig: 0,
                b_pos: 1800,
                shared: 10,
                strand: 0,
            },
        ];
        let chains = chain_hits(&shifted, 10, 5, 1000, 8, 2000);
        assert_eq!(chains.len(), 2, "diag shift beyond band must not merge");
    }

    #[test]
    fn psl_block_coordinates() {
        let seq: Vec<u8> = (0..100u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let (ia, ib) = (build(&seq), build(&seq));
        let fwd = Chain {
            a_contig: 0,
            b_contig: 0,
            strand: 0,
            a_start: 10,
            a_end: 50,
            b_start: 20,
            b_end: 60,
            seeds: 1,
            diag: -10,
            diag_span: 0,
        };
        let p = chain_to_psl(&fwd, &ia, &ib);
        assert_eq!(p.strand, "+");
        assert_eq!((p.q_start, p.q_end), (20, 60));
        assert_eq!((p.t_start, p.t_end), (10, 50));

        let rev = Chain {
            strand: 1,
            b_start: 20,
            b_end: 60,
            ..fwd
        };
        let p = chain_to_psl(&rev, &ia, &ib);
        assert_eq!(p.strand, "-");
        // original query coordinates: reverse_range(20..60, len 100) = (40, 80)
        assert_eq!((p.q_start, p.q_end), (40, 80));
    }

    #[test]
    fn extend_chain_produces_scoring_psl() {
        let seq = pseudo_random_seq(500, 3);
        let (ia, ib) = (build(&seq), build(&seq));
        let a_seqs = vec![(String::from("c"), seq.clone())];
        let b_seqs = vec![(String::from("c"), seq)];
        let params = AlignParams::default();
        let psls = align_to_psl_ext(&ia, &ib, &params, &a_seqs, &b_seqs).unwrap();
        assert!(!psls.is_empty());
        assert!(
            psls.iter().all(|p| p.match_count + p.mismatch_count > 0),
            "extended blocks must carry alignment counts"
        );
        assert!(psls.iter().any(|p| p.strand == "+" && p.match_count > 0));
        assert!(psls
            .iter()
            .all(|p| p.q_start >= 0 && p.q_end <= p.q_size as i32));
        assert!(psls
            .iter()
            .all(|p| p.t_start >= 0 && p.t_end <= p.t_size as i32));
    }

    #[test]
    fn extend_chain_rc_query() {
        let seq = pseudo_random_seq(500, 9);
        let (ia, ir) = (build(&seq), build(&rev_comp(&seq).collect::<Vec<u8>>()));
        let a_seqs = vec![(String::from("c"), seq.clone())];
        let rc: Vec<u8> = rev_comp(&seq).collect();
        let b_seqs = vec![(String::from("c"), rc)];
        // Exact-only seeds: partial matches could pair unrelated same-strand
        // positions in this tiny random sequence and break the pure-RC check.
        let params = AlignParams {
            min_shared: Some(10),
            ..AlignParams::default()
        };
        let psls = align_to_psl_ext(&ia, &ir, &params, &a_seqs, &b_seqs).unwrap();
        assert!(!psls.is_empty());
        assert!(
            psls.iter().all(|p| p.strand == "-" && p.match_count > 0),
            "RC query must extend to minus-strand PSL"
        );
        assert!(psls
            .iter()
            .all(|p| p.q_start >= 0 && p.q_end <= p.q_size as i32));
    }

    #[test]
    fn extend_chain_windows_long_interval() {
        let seq = pseudo_random_seq(400, 21);
        let (ia, ib) = (build(&seq), build(&seq));
        let a_seqs = vec![(String::from("c"), seq.clone())];
        let b_seqs = vec![(String::from("c"), seq)];
        let params = AlignParams::default();
        let hits = merge_seed_hits(&ia, &ib, params.freq, 10).unwrap();
        let chains = chain_hits(
            &hits,
            ia.k as u32,
            params.min_span,
            params.max_gap,
            params.band,
            params.merge_gap,
        );
        assert_eq!(chains.len(), 1);
        // Tiny windows force the interval to be split into multiple records.
        let psls = extend_chain(&chains[0], &ia, &ib, &a_seqs, &b_seqs, 100, 90);
        assert!(
            psls.len() >= 3,
            "windowed extension must split: {}",
            psls.len()
        );
        let covered: u32 = psls.iter().map(|p| (p.q_end - p.q_start) as u32).sum();
        assert!(covered >= 300, "windowed coverage too low: {covered}");
        assert!(psls.iter().all(|p| p.match_count > 0));
    }
}
