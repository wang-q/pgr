//! Two-index merge alignment: seed hits -> anti-diagonal chains -> PSL blocks.

use super::dist::validate_compatible;
use super::PgiIndex;
use crate::libs::alignment::align_banded_local;
use crate::libs::alignment::coords::reverse_range_pair;
use crate::libs::alignment::wave::local_alignment;
use crate::libs::fmt::psl::Psl;
use crate::libs::nt::{complement, rev_comp};
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
    /// Chaining workflow: FastGA tube chaining or the default greedy chains.
    pub workflow: Workflow,
}

/// Chaining workflow.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Workflow {
    /// FastGA `align_contigs` tube chaining.
    Tube,
    /// Default greedy anti-diagonal chains.
    #[default]
    Greedy,
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
            workflow: Workflow::Greedy,
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

/// One FastGA-style tube: a seed chain with its anti-diagonal range
/// (`anti = a_pos + b_pos`) and diagonal band (`diag = a_pos - b_pos`).
#[derive(Debug, Clone, Copy)]
pub struct Tube {
    pub a_contig: u32,
    pub b_contig: u32,
    pub strand: u8,
    /// Anti-diagonal low end (bp).
    pub anti_low: i64,
    /// Anti-diagonal high end (bp, exclusive).
    pub anti_high: i64,
    /// Diagonal band low (bp).
    pub diag_min: i64,
    /// Diagonal band high (bp).
    pub diag_max: i64,
}

/// FastGA `align_contigs` tube chaining over seed hits.
///
/// Seeds are bucketed by diagonal (`diag >> 6`, width 64), and adjacent
/// buckets are merged in a-position order into tubes. A tube extends while
/// the next seed's anti-diagonal stays within `CHAIN_BREAK` (1000 bp) of the
/// current high; its anti coverage is the union of seed extents (`shared`),
/// and a tube with coverage at least `CHAIN_MIN` (85 bp) is emitted.
pub fn chain_tubes(hits: &[SeedHit], k: u32) -> Vec<Tube> {
    const BUCK: i64 = 64;
    const BREAK: i64 = 1000;
    const MIN_COV: u64 = 85;

    // Group by (a_contig, b_contig, strand); within a group sort by
    // (diagonal bucket, anti) so adjacent buckets are contiguous and each
    // bucket is anti-ordered (FastGA merges its stream by ipost = anti).
    let mut sorted: Vec<&SeedHit> = hits.iter().collect();
    // Pack the sort key into one u128: (a_contig, b_contig, strand,
    // diagonal bucket + offset, anti). The bucket offset keeps negative
    // diagonals orderable under unsigned packing.
    const BUCK_OFF: i64 = 1_000_000;
    sorted.par_sort_unstable_by_key(|h| {
        let diag = h.a_pos as i64 - h.b_pos as i64;
        let bucket = (diag.div_euclid(64) + BUCK_OFF) as u64;
        let anti = (h.a_pos as i64 + h.b_pos as i64) as u64;
        ((h.a_contig as u128) << 88)
            | ((h.b_contig as u128) << 72)
            | ((h.strand as u128) << 71)
            | ((bucket as u128) << 24)
            | anti as u128
    });

    let mut tubes = Vec::new();
    let mut start = 0usize;
    while start < sorted.len() {
        let g = sorted[start];
        let mut end = start + 1;
        while end < sorted.len()
            && sorted[end].a_contig == g.a_contig
            && sorted[end].b_contig == g.b_contig
            && sorted[end].strand == g.strand
        {
            end += 1;
        }
        tubes_for_group(
            &sorted[start..end],
            k,
            g.a_contig,
            g.b_contig,
            g.strand,
            BUCK,
            BREAK,
            MIN_COV,
            &mut tubes,
        );
        start = end;
    }
    tubes
}

#[allow(clippy::too_many_arguments)]
fn tubes_for_group(
    seeds: &[&SeedHit],
    k: u32,
    a_contig: u32,
    b_contig: u32,
    strand: u8,
    buck: i64,
    brk: i64,
    min_cov: u64,
    tubes: &mut Vec<Tube>,
) {
    // Bucket the seeds by diagonal.
    let mut buckets: Vec<(i64, Vec<&SeedHit>)> = Vec::new();
    for h in seeds {
        let diag = h.a_pos as i64 - h.b_pos as i64;
        let b = diag.div_euclid(buck);
        match buckets.last_mut() {
            Some((last, v)) if *last == b => v.push(h),
            _ => buckets.push((b, vec![*h])),
        }
    }

    // Process each adjacent bucket pair (c, c+1), merging by a_pos.
    for w in 0..buckets.len() {
        let (cb, b_seeds) = &buckets[w];
        let m_seeds: &[&SeedHit] = if w + 1 < buckets.len() && buckets[w + 1].0 == cb + 1 {
            &buckets[w + 1].1
        } else {
            &[]
        };
        // Merge by anti (ties prefer the lower bucket, like FastGA's
        // ipost-ordered merge; a_pos order breaks when the diagonal drifts
        // inside a bucket and starves the coverage accounting).
        let mut bi = 0usize;
        let mut mi = 0usize;
        let mut alow: i64 = i64::MAX;
        let mut ahgh: i64 = -brk;
        let mut cov: u64 = 0;
        let mut dgmin = 2 * buck;
        let mut dgmax = 0i64;
        while bi < b_seeds.len() || mi < m_seeds.len() {
            let anti_b = b_seeds.get(bi).map(|h| h.a_pos as i64 + h.b_pos as i64);
            let anti_m = m_seeds.get(mi).map(|h| h.a_pos as i64 + h.b_pos as i64);
            let (h, side) = match (anti_b, anti_m) {
                (Some(ab), Some(am)) if am < ab => (m_seeds[mi], 2u8),
                (Some(_), _) => (b_seeds[bi], 1u8),
                (None, Some(_)) => (m_seeds[mi], 2u8),
                (None, None) => unreachable!(),
            };
            if side == 1 {
                bi += 1;
            } else {
                mi += 1;
            }
            let anti = h.a_pos as i64 + h.b_pos as i64;
            let ext = (h.shared as i64).clamp(1, k as i64);
            let dg = h.a_pos as i64 - h.b_pos as i64 - cb * buck;
            if anti < ahgh + brk {
                let cps = anti + ext;
                if cps > ahgh {
                    if anti >= ahgh {
                        cov += ext as u64;
                    } else {
                        cov += (cps - ahgh) as u64;
                    }
                    ahgh = cps;
                }
                dgmin = dgmin.min(dg);
                dgmax = dgmax.max(dg);
                alow = alow.min(anti);
            } else {
                if cov >= min_cov {
                    tubes.push(Tube {
                        a_contig,
                        b_contig,
                        strand,
                        anti_low: alow,
                        anti_high: ahgh,
                        diag_min: cb * buck + dgmin,
                        diag_max: cb * buck + dgmax,
                    });
                }
                alow = anti;
                ahgh = anti + ext;
                cov = ext as u64;
                dgmin = dg;
                dgmax = dg;
            }
        }
        if cov >= min_cov {
            tubes.push(Tube {
                a_contig,
                b_contig,
                strand,
                anti_low: alow,
                anti_high: ahgh,
                diag_min: cb * buck + dgmin,
                diag_max: cb * buck + dgmax,
            });
        }
    }
}

/// FastGA `BUCK_ANTI`: tube processing slides the mid-line by 128 bp.
const BUCK_ANTI: i64 = 128;
/// FastGA `alnMin = ALIGN_MIN - 50` (default `-l 100`).
const TUBE_MIN_LEN: i64 = 50;
/// FastGA `alnRate = ALIGN_RATE + 0.05` (default `-i 0.7`, mismatch = 1).
const TUBE_MIN_RATE: f64 = 0.35;

/// Extend one tube with FastGA `align_contigs` semantics: slide the mid-line
/// anti-diagonal by `BUCK_ANTI` through the tube box and call the mid-line
/// wave `Local_Alignment` (full contig sequences, tube diagonal band) at each
/// position. `alast` carries the furthest anti already aligned by an earlier
/// tube of the same (contig pair, strand) group.
#[allow(clippy::too_many_arguments)]
fn extend_tube(
    tube: &Tube,
    alast: &mut i64,
    a: &PgiIndex,
    b: &PgiIndex,
    a_seqs: &[(String, Vec<u8>)],
    b_seqs: &[(String, Vec<u8>)],
    a_revs: &[Vec<u8>],
    b_revs: &[Vec<u8>],
    b_comps: &[Vec<u8>],
) -> Vec<Psl> {
    let Some((_, a_int)) = a_seqs.get(tube.a_contig as usize) else {
        return Vec::new();
    };
    let Some((_, b_int)) = b_seqs.get(tube.b_contig as usize) else {
        return Vec::new();
    };
    // Query contig in orientation space (RC for minus strand), once.
    let q: Vec<u8> = if tube.strand == 0 {
        b_int.to_vec()
    } else {
        rev_comp(b_int).collect()
    };
    let rt = &a_revs[tube.a_contig as usize];
    let rq = if tube.strand == 0 {
        &b_revs[tube.b_contig as usize]
    } else {
        &b_comps[tube.b_contig as usize]
    };
    let mut alow = tube.anti_low.max(*alast);
    let ahgh = tube.anti_high - BUCK_ANTI;
    let mut dgmin = tube.diag_min;
    let dgmax = tube.diag_max;
    if ahgh <= *alast {
        return Vec::new(); // already covered by an earlier tube
    }
    // Large-tube homology gate: tubes whose region never reaches ~55%
    // identity on any sampled diagonal of the band can only produce rejected
    // wave calls (91% of all calls come from such tubes), so skip them.
    // Only large tubes are checked (small tubes are cheap and thin conserved
    // regions there are caught by the waves).
    if ahgh - alow > 10_000 {
        let a_len = a_int.len() as i64;
        let b_len = q.len() as i64;
        let mut best = 0i64;
        for dgi in 0..9 {
            let dg = dgmin + (dgmax - dgmin) * dgi / 8;
            let mut m = 0i64;
            let mut anti = alow;
            while anti < ahgh + 128 {
                let a = (anti + dg) / 2;
                let b = (anti - dg) / 2;
                if a >= 0 && b >= 0 && a < a_len && b < b_len && q[b as usize] == a_int[a as usize]
                {
                    m += 1;
                }
                if anti >= alow + 128 {
                    let a0 = (anti - 128 + dg) / 2;
                    let b0 = (anti - 128 - dg) / 2;
                    if a0 >= 0
                        && b0 >= 0
                        && a0 < a_len
                        && b0 < b_len
                        && q[b0 as usize] == a_int[a0 as usize]
                    {
                        m -= 1;
                    }
                    best = best.max(m);
                }
                anti += 2;
            }
        }
        if best < 32 {
            return Vec::new();
        }
    }
    let mut out = Vec::new();
    while alow < ahgh {
        let mut amid = alow + BUCK_ANTI;
        if amid > ahgh {
            amid = ahgh;
        }
        if amid + dgmin < 0 {
            dgmin = -amid;
            if dgmin > dgmax {
                break;
            }
        }
        if let Some(aln) = local_alignment(&q, a_int, rt, rq, dgmin, dgmax, amid) {
            let rlen = (aln.t_end - aln.t_start) as i64;
            if rlen >= TUBE_MIN_LEN && TUBE_MIN_RATE * rlen as f64 >= aln.diffs as f64 {
                let strand = if tube.strand == 0 { "+" } else { "-" };
                if let Some(psl) = Psl::from_align(
                    &b.contigs[tube.b_contig as usize].0,
                    b.contigs[tube.b_contig as usize].1 as u32,
                    aln.q_start as i32,
                    aln.q_end as i32,
                    &String::from_utf8_lossy(&aln.q_aln),
                    &a.contigs[tube.a_contig as usize].0,
                    a.contigs[tube.a_contig as usize].1 as u32,
                    aln.t_start as i32,
                    aln.t_end as i32,
                    &String::from_utf8_lossy(&aln.t_aln),
                    strand,
                ) {
                    out.push(psl);
                }
            }
            let eant = (aln.t_end + aln.q_end) as i64;
            alow = if eant <= alow { amid } else { eant };
        } else {
            alow = amid;
        }
    }
    *alast = alow;
    out
}

/// Drop blocks that overlap an earlier block of the same (contig pair,
/// strand) on both axes by at least 95% of their own span.
///
/// Adjacent diagonal-bucket tubes align the same region twice; FastGA's
/// sequential `alast` skip prevents that, the parallel tube pass needs this
/// post-filter instead. The threshold is high on purpose: blocks from the
/// same tube's successive calls can overlap ~87% while one extends the
/// other (dropping the extension would lose real coverage).
fn dedupe_contained(blocks: &mut Vec<Psl>) {
    if blocks.len() < 2 {
        return;
    }
    blocks.sort_by_key(|p| {
        (
            p.t_name.clone(),
            p.q_name.clone(),
            p.strand.clone(),
            p.t_start,
            p.q_start,
        )
    });
    let mut kept: Vec<Psl> = Vec::with_capacity(blocks.len());
    for b in blocks.drain(..) {
        let dup = kept.iter().rev().take(64).any(|k| {
            k.t_name == b.t_name
                && k.q_name == b.q_name
                && k.strand == b.strand
                && overlap_frac(k.t_start, k.t_end, b.t_start, b.t_end) >= 0.95
                && overlap_frac(k.q_start, k.q_end, b.q_start, b.q_end) >= 0.95
        });
        if !dup {
            kept.push(b);
        }
    }
    *blocks = kept;
}

/// Fraction of `[b1, b2)` covered by `[a1, a2)`.
fn overlap_frac(a1: i32, a2: i32, b1: i32, b2: i32) -> f64 {
    let own = b2 - b1;
    if own <= 0 {
        return 0.0;
    }
    let ov = a2.min(b2) - a1.max(b1);
    (ov.max(0) as f64) / own as f64
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
    match params.workflow {
        Workflow::Tube => anyhow::bail!(
            "the tube workflow needs --ref-seq/--query-seq (wave extension requires sequences)"
        ),
        Workflow::Greedy => {
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
    }
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
    if params.workflow == Workflow::Tube {
        // Tubes are independent alignment tasks; the `alast` overlap skip of
        // FastGA's sequential scan is replaced by a parallel pass plus a
        // containment dedup on the emitted blocks.
        let tubes = chain_tubes(&hits, a.k as u32);
        // FastGA's chains break at ~100 kb (its adaptive seed stream is
        // sparser); cap oversized tubes so one giant tube cannot serialize
        // the mid-line slides.
        const TUBE_ANTI_CAP: i64 = 40_000;
        let tubes: Vec<Tube> = tubes
            .into_iter()
            .flat_map(|t| {
                let span = t.anti_high - t.anti_low;
                if span <= TUBE_ANTI_CAP {
                    vec![t]
                } else {
                    let n = (span + TUBE_ANTI_CAP - 1) / TUBE_ANTI_CAP;
                    (0..n)
                        .map(|i| {
                            let mut piece = t;
                            piece.anti_low = t.anti_low + i * TUBE_ANTI_CAP;
                            piece.anti_high =
                                (t.anti_low + (i + 1) * TUBE_ANTI_CAP).min(t.anti_high);
                            piece
                        })
                        .collect()
                }
            })
            .collect();
        let a_revs: Vec<Vec<u8>> = a_seqs
            .iter()
            .map(|(_, s)| s.iter().rev().copied().collect())
            .collect();
        let b_revs: Vec<Vec<u8>> = b_seqs
            .iter()
            .map(|(_, s)| s.iter().rev().copied().collect())
            .collect();
        let b_comps: Vec<Vec<u8>> = b_seqs
            .iter()
            .map(|(_, s)| complement(s).collect())
            .collect();
        let records: Vec<Vec<Psl>> = tubes
            .par_iter()
            .map(|t| {
                let mut alast = 0i64;
                extend_tube(
                    t, &mut alast, a, b, a_seqs, b_seqs, &a_revs, &b_revs, &b_comps,
                )
            })
            .collect();
        let mut out: Vec<Psl> = records.into_iter().flatten().collect();
        dedupe_contained(&mut out);
        return Ok(out);
    }
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
    fn tubes_join_colinear_seeds() {
        let hits = vec![
            SeedHit {
                a_contig: 0,
                a_pos: 100,
                b_contig: 0,
                b_pos: 100,
                strand: 0,
                shared: 40,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 130,
                b_contig: 0,
                b_pos: 130,
                strand: 0,
                shared: 40,
            },
            SeedHit {
                a_contig: 0,
                a_pos: 160,
                b_contig: 0,
                b_pos: 160,
                strand: 0,
                shared: 40,
            },
        ];
        let tubes = chain_tubes(&hits, 40);
        assert_eq!(tubes.len(), 1);
        assert_eq!((tubes[0].anti_low, tubes[0].anti_high), (200, 360));
        assert_eq!((tubes[0].diag_min, tubes[0].diag_max), (0, 0));
    }

    #[test]
    fn tube_merge_uses_anti_order_when_diagonal_drifts() {
        // Seeds in one diagonal bucket whose anti order differs from a_pos
        // order (the diagonal drifts inside the bucket): the tube must start
        // at the smallest anti and accumulate full coverage. Merging by a_pos
        // starves the coverage (cov < CHAIN_MIN) and drops the tube.
        let mut hits = Vec::new();
        // diag -4492 at a_pos 112701 (anti 229894), sorted before -4490.
        hits.push(SeedHit {
            a_contig: 0,
            a_pos: 112_701,
            b_contig: 0,
            b_pos: 117_193,
            strand: 0,
            shared: 40,
        });
        // diag -4490, dense run starting at a_pos 111544 (anti 227578).
        for i in (0..200).step_by(2) {
            hits.push(SeedHit {
                a_contig: 0,
                a_pos: 111_544 + i,
                b_contig: 0,
                b_pos: 116_034 + i,
                strand: 0,
                shared: 40,
            });
        }
        let tubes = chain_tubes(&hits, 40);
        assert!(!tubes.is_empty(), "dense collinear seeds must form a tube");
        let t = &tubes[0];
        assert!(
            t.anti_low <= 227_578,
            "tube must start at the smallest anti, got {}",
            t.anti_low
        );
        assert!(
            t.anti_high >= 227_578 + 200 + 40,
            "tube must cover the dense run, got {}",
            t.anti_high
        );
    }

    #[test]
    fn dedupe_keeps_blocks_that_extend_earlier_ones() {
        // Two blocks from successive calls of one tube: ~87% overlap on both
        // axes while the later one extends further (the extension carries
        // real coverage). The dedupe must keep both.
        let mk = |t_start: i32, t_end: i32, q_start: i32, q_end: i32| Psl {
            q_name: String::from("q"),
            q_size: 5_000_000,
            q_start,
            q_end,
            t_name: String::from("t"),
            t_size: 5_000_000,
            t_start,
            t_end,
            strand: String::from("+"),
            block_count: 1,
            block_sizes: vec![(t_end - t_start) as u32],
            q_starts: vec![q_start as u32],
            t_starts: vec![t_start as u32],
            ..Default::default()
        };
        let short = mk(4_094_464, 4_115_709, 4_894_644, 4_915_804);
        let extended = mk(4_094_559, 4_118_787, 4_894_644, 4_918_971);
        let mut blocks = vec![extended, short];
        dedupe_contained(&mut blocks);
        assert_eq!(blocks.len(), 2, "the extended block must not be dropped");
    }

    #[test]
    fn tubes_break_on_anti_gap() {
        let mk = |a: u32, b: u32| SeedHit {
            a_contig: 0,
            a_pos: a,
            b_contig: 0,
            b_pos: b,
            strand: 0,
            shared: 40,
        };
        let hits = vec![
            mk(100, 100),
            mk(130, 130),
            mk(160, 160),
            mk(2000, 2000),
            mk(2030, 2030),
            mk(2060, 2060),
        ];
        let tubes = chain_tubes(&hits, 40);
        assert_eq!(tubes.len(), 2, "anti gap must split tubes");
        assert_eq!(tubes[0].anti_high, 360);
        assert_eq!(tubes[1].anti_low, 4000);
    }

    #[test]
    fn isolated_seed_produces_no_tube() {
        let hits = vec![SeedHit {
            a_contig: 0,
            a_pos: 100,
            b_contig: 0,
            b_pos: 100,
            strand: 0,
            shared: 40,
        }];
        assert!(chain_tubes(&hits, 40).is_empty(), "coverage 40 < CHAIN_MIN");
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
