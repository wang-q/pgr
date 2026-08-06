//! Two-index merge alignment: seed hits -> anti-diagonal chains -> PSL blocks.

use super::dist::validate_compatible;
use super::wave::{local_alignment, TrimSpec};
use super::{unpack_position, PgiIndex, PgiQuery, PgiStream};
use crate::libs::alignment::coords::{reverse_range, reverse_range_pair};
use crate::libs::fmt::psl::Psl;
use crate::libs::nt::{complement, rc_key, rev_comp};
use rayon::prelude::*;
use std::io::Read;

/// Alignment parameters; defaults follow FastGA (`-f 10`, `-c 85`, `-s 1000`).
#[derive(Debug, Clone, Copy)]
pub struct AlignParams {
    /// K-mers occurring at least this often on either side are skipped as
    /// seeds (FastGA's frequency cutoff).
    pub freq: u32,
    /// Minimum shared seed length (bp); `None` selects the workflow default
    /// (FastGA's plen floor of 12).
    pub min_shared: Option<usize>,
}

impl Default for AlignParams {
    fn default() -> Self {
        Self {
            freq: 10,
            min_shared: None,
        }
    }
}

/// A seed hit between two indexes, orientation resolved into a shared
/// (contig_a, contig_b, strand) coordinate space.
#[derive(Debug, Clone, Copy)]
pub struct SeedHit {
    /// Reference (a) contig id.
    pub a_contig: u16,
    /// Reference position, forward coordinates.
    pub a_pos: u32,
    /// Query (b) contig id.
    pub b_contig: u16,
    /// Query position in orientation space (RC space when `strand` is 1).
    pub b_pos: u32,
    /// Shared prefix length (bp) of the two k-mers.
    pub shared: u16,
    /// 0 = forward, 1 = reverse (query window is the RC of the reference).
    pub strand: u8,
}

/// Cursor for the 归并式 (sequential) merge over `b` within one ascending
/// batch of `a` entries.
///
/// The floor window (entries sharing the `min_shared`-base prefix, the
/// widest window used) moves monotonically as `a` ascends, so its bounds
/// advance by [`PgiQuery::entry_lower_bound_ge`] scanning instead of a
/// per-entry binary search. All narrower windows (start / maximal-prefix)
/// are sub-windows of the floor window, so they are located by advancing
/// within it. `first` forces a batch-header binary search.
#[derive(Debug, Clone, Copy)]
struct MergeCursor {
    /// Previous floor-window lower bound (a group start in `b`).
    f0: usize,
    /// Previous floor-window upper bound.
    f1: usize,
    /// Whether the next surviving entry is the batch's first (binary search).
    first: bool,
}

impl Default for MergeCursor {
    fn default() -> Self {
        Self {
            f0: 0,
            f1: 0,
            first: true,
        }
    }
}

/// Merge two sorted k-mer indexes, emitting adaptamer-style seed hits for
/// every k-mer pair sharing at least `min_shared` bases (mirroring FastGA's
/// lcp-driven merge with FastGA's adaptamer seed selection: each `a` entry
/// seeds only its *longest* shared prefix against `b`, and the frequency
/// filter applies to the extended range at that length (not the fixed
/// `min_shared` floor window). Exact k-mers are the `min_shared == k` case.
///
/// Both-strand entries are emitted by `pgi build`, so shared keys with equal
/// strand flags mean a forward hit, and differing flags a reverse hit.
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
    // The per-entry prefix query against `b` is independent, so split the
    // `a` entries into chunks and merge them in parallel (the chaining passes
    // re-sort the hits, so the output order is free).
    let hits: Vec<SeedHit> = a
        .entries
        .par_chunks(4096)
        .map(|ents| -> anyhow::Result<Vec<SeedHit>> {
            let mut hits = Vec::new();
            let mut prev_kmer = None;
            let mut cur = MergeCursor::default();
            for ea in ents {
                let ap = &a.positions[ea.pos_start as usize..(ea.pos_start + ea.freq) as usize];
                hits.extend(emit_entry_hits(
                    ea.kmer, ea.freq, ap, b, freq, min_shared, k, prev_kmer, &mut cur,
                )?);
                prev_kmer = Some(ea.kmer);
            }
            Ok(hits)
        })
        .collect::<anyhow::Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect();
    Ok(hits)
}

/// FastGA-style seed merge with the reference index streamed from disk: `a`
/// entries are read in batches and merged in parallel against the query
/// index `b` (resident or memory-mapped), so neither index is materialized
/// in full.
pub fn merge_seed_hits_from_stream<R: Read + Send, B: PgiQuery + Sync>(
    a: &mut PgiStream<R>,
    b: &B,
    freq: u32,
    min_shared: usize,
) -> anyhow::Result<Vec<SeedHit>> {
    validate_compatible(
        &PgiIndex {
            k: a.header().k,
            smer: a.header().smer,
            window: a.header().window,
            contigs: Vec::new(),
            entries: Vec::new(),
            positions: Vec::new(),
        },
        b,
    )?;
    let k = a.header().k;
    anyhow::ensure!(
        min_shared >= 1 && min_shared <= k,
        "min_shared must be in 1..={k}"
    );
    // One rayon task per batch keeps a handful of batches in flight instead
    // of the whole reference index.
    const BATCH_ENTRIES: usize = 8192;
    let batches =
        std::iter::from_fn(|| Some(a.next_batch(BATCH_ENTRIES))).take_while(|r| match r {
            Ok(v) => !v.is_empty(),
            Err(_) => true,
        });
    let parts: Vec<Vec<SeedHit>> = batches
        .par_bridge()
        .map(|batch| {
            let mut hits = Vec::new();
            let mut prev_kmer = None;
            let mut cur = MergeCursor::default();
            for (ea, poss) in batch? {
                hits.extend(emit_entry_hits(
                    ea.kmer, ea.freq, &poss, b, freq, min_shared, k, prev_kmer, &mut cur,
                )?);
                prev_kmer = Some(ea.kmer);
            }
            Ok(hits)
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    Ok(parts.into_iter().flatten().collect())
}

/// Emit FastGA-style seed hits for one `a` entry against the full `b` index.
///
/// `a_positions` is the entry's packed positions (from the full index or a
/// streamed batch). Each entry seeds only its longest shared prefix, and the
/// frequency filter applies to the extended range at that length.
#[allow(clippy::too_many_arguments)]
fn emit_entry_hits<B: PgiQuery>(
    ea_kmer: u128,
    ea_freq: u32,
    a_positions: &[u64],
    b: &B,
    freq: u32,
    min_shared: usize,
    k: usize,
    prev_kmer: Option<u128>,
    cur: &mut MergeCursor,
) -> anyhow::Result<Vec<SeedHit>> {
    let mut hits = Vec::new();
    // FastGA excludes k-mers whose count is >= the frequency cutoff on both
    // sides at index-build time (GIXmake: "only k-mers whose count is less
    // than the adaptamer frequency cutoff are in the index"); the query-side
    // extended-range filter below drops `occ >= freq` for the same reason.
    // Dropping only `> freq` here would let a k-mer occurring exactly `freq`
    // times on the reference side seed hits that FastGA never emits.
    if ea_freq >= freq {
        return Ok(hits);
    }
    // FastGA stores each position under its canonical orientation only; the
    // reverse-complement key would duplicate every physical hit (both strands
    // are emitted by `pgi build`).
    if ea_kmer > rc_key(ea_kmer, k) {
        return Ok(hits);
    }
    let k_bits = 2 * k;
    let mask = if k_bits >= 128 {
        u128::MAX
    } else {
        (1u128 << k_bits) - 1
    };
    // Prefix range of `b` entries sharing the first `len` bases with the key.
    let window_bounds = |len: usize| -> (u128, u128) {
        let r = 1u128 << (k_bits - 2 * len);
        let lo = ea_kmer & !(r - 1) & mask;
        // The prefix boundary `lo + r` can equal `2^(2k)`; for k=64 that is
        // `2^128`, which `u128` cannot represent. Saturate instead of
        // overflowing (a debug panic / release wrap): the top-boundary key
        // `u128::MAX` is the only one the saturated range excludes, an
        // extreme all-T 64-mer edge that is not a real seed.
        let hi = lo.saturating_add(r);
        (lo, hi)
    };
    // Floor window (widest, at `min_shared`): monotonic across ascending `a`
    // entries, so its bounds advance sequentially from the cursor except at
    // the batch's first surviving entry (binary search).
    let (flo, fhi) = window_bounds(min_shared);
    let (f0, f1) = if cur.first {
        let (f0, f1) = b.entry_range(flo, fhi);
        cur.f0 = f0;
        cur.f1 = f1;
        cur.first = false;
        (f0, f1)
    } else {
        let f0 = b.entry_lower_bound_ge(flo, cur.f0);
        let f1 = b.entry_lower_bound_ge(fhi, cur.f0.max(cur.f1));
        cur.f0 = f0;
        cur.f1 = f1;
        (f0, f1)
    };
    // Lcp propagation (FastGA `new_merge_thread`'s `vlcp` table): an entry
    // shares at least `lcp(prev, cur)` bases with every match its predecessor
    // had, so the scan can start from that (usually much narrower) prefix
    // window instead of the `min_shared` floor. When the predecessor's prefix
    // does not carry over (empty narrow window), fall back to the floor
    // window: the longest match may still reach `min_shared`.
    let start = prev_kmer
        .map(|pk| shared_prefix(pk, ea_kmer, k).max(min_shared as u32))
        .unwrap_or(min_shared as u32);
    let (mut j0, mut j) = if start as usize == min_shared {
        (f0, f1)
    } else {
        // A narrower prefix window is a sub-window of the floor window; find
        // it by advancing within `[f0, f1)` instead of a full binary search.
        let (slo, shi) = window_bounds(start as usize);
        let s0 = b.entry_lower_bound_ge(slo, f0);
        let s1 = b.entry_lower_bound_ge(shi, s0);
        (s0, s1)
    };
    if j0 == j && start as usize > min_shared {
        (j0, j) = (f0, f1);
    }
    if j == j0 {
        return Ok(hits);
    }
    // Maximal shared prefix over the scan window (FastGA extends each entry
    // to its longest match before filtering).
    let max_shared_over = |mut i: usize, j: usize| -> u32 {
        let mut m = 0u32;
        while i < j {
            // Entries with `freq >= cutoff` are absent from FastGA's GIX
            // index; treat them as absent here too (they must not influence
            // the maximal shared prefix).
            if b.entry_freq(i) >= freq {
                i = b.entry_next(i);
                continue;
            }
            m = m.max(shared_prefix(ea_kmer, b.entry_kmer(i), k));
            i = b.entry_next(i);
        }
        m
    };
    let mut m = max_shared_over(j0, j);
    if m < min_shared as u32 && start as usize > min_shared {
        // FastGA's GIX index omits `>= freq` k-mers, so its lcp-narrowed
        // window is empty whenever no under-frequency entry shares the lcp
        // prefix, and it falls back to the floor window. pgr keeps those
        // entries in the index, so the narrowed window can be non-empty yet
        // hold only high-frequency k-mers, leaving `m < min_shared` even
        // though an under-frequency match exists in the floor window below
        // the lcp. Recover the same way the empty-window path does.
        (j0, j) = (f0, f1);
        m = max_shared_over(j0, j);
    }
    if m < min_shared as u32 {
        return Ok(hits);
    }
    // Restrict to the maximal-prefix range (its occurrence count must stay
    // under `freq` -- extended-range filter, not the floor window). It is a
    // sub-window of the scan range, so locate it by advancing within it.
    let (mlo, mhi) = window_bounds(m as usize);
    let m0 = b.entry_lower_bound_ge(mlo, j0);
    let m1 = b.entry_lower_bound_ge(mhi, m0);
    let mut occ: u32 = 0;
    let mut i = m0;
    while i < m1 {
        let f = b.entry_freq(i);
        if f < freq {
            occ += f;
        }
        i = b.entry_next(i);
    }
    if occ == 0 || occ >= freq {
        return Ok(hits);
    }
    let mut i = m0;
    while i < m1 {
        if b.entry_freq(i) < freq {
            for &arec in a_positions {
                let (ac, apos, astrand) = unpack_position(arec);
                for brec in b.entry_positions(i) {
                    let (bc, bpos, bstrand) = unpack_position(brec);
                    crate::libs::pgi::validate_record(bc, bpos, k, b.contigs())?;
                    let fwd = astrand == bstrand;
                    let b_len = b.contigs()[bc as usize].1;
                    let oriented = if fwd {
                        bpos as u64
                    } else {
                        b_len - k as u64 - bpos as u64
                    };
                    hits.push(SeedHit {
                        a_contig: ac as u16,
                        a_pos: apos,
                        b_contig: bc as u16,
                        b_pos: oriented as u32,
                        shared: m as u16,
                        strand: u8::from(!fwd),
                    });
                }
            }
        }
        i = b.entry_next(i);
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

/// Effective minimum shared seed length: `params.min_shared` or FastGA's
/// adaptamer plen floor (12), capped at `k`.
fn effective_min_shared(k: usize, params: &AlignParams) -> usize {
    match params.min_shared {
        Some(v) => v.min(k),
        // FastGA's adaptamer merge emits seeds with plen in 12..=k; the
        // floor recovers indel-shifted regions.
        None => 12.min(k),
    }
}

/// Current peak RSS in MB, read from `/proc/self/status` (Linux only).
fn peak_rss_mb() -> u64 {
    let Ok(s) = std::fs::read_to_string("/proc/self/status") else {
        return 0;
    };
    s.lines()
        .find(|l| l.starts_with("VmHWM:"))
        .and_then(|l| l.split_whitespace().nth(1))
        .and_then(|v| v.parse::<u64>().ok())
        .map(|kb| kb / 1024)
        .unwrap_or(0)
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
    /// Seed span start on the reference.
    pub a_start: u32,
    /// Seed span end (exclusive, includes `k`) on the reference.
    pub a_end: u32,
    /// Seed span start on the query, in orientation space.
    pub b_start: u32,
    /// Seed span end (exclusive, includes `k`) on the query, orientation space.
    pub b_end: u32,
}

/// FastGA `align_contigs` tube chaining over seed hits.
///
/// Seeds are bucketed by diagonal (`diag >> 6`, width 64), and adjacent
/// buckets are merged in a-position order into tubes. A tube extends while
/// the next seed's anti-diagonal stays within `CHAIN_BREAK` (2000 bp, FastGA's
/// internal value: `-s 1000` is doubled into anti-diagonal space) of the
/// current high; its anti coverage is the union of seed extents (`shared`,
/// single-axis, so `CHAIN_MIN` 85 == FastGA's 170 in anti space), and a tube
/// with coverage at least `CHAIN_MIN` (85 bp) is emitted.
pub fn chain_tubes(hits: &[SeedHit], k: u32) -> Vec<Tube> {
    const BUCK: i64 = 64;
    const BREAK: i64 = 2000;
    const MIN_COV: u64 = 85;

    // Group by (a_contig, b_contig, strand); within a group sort by
    // (diagonal bucket, anti) so adjacent buckets are contiguous and each
    // bucket is anti-ordered (FastGA merges its stream by ipost = anti).
    // Pack the sort key into one u128: (a_contig, b_contig, strand,
    // diagonal bucket + offset, anti). The bucket offset keeps negative
    // diagonals orderable under unsigned packing: diag = a_pos - b_pos
    // spans -2^32+1..2^32-1, so diag/64 spans -2^26..2^26-1 and the offset
    // must be at least 2^26 — 1,000,000 only covered diagonals down to
    // ~-64 Mb, and deeper negatives wrapped through `as u64`, spilling
    // 0xFFFF into the strand/contig fields and interleaving contigs.
    // `anti` needs 40 bits (a_pos + b_pos exceeds 2^24 on any genome larger
    // than ~8 Mb and 2^32 on >2.1 Gb contig pairs) and `bucket` 32;
    // overflowing into the bucket field interleaves diagonal buckets and
    // fragments tubes.
    const BUCK_OFF: i64 = 1 << 26;
    // Precompute the packed keys into a flat key array and radix-sort the
    // index permutation alongside (the MSD radix keeps a `u128` key array
    // plus a `u32` index array, ~40% smaller than `(u128, u32)` tuples).
    let mut keys: Vec<u128> = Vec::with_capacity(hits.len());
    for h in hits {
        let diag = h.a_pos as i64 - h.b_pos as i64;
        let bucket = (diag.div_euclid(64) + BUCK_OFF) as u64;
        let anti = (h.a_pos as i64 + h.b_pos as i64) as u64;
        keys.push(
            ((h.a_contig as u128) << 89)
                | ((h.b_contig as u128) << 73)
                | ((h.strand as u128) << 72)
                | ((bucket as u128) << 40)
                | anti as u128,
        );
    }
    let mut order: Vec<u32> = (0..hits.len() as u32).collect();
    crate::libs::ds::radix_sort::radix_sort_u128_par(&mut keys, &mut order, 112);

    let mut tubes = Vec::new();
    let mut start = 0usize;
    while start < order.len() {
        let g = &hits[order[start] as usize];
        let mut end = start + 1;
        while end < order.len()
            && hits[order[end] as usize].a_contig == g.a_contig
            && hits[order[end] as usize].b_contig == g.b_contig
            && hits[order[end] as usize].strand == g.strand
        {
            end += 1;
        }
        tubes_for_group(
            &order[start..end],
            hits,
            k,
            g.a_contig as u32,
            g.b_contig as u32,
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
    seeds: &[u32],
    hits: &[SeedHit],
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
    let mut buckets: Vec<(i64, Vec<u32>)> = Vec::new();
    for &s in seeds {
        let h = &hits[s as usize];
        let diag = h.a_pos as i64 - h.b_pos as i64;
        let b = diag.div_euclid(buck);
        match buckets.last_mut() {
            Some((last, v)) if *last == b => v.push(s),
            _ => buckets.push((b, vec![s])),
        }
    }

    // Process each adjacent bucket pair (c, c+1), merging by a_pos.
    // Each adjacent bucket pair is an independent merge; process them in
    // parallel and concatenate (rayon preserves the pair order).
    let pair_tubes: Vec<Vec<Tube>> = (0..buckets.len())
        .into_par_iter()
        .map(|w| {
            let (cb, b_seeds) = &buckets[w];
            let m_seeds: &[u32] = if w + 1 < buckets.len() && buckets[w + 1].0 == cb + 1 {
                &buckets[w + 1].1
            } else {
                &[]
            };
            let mut local = Vec::new();
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
            let mut a_lo = i64::MAX;
            let mut a_hi = 0i64;
            let mut b_lo = i64::MAX;
            let mut b_hi = 0i64;
            while bi < b_seeds.len() || mi < m_seeds.len() {
                let anti_b = b_seeds
                    .get(bi)
                    .map(|&s| hits[s as usize].a_pos as i64 + hits[s as usize].b_pos as i64);
                let anti_m = m_seeds
                    .get(mi)
                    .map(|&s| hits[s as usize].a_pos as i64 + hits[s as usize].b_pos as i64);
                let (h, side) = match (anti_b, anti_m) {
                    (Some(ab), Some(am)) if am < ab => (&hits[m_seeds[mi] as usize], 2u8),
                    (Some(_), _) => (&hits[b_seeds[bi] as usize], 1u8),
                    (None, Some(_)) => (&hits[m_seeds[mi] as usize], 2u8),
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
                    a_lo = a_lo.min(h.a_pos as i64);
                    a_hi = a_hi.max(h.a_pos as i64 + k as i64);
                    b_lo = b_lo.min(h.b_pos as i64);
                    b_hi = b_hi.max(h.b_pos as i64 + k as i64);
                } else {
                    if cov >= min_cov {
                        local.push(Tube {
                            a_contig,
                            b_contig,
                            strand,
                            anti_low: alow,
                            anti_high: ahgh,
                            diag_min: cb * buck + dgmin,
                            diag_max: cb * buck + dgmax,
                            a_start: a_lo as u32,
                            a_end: a_hi as u32,
                            b_start: b_lo as u32,
                            b_end: b_hi as u32,
                        });
                    }
                    alow = anti;
                    ahgh = anti + ext;
                    cov = ext as u64;
                    dgmin = dg;
                    dgmax = dg;
                    a_lo = h.a_pos as i64;
                    a_hi = h.a_pos as i64 + k as i64;
                    b_lo = h.b_pos as i64;
                    b_hi = h.b_pos as i64 + k as i64;
                }
            }
            if cov >= min_cov {
                local.push(Tube {
                    a_contig,
                    b_contig,
                    strand,
                    anti_low: alow,
                    anti_high: ahgh,
                    diag_min: cb * buck + dgmin,
                    diag_max: cb * buck + dgmax,
                    a_start: a_lo as u32,
                    a_end: a_hi as u32,
                    b_start: b_lo as u32,
                    b_end: b_hi as u32,
                });
            }
            local
        })
        .collect();
    tubes.extend(pair_tubes.into_iter().flatten());
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
    b: &impl PgiQuery,
    a_seqs: &[(String, Vec<u8>)],
    b_seqs: &[(String, Vec<u8>)],
    a_revs: &[Vec<u8>],
    b_revs: &[Vec<u8>],
    b_comps: &[Vec<u8>],
    b_rcs: &[Vec<u8>],
    self_mode: bool,
    spec: &TrimSpec,
) -> Vec<Psl> {
    let Some((_, a_int)) = a_seqs.get(tube.a_contig as usize) else {
        return Vec::new();
    };
    let Some((_, b_int)) = b_seqs.get(tube.b_contig as usize) else {
        return Vec::new();
    };
    // Query contig in orientation space (RC for minus strand). The RC copies
    // are precomputed once and shared: per-tube copies of the ~5.5 MB main
    // chromosome used to serialize the allocator and inflate the extend peak.
    let q = if tube.strand == 0 {
        std::borrow::Cow::Borrowed(b_int)
    } else {
        std::borrow::Cow::Borrowed(&b_rcs[tube.b_contig as usize])
    };
    let q: &[u8] = &q;
    let rt = &a_revs[tube.a_contig as usize];
    let rq = if tube.strand == 0 {
        &b_revs[tube.b_contig as usize]
    } else {
        &b_comps[tube.b_contig as usize]
    };
    // FastGA self mode: only a same-contig forward tube is a self-alignment;
    // its wave must not cross diagonal 0 (the exact self-identity line).
    let selfie = self_mode && tube.a_contig == tube.b_contig && tube.strand == 0;
    let mut alow = tube.anti_low.max(*alast);
    let ahgh = tube.anti_high - BUCK_ANTI;
    let mut dgmin = tube.diag_min;
    let dgmax = tube.diag_max;
    if ahgh <= *alast {
        return Vec::new(); // already covered by an earlier tube
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
        if let Some(aln) = local_alignment(q, a_int, rt, rq, dgmin, dgmax, amid, selfie, spec) {
            let rlen = (aln.t_end - aln.t_start) as i64;
            if rlen >= TUBE_MIN_LEN && TUBE_MIN_RATE * rlen as f64 >= aln.diffs as f64 {
                let strand = if tube.strand == 0 { "+" } else { "-" };
                let mut q_start = aln.q_start as i32;
                let mut q_end = aln.q_end as i32;
                if tube.strand == 1 {
                    // `Psl::from_align` expects plus-strand coordinates; wave
                    // alignments are in RC space for minus-strand tubes.
                    reverse_range(
                        &mut q_start,
                        &mut q_end,
                        b.contigs()[tube.b_contig as usize].1 as i32,
                    );
                }
                if let Some(psl) = Psl::from_align(
                    &b.contigs()[tube.b_contig as usize].0,
                    b.contigs()[tube.b_contig as usize].1 as u32,
                    q_start,
                    q_end,
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
                // A genuinely smaller block contained in a much larger one
                // (e.g. a copy-pair hit inside the exact self-identity
                // diagonal of an explicit same-file pair) is a distinct
                // alignment, not a duplicate; the adjacent-tube duplicates
                // this filter targets have similar spans.
                && (k.t_end - k.t_start) <= 4 * (b.t_end - b.t_start)
                && (k.q_end - k.q_start) <= 4 * (b.q_end - b.q_start)
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

/// Convert one tube into a single-block PSL record (q = query, t = reference)
/// from its seed span; used when no extension sequences are available.
///
/// Reverse-strand blocks carry `q_start`/`q_end` in original query coordinates
/// (converted from orientation space), matching the pgr/UCSC PSL convention.
pub fn tube_to_psl(tube: &Tube, a: &PgiIndex, b: &impl PgiQuery) -> Psl {
    let (a_name, a_len) = {
        let (name, len) = &a.contigs[tube.a_contig as usize];
        (name, *len as u32)
    };
    let (b_name, b_len) = {
        let (name, len) = &b.contigs()[tube.b_contig as usize];
        (name, *len as u32)
    };
    let (q_start, q_end, strand) = if tube.strand == 0 {
        (tube.b_start, tube.b_end, "+")
    } else {
        let (s, e) = reverse_range_pair(tube.b_start, tube.b_end, b_len);
        (s, e, "-")
    };
    let mut psl = Psl::new();
    psl.q_name = b_name.clone();
    psl.q_size = b_len;
    psl.q_start = q_start as i32;
    psl.q_end = q_end as i32;
    psl.t_name = a_name.clone();
    psl.t_size = a_len;
    psl.t_start = tube.a_start as i32;
    psl.t_end = tube.a_end as i32;
    psl.strand = strand.to_string();
    psl.block_count = 1;
    psl.block_sizes.push(q_end - q_start);
    psl.q_starts.push(if tube.strand == 0 {
        q_start
    } else {
        b_len - q_end
    });
    psl.t_starts.push(tube.a_start);
    psl
}

/// Align two compatible indexes: merge seeds, chain tubes, and emit one
/// geometric PSL block per tube.
pub fn align_to_psl(a: &PgiIndex, b: &PgiIndex, params: &AlignParams) -> anyhow::Result<Vec<Psl>> {
    let hits = merge_seed_hits(a, b, params.freq, effective_min_shared(a.k, params))?;
    let tubes = chain_tubes(&hits, a.k as u32);
    Ok(tubes.iter().map(|t| tube_to_psl(t, a, b)).collect())
}

/// Same as [`align_to_psl`] with the reference index streamed from disk.
/// Drop hits that coincide exactly with the query's own position (same contig,
/// same coordinate, forward strand): a self-alignment must not report a
/// segment as its own copy (FastGA's self mode skips identical diagonals).
fn drop_self_hits(hits: &mut Vec<SeedHit>) {
    hits.retain(|h| !(h.a_contig == h.b_contig && h.a_pos == h.b_pos && h.strand == 0));
}

/// Same as [`align_to_psl`] with the reference index streamed from disk;
/// `is_self` drops the exact self-identity hits of a single-genome alignment.
pub fn align_to_psl_streaming<R: Read + Send, B: PgiQuery + Sync>(
    a: &mut PgiStream<R>,
    b: &B,
    params: &AlignParams,
    is_self: bool,
) -> anyhow::Result<Vec<Psl>> {
    let k = a.header().k;
    let mut hits = merge_seed_hits_from_stream(a, b, params.freq, effective_min_shared(k, params))?;
    if is_self {
        drop_self_hits(&mut hits);
    }
    let tubes = chain_tubes(&hits, k as u32);
    let a = PgiIndex {
        k,
        smer: a.header().smer,
        window: a.header().window,
        contigs: a.header().contigs.clone(),
        entries: Vec::new(),
        positions: Vec::new(),
    };
    Ok(tubes.iter().map(|t| tube_to_psl(t, &a, b)).collect())
}

/// Align two indexes, extending chains when sequences are provided and
/// falling back to plain blocks otherwise.
pub fn align_to_psl_ext(
    a: PgiIndex,
    mut b: PgiIndex,
    params: &AlignParams,
    a_seqs: &[(String, Vec<u8>)],
    b_seqs: &[(String, Vec<u8>)],
    is_self: bool,
) -> anyhow::Result<Vec<Psl>> {
    let min_shared = effective_min_shared(a.k, params);
    let mut hits = merge_seed_hits(&a, &b, params.freq, min_shared)?;
    if is_self {
        drop_self_hits(&mut hits);
    }
    log::info!(
        "merge: {} seed hits (min-shared={min_shared}, freq={})",
        hits.len(),
        params.freq
    );
    log::debug!("peak RSS after merge: {} MB", peak_rss_mb());
    // The extension phase only needs the contig tables; free the resident
    // query tables before the parallel pass (the mmap path never allocates
    // them in the first place).
    b.entries = Vec::new();
    b.positions = Vec::new();
    psls_from_hits(a, b, hits, a_seqs, b_seqs, is_self)
}

/// Align two indexes with the reference (`a`) streamed from disk and the
/// query (`b`) read through a [`PgiQuery`] view (resident or mmap'd): the
/// merge reads `a` in batches, so no index is materialized in full.
/// Same as [`align_to_psl_ext`] with the reference streamed and the query read
/// through a [`PgiQuery`] view; `is_self` drops exact self-identity hits.
pub fn align_to_psl_ext_streaming<R: Read + Send, B: PgiQuery + Sync>(
    mut a: PgiStream<R>,
    b: B,
    params: &AlignParams,
    a_seqs: &[(String, Vec<u8>)],
    b_seqs: &[(String, Vec<u8>)],
    is_self: bool,
) -> anyhow::Result<Vec<Psl>> {
    let header = a.header().clone();
    let min_shared = effective_min_shared(header.k, params);
    let t0 = std::time::Instant::now();
    let mut hits = merge_seed_hits_from_stream(&mut a, &b, params.freq, min_shared)?;
    if is_self {
        drop_self_hits(&mut hits);
    }
    log::debug!("merge: {} ms", t0.elapsed().as_millis());
    log::info!(
        "merge: {} seed hits (min-shared={min_shared}, freq={})",
        hits.len(),
        params.freq
    );
    log::debug!("peak RSS after merge: {} MB", peak_rss_mb());
    let a = PgiIndex {
        k: header.k,
        smer: header.smer,
        window: header.window,
        contigs: header.contigs,
        entries: Vec::new(),
        positions: Vec::new(),
    };
    psls_from_hits(a, b, hits, a_seqs, b_seqs, is_self)
}

/// Chain and extend a seed-hit list into PSL blocks.
///
/// `a` only contributes its contig table after the merge (its entries and
/// positions are dropped before the extension phase); callers may pass an
/// index with empty tables when the reference was streamed.
fn psls_from_hits(
    mut a: PgiIndex,
    b: impl PgiQuery + Sync,
    hits: Vec<SeedHit>,
    a_seqs: &[(String, Vec<u8>)],
    b_seqs: &[(String, Vec<u8>)],
    is_self: bool,
) -> anyhow::Result<Vec<Psl>> {
    // Tubes are independent alignment tasks; the `alast` overlap skip of
    // FastGA's sequential scan is replaced by a parallel pass plus a
    // containment dedup on the emitted blocks.
    let t0 = std::time::Instant::now();
    let tubes = chain_tubes(&hits, a.k as u32);
    log::debug!("chain_tubes: {} ms", t0.elapsed().as_millis());
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
                        piece.anti_high = (t.anti_low + (i + 1) * TUBE_ANTI_CAP).min(t.anti_high);
                        piece
                    })
                    .collect()
            }
        })
        .collect();
    log::debug!("peak RSS after chain_tubes: {} MB", peak_rss_mb());
    drop(hits); // the tubes and sequences are enough for extension
                // The k-mer entries/positions are only needed for the merge; drop the
                // reference index's tables before the parallel extension phase.
    drop(std::mem::take(&mut a.entries));
    drop(std::mem::take(&mut a.positions));
    log::debug!("peak RSS after dropping indexes: {} MB", peak_rss_mb());
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
    let b_rcs: Vec<Vec<u8>> = b_seqs.iter().map(|(_, s)| rev_comp(s).collect()).collect();
    let t_ext = std::time::Instant::now();
    let spec = TrimSpec::for_seqs(a_seqs);
    let records: Vec<Vec<Psl>> = tubes
        .par_iter()
        .map(|t| {
            let mut alast = 0i64;
            extend_tube(
                t, &mut alast, &a, &b, a_seqs, b_seqs, &a_revs, &b_revs, &b_comps, &b_rcs, is_self,
                &spec,
            )
        })
        .collect();
    log::debug!("extend: {} ms", t_ext.elapsed().as_millis());
    let mut out: Vec<Psl> = records.into_iter().flatten().collect();
    dedupe_contained(&mut out);
    log::debug!("peak RSS after extend: {} MB", peak_rss_mb());
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::pgi::build::build_from_seqs;
    use crate::libs::pgi::PgiEntry;

    fn pseudo_random_seq(len: usize, seed: u64) -> Vec<u8> {
        let bases = *b"ACGT";
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
        build_from_seqs(
            vec![(String::from("c"), seq.to_vec())],
            10,
            4,
            2,
            false,
            false,
        )
        .unwrap()
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

        // Random sequence: k-mers are mostly unique, so freq=2 keeps hits
        // (FastGA's filter skips ranges with `>= freq` occurrences).
        let rnd = pseudo_random_seq(300, 7);
        let rnd_idx = build(&rnd);
        let hits = merge_seed_hits(&rnd_idx, &rnd_idx, 2, 10).unwrap();
        assert!(!hits.is_empty(), "unique k-mers pass freq=2");
    }

    /// Binary-search-only reference of the sequential merge (the pre-方案 A
    /// `emit_entry_hits`): every window is located by `entry_range`, so the
    /// sequential-advance path must produce byte-identical hits.
    #[allow(clippy::too_many_arguments)]
    fn emit_entry_hits_ref<B: PgiQuery>(
        ea_kmer: u128,
        ea_freq: u32,
        a_positions: &[u64],
        b: &B,
        freq: u32,
        min_shared: usize,
        k: usize,
        prev_kmer: Option<u128>,
    ) -> anyhow::Result<Vec<SeedHit>> {
        let mut hits = Vec::new();
        if ea_freq >= freq || ea_kmer > rc_key(ea_kmer, k) {
            return Ok(hits);
        }
        let k_bits = 2 * k;
        let mask = if k_bits >= 128 {
            u128::MAX
        } else {
            (1u128 << k_bits) - 1
        };
        let window = |len: usize| {
            let r = 1u128 << (k_bits - 2 * len);
            let lo = ea_kmer & !(r - 1) & mask;
            b.entry_range(lo, lo.saturating_add(r))
        };
        let start = prev_kmer
            .map(|pk| shared_prefix(pk, ea_kmer, k).max(min_shared as u32))
            .unwrap_or(min_shared as u32);
        let (mut j0, mut j) = window(start as usize);
        if j0 == j && start as usize > min_shared {
            (j0, j) = window(min_shared);
        }
        if j == j0 {
            return Ok(hits);
        }
        let max_shared_over = |mut i: usize, j: usize| -> u32 {
            let mut m = 0u32;
            while i < j {
                if b.entry_freq(i) >= freq {
                    i = b.entry_next(i);
                    continue;
                }
                m = m.max(shared_prefix(ea_kmer, b.entry_kmer(i), k));
                i = b.entry_next(i);
            }
            m
        };
        let mut m = max_shared_over(j0, j);
        if m < min_shared as u32 && start as usize > min_shared {
            (j0, j) = window(min_shared);
            m = max_shared_over(j0, j);
        }
        if m < min_shared as u32 {
            return Ok(hits);
        }
        let (m0, m1) = window(m as usize);
        let mut occ: u32 = 0;
        let mut i = m0;
        while i < m1 {
            let f = b.entry_freq(i);
            if f < freq {
                occ += f;
            }
            i = b.entry_next(i);
        }
        if occ == 0 || occ >= freq {
            return Ok(hits);
        }
        let mut i = m0;
        while i < m1 {
            if b.entry_freq(i) < freq {
                for &arec in a_positions {
                    let (ac, apos, astrand) = unpack_position(arec);
                    for brec in b.entry_positions(i) {
                        let (bc, bpos, bstrand) = unpack_position(brec);
                        let fwd = astrand == bstrand;
                        let b_len = b.contigs()[bc as usize].1;
                        let oriented = if fwd {
                            bpos as u64
                        } else {
                            b_len - k as u64 - bpos as u64
                        };
                        hits.push(SeedHit {
                            a_contig: ac as u16,
                            a_pos: apos,
                            b_contig: bc as u16,
                            b_pos: oriented as u32,
                            shared: m as u16,
                            strand: u8::from(!fwd),
                        });
                    }
                }
            }
            i = b.entry_next(i);
        }
        Ok(hits)
    }

    fn merge_ref<B: PgiQuery>(a: &PgiIndex, b: &B, freq: u32, min_shared: usize) -> Vec<SeedHit> {
        let k = a.k;
        let mut hits = Vec::new();
        let mut prev = None;
        for ea in &a.entries {
            let ap = &a.positions[ea.pos_start as usize..(ea.pos_start + ea.freq) as usize];
            hits.extend(
                emit_entry_hits_ref(ea.kmer, ea.freq, ap, b, freq, min_shared, k, prev).unwrap(),
            );
            prev = Some(ea.kmer);
        }
        hits
    }

    #[test]
    fn sequential_merge_matches_binary_reference() {
        // The 归并式 (sequential-advance) merge must produce byte-identical
        // seed hits to the pre-optimization binary-search reference, on
        // indexes with varied lcp propagation (long A-runs, RC query) and
        // both resident and mmap query views.
        let a = build_from_seqs(
            vec![
                (String::from("c1"), pseudo_random_seq(3000, 21)),
                (String::from("c2"), pseudo_random_seq(2000, 22)),
            ],
            40,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let b = build_from_seqs(
            vec![
                (String::from("c1"), pseudo_random_seq(3000, 23)),
                (String::from("c2"), pseudo_random_seq(2000, 24)),
            ],
            40,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let key = |h: &SeedHit| (h.a_contig, h.a_pos, h.b_contig, h.b_pos, h.shared, h.strand);

        let mut ref_resident: Vec<_> = merge_ref(&a, &b, 10, 12).iter().map(key).collect();
        let mut seq_resident: Vec<_> = merge_seed_hits(&a, &b, 10, 12)
            .unwrap()
            .iter()
            .map(key)
            .collect();
        ref_resident.sort_unstable();
        seq_resident.sort_unstable();
        assert_eq!(ref_resident, seq_resident, "resident sequential mismatch");

        // Mmap query view: the generic merge is the streaming one (it is the
        // only path that accepts a `B: PgiQuery` query side).
        let mut a_bytes = Vec::new();
        a.write(&mut a_bytes).unwrap();
        let mut b_bytes = Vec::new();
        b.write(&mut b_bytes).unwrap();
        let dir = tempfile::TempDir::new().unwrap();
        let q_path = dir.path().join("q.pgi");
        std::fs::write(&q_path, &b_bytes).unwrap();
        let mapped = crate::libs::pgi::PgiMmap::open(&q_path).unwrap();
        let mut ref_mmap: Vec<_> = merge_ref(&a, &mapped, 10, 12).iter().map(key).collect();
        let seq_mmap = {
            let mut s = PgiStream::open(std::io::Cursor::new(&a_bytes)).unwrap();
            merge_seed_hits_from_stream(&mut s, &mapped, 10, 12).unwrap()
        };
        let mut seq_mmap: Vec<_> = seq_mmap.iter().map(key).collect();
        ref_mmap.sort_unstable();
        seq_mmap.sort_unstable();
        assert_eq!(ref_mmap, seq_mmap, "mmap sequential mismatch");
        assert!(!ref_mmap.is_empty(), "expected shared seeds");
    }

    #[test]
    fn streamed_merge_matches_resident_query() {
        // The memory-mapped query view must decode the same seed hits as the
        // fully loaded index (entry ranges, frequencies, positions).
        let a = build_from_seqs(
            vec![
                (String::from("c1"), pseudo_random_seq(2000, 11)),
                (String::from("c2"), pseudo_random_seq(1500, 5)),
            ],
            10,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let b = build_from_seqs(
            vec![
                (String::from("c1"), pseudo_random_seq(2000, 12)),
                (String::from("c2"), pseudo_random_seq(1500, 6)),
            ],
            10,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let mut ref_bytes = Vec::new();
        a.write(&mut ref_bytes).unwrap();
        let mut q_bytes = Vec::new();
        b.write(&mut q_bytes).unwrap();
        let dir = tempfile::TempDir::new().unwrap();
        let q_path = dir.path().join("query.pgi");
        std::fs::write(&q_path, &q_bytes).unwrap();
        let mapped = crate::libs::pgi::PgiMmap::open(&q_path).unwrap();

        let resident = {
            let mut s = PgiStream::open(std::io::Cursor::new(&ref_bytes)).unwrap();
            merge_seed_hits_from_stream(&mut s, &b, 10, 10).unwrap()
        };
        let mmap = {
            let mut s = PgiStream::open(std::io::Cursor::new(&ref_bytes)).unwrap();
            merge_seed_hits_from_stream(&mut s, &mapped, 10, 10).unwrap()
        };
        let key = |h: &SeedHit| (h.a_contig, h.a_pos, h.b_contig, h.b_pos, h.shared, h.strand);
        let mut rk: Vec<_> = resident.iter().map(key).collect();
        let mut mk: Vec<_> = mmap.iter().map(key).collect();
        rk.sort_unstable();
        mk.sort_unstable();
        assert_eq!(rk, mk, "mmap query must produce identical seed hits");
        assert!(!rk.is_empty(), "expected shared seeds");
    }

    #[test]
    fn mmap_merge_rejects_out_of_range_contig() {
        // Regression: a crafted query index whose occurrence record carries a
        // contig id beyond the contig table used to panic in `emit_entry_hits`
        // (`b.contigs()[bc]`). The lazy mmap decode must surface a friendly
        // error when the corrupted record is actually merged.
        let seq = pseudo_random_seq(2000, 77);
        let idx = build_from_seqs(
            vec![(String::from("c"), seq.clone())],
            10,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let mut bytes = Vec::new();
        idx.write(&mut bytes).unwrap();
        // Corrupt the first *canonical* record: non-canonical keys are never
        // merged (they are filtered in `emit_entry_hits`), so the corrupted
        // record must be one the merge actually decodes.
        let (header, _n, layout, records_off) =
            crate::libs::pgi::parse_header_bytes(&bytes).unwrap();
        let rec_size = layout.size();
        for i in 0..(bytes.len() - records_off) / rec_size {
            let off = records_off + i * rec_size;
            let kmer =
                crate::libs::pgi::unpack_kmer(&bytes[off..off + layout.kmer_bytes], header.k);
            if kmer <= rc_key(kmer, header.k) {
                bytes[off + layout.kmer_bytes + layout.pos_bytes] = 0x7f;
                break;
            }
        }
        let dir = tempfile::TempDir::new().unwrap();
        let q_path = dir.path().join("bad.pgi");
        std::fs::write(&q_path, &bytes).unwrap();
        let mapped = crate::libs::pgi::PgiMmap::open(&q_path).unwrap();

        let mut ref_bytes = Vec::new();
        idx.write(&mut ref_bytes).unwrap();
        let mut s = PgiStream::open(std::io::Cursor::new(&ref_bytes)).unwrap();
        let err = merge_seed_hits_from_stream(&mut s, &mapped, 10, 10).unwrap_err();
        assert!(err.to_string().contains("out of range"), "got: {err}");
    }

    #[test]
    fn drop_self_hits_filters_exact_identity() {
        let mk = |a_contig: u16, a_pos: u32, b_contig: u16, b_pos: u32, strand: u8| SeedHit {
            a_contig,
            a_pos,
            b_contig,
            b_pos,
            shared: 40,
            strand,
        };
        let mut hits = vec![
            mk(0, 100, 0, 100, 0), // exact self-identity -> dropped
            mk(0, 100, 0, 500, 0), // intra-genome repeat -> kept
            mk(0, 100, 0, 100, 1), // reverse copy at the same spot -> kept
            mk(0, 100, 1, 100, 0), // cross-contig homology -> kept
        ];
        drop_self_hits(&mut hits);
        assert_eq!(hits.len(), 3);
        assert!(hits
            .iter()
            .all(|h| { !(h.a_contig == h.b_contig && h.a_pos == h.b_pos && h.strand == 0) }));
        assert!(hits.iter().any(|h| h.b_pos == 500 && h.strand == 0));
        assert!(hits.iter().any(|h| h.strand == 1));
        assert!(hits.iter().any(|h| h.a_contig != h.b_contig));
    }

    #[test]
    fn merge_k64_high_key_no_prefix_overflow() {
        // Regression: at k=64 the prefix boundary `lo + r` reaches `2^128`,
        // which `u128` cannot represent. An `a` entry whose leading 12 bases
        // are T (bits 104..128 set) drove `window(12)` to `hi = lo + r =
        // 2^128`, overflowing `u128` (a debug panic / release wrap). The
        // merge must not panic.
        let k64 = PgiIndex {
            k: 64,
            smer: 8,
            window: 5,
            contigs: vec![(String::from("c"), 1000)],
            entries: Vec::new(),
            positions: Vec::new(),
        };
        // K = 12×T then 52×A: canonical (K <= rc_key(K), since its RC is
        // 52×T then 12×A), and window(12) gives lo = 0xFFFFFF << 104,
        // hi = lo + 2^104 = 2^128 (the overflow).
        let k = 0xFFFF_FF00_0000_0000_0000_0000_0000_0000u128;
        let a = PgiIndex {
            entries: vec![PgiEntry {
                kmer: k,
                pos_start: 0,
                freq: 1,
            }],
            positions: vec![crate::libs::pgi::pack_position(0, 0, 0)],
            ..k64.clone()
        };
        let b = k64;
        let hits = merge_seed_hits(&a, &b, 10, 12).unwrap();
        assert!(hits.is_empty(), "no b entries, so no hits: {hits:?}");
    }

    #[test]
    fn merge_keeps_only_maximal_shared_prefix() {
        // k=32, all-A key against a 31-base and a 25-base match: only the
        // longest match is emitted (FastGA adaptamer semantics).
        let mk = |entries: Vec<PgiEntry>, positions: Vec<u64>| PgiIndex {
            k: 32,
            smer: 8,
            window: 5,
            contigs: vec![(String::from("c"), 1000)],
            entries,
            positions,
        };
        let a = mk(
            vec![PgiEntry {
                kmer: 0,
                pos_start: 0,
                freq: 1,
            }],
            vec![crate::libs::pgi::pack_position(0, 10, 0)],
        );
        // First 31 bases shared (A..AT), then first 25 bases shared (A..AT);
        // entries must stay ascending by k-mer.
        let b = mk(
            vec![
                PgiEntry {
                    kmer: 3,
                    pos_start: 0,
                    freq: 1,
                },
                PgiEntry {
                    kmer: 3 << 12,
                    pos_start: 1,
                    freq: 1,
                },
            ],
            vec![
                crate::libs::pgi::pack_position(0, 50, 0),
                crate::libs::pgi::pack_position(0, 60, 0),
            ],
        );
        let hits = merge_seed_hits(&a, &b, 10, 20).unwrap();
        assert_eq!(hits.len(), 1, "only the maximal match is emitted");
        assert_eq!(hits[0].shared, 31);
        assert_eq!((hits[0].a_pos, hits[0].b_pos), (10, 50));
    }

    #[test]
    fn merge_filters_extended_range_by_occurrences() {
        // The extended range at the maximal shared prefix holds `freq` or
        // more total occurrences of *under-frequency* entries: the whole
        // entry is skipped (FastGA's extended-range filter sums the index
        // counts over the range and drops it at `>= FREQ`).
        let mk_a = PgiIndex {
            k: 32,
            smer: 8,
            window: 5,
            contigs: vec![(String::from("c"), 1000)],
            entries: vec![PgiEntry {
                kmer: 0,
                pos_start: 0,
                freq: 1,
            }],
            positions: vec![crate::libs::pgi::pack_position(0, 10, 0)],
        };
        // Three distinct 25-base-matching k-mers, 4 occurrences each (each
        // entry is under the cutoff, but 12 >= freq=10 in total).
        let b = PgiIndex {
            k: 32,
            smer: 8,
            window: 5,
            contigs: vec![(String::from("c"), 1000)],
            entries: vec![
                PgiEntry {
                    kmer: 1 << 12,
                    pos_start: 0,
                    freq: 4,
                },
                PgiEntry {
                    kmer: 1 << 13,
                    pos_start: 4,
                    freq: 4,
                },
                PgiEntry {
                    kmer: 3 << 12,
                    pos_start: 8,
                    freq: 4,
                },
            ],
            positions: (0..10)
                .map(|i| crate::libs::pgi::pack_position(0, 50 + i, 0))
                .collect(),
        };
        let hits = merge_seed_hits(&mk_a, &b, 10, 20).unwrap();
        assert!(hits.is_empty(), "extended range is too frequent: {hits:?}");
    }

    #[test]
    fn freq_boundary_drops_exact_freq_on_reference_side() {
        // Regression: the reference-side frequency check used `> freq`
        // (keeping `== freq`) while the query side (and FastGA's GIX build)
        // drop `>= freq`. A k-mer occurring exactly `freq` times on the
        // reference but rarely on the query must not seed either.
        let mk = |freq: u32, positions: Vec<u64>| PgiIndex {
            k: 32,
            smer: 8,
            window: 5,
            contigs: vec![(String::from("c"), 1000)],
            entries: vec![PgiEntry {
                kmer: 0,
                pos_start: 0,
                freq,
            }],
            positions,
        };
        let a = mk(
            2,
            vec![
                crate::libs::pgi::pack_position(0, 10, 0),
                crate::libs::pgi::pack_position(0, 20, 0),
            ],
        );
        let b = mk(1, vec![crate::libs::pgi::pack_position(0, 50, 0)]);
        let hits = merge_seed_hits(&a, &b, 2, 32).unwrap();
        assert!(hits.is_empty(), "exact-freq k-mers must not seed: {hits:?}");
    }

    #[test]
    fn exact_freq_query_entries_are_absent_not_range_killers() {
        // FastGA's GIX index excludes k-mers with count >= FREQ at build
        // time, so an `== freq` entry must behave as absent: it neither
        // raises the maximal shared prefix nor drops the extended range. A
        // rare (freq 1) 25-base match next to an `== freq` 31-base entry must
        // still seed at the rare entry's prefix length.
        let mk_a = PgiIndex {
            k: 32,
            smer: 8,
            window: 5,
            contigs: vec![(String::from("c"), 1000)],
            entries: vec![PgiEntry {
                kmer: 0,
                pos_start: 0,
                freq: 1,
            }],
            positions: vec![crate::libs::pgi::pack_position(0, 10, 0)],
        };
        let b = PgiIndex {
            k: 32,
            smer: 8,
            window: 5,
            contigs: vec![(String::from("c"), 1000)],
            entries: vec![
                PgiEntry {
                    kmer: 3,
                    pos_start: 0,
                    freq: 2,
                },
                PgiEntry {
                    kmer: 3 << 12,
                    pos_start: 2,
                    freq: 1,
                },
            ],
            positions: vec![
                crate::libs::pgi::pack_position(0, 50, 0),
                crate::libs::pgi::pack_position(0, 51, 0),
                crate::libs::pgi::pack_position(0, 60, 0),
            ],
        };
        let hits = merge_seed_hits(&mk_a, &b, 2, 20).unwrap();
        assert_eq!(hits.len(), 1, "rare entry must seed: {hits:?}");
        assert_eq!(hits[0].shared, 25);
        assert_eq!(hits[0].b_pos, 60);
    }

    #[test]
    fn lcp_narrowed_window_all_high_freq_falls_back_to_floor() {
        // Regression: the lcp-narrowed window for an `a` entry can be
        // non-empty yet hold only `>= freq` (high-frequency) k-mers, which
        // the merge treats as absent. FastGA's GIX index omits those k-mers,
        // so its narrowed window would be empty and it falls back to the
        // floor window; pgr keeps them in the index, so the narrowed window
        // looked non-empty and the merge returned no hits, missing an
        // under-frequency seed in the floor window below the lcp.
        let mk = |entries: Vec<PgiEntry>, positions: Vec<u64>| PgiIndex {
            k: 8,
            smer: 4,
            window: 2,
            contigs: vec![(String::from("c"), 1000)],
            entries,
            positions,
        };
        // a1 = AAAAAAAA (0x0000), a0 = AAAAAATT (0x000F); lcp(a1, a0) = 6.
        let a = mk(
            vec![
                PgiEntry {
                    kmer: 0x0000,
                    pos_start: 0,
                    freq: 1,
                },
                PgiEntry {
                    kmer: 0x000F,
                    pos_start: 1,
                    freq: 1,
                },
            ],
            vec![
                crate::libs::pgi::pack_position(0, 0, 0),
                crate::libs::pgi::pack_position(0, 10, 0),
            ],
        );
        // b' = AAAAAACA (0x0004) shares 6 bases -> in the narrowed window,
        // but is high-frequency (freq 2 == cutoff). b = AAAAACAA (0x0040)
        // shares 4 bases (< lcp 6, >= min-shared 4) -> in the floor window
        // only, under-frequency.
        let b = mk(
            vec![
                PgiEntry {
                    kmer: 0x0004,
                    pos_start: 0,
                    freq: 2,
                },
                PgiEntry {
                    kmer: 0x0040,
                    pos_start: 1,
                    freq: 1,
                },
            ],
            vec![
                crate::libs::pgi::pack_position(0, 60, 0),
                crate::libs::pgi::pack_position(0, 50, 0),
            ],
        );
        let hits = merge_seed_hits(&a, &b, 2, 4).unwrap();
        // a0 (AAAAAATT) must still seed the under-frequency 4-base match even
        // though its lcp-narrowed window held only a high-frequency entry.
        assert!(
            hits.iter().any(|h| h.a_pos == 10 && h.shared == 4),
            "a0 must fall back to the floor window and seed its 4-base match: {hits:?}"
        );
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
    fn dedupe_keeps_small_block_inside_large_one() {
        // Regression: a copy-pair hit (small block) contained on both axes in
        // the exact self-identity diagonal (large block) of an explicit
        // same-file pair used to be dropped as "duplicate", losing the real
        // homology. A genuinely distinct alignment has a very different span.
        let mk = |t_start: i32, t_end: i32, q_start: i32, q_end: i32| Psl {
            q_name: String::from("q"),
            q_size: 100_000,
            q_start,
            q_end,
            t_name: String::from("t"),
            t_size: 100_000,
            t_start,
            t_end,
            strand: String::from("+"),
            block_count: 1,
            block_sizes: vec![(t_end - t_start) as u32],
            q_starts: vec![q_start as u32],
            t_starts: vec![t_start as u32],
            ..Default::default()
        };
        let large = mk(0, 30_600, 0, 30_600);
        let small = mk(10_004, 11_201, 13_204, 14_401);
        let mut blocks = vec![large, small];
        dedupe_contained(&mut blocks);
        assert_eq!(blocks.len(), 2, "the pair hit must survive: {blocks:?}");
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
    fn tube_sort_key_supports_large_anti_coordinates() {
        // anti = a_pos + b_pos > 2^24 must not corrupt the u128 sort key
        // (the 24-bit anti field overflowed into the bucket field on
        // genomes larger than ~8 Mb, interleaving diagonal buckets).
        let mk = |a: u32, b: u32| SeedHit {
            a_contig: 0,
            a_pos: a,
            b_contig: 0,
            b_pos: b,
            strand: 0,
            shared: 40,
        };
        // Three collinear seeds on diagonal 0 with anti ~25,000,000 (> 2^24).
        let mut hits = vec![
            mk(12_500_000, 12_500_000),
            mk(12_500_020, 12_500_020),
            mk(12_500_040, 12_500_040),
        ];
        // One seed on diagonal 64 whose anti (~8.22M, < 2^24) sorts between
        // the big seeds' anti_low values; with the overflow it interleaves
        // and fragments the diagonal-0 run.
        hits.push(mk(4_111_444, 4_111_380));
        let tubes = chain_tubes(&hits, 40);
        assert_eq!(tubes.len(), 1, "diag-0 run must stay one tube: {tubes:?}");
        assert_eq!((tubes[0].diag_min, tubes[0].diag_max), (0, 0));
        assert_eq!(
            (tubes[0].anti_low, tubes[0].anti_high),
            (25_000_000, 25_000_120)
        );
    }

    #[test]
    fn tube_sort_key_supports_deeply_negative_diagonals() {
        // Regression: the diagonal bucket offset (1,000,000) only covered
        // diagonals down to ~-64 Mb; a seed pair 80 Mb apart (diag -80M)
        // made the bucket negative, which wrapped through `as u64` and
        // spilled into the strand/contig fields of the packed sort key.
        // Long-range repeats on large contigs must still form tubes.
        let mk = |a: u32, b: u32| SeedHit {
            a_contig: 0,
            a_pos: a,
            b_contig: 0,
            b_pos: b,
            strand: 0,
            shared: 40,
        };
        // Dense run on diagonal -80,000,000 (a=0..200, b=80,000,000..).
        let mut hits: Vec<SeedHit> = (0..200u32).map(|i| mk(i, 80_000_000 + i)).collect();
        // A second dense run on diagonal 0.
        hits.extend((0..200u32).map(|i| mk(1_000_000 + i, 1_000_000 + i)));
        let tubes = chain_tubes(&hits, 40);
        assert_eq!(
            tubes.len(),
            2,
            "each diagonal run must form one tube: {tubes:?}"
        );
        let d0 = tubes.iter().find(|t| t.diag_min == -80_000_000).unwrap();
        assert_eq!(d0.diag_max, -80_000_000);
        assert_eq!((d0.anti_low, d0.anti_high), (80_000_000, 80_000_438));
        let main = tubes.iter().find(|t| t.diag_min == 0).unwrap();
        assert_eq!(main.diag_max, 0);
        assert_eq!((main.anti_low, main.anti_high), (2_000_000, 2_000_438));
    }

    #[test]
    fn tube_sort_key_does_not_mix_contigs_at_negative_diagonals() {
        // Regression: with the 1,000,000 bucket offset, seeds on diagonals
        // below ~-64 Mb wrapped their bucket through `as u64`, spilling
        // 0xFFFF into the strand/contig key fields. Every wrapped seed then
        // shared the same corrupted (contig, strand) fields, so seeds from
        // different contigs interleaved by anti and their tubes fragmented.
        let mk = |a_contig: u16, a: u32, b: u32| SeedHit {
            a_contig,
            a_pos: a,
            b_contig: 0,
            b_pos: b,
            strand: 0,
            shared: 40,
        };
        // Two contigs, each with a dense run on diagonal -80,000,000.
        let hits: Vec<SeedHit> = (0..200u32)
            .flat_map(|i| vec![mk(0, i, 80_000_000 + i), mk(1, i, 80_000_000 + i)])
            .collect();
        let tubes = chain_tubes(&hits, 40);
        assert_eq!(
            tubes.len(),
            2,
            "one tube per contig expected, got {tubes:?}"
        );
        let t0 = tubes.iter().find(|t| t.a_contig == 0).unwrap();
        let t1 = tubes.iter().find(|t| t.a_contig == 1).unwrap();
        assert_eq!((t0.anti_low, t0.anti_high), (80_000_000, 80_000_438));
        assert_eq!((t1.anti_low, t1.anti_high), (80_000_000, 80_000_438));
    }

    #[test]
    fn psl_block_coordinates() {
        let seq: Vec<u8> = (0..100u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let (ia, ib) = (build(&seq), build(&seq));
        let fwd = Tube {
            a_contig: 0,
            b_contig: 0,
            strand: 0,
            anti_low: 30,
            anti_high: 110,
            diag_min: -10,
            diag_max: -10,
            a_start: 10,
            a_end: 50,
            b_start: 20,
            b_end: 60,
        };
        let p = tube_to_psl(&fwd, &ia, &ib);
        assert_eq!(p.strand, "+");
        assert_eq!((p.q_start, p.q_end), (20, 60));
        assert_eq!((p.t_start, p.t_end), (10, 50));

        let rev = Tube {
            strand: 1,
            b_start: 20,
            b_end: 60,
            ..fwd
        };
        let p = tube_to_psl(&rev, &ia, &ib);
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
        let psls = align_to_psl_ext(ia, ib, &params, &a_seqs, &b_seqs, false).unwrap();
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
        let b_seqs = vec![(String::from("c"), rc.clone())];
        // Exact-only seeds: partial matches could pair unrelated same-strand
        // positions in this tiny random sequence and break the pure-RC check.
        let params = AlignParams {
            min_shared: Some(10),
            ..AlignParams::default()
        };
        let psls = align_to_psl_ext(ia, ir, &params, &a_seqs, &b_seqs, false).unwrap();
        assert!(!psls.is_empty());
        assert!(
            psls.iter().all(|p| p.strand == "-" && p.match_count > 0),
            "RC query must extend to minus-strand PSL"
        );
        for p in &psls {
            assert!(p.q_start >= 0 && p.q_end <= p.q_size as i32);
            // UCSC/`psl chain` minus-strand convention: qStart/qEnd are on the
            // plus strand and the internal qStarts are in RC frame (the first
            // internal start mirrors qEnd; `calc_block_score` reads them as
            // RC coordinates). Regression: a swapped frame scored negative in
            // `psl chain` and silently dropped every minus-strand block.
            assert_eq!(
                p.q_starts[0] as i32,
                p.q_size as i32 - p.q_end,
                "minus-strand qStarts must be in RC frame"
            );
            assert!(
                p.q_starts.windows(2).all(|w| w[1] > w[0]),
                "minus-strand qStarts must ascend"
            );
            for (i, ((&sz, &qs), &ts)) in p
                .block_sizes
                .iter()
                .zip(&p.q_starts)
                .zip(&p.t_starts)
                .enumerate()
            {
                let q_plus = (p.q_size as i32 - (qs + sz) as i32) as usize;
                let q_seg = &rc[q_plus..q_plus + sz as usize];
                let t_seg = &seq[ts as usize..(ts + sz) as usize];
                let rev_comp_seg: Vec<u8> = rev_comp(q_seg).collect();
                let mm = rev_comp_seg
                    .iter()
                    .zip(t_seg)
                    .filter(|(a, b)| a != b)
                    .count();
                assert!(
                    mm <= (sz / 10) as usize,
                    "segment {i} must match the RC query ({mm} mismatches)"
                );
            }
        }
    }

    #[test]
    fn tube_self_keeps_diagonal_away_from_zero() {
        // Tandem repeat of two 400 bp copies: self-alignment must find the
        // copy pairs and must not emit a same-contig forward block whose
        // target/query intervals overlap (FastGA self mode forbids crossing
        // diagonal 0, the exact self-identity line).
        let copy = pseudo_random_seq(400, 7);
        let seq: Vec<u8> = copy.iter().chain(copy.iter()).copied().collect();
        let (ia, ib) = (build(&seq), build(&seq));
        let a_seqs = vec![(String::from("c"), seq.clone())];
        let b_seqs = vec![(String::from("c"), seq)];
        let params = AlignParams::default();
        let psls = align_to_psl_ext(ia, ib, &params, &a_seqs, &b_seqs, true).unwrap();
        assert!(!psls.is_empty(), "copy pair must be found");
        for p in &psls {
            if p.q_name == p.t_name && p.strand == "+" {
                // Each block is a gapless run with a constant diagonal; self
                // mode must never emit a block on diagonal 0 (same-position
                // self-identity).
                for (i, ((&sz, &qs), &ts)) in p
                    .block_sizes
                    .iter()
                    .zip(&p.q_starts)
                    .zip(&p.t_starts)
                    .enumerate()
                {
                    assert!(
                        ts as i64 != qs as i64,
                        "self block {i} sits on diagonal 0: q {}..{} t {}..{}",
                        qs,
                        qs + sz,
                        ts,
                        ts + sz
                    );
                }
            }
        }
    }

    #[test]
    fn multi_contig_align_with_rc_query() {
        // Three contigs; the query's second contig is the reverse complement
        // and the third carries scattered point mutations. Regression:
        // contig grouping, the streamed merge, and the minus-strand PSL
        // convention must all hold on multi-contig inputs.
        let c1 = pseudo_random_seq(2000, 1);
        let c2 = pseudo_random_seq(1500, 2);
        let c3 = pseudo_random_seq(1000, 3);
        let qc2: Vec<u8> = rev_comp(&c2).collect();
        let mut qc3 = c3.clone();
        for i in (0..qc3.len()).step_by(50) {
            qc3[i] = match qc3[i] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        let ia = build_from_seqs(
            vec![
                (String::from("c1"), c1.clone()),
                (String::from("c2"), c2.clone()),
                (String::from("c3"), c3.clone()),
            ],
            10,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let ib = build_from_seqs(
            vec![
                (String::from("c1"), c1.clone()),
                (String::from("c2"), qc2.clone()),
                (String::from("c3"), qc3.clone()),
            ],
            10,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let a_seqs = vec![
            (String::from("c1"), c1.clone()),
            (String::from("c2"), c2),
            (String::from("c3"), c3),
        ];
        let b_seqs = vec![
            (String::from("c1"), c1),
            (String::from("c2"), qc2),
            (String::from("c3"), qc3),
        ];
        let params = AlignParams {
            min_shared: Some(10),
            ..AlignParams::default()
        };
        let psls = align_to_psl_ext(ia, ib, &params, &a_seqs, &b_seqs, false).unwrap();
        assert!(psls.len() >= 3, "expected one block per contig pair");
        assert!(
            psls.iter()
                .any(|p| p.q_name == "c1" && p.strand == "+" && p.match_count >= 1900),
            "c1 must align forward near full length"
        );
        assert!(
            psls.iter()
                .any(|p| p.q_name == "c2" && p.strand == "-" && p.match_count >= 1400),
            "the RC query contig must align on the minus strand"
        );
        assert!(
            psls.iter()
                .any(|p| p.q_name == "c3" && p.strand == "+" && p.mismatch_count > 0),
            "the mutated contig must align with mismatches"
        );
        assert!(psls.iter().all(|p| {
            p.q_start >= 0
                && p.q_end <= p.q_size as i32
                && p.t_start >= 0
                && p.t_end <= p.t_size as i32
        }));
    }
}
