//! In-memory `merge_seed_hits` timing vs the streaming command path, plus a
//! skip-scan variant to isolate the cost of the max-shared-prefix window
//! scan (the part the LCP optimizations touch). Includes a split-profile
//! variant that times the four sub-components of `emit_entry_hits`:
//! (a) `entry_range` binary search, (b) the max-m window scan, (c)
//! `shared_prefix` calls, (d) position emission (design: pgi-query-layer.md
//! 阶段 0).
//!
//! Usage: cargo run --release --example merge_mem_bench -- ref.pgi query.pgi

use pgr::libs::pgi::align::SeedHit;
use pgr::libs::pgi::{PgiIndex, PgiQuery};
use rayon::prelude::*;
use std::time::Instant;

const CID_MASK: u64 = (1 << 20) - 1;
const STRAND_OFF: u32 = 52;

fn shared_prefix(a: u128, b: u128, k: usize) -> u32 {
    let x = a ^ b;
    if x == 0 {
        return k as u32;
    }
    let bitlen = 128 - x.leading_zeros();
    ((2 * k).saturating_sub(bitlen as usize) / 2) as u32
}

/// Accumulated split timers for one merge pass.
#[derive(Default)]
struct Probe {
    /// Time inside `entry_range` (binary search) calls.
    range_ns: u128,
    /// Time inside the max-shared-prefix window scan.
    scan_ns: u128,
    /// Time inside the final position-emission loop.
    emit_ns: u128,
    /// Number of `entry_range` calls.
    range_calls: u64,
    /// Number of b entries visited in the max-m scan.
    scan_entries: u64,
    /// Number of `shared_prefix` calls.
    shared_prefix_calls: u64,
    /// Number of a entries processed.
    entries: u64,
    /// Number of empty-start-window fallbacks to the floor window.
    fallbacks: u64,
}

impl Probe {
    fn merge(&mut self, o: &Probe) {
        self.range_ns += o.range_ns;
        self.scan_ns += o.scan_ns;
        self.emit_ns += o.emit_ns;
        self.range_calls += o.range_calls;
        self.scan_entries += o.scan_entries;
        self.shared_prefix_calls += o.shared_prefix_calls;
        self.entries += o.entries;
        self.fallbacks += o.fallbacks;
    }
}

/// Same semantics as the library's `emit_entry_hits`, with optional skip of
/// the max-shared-prefix window scan (`m = start`) and a `Probe` accumulator.
#[allow(clippy::too_many_arguments)]
fn emit_hits<B: PgiQuery>(
    ea_kmer: u128,
    ea_freq: u32,
    a_positions: &[u64],
    b: &B,
    freq: u32,
    min_shared: usize,
    k: usize,
    prev_kmer: Option<u128>,
    skip_scan: bool,
    probe: &mut Probe,
) -> anyhow::Result<Vec<SeedHit>> {
    let mut hits = Vec::new();
    probe.entries += 1;
    if ea_freq >= freq || ea_kmer > pgr::libs::nt::rc_key(ea_kmer, k) {
        return Ok(hits);
    }
    let k_bits = 2 * k;
    let mask = if k_bits >= 128 {
        u128::MAX
    } else {
        (1u128 << k_bits) - 1
    };
    let window = |len: usize, probe: &mut Probe| -> (usize, usize) {
        let r = 1u128 << (k_bits - 2 * len);
        let lo = ea_kmer & !(r - 1) & mask;
        let hi = lo.saturating_add(r);
        let t = Instant::now();
        let r = b.entry_range(lo, hi);
        probe.range_ns += t.elapsed().as_nanos();
        probe.range_calls += 1;
        r
    };
    let start = prev_kmer
        .map(|pk| shared_prefix(pk, ea_kmer, k).max(min_shared as u32))
        .unwrap_or(min_shared as u32);
    let (mut j0, mut j) = window(start as usize, probe);
    if j0 == j && start as usize > min_shared {
        probe.fallbacks += 1;
        (j0, j) = window(min_shared, probe);
    }
    if j == j0 {
        return Ok(hits);
    }
    let m: u32 = if skip_scan {
        start
    } else {
        let mut m: u32 = 0;
        let mut i = j0;
        while i < j {
            if b.entry_freq(i) >= freq {
                i = b.entry_next(i);
                continue;
            }
            probe.scan_entries += 1;
            probe.shared_prefix_calls += 1;
            let t = Instant::now();
            m = m.max(shared_prefix(ea_kmer, b.entry_kmer(i), k));
            probe.scan_ns += t.elapsed().as_nanos();
            i = b.entry_next(i);
        }
        m
    };
    if m < min_shared as u32 {
        return Ok(hits);
    }
    let (m0, m1) = window(m as usize, probe);
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
    let t = Instant::now();
    let mut i = m0;
    while i < m1 {
        if b.entry_freq(i) < freq {
            for &arec in a_positions {
                let ac = ((arec >> 32) & CID_MASK) as u32;
                let apos = (arec & 0xffff_ffff) as u32;
                let astrand = ((arec >> STRAND_OFF) & 1) as u8;
                for brec in b.entry_positions(i) {
                    let bc = ((brec >> 32) & CID_MASK) as u32;
                    let bpos = (brec & 0xffff_ffff) as u32;
                    let bstrand = ((brec >> STRAND_OFF) & 1) as u8;
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
    probe.emit_ns += t.elapsed().as_nanos();
    Ok(hits)
}

fn merge_all(
    a: &PgiIndex,
    b: &PgiIndex,
    use_lcp: bool,
    skip_scan: bool,
) -> anyhow::Result<Vec<SeedHit>> {
    let kb = a.key_bytes();
    let hits: Vec<SeedHit> = a
        .keys
        .par_chunks(kb * 4096)
        .zip(a.entries.par_chunks(4096))
        .map(|(keys_chunk, ents)| -> anyhow::Result<Vec<SeedHit>> {
            let mut hits = Vec::new();
            let mut prev = None;
            let mut probe = Probe::default();
            for (ei, ea) in ents.iter().enumerate() {
                let ea_kmer = pgr::libs::kmer::key::Kmer::from_bytes(
                    a.k,
                    &keys_chunk[ei * kb..(ei + 1) * kb],
                )
                .to_u128();
                let ap = &a.positions[ea.pos_start as usize..(ea.pos_start + ea.freq) as usize];
                hits.extend(emit_hits(
                    ea_kmer,
                    ea.freq,
                    ap,
                    b,
                    10,
                    12,
                    a.k,
                    if use_lcp { prev } else { None },
                    skip_scan,
                    &mut probe,
                )?);
                prev = Some(ea_kmer);
            }
            let _ = probe;
            Ok(hits)
        })
        .collect::<anyhow::Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect();
    Ok(hits)
}

/// Split-profile a full merge pass: times (a) entry_range, (b) max-m scan,
/// (c) shared_prefix calls, (d) emission, aggregated across rayon chunks.
fn profile_merge(
    a: &PgiIndex,
    b: &PgiIndex,
    use_lcp: bool,
) -> anyhow::Result<(Vec<SeedHit>, Probe)> {
    let kb = a.key_bytes();
    let probes: Vec<Probe> = a
        .keys
        .par_chunks(kb * 4096)
        .zip(a.entries.par_chunks(4096))
        .map(|(keys_chunk, ents)| -> anyhow::Result<Probe> {
            let mut hits = Vec::new();
            let mut prev = None;
            let mut probe = Probe::default();
            for (ei, ea) in ents.iter().enumerate() {
                let ea_kmer = pgr::libs::kmer::key::Kmer::from_bytes(
                    a.k,
                    &keys_chunk[ei * kb..(ei + 1) * kb],
                )
                .to_u128();
                let ap = &a.positions[ea.pos_start as usize..(ea.pos_start + ea.freq) as usize];
                hits.extend(emit_hits(
                    ea_kmer,
                    ea.freq,
                    ap,
                    b,
                    10,
                    12,
                    a.k,
                    if use_lcp { prev } else { None },
                    false,
                    &mut probe,
                )?);
                prev = Some(ea_kmer);
            }
            let _ = hits;
            Ok(probe)
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let mut all = Probe::default();
    for p in probes {
        all.merge(&p);
    }
    // Re-run for the hit count (the probe pass discards hits to keep the
    // timing loop lean).
    let hits = merge_all(a, b, use_lcp, false)?;
    Ok((hits, all))
}

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let ref_path = args.next().expect("ref.pgi");
    let query_path = args.next().expect("query.pgi");

    let t = Instant::now();
    let mut rf = std::fs::File::open(&ref_path)?;
    let a = PgiIndex::read(&mut rf)?;
    let mut qf = std::fs::File::open(&query_path)?;
    let b = PgiIndex::read(&mut qf)?;
    eprintln!("load both indexes: {:?}", t.elapsed());

    let t = Instant::now();
    let hits = pgr::libs::pgi::align::merge_seed_hits(&a, &b, 10, 12)?;
    eprintln!(
        "merge_seed_hits (lib, full scan): {:?}, {} seed hits",
        t.elapsed(),
        hits.len()
    );
    let t = Instant::now();
    let hits_nolcp = merge_all(&a, &b, false, false)?;
    eprintln!(
        "merge_all (no lcp, full scan): {:?}, {} seed hits",
        t.elapsed(),
        hits_nolcp.len()
    );
    let t = Instant::now();
    let hits_full = merge_all(&a, &b, true, false)?;
    eprintln!(
        "merge_all (lcp, full scan):    {:?}, {} seed hits",
        t.elapsed(),
        hits_full.len()
    );
    let t = Instant::now();
    let hits_skip = merge_all(&a, &b, true, true)?;
    eprintln!(
        "merge_all (lcp, skip scan):    {:?}, {} seed hits",
        t.elapsed(),
        hits_skip.len()
    );

    // Split profile of the production path (lcp on).
    let (p_hits, p) = profile_merge(&a, &b, true)?;
    eprintln!(
        "\nsplit profile (lcp, {:?}): {} seed hits",
        a.entries.len(),
        p_hits.len()
    );
    eprintln!(
        "  a entries: {}\n  entry_range calls: {} ({} ns/entry {:.2} call/entry)",
        p.entries,
        p.range_calls,
        p.range_ns / p.entries.max(1) as u128,
        p.range_calls as f64 / p.entries.max(1) as f64
    );
    eprintln!(
        "  range: {:>8.3} ms ({:>5.1}%)  scan: {:>8.3} ms ({:>5.1}%)  emit: {:>8.3} ms ({:>5.1}%)",
        p.range_ns as f64 / 1e6,
        p.range_ns as f64 / (p.range_ns + p.scan_ns + p.emit_ns).max(1) as f64 * 100.0,
        p.scan_ns as f64 / 1e6,
        p.scan_ns as f64 / (p.range_ns + p.scan_ns + p.emit_ns).max(1) as f64 * 100.0,
        p.emit_ns as f64 / 1e6,
        p.emit_ns as f64 / (p.range_ns + p.scan_ns + p.emit_ns).max(1) as f64 * 100.0,
    );
    eprintln!(
        "  scan_entries: {} ({:.2}/entry)  shared_prefix calls: {} ({:.2}/entry)  fallbacks: {} ({:.2}%)",
        p.scan_entries,
        p.scan_entries as f64 / p.entries.max(1) as f64,
        p.shared_prefix_calls,
        p.shared_prefix_calls as f64 / p.entries.max(1) as f64,
        p.fallbacks,
        p.fallbacks as f64 / p.entries.max(1) as f64 * 100.0,
    );
    Ok(())
}
