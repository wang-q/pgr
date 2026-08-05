//! In-memory `merge_seed_hits` timing vs the streaming command path, plus a
//! skip-scan variant to isolate the cost of the max-shared-prefix window
//! scan (the part the LCP optimizations touch).
//!
//! Usage: cargo run --release --example merge_mem_bench -- ref.pgi query.pgi

use pgr::libs::pgi::align::merge_seed_hits;
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

/// Simplified `emit_entry_hits` (same semantics as the library), with an
/// optional skip of the max-shared-prefix window scan (`m = start`).
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
) -> anyhow::Result<Vec<SeedHit>> {
    let mut hits = Vec::new();
    if ea_freq >= freq || ea_kmer > pgr::libs::nt::rc_key(ea_kmer, k) {
        return Ok(hits);
    }
    let k_bits = 2 * k;
    let mask = if k_bits >= 128 {
        u128::MAX
    } else {
        (1u128 << k_bits) - 1
    };
    let window = |len: usize| -> (usize, usize) {
        let r = 1u128 << (k_bits - 2 * len);
        let lo = ea_kmer & !(r - 1) & mask;
        let hi = lo.saturating_add(r);
        b.entry_range(lo, hi)
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
            m = m.max(shared_prefix(ea_kmer, b.entry_kmer(i), k));
            i = b.entry_next(i);
        }
        m
    };
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
    Ok(hits)
}

fn merge_all(
    a: &PgiIndex,
    b: &PgiIndex,
    use_lcp: bool,
    skip_scan: bool,
) -> anyhow::Result<Vec<SeedHit>> {
    let hits: Vec<SeedHit> = a
        .entries
        .par_chunks(4096)
        .map(|ents| -> anyhow::Result<Vec<SeedHit>> {
            let mut hits = Vec::new();
            let mut prev = None;
            let mut fallbacks = 0usize;
            for ea in ents {
                let ap = &a.positions[ea.pos_start as usize..(ea.pos_start + ea.freq) as usize];
                hits.extend(emit_hits(
                    ea.kmer,
                    ea.freq,
                    ap,
                    b,
                    10,
                    12,
                    a.k,
                    if use_lcp { prev } else { None },
                    skip_scan,
                )?);
                if use_lcp {
                    let k = a.k;
                    let start = prev
                        .map(|pk| shared_prefix(pk, ea.kmer, k).max(12u32))
                        .unwrap_or(12u32) as usize;
                    let k_bits = 2 * k;
                    let mask = if k_bits >= 128 {
                        u128::MAX
                    } else {
                        (1u128 << k_bits) - 1
                    };
                    let r = 1u128 << (k_bits - 2 * start);
                    let lo = ea.kmer & !(r - 1) & mask;
                    let hi = lo.saturating_add(r);
                    let (w0, w1) = b.entry_range(lo, hi);
                    if w0 == w1 {
                        fallbacks += 1;
                    }
                }
                prev = Some(ea.kmer);
            }
            if fallbacks > 0 {
                eprintln!("  chunk fallbacks: {fallbacks}");
            }
            Ok(hits)
        })
        .collect::<anyhow::Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect();
    Ok(hits)
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
    let hits = merge_seed_hits(&a, &b, 10, 12)?;
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
    Ok(())
}
