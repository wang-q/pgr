//! Phase-0 PoC: does the cross-genome *shared syncmer* signal from two `.pgi`
//! locate collinear regions that the pairwise `align pgi` PSL anchors miss?
//!
//! Mode A (3 args): shared-syncmer coverage vs anchored PSL coverage per genome.
//! Mode B (4 args, + lastz.psl): classify the PAF-uncovered region on the
//!   target (genome1) side against lastz ground truth:
//!     uncovered = genome1 − pgi anchors
//!     capturable  = uncovered ∩ lastz            (high-divergence collinear)
//!     nonortholog = uncovered − lastz            (strain-specific / no ortholog)
//!
//! Usage:
//!   synteny_poc <pgi1> <pgi2> <anchors.psl>
//!   synteny_poc <pgi1> <pgi2> <anchors.psl> <lastz.psl>
//!
//! pgi1 = reference (PSL target side), pgi2 = query (PSL query side).

use pgr::libs::pgi::mmap::PgiMmap;
use pgr::libs::pgi::PgiQuery;
use std::collections::HashMap;
use std::path::Path;

const CID_MASK: u64 = (1 << 20) - 1;
const POS_MASK: u64 = (1 << 32) - 1;
const STRAND_OFF: u32 = 52;

fn unpack(rec: u64) -> (u32, u32, u8) {
    (
        ((rec >> 32) & CID_MASK) as u32,
        (rec & POS_MASK) as u32,
        ((rec >> STRAND_OFF) & 1) as u8,
    )
}

/// Per-contig sorted non-overlapping [start,end) intervals.
#[derive(Default, Clone)]
struct ContigCover {
    ints: Vec<(u32, u32)>,
}

impl ContigCover {
    fn covered(&self, pos: u32) -> bool {
        self.ints
            .binary_search_by(|&(s, e)| {
                if pos < s {
                    std::cmp::Ordering::Greater
                } else if pos >= e {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .is_ok()
    }
    fn bp(&self) -> u64 {
        self.ints.iter().map(|&(s, e)| (e - s) as u64).sum()
    }
}

/// Build per-contig coverage from a PSL. `name_i`/`start_i`/`end_i` are the
/// whitespace field indices for the side of interest.
fn psl_coverage(
    path: &Path,
    contig_names: &[String],
    name_i: usize,
    start_i: usize,
    end_i: usize,
) -> Vec<ContigCover> {
    let name_to_id: HashMap<&str, usize> = contig_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let mut out: Vec<ContigCover> = (0..contig_names.len())
        .map(|_| ContigCover::default())
        .collect();
    let mut raw: Vec<(usize, u32, u32)> = Vec::new();
    let data = std::fs::read_to_string(path).expect("read psl");
    for line in data.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() <= end_i {
            continue;
        }
        let Some(&id) = name_to_id.get(f[name_i]) else {
            continue;
        };
        let start: u32 = f[start_i].parse().unwrap();
        let end: u32 = f[end_i].parse().unwrap();
        raw.push((id, start, end));
    }
    for &(id, s, e) in &raw {
        out[id].ints.push((s, e));
    }
    for cov in out.iter_mut() {
        cov.ints.sort_unstable();
        let mut merged: Vec<(u32, u32)> = Vec::new();
        for (s, e) in cov.ints.drain(..) {
            if let Some(last) = merged.last_mut() {
                if s <= last.1 {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            merged.push((s, e));
        }
        cov.ints = merged;
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert!(
        args.len() == 4 || args.len() == 5,
        "usage: synteny_poc <pgi1> <pgi2> <anchors.psl> [lastz.psl]"
    );
    let pgi1 = PgiMmap::open(Path::new(&args[1])).expect("open pgi1");
    let pgi2 = PgiMmap::open(Path::new(&args[2])).expect("open pgi2");

    let contigs1: Vec<String> = pgi1.contigs().iter().map(|(n, _)| n.clone()).collect();
    let contigs2: Vec<String> = pgi2.contigs().iter().map(|(n, _)| n.clone()).collect();
    let len1: Vec<u64> = pgi1.contigs().iter().map(|(_, l)| *l).collect();

    let cov1 = psl_coverage(Path::new(&args[3]), &contigs1, 13, 15, 16);
    let cov2 = psl_coverage(Path::new(&args[3]), &contigs2, 9, 11, 12);

    // --- Mode A: shared-syncmer coverage ---
    let (i0, i1) = pgi1.entry_range(0, u128::MAX);
    let mut n_shared_kmers = 0u64;
    let mut g1_shared_pos = 0u64;
    let mut g1_shared_cov = 0u64;
    let mut g1_shared_uncov = 0u64;
    let mut i = i0;
    while i < i1 {
        let kmer = pgi1.entry_kmer(i);
        let (j0, j1) = pgi2.entry_range(kmer, kmer + 1);
        if j0 < j1 {
            n_shared_kmers += 1;
            let npos = pgi1.entry_positions(i).count() as u64;
            g1_shared_pos += npos;
            for rec in pgi1.entry_positions(i) {
                let (cid, pos, _) = unpack(rec);
                if (cid as usize) < cov1.len() && cov1[cid as usize].covered(pos) {
                    g1_shared_cov += 1;
                } else {
                    g1_shared_uncov += 1;
                }
            }
        }
        i = pgi1.entry_next(i);
    }

    let total1: u64 = len1.iter().sum();
    let cov_len1: u64 = cov1.iter().map(|c| c.bp()).sum();
    let total2: u64 = pgi2.contigs().iter().map(|(_, l)| *l).sum();
    let cov_len2: u64 = cov2.iter().map(|c| c.bp()).sum();
    println!(
        "genome1 total_len={} anchor_cov={} ({:.2}%)",
        total1,
        cov_len1,
        100.0 * cov_len1 as f64 / total1 as f64
    );
    println!(
        "genome2 total_len={} anchor_cov={} ({:.2}%)",
        total2,
        cov_len2,
        100.0 * cov_len2 as f64 / total2 as f64
    );
    println!("shared kmers (present in both): {}", n_shared_kmers);
    println!(
        "genome1 shared-pos: {} | covered {} ({:.2}%) | UNCOVERED {} ({:.2}%)",
        g1_shared_pos,
        g1_shared_cov,
        100.0 * g1_shared_cov as f64 / g1_shared_pos as f64,
        g1_shared_uncov,
        100.0 * g1_shared_uncov as f64 / g1_shared_pos as f64
    );

    // --- Mode B: uncovered-region composition (genome1 / target side) ---
    if args.len() == 5 {
        let lz = psl_coverage(Path::new(&args[4]), &contigs1, 13, 15, 16);
        let mut uncov = 0u64; // not covered by pgi anchors
        let mut capturable = 0u64; // uncovered & lastz
        let mut nonortholog = 0u64; // uncovered & !lastz
        for (cid, l) in len1.iter().enumerate() {
            for pos in 0..*l {
                if cov1[cid].covered(pos as u32) {
                    continue;
                }
                uncov += 1;
                if lz[cid].covered(pos as u32) {
                    capturable += 1;
                } else {
                    nonortholog += 1;
                }
            }
        }
        let lz_bp: u64 = lz.iter().map(|c| c.bp()).sum();
        println!("---- genome1 uncovered-region composition (vs lastz) ----");
        println!(
            "lastz cov: {} ({:.2}%)",
            lz_bp,
            100.0 * lz_bp as f64 / total1 as f64
        );
        println!(
            "uncovered by pgi anchors: {} ({:.2}%)",
            uncov,
            100.0 * uncov as f64 / total1 as f64
        );
        println!(
            "  capturable (uncovered & lastz): {} ({:.2}% of genome; {:.2}% of uncovered)",
            capturable,
            100.0 * capturable as f64 / total1 as f64,
            100.0 * capturable as f64 / uncov as f64
        );
        println!(
            "  non-ortholog (uncovered & !lastz): {} ({:.2}% of genome; {:.2}% of uncovered)",
            nonortholog,
            100.0 * nonortholog as f64 / total1 as f64,
            100.0 * nonortholog as f64 / uncov as f64
        );
    }
}
