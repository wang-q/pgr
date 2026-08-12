//! Exact overlap detection between unitigs (OLC stage S1).
//!
//! Unitigs are error-free maximal paths from the de Bruijn graph, so
//! overlaps are detected with an exact seed-and-verify scheme: a canonical
//! k-mer index (same packed FastK layout as `libs::map`) is built over all
//! unitigs, the boundary k-mers of each unitig are looked up, and every
//! candidate is extended to the maximal exact overlap in both orientations.

use crate::libs::ds::radix_sort::radix_sort_bytes;
use crate::libs::kmer::canonical_keys;
use crate::libs::kmer::key::Kmer;
use anyhow::Result;
use rayon::prelude::*;

/// One pseudo-read (unitig) with its sequence.
pub struct Unitig {
    /// Unique name (callers must disambiguate across input files).
    pub name: String,
    /// Sequence bases (ACGT, case-insensitive; N windows are skipped).
    pub seq: Vec<u8>,
}

/// Overlap type: end-to-end (dovetail) or one sequence inside the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OverlapType {
    Dovetail,
    Contain,
}

/// A verified exact overlap between two unitigs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Overlap {
    /// Query unitig index.
    pub qid: usize,
    /// Target unitig index.
    pub tid: usize,
    /// Strand of the target relative to the query (`+` same, `-` reverse).
    pub strand: char,
    /// Query interval covered by the overlap (0-based, half-open).
    pub q_start: usize,
    pub q_end: usize,
    /// Target interval covered by the overlap.
    pub t_start: usize,
    pub t_end: usize,
    /// Overlap length in bases (`q_end - q_start`).
    pub length: usize,
    pub otype: OverlapType,
}

/// Options for [`find_overlaps`].
pub struct OverlapOptions {
    /// Seed k-mer length (clamped to the shortest unitig).
    pub seed_k: usize,
    /// Minimum accepted overlap length in bases.
    pub min_overlap: usize,
}

/// Finds all exact overlaps between `unitigs` (excluding self).
///
/// The seed k is clamped to `min(seed_k, shortest unitig length)` so every
/// unitig remains queryable. Output is sorted by (qid, tid, strand,
/// q_start, t_start) and deduplicated, so results are deterministic.
pub fn find_overlaps(unitigs: &[Unitig], opts: &OverlapOptions) -> Result<Vec<Overlap>> {
    let seed_k = opts
        .seed_k
        .min(unitigs.iter().map(|u| u.seq.len()).min().unwrap_or(0));
    anyhow::ensure!(
        seed_k >= 1,
        "cannot overlap: unitigs are empty or shorter than 1 bp"
    );

    let key_bytes = seed_k.div_ceil(4);
    let mut keys = Vec::new();
    let mut payloads = Vec::new();
    for (cid, u) in unitigs.iter().enumerate() {
        canonical_keys(&u.seq, seed_k, |pos, key| {
            keys.extend_from_slice(key.to_bytes());
            payloads.push(((cid as u64) << 32) | pos as u64);
        });
    }
    radix_sort_bytes(&mut keys, key_bytes, &mut payloads);

    let mut overlaps: Vec<Overlap> = unitigs
        .par_iter()
        .enumerate()
        .flat_map(|(qid, u)| query_unitig(qid, &u.seq, seed_k, &keys, &payloads, unitigs, opts))
        .collect();
    overlaps.sort_by(|a, b| {
        (a.qid, a.tid, a.strand, a.q_start, a.t_start, a.length)
            .cmp(&(b.qid, b.tid, b.strand, b.q_start, b.t_start, b.length))
    });
    overlaps.dedup_by(|a, b| {
        a.qid == b.qid
            && a.tid == b.tid
            && a.strand == b.strand
            && a.q_start == b.q_start
            && a.q_end == b.q_end
            && a.t_start == b.t_start
            && a.t_end == b.t_end
    });
    Ok(overlaps)
}

/// Queries the boundary k-mers of unitig `qid` and verifies every candidate.
fn query_unitig(
    qid: usize,
    seq: &[u8],
    seed_k: usize,
    keys: &[u8],
    payloads: &[u64],
    unitigs: &[Unitig],
    opts: &OverlapOptions,
) -> Vec<Overlap> {
    let n = seq.len();
    let mut out = Vec::new();
    let mut positions = vec![0usize];
    if n > seed_k {
        positions.push(n - seed_k);
    }
    for &p in &positions {
        let mut key: Option<Kmer> = None;
        canonical_keys(seq, seed_k, |pos, kmer| {
            if pos == p {
                key = Some(*kmer);
            }
        });
        let Some(key) = key else { continue };
        let mut idx = lower_bound(keys, seed_k.div_ceil(4), &key);
        let end = upper_bound(keys, seed_k.div_ceil(4), &key);
        while idx < end {
            let payload = payloads[idx];
            let tid = (payload >> 32) as usize;
            let tpos = (payload & 0xffff_ffff) as usize;
            if tid != qid {
                if let Some(ov) = verify_seed(seq, p, &unitigs[tid].seq, tpos, seed_k) {
                    if ov.length >= opts.min_overlap {
                        out.push(ov.wrap(qid, tid));
                    }
                }
            }
            idx += 1;
        }
    }
    out
}

/// Verifies the seed window `q[p..p+k]` against `t[tpos..tpos+k]` and
/// extends to the maximal exact overlap; returns it or `None`.
fn verify_seed(q: &[u8], p: usize, t: &[u8], tpos: usize, k: usize) -> Option<Overlap> {
    let n = q.len();
    let m = t.len();
    let mut best: Option<(usize, usize, usize, usize, char)> = None;
    for (plus, strand) in [(true, '+'), (false, '-')] {
        if let Some((qs, qe, ts, te)) = extend(q, p, t, tpos, k, plus) {
            let cand = (qs, qe, ts, te, strand);
            let keep = match best {
                None => true,
                Some(b) => {
                    let l = qe - qs;
                    let bl = b.1 - b.0;
                    l > bl || (l == bl && strand == '+')
                }
            };
            if keep {
                best = Some(cand);
            }
        }
    }
    let (qs, qe, ts, te, strand) = best?;
    let otype = if (qs == 0 && qe == n) || (ts == 0 && te == m) {
        OverlapType::Contain
    } else {
        OverlapType::Dovetail
    };
    Some(Overlap {
        qid: 0,
        tid: 0,
        strand,
        q_start: qs,
        q_end: qe,
        t_start: ts,
        t_end: te,
        length: qe - qs,
        otype,
    })
}

/// Extends the seed in one orientation; `plus` = same-strand alignment.
fn extend(
    q: &[u8],
    p: usize,
    t: &[u8],
    tpos: usize,
    k: usize,
    plus: bool,
) -> Option<(usize, usize, usize, usize)> {
    let n = q.len();
    let m = t.len();
    if plus {
        if q[p..p + k] != t[tpos..tpos + k] {
            return None;
        }
    } else {
        for i in 0..k {
            if q[p + i] != complement(t[tpos + k - 1 - i]) {
                return None;
            }
        }
    }
    let mut qs = p;
    let mut qe = p + k;
    let mut ts = tpos;
    let mut te = tpos + k;
    if plus {
        while qs > 0 && ts > 0 && q[qs - 1] == t[ts - 1] {
            qs -= 1;
            ts -= 1;
        }
        while qe < n && te < m && q[qe] == t[te] {
            qe += 1;
            te += 1;
        }
    } else {
        while qs > 0 && te < m && q[qs - 1] == complement(t[te]) {
            qs -= 1;
            te += 1;
        }
        while qe < n && ts > 0 && q[qe] == complement(t[ts - 1]) {
            qe += 1;
            ts -= 1;
        }
    }
    Some((qs, qe, ts, te))
}

/// Complement of an ASCII base (both cases).
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        x => x,
    }
}

/// First index whose packed key is `>= key`.
fn lower_bound(keys: &[u8], key_bytes: usize, key: &Kmer) -> usize {
    let mut lo = 0usize;
    let mut hi = keys.len() / key_bytes;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if &keys[mid * key_bytes..(mid + 1) * key_bytes] < key.to_bytes() {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

/// First index whose packed key is `> key`.
fn upper_bound(keys: &[u8], key_bytes: usize, key: &Kmer) -> usize {
    let mut lo = 0usize;
    let mut hi = keys.len() / key_bytes;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if &keys[mid * key_bytes..(mid + 1) * key_bytes] <= key.to_bytes() {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

impl Overlap {
    /// Wraps a verified interval into a full record with q/t ids.
    fn wrap(mut self, qid: usize, tid: usize) -> Self {
        self.qid = qid;
        self.tid = tid;
        self
    }
}

/// Drops unitigs fully contained in a longer (or equal-length, lower-id)
/// unitig and remaps the surviving ids.
///
/// Containment is read from `ov:A:C` overlaps that cover the whole shorter
/// sequence. A contained unitig's sequence and every exact overlap it has
/// are also present in its container, so no sequence content is lost; the
/// overlap graph shrinks and greedy path selection changes (multi-k
/// redundancy is the main thing that used to break chains — Lambda: 90 ->
/// 22 unitigs, 16 fragments -> 1 full-genome contig containing all 16).
pub fn filter_contained(unitigs: &[Unitig], overlaps: &[Overlap]) -> (Vec<Unitig>, Vec<Overlap>) {
    let n = unitigs.len();
    let mut contained = vec![false; n];
    for ov in overlaps {
        if ov.otype != OverlapType::Contain {
            continue;
        }
        let qlen = unitigs[ov.qid].seq.len();
        let tlen = unitigs[ov.tid].seq.len();
        if ov.q_start == 0 && ov.q_end == qlen {
            mark_contained(&mut contained, ov.qid, ov.tid, qlen, tlen);
        }
        if ov.t_start == 0 && ov.t_end == tlen {
            mark_contained(&mut contained, ov.tid, ov.qid, tlen, qlen);
        }
    }

    let mut new_id = vec![usize::MAX; n];
    let mut filtered = Vec::with_capacity(n);
    for (i, u) in unitigs.iter().enumerate() {
        if !contained[i] {
            new_id[i] = filtered.len();
            filtered.push(Unitig {
                name: u.name.clone(),
                seq: u.seq.clone(),
            });
        }
    }
    let filtered_overlaps = overlaps
        .iter()
        .filter_map(|ov| {
            let qid = new_id[ov.qid];
            let tid = new_id[ov.tid];
            (qid != usize::MAX && tid != usize::MAX).then_some(Overlap {
                qid,
                tid,
                strand: ov.strand,
                q_start: ov.q_start,
                q_end: ov.q_end,
                t_start: ov.t_start,
                t_end: ov.t_end,
                length: ov.length,
                otype: ov.otype,
            })
        })
        .collect();
    (filtered, filtered_overlaps)
}

/// Marks `c` as contained in `d` when `d` is longer (or equal with a lower
/// id, so identical duplicates keep the first one).
fn mark_contained(contained: &mut [bool], c: usize, d: usize, clen: usize, dlen: usize) {
    if !contained[c] && (dlen > clen || (dlen == clen && d < c)) {
        contained[c] = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unitigs(names: &[&str], seqs: &[&str]) -> Vec<Unitig> {
        names
            .iter()
            .zip(seqs)
            .map(|(n, s)| Unitig {
                name: (*n).to_string(),
                seq: s.as_bytes().to_vec(),
            })
            .collect()
    }

    fn find(us: &[Unitig], seed_k: usize, min_overlap: usize) -> Vec<Overlap> {
        find_overlaps(
            us,
            &OverlapOptions {
                seed_k,
                min_overlap,
            },
        )
        .unwrap()
    }

    /// Two unitigs sharing a 10 bp suffix/prefix overlap.
    #[test]
    fn dovetail_suffix_prefix() {
        let a = "TTTTTTTTTTACGTACGTAC"; // suffix 10 bp
        let b = "ACGTACGTACGGGGGGGGGG"; // prefix 10 == a suffix 10
        let us = unitigs(&["a", "b"], &[a, b]);
        let ovs = find(&us, 5, 8);
        let ov = ovs
            .iter()
            .find(|o| o.qid == 0 && o.tid == 1)
            .expect("a->b overlap");
        // a suffix (10..20) == b prefix (0..10) == ACGTACGTAC
        assert_eq!(ov.strand, '+');
        assert_eq!(ov.q_start, 10);
        assert_eq!(ov.q_end, 20);
        assert_eq!(ov.t_start, 0);
        assert_eq!(ov.t_end, 10);
        assert_eq!(ov.length, 10);
        assert_eq!(ov.otype, OverlapType::Dovetail);
    }

    /// Reverse-complement orientation is found and verified.
    #[test]
    fn dovetail_reverse_strand() {
        // a's suffix "ACGTACGTAC" matches rc(b) where b starts with that rc.
        let a = "TTTTTTTTTTACGTACGTAC"; // suffix 10 bp
        let b = "GTACGTACGTCCCCCCCCCC"; // rc of "ACGTACGTAC" = "GTACGTACGT"
        let us = unitigs(&["a", "b"], &[a, b]);
        let ovs = find(&us, 5, 8);
        let ov = ovs
            .iter()
            .find(|o| o.qid == 0 && o.tid == 1 && o.length == 10)
            .expect("a->b 10 bp overlap");
        // a's suffix rc-matches b's prefix; a shorter same-strand overlap
        // also exists, but the longest one is the reverse-strand dovetail.
        assert_eq!(ov.strand, '-');
        assert_eq!(ov.q_start, 10);
        assert_eq!(ov.q_end, 20);
        assert_eq!(ov.t_start, 0);
        assert_eq!(ov.t_end, 10);
        assert_eq!(ov.length, 10);
    }

    /// A short unitig contained in a longer one.
    #[test]
    fn contain_classified() {
        let a = "AAAAGGTTAACCGGTTCCCC"; // b = a[4..16]
        let b = "GGTTAACCGGTT";
        let us = unitigs(&["a", "b"], &[a, b]);
        let ovs = find(&us, 5, 10);
        let ov = ovs
            .iter()
            .find(|o| o.qid == 1 && o.tid == 0)
            .expect("b inside a");
        assert_eq!(ov.otype, OverlapType::Contain);
        assert_eq!(ov.q_start, 0);
        assert_eq!(ov.q_end, 12);
        assert_eq!(ov.t_start, 4);
        assert_eq!(ov.t_end, 16);
    }

    /// No self overlaps and nothing below the minimum length.
    #[test]
    fn filters_self_and_short() {
        let a = "TTTTTTTTTTACGTACGTAC";
        let b = "ACGTACGTACGGGGGGGGGG";
        let us = unitigs(&["a", "b"], &[a, b]);
        let ovs = find(&us, 5, 10);
        assert!(!ovs.iter().any(|o| o.qid == o.tid));
        assert!(ovs.iter().all(|o| o.length >= 10));
    }

    /// Repetitive seeds still need full extension to the true boundaries.
    #[test]
    fn repeat_kmers_need_full_verification() {
        let a = "AAAAAAAAAAAAAAAAAAAA"; // 20 A
        let b = "AAAAAAAAAAAAAAAAAAAAA"; // 21 A
        let us = unitigs(&["a", "b"], &[a, b]);
        let ovs = find(&us, 5, 5);
        let ov = ovs
            .iter()
            .find(|o| o.qid == 0 && o.tid == 1)
            .expect("a->b overlap");
        // b contains a: a[0..20] == b[0..20]
        assert_eq!(ov.q_start, 0);
        assert_eq!(ov.q_end, 20);
        assert_eq!(ov.t_start, 0);
        assert_eq!(ov.t_end, 20);
        assert_eq!(ov.otype, OverlapType::Contain);
    }

    /// Contained unitigs are dropped, ids remapped, overlaps filtered.
    #[test]
    fn filter_contained_drops_and_remaps() {
        let us = unitigs(
            &["u0", "u1", "u2"],
            &["AAAACCCCGGGGTTTT", "CCCCGGGG", "AAAACCCCAAAA"],
        );
        let ovs = find(&us, 5, 8);
        let (fu, fo) = filter_contained(&us, &ovs);
        assert_eq!(fu.len(), 2);
        assert_eq!(fu[0].name, "u0");
        assert_eq!(fu[1].name, "u2");
        assert!(fo.iter().all(|o| o.qid < 2 && o.tid < 2));
        assert!(fo.iter().all(|o| o.otype != OverlapType::Contain));
    }

    /// Identical unitigs: the lower id wins.
    #[test]
    fn filter_contained_equal_duplicates() {
        let seq = "AAAACCCCGGGGTTTT";
        let us = unitigs(&["a", "b", "c"], &[seq, seq, "TTTTGGGGCCCCAAAA"]);
        let ovs = find(&us, 5, 8);
        let (fu, _) = filter_contained(&us, &ovs);
        assert_eq!(fu.len(), 2);
        assert_eq!(fu[0].name, "a");
        assert_eq!(fu[1].name, "c");
    }
}
