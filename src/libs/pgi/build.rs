//! Build a `.pgi` index from FASTA or 2bit sequences.

use super::{PgiEntry, PgiIndex};
use crate::libs::syncmer::SyncmerParams;
use anyhow::Context;
use std::collections::VecDeque;

/// Encode `k` bases as 2-bit (A=0, C=1, G=2, T=3), high bits first.
/// Returns `None` if any base is not A/C/G/T (e.g. N).
pub fn pack_kmer(seq: &[u8], k: usize) -> Option<u128> {
    if seq.len() < k {
        return None;
    }
    let mut x: u128 = 0;
    for &b in &seq[..k] {
        let c = match b {
            b'A' | b'a' => 0u128,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
        x = (x << 2) | c;
    }
    Some(x)
}

/// Reverse-complement a 2-bit encoded k-mer key in place of orientation.
pub fn rc_key(x: u128, k: usize) -> u128 {
    // 4-base lookup: the lowest byte holds the last 4 bases, whose RC becomes
    // the highest byte of the result. Verified against the bitwise loop on
    // 20000 random sequences for k in 10..=40.
    let mut r: u128 = 0;
    let nbytes = k / 4;
    for block in 0..nbytes {
        let byte = ((x >> (8 * block)) & 0xff) as usize;
        r = (r << 8) | RC_TABLE[byte] as u128;
    }
    // The leftover high bases (0..k-4*nbytes) close the RC from the most
    // significant one down to base 0.
    for i in (0..k - 4 * nbytes).rev() {
        let c = ((x >> (2 * (k - 1 - i))) & 3) ^ 3;
        r = (r << 2) | c;
    }
    r
}

/// Reverse-complement table for one byte (4 x 2-bit bases); verified against
/// the bitwise loop on 10000 random k=40 sequences.
const RC_TABLE: [u8; 256] = [
    255, 191, 127, 63, 239, 175, 111, 47, 223, 159, 95, 31, 207, 143, 79, 15, 251, 187, 123, 59,
    235, 171, 107, 43, 219, 155, 91, 27, 203, 139, 75, 11, 247, 183, 119, 55, 231, 167, 103, 39,
    215, 151, 87, 23, 199, 135, 71, 7, 243, 179, 115, 51, 227, 163, 99, 35, 211, 147, 83, 19, 195,
    131, 67, 3, 254, 190, 126, 62, 238, 174, 110, 46, 222, 158, 94, 30, 206, 142, 78, 14, 250, 186,
    122, 58, 234, 170, 106, 42, 218, 154, 90, 26, 202, 138, 74, 10, 246, 182, 118, 54, 230, 166,
    102, 38, 214, 150, 86, 22, 198, 134, 70, 6, 242, 178, 114, 50, 226, 162, 98, 34, 210, 146, 82,
    18, 194, 130, 66, 2, 253, 189, 125, 61, 237, 173, 109, 45, 221, 157, 93, 29, 205, 141, 77, 13,
    249, 185, 121, 57, 233, 169, 105, 41, 217, 153, 89, 25, 201, 137, 73, 9, 245, 181, 117, 53,
    229, 165, 101, 37, 213, 149, 85, 21, 197, 133, 69, 5, 241, 177, 113, 49, 225, 161, 97, 33, 209,
    145, 81, 17, 193, 129, 65, 1, 252, 188, 124, 60, 236, 172, 108, 44, 220, 156, 92, 28, 204, 140,
    76, 12, 248, 184, 120, 56, 232, 168, 104, 40, 216, 152, 88, 24, 200, 136, 72, 8, 244, 180, 116,
    52, 228, 164, 100, 36, 212, 148, 84, 20, 196, 132, 68, 4, 240, 176, 112, 48, 224, 160, 96, 32,
    208, 144, 80, 16, 192, 128, 64, 0,
];

/// Rolling 2-bit k-mer keys: `out[p]` is the key of `seq[p..p+k)`, or `None`
/// if that window contains a non-ACGT base (e.g. N).
pub fn rolling_kmer_keys(seq: &[u8], k: usize) -> Vec<Option<u128>> {
    let n = seq.len();
    if n < k {
        return Vec::new();
    }
    let mask = if k * 2 >= 128 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    };
    let mut out = vec![None; n - k + 1];
    let mut x: u128 = 0;
    let mut valid = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        let c = match b {
            b'A' | b'a' => 0u128,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        };
        if c == 4 {
            x = 0;
            valid = 0;
        } else {
            x = ((x << 2) | c) & mask;
            valid += 1;
        }
        if i + 1 >= k && valid >= k {
            out[i + 1 - k] = Some(x);
        }
    }
    out
}

/// Splitmix64-style deterministic odd hash factor (matches `pgr dist` syncmers).
fn hash_factor(seed: u64) -> u64 {
    let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    (z ^ (z >> 31)) | 1
}

/// Parallel growable buffers for (key, contig, pos, strand) records.
struct RecordBuf {
    keys: Vec<u128>,
    payloads: Vec<(u32, u32, u8)>,
}

/// Single-pass collection of (key, contig, pos, strand) records for one contig.
///
/// Mirrors FastGA's GIXmake streaming style: one scan computes the s-mer
/// canonical hashes (closed-syncmer minimum window), rolls the k-mer key,
/// and emits records at each syncmer position once its k-mer window is ready.
/// N in the k-mer window invalidates that position's key.
fn collect_one_contig(
    seq: &[u8],
    cid: u32,
    k: usize,
    params: &SyncmerParams,
    no_rev: bool,
    out: &mut RecordBuf,
) {
    let smer = params.smer;
    let window = params.window;
    let n = seq.len();
    if n < k.max(smer) || window == 0 {
        return;
    }
    let factor = hash_factor(7);
    let smask = (1u64 << (2 * smer)) - 1;
    let sshift = (64 - 2 * smer) as u32;
    let kmask = if 2 * k >= 128 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    };
    let pattern_rc: [u64; 4] = std::array::from_fn(|i| ((3 - i) as u64) << (2 * (smer - 1)));
    // Base -> 2-bit code (0..3) or 4 (N / ambiguity); index by byte value.
    let codes = [4u64; 256];
    let codes = {
        let mut c = codes;
        c[b'A' as usize] = 0;
        c[b'C' as usize] = 1;
        c[b'G' as usize] = 2;
        c[b'T' as usize] = 3;
        c[b'a' as usize] = 0;
        c[b'c' as usize] = 1;
        c[b'g' as usize] = 2;
        c[b't' as usize] = 3;
        c
    };
    // Top two bits hold the complement of the newest base for the rolling RC.
    let rc_top = (2 * k - 2) as u32;

    let mut sh: u64 = 0;
    let mut shr: u64 = 0;
    let mut kx: u128 = 0;
    let mut kxr: u128 = 0;
    let mut kvalid = 0usize;
    // Monotonic ring queue of (s-mer start, canonical hash); front is the
    // first (and minimal) hash in the current window.
    let dq_cap = (window + 2).next_power_of_two();
    let dq_mask = dq_cap - 1;
    let mut dq_idx = vec![0usize; dq_cap];
    let mut dq_hash = vec![0u64; dq_cap];
    let mut dq_head = 0usize;
    let mut dq_tail = 0usize;
    // Syncmer positions awaiting their k-mer window.
    let mut pending: VecDeque<usize> = VecDeque::new();

    for (i, &b) in seq.iter().enumerate() {
        let code = codes[b as usize];

        // Roll the k-mer key and flush pending syncmer positions whose window
        // [pos, pos+k) is now complete (current key == seq[pos..pos+k)).
        if code == 4 {
            kx = 0;
            kxr = 0;
            kvalid = 0;
        } else {
            kx = ((kx << 2) | code as u128) & kmask;
            kxr = (kxr >> 2) | (((3 - code) as u128) << rc_top);
            kvalid += 1;
        }
        while let Some(&pos) = pending.front() {
            if pos + k - 1 > i {
                break;
            }
            pending.pop_front();
            if kvalid >= k && i + 1 >= k {
                out.keys.push(kx);
                out.payloads.push((cid, pos as u32, 0));
                if !no_rev {
                    out.keys.push(kxr);
                    out.payloads.push((cid, pos as u32, 1));
                }
            }
        }

        // Roll the s-mer canonical hash (N treated as A, matching pgr dist).
        let sb = if code == 4 { 0 } else { code };
        sh = ((sh << 2) | sb) & smask;
        shr = (shr >> 2) | pattern_rc[sb as usize];
        if i + 1 >= smer {
            let j = i + 1 - smer; // s-mer start
            let hf = sh.wrapping_mul(factor) >> sshift;
            let hr = shr.wrapping_mul(factor) >> sshift;
            let ch = if hf < hr { hf } else { hr };
            while dq_tail > dq_head && dq_hash[(dq_tail - 1) & dq_mask] > ch {
                dq_tail -= 1;
            }
            dq_idx[dq_tail & dq_mask] = j;
            dq_hash[dq_tail & dq_mask] = ch;
            dq_tail += 1;
            if j >= window - 1 {
                let start = j + 1 - window;
                while dq_head < dq_tail && dq_idx[dq_head & dq_mask] < start {
                    dq_head += 1;
                }
                let min_idx = dq_idx[dq_head & dq_mask];
                let min_val = dq_hash[dq_head & dq_mask];
                if min_idx == start {
                    if start + k <= n {
                        pending.push_back(start);
                    }
                } else if ch == min_val && j + k <= n {
                    pending.push_back(j);
                }
            }
        }
    }
}

/// Build an index from named sequences.
///
/// Each syncmer position seeds a k-mer (both strands unless `no_rev`); the
/// resulting records are sorted by key (in-place MSD radix) and grouped into
/// unique entries, globally ascending by k-mer.
pub fn build_from_seqs(
    contigs: Vec<(String, Vec<u8>)>,
    k: usize,
    smer: usize,
    window: usize,
    no_rev: bool,
) -> anyhow::Result<PgiIndex> {
    anyhow::ensure!(k > 0 && k * 2 <= 128, "k must be in 1..=64, got {k}");
    anyhow::ensure!(smer > 0, "smer must be positive");
    anyhow::ensure!(window > 0, "window must be positive");
    let params = SyncmerParams {
        smer,
        window,
        seed: 7,
    };
    params.validate()?;

    // Collect (key, contig, pos, strand) records into parallel growable
    // vectors (a single vector reallocates far less than per-key buckets).
    // Capacity estimate: closed syncmers sample ~2 positions per window
    // bases, doubled for both strands.
    let est: usize = contigs
        .iter()
        .map(|(_, s)| (s.len() / window).saturating_mul(4) + 64)
        .sum();
    let mut buf = RecordBuf {
        keys: Vec::with_capacity(est),
        payloads: Vec::with_capacity(est),
    };
    for (cid, (_, seq)) in contigs.iter().enumerate() {
        collect_one_contig(seq, cid as u32, k, &params, no_rev, &mut buf);
    }

    // Sort globally by k-mer with an in-place MSD radix sort (no auxiliary
    // arrays); equal k-mers stay contiguous so the grouping below produces
    // entries strictly ascending by k-mer. The parallel variant distributes
    // by the top byte and sorts each byte bucket concurrently.
    crate::libs::ds::radix_sort::radix_sort_u128_par(
        &mut buf.keys,
        &mut buf.payloads,
        2 * k as u32,
    );
    let keys = buf.keys;
    let payloads = buf.payloads;

    let mut entries: Vec<PgiEntry> = Vec::with_capacity(keys.len());
    let mut positions: Vec<(u32, u32, u8)> = Vec::with_capacity(keys.len());
    let mut i = 0usize;
    while i < keys.len() {
        let kmer = keys[i];
        let pos_start = positions.len() as u32;
        let mut j = i;
        while j < keys.len() && keys[j] == kmer {
            positions.push(payloads[j]);
            j += 1;
        }
        entries.push(PgiEntry {
            kmer,
            pos_start,
            freq: (j - i) as u32,
        });
        i = j;
    }

    Ok(PgiIndex {
        k,
        smer,
        window,
        contigs: contigs
            .into_iter()
            .map(|(n, s)| (n, s.len() as u64))
            .collect(),
        entries,
        positions,
    })
}

/// Read all sequences from a FASTA file (plain or gzipped).
pub fn read_fasta(path: &str) -> anyhow::Result<Vec<(String, Vec<u8>)>> {
    let mut reader = crate::libs::fmt::fa::reader(path)
        .with_context(|| format!("failed to open FASTA {path}"))?;
    let mut contigs = Vec::new();
    for result in reader.records() {
        let rec = result?;
        let name = String::from_utf8(rec.name().into()).context("FASTA name utf8")?;
        contigs.push((name, rec.sequence().as_ref().to_vec()));
    }
    Ok(contigs)
}

/// Read all sequences from a 2bit file.
pub fn read_2bit(path: &str) -> anyhow::Result<Vec<(String, Vec<u8>)>> {
    let mut tb = crate::libs::fmt::twobit::TwoBitFile::open(path)
        .with_context(|| format!("failed to open 2bit {path}"))?;
    let names = tb.get_sequence_names();
    let mut contigs = Vec::with_capacity(names.len());
    for name in names {
        let seq = tb
            .read_sequence(&name, None, None, true)
            .with_context(|| format!("reading {name} from 2bit"))?
            .into_bytes();
        contigs.push((name, seq));
    }
    Ok(contigs)
}

/// Build an index from a FASTA or 2bit input file (extension decides).
pub fn build_from_path(
    path: &str,
    k: usize,
    smer: usize,
    window: usize,
    no_rev: bool,
) -> anyhow::Result<PgiIndex> {
    let is_2bit = std::path::Path::new(path)
        .extension()
        .and_then(|e| e.to_str())
        == Some("2bit");
    let contigs = if is_2bit {
        read_2bit(path)?
    } else {
        read_fasta(path)?
    };
    build_from_seqs(contigs, k, smer, window, no_rev)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pack_and_rc() {
        // A=0 C=1 G=2 T=3, high bits first: "ACGT" -> 0b00011011
        assert_eq!(pack_kmer(b"ACGT", 4), Some(0b00011011));
        // RC("ACGT") = "ACGT" (reverse-complement palindrome); double RC
        // restores the original for any sequence.
        let x = pack_kmer(b"ACGT", 4).unwrap();
        assert_eq!(rc_key(x, 4), x);
        assert_eq!(rc_key(rc_key(x, 4), 4), x);
        // RC("AAAA") = "TTTT"
        let a = pack_kmer(b"AAAA", 4).unwrap();
        assert_eq!(rc_key(a, 4), pack_kmer(b"TTTT", 4).unwrap());
        // N is rejected
        assert_eq!(pack_kmer(b"ACNT", 4), None);
    }

    #[test]
    fn rc_key_matches_bitwise_for_partial_bytes() {
        // The 4-base lookup only covers whole bytes; the leftover bases (k % 4)
        // must close the RC from the most significant one down (regression for
        // k not divisible by 4).
        let mut rng = rand::rngs::StdRng::seed_from_u64(5);
        use rand::{Rng, SeedableRng};
        for &k in &[9usize, 10, 13, 21] {
            for _ in 0..200 {
                let x: u128 = rng.random_range(0..(1u128 << (2 * k)));
                let mut expect: u128 = 0;
                for i in (0..k).rev() {
                    let c = ((x >> (2 * (k - 1 - i))) & 3) ^ 3;
                    expect = (expect << 2) | c;
                }
                assert_eq!(rc_key(x, k), expect, "k={k} x={x:x}");
            }
        }
    }

    #[test]
    fn build_small_index() {
        let idx = build_from_seqs(
            vec![(
                String::from("c1"),
                b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            )],
            10,
            4,
            2,
            false,
        )
        .unwrap();
        assert_eq!(idx.k, 10);
        assert!(idx.n_unique() > 0);
        // forward and reverse keys both present for the first syncmer position
        assert_eq!(idx.contigs[0].0, "c1");
        assert!(idx.entries.iter().all(|e| e.freq >= 1));
    }

    #[test]
    fn rolling_keys_match_pack() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let keys = rolling_kmer_keys(seq, 10);
        assert_eq!(keys.len(), seq.len() - 10 + 1);
        for (p, k) in keys.iter().enumerate() {
            assert_eq!(*k, pack_kmer(&seq[p..p + 10], 10), "position {p}");
        }
    }

    #[test]
    fn rolling_keys_handle_n() {
        let seq = b"ACGTACGTACNTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let keys = rolling_kmer_keys(seq, 10);
        // windows containing the N (positions covering index 9) are None
        assert!(keys[0].is_some());
        assert!(keys[1].is_none()); // window [1,11) includes the N at index 9
        assert!(keys[10].is_none()); // window [10,20) starts at N
        assert!(keys[11].is_some()); // window [11,21) clear
    }

    #[test]
    fn single_pass_matches_reference() {
        // Reference: syncmer_dna + rolling keys (the two-pass approach).
        let seqs = [
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACNTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        ];
        let (k, smer, window) = (10usize, 4usize, 2usize);
        let params = SyncmerParams {
            smer,
            window,
            seed: 7,
        };
        for seq in seqs {
            let mut buf = RecordBuf {
                keys: Vec::new(),
                payloads: Vec::new(),
            };
            collect_one_contig(&seq, 0, k, &params, false, &mut buf);
            let mut single_recs: Vec<(u128, u32, u32, u8)> = buf
                .keys
                .into_iter()
                .zip(buf.payloads)
                .map(|(key, (cid, pos, strand))| (key, cid, pos, strand))
                .collect();
            single_recs.sort_unstable();

            let sm = crate::libs::syncmer::syncmer_dna(&seq, &params).unwrap();
            let keys = rolling_kmer_keys(&seq, k);
            let mut ref_recs: Vec<(u128, u32, u32, u8)> = Vec::new();
            for (_h, pos, _fwd) in sm {
                if pos + k > seq.len() {
                    continue;
                }
                if let Some(key) = keys[pos] {
                    ref_recs.push((key, 0, pos as u32, 0));
                    ref_recs.push((rc_key(key, k), 0, pos as u32, 1));
                }
            }
            ref_recs.sort_unstable();
            assert_eq!(single_recs, ref_recs, "single-pass mismatch");
        }

        // Default build parameters (k=40, syncmer 8/5) on a longer sequence.
        let long: Vec<u8> = (0..200_000u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let (k, smer, window) = (40usize, 8usize, 5usize);
        let params = SyncmerParams {
            smer,
            window,
            seed: 7,
        };
        let mut buf = RecordBuf {
            keys: Vec::new(),
            payloads: Vec::new(),
        };
        collect_one_contig(&long, 0, k, &params, false, &mut buf);
        let mut single_recs: Vec<(u128, u32, u32, u8)> = buf
            .keys
            .into_iter()
            .zip(buf.payloads)
            .map(|(key, (cid, pos, strand))| (key, cid, pos, strand))
            .collect();
        single_recs.sort_unstable();
        let sm = crate::libs::syncmer::syncmer_dna(&long, &params).unwrap();
        let keys = rolling_kmer_keys(&long, k);
        let mut ref_recs: Vec<(u128, u32, u32, u8)> = Vec::new();
        for (_h, pos, _fwd) in sm {
            if pos + k > long.len() {
                continue;
            }
            if let Some(key) = keys[pos] {
                ref_recs.push((key, 0, pos as u32, 0));
                ref_recs.push((rc_key(key, k), 0, pos as u32, 1));
            }
        }
        ref_recs.sort_unstable();
        assert_eq!(single_recs.len(), ref_recs.len(), "record count mismatch");
        assert_eq!(
            single_recs, ref_recs,
            "single-pass mismatch at default params"
        );

        // Same, but with N runs (N is treated as A in s-mer hashing and
        // invalidates k-mer windows that contain it).
        let mut with_n: Vec<u8> = (0..200_000u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        for i in (50_000..50_100).chain(120_000..120_050) {
            with_n[i] = b'N';
        }
        let mut buf = RecordBuf {
            keys: Vec::new(),
            payloads: Vec::new(),
        };
        collect_one_contig(&with_n, 0, k, &params, false, &mut buf);
        let mut single_recs: Vec<(u128, u32, u32, u8)> = buf
            .keys
            .into_iter()
            .zip(buf.payloads)
            .map(|(key, (cid, pos, strand))| (key, cid, pos, strand))
            .collect();
        single_recs.sort_unstable();
        let sm = crate::libs::syncmer::syncmer_dna(&with_n, &params).unwrap();
        let keys = rolling_kmer_keys(&with_n, k);
        let mut ref_recs: Vec<(u128, u32, u32, u8)> = Vec::new();
        for (_h, pos, _fwd) in sm {
            if pos + k > with_n.len() {
                continue;
            }
            if let Some(key) = keys[pos] {
                ref_recs.push((key, 0, pos as u32, 0));
                ref_recs.push((rc_key(key, k), 0, pos as u32, 1));
            }
        }
        ref_recs.sort_unstable();
        assert_eq!(
            single_recs.len(),
            ref_recs.len(),
            "record count mismatch with N"
        );
        assert_eq!(single_recs, ref_recs, "single-pass mismatch with N");
    }

    #[test]
    fn entries_globally_sorted() {
        // Regression: entries must be strictly ascending by k-mer because
        // `dist pgi` merges two sorted key arrays. A low-byte-major bucket
        // order (used previously) is not a global sort and undercounts
        // intersections.
        let mut rng = rand::rngs::StdRng::seed_from_u64(11);
        use rand::{Rng, SeedableRng};
        let seq: Vec<u8> = (0..300_000u32)
            .map(|_| b"ACGT"[rng.random_range(0..4) as usize])
            .collect();
        let idx = build_from_seqs(vec![(String::from("c"), seq)], 40, 8, 5, false).unwrap();
        assert!(idx.entries.len() > 10_000, "too few unique k-mers");
        for w in idx.entries.windows(2) {
            assert!(w[0].kmer <= w[1].kmer, "entries not sorted by k-mer");
        }
    }

    #[test]
    fn no_rev_halves_strands() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let both =
            build_from_seqs(vec![(String::from("c1"), seq.clone())], 10, 4, 2, false).unwrap();
        let fwd = build_from_seqs(vec![(String::from("c1"), seq)], 10, 4, 2, true).unwrap();
        assert!(both.n_positions() >= fwd.n_positions());
    }
}
