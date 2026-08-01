//! Build a `.pgi` index from FASTA or 2bit sequences.

use super::{PgiEntry, PgiIndex};
use crate::libs::syncmer::{syncmer_dna, SyncmerParams};
use anyhow::Context;

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
    // x's lowest 2-bit group is the last base; the RC key's highest group is
    // the complement of that base, so iterate the groups low -> high.
    let mut r: u128 = 0;
    for i in 0..k {
        let c = ((x >> (2 * i)) & 3) ^ 3;
        r = (r << 2) | c;
    }
    r
}

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

/// Build an index from named sequences.
///
/// Each syncmer position seeds a k-mer (both strands unless `no_rev`); the
/// resulting records are sorted by key and grouped into unique entries.
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

    // Collect raw (key, contig, pos, strand) records into one growable vector
    // (a single vector reallocates far less than 256 per-key buckets).
    // Capacity estimate: ~2 syncmer-sampled positions per window of
    // `window+1` bases, doubled for both strands.
    let est: usize = contigs
        .iter()
        .map(|(_, s)| (s.len() / (window + 1)).saturating_mul(2) * 2)
        .sum();
    let mut records: Vec<(u128, u32, u32, u8)> = Vec::with_capacity(est);
    for (cid, (_, seq)) in contigs.iter().enumerate() {
        if seq.len() < k {
            continue;
        }
        let sm = syncmer_dna(seq, &params)?;
        let keys = rolling_kmer_keys(seq, k);
        for (_h, pos, _is_fwd) in sm {
            let p = pos;
            if p + k > seq.len() {
                continue;
            }
            let Some(key_fwd) = keys[p] else {
                continue; // k-mer contains N or ambiguity; skip
            };
            records.push((key_fwd, cid as u32, p as u32, 0));
            if !no_rev {
                let rev = rc_key(key_fwd, k);
                records.push((rev, cid as u32, p as u32, 1));
            }
        }
    }

    // Counting-sort bucketing by the key's lowest byte; bucket-internal sort
    // keeps the merged order equal to a full sort.
    const NBUCKETS: usize = 256;
    let mut counts = [0usize; NBUCKETS];
    for r in &records {
        counts[(r.0 & 0xff) as usize] += 1;
    }
    let mut offsets = [0usize; NBUCKETS];
    let mut cum = 0usize;
    for b in 0..NBUCKETS {
        offsets[b] = cum;
        cum += counts[b];
    }
    let mut bucketed = vec![(0u128, 0u32, 0u32, 0u8); records.len()];
    let mut next = offsets;
    for r in records {
        let b = (r.0 & 0xff) as usize;
        bucketed[next[b]] = r;
        next[b] += 1;
    }
    let mut entries: Vec<PgiEntry> = Vec::new();
    let mut positions: Vec<(u32, u32, u8)> = Vec::with_capacity(bucketed.len());
    for b in 0..NBUCKETS {
        let (s, e) = (offsets[b], offsets[b] + counts[b]);
        bucketed[s..e].sort_unstable();
        let mut i = 0usize;
        while i < e - s {
            let key = bucketed[s + i].0;
            let pos_start = positions.len() as u32;
            let mut j = i;
            while j < e - s && bucketed[s + j].0 == key {
                positions.push((bucketed[s + j].1, bucketed[s + j].2, bucketed[s + j].3));
                j += 1;
            }
            entries.push(PgiEntry {
                kmer: key,
                pos_start,
                freq: (j - i) as u32,
            });
            i = j;
        }
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
    fn no_rev_halves_strands() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let both =
            build_from_seqs(vec![(String::from("c1"), seq.clone())], 10, 4, 2, false).unwrap();
        let fwd = build_from_seqs(vec![(String::from("c1"), seq)], 10, 4, 2, true).unwrap();
        assert!(both.n_positions() >= fwd.n_positions());
    }
}
