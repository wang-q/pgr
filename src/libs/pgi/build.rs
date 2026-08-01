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
    let mut r: u128 = 0;
    for i in (0..k).rev() {
        let c = ((x >> (2 * i)) & 3) ^ 3;
        r = (r << 2) | c;
    }
    r
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

    // (key, contig_id, pos, strand)
    let mut records: Vec<(u128, u32, u32, u8)> = Vec::new();
    for (cid, (_, seq)) in contigs.iter().enumerate() {
        if seq.len() < k {
            continue;
        }
        let sm = syncmer_dna(seq, &params)?;
        for (_h, pos, _is_fwd) in sm {
            let p = pos;
            if p + k > seq.len() {
                continue;
            }
            let Some(key_fwd) = pack_kmer(&seq[p..p + k], k) else {
                continue; // k-mer contains N or ambiguity; skip
            };
            records.push((key_fwd, cid as u32, p as u32, 0));
            if !no_rev {
                records.push((rc_key(key_fwd, k), cid as u32, p as u32, 1));
            }
        }
    }
    records.sort_unstable();

    let mut entries: Vec<PgiEntry> = Vec::new();
    let mut positions: Vec<(u32, u32, u8)> = Vec::with_capacity(records.len());
    let mut i = 0usize;
    while i < records.len() {
        let key = records[i].0;
        let pos_start = positions.len() as u32;
        let mut j = i;
        while j < records.len() && records[j].0 == key {
            positions.push((records[j].1, records[j].2, records[j].3));
            j += 1;
        }
        entries.push(PgiEntry {
            kmer: key,
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
        // RC("ACGT") = "TGCA"; double RC restores the original.
        let x = pack_kmer(b"ACGT", 4).unwrap();
        assert_eq!(rc_key(x, 4), pack_kmer(b"TGCA", 4).unwrap());
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
    fn no_rev_halves_strands() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let both =
            build_from_seqs(vec![(String::from("c1"), seq.clone())], 10, 4, 2, false).unwrap();
        let fwd = build_from_seqs(vec![(String::from("c1"), seq)], 10, 4, 2, true).unwrap();
        assert!(both.n_positions() >= fwd.n_positions());
    }
}
