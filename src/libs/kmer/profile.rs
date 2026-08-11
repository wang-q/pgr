//! Per-sequence k-mer profiles (FastK `-p` / `-p:<table>` equivalents).

use super::KmerTable;
use crate::libs::ds::radix_sort::radix_sort_bytes_par;
use anyhow::Context;
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;

/// FastK caps per-k-mer counts at 32767 (`0x7fff`); profile values are `u16`.
const PROFILE_CAP: u16 = 0x7fff;

/// Profile value per k-mer position from the dataset-wide count (self).
///
/// `profiles[i][p]` is the count of the canonical k-mer at position `p` of
/// sequence `i`, or 0 when the window contains N (FastK splits on gaps).
pub fn self_profiles(seqs: &[Vec<u8>], k: usize, table: &KmerTable) -> Vec<Vec<u16>> {
    table_profiles(seqs, k, table)
}

/// Profile value per k-mer position from the repeat-table count (relative).
///
/// Values are the table count of the k-mer, or 0 when the window contains N
/// or the k-mer is absent from the table (FastK `-p:<table>` semantics).
pub fn relative_profiles(seqs: &[Vec<u8>], k: usize, table: &KmerTable) -> Vec<Vec<u16>> {
    table_profiles(seqs, k, table)
}

/// Shared scan: each valid window's canonical key is looked up in `table` and
/// its count (capped to the FastK 32767 limit) becomes the profile value;
/// N-containing windows and missing keys profile as 0.
///
/// Implementation follows FastK's merge strategy instead of per-window
/// table search: collect all window keys, sort them (parallel radix), and
/// merge against the sorted table once. Per-window random lookups on the
/// multi-MB key table are latency-bound regardless of search structure
/// (measured 2026-08-09, see notes/benchmarks/bench-profile-hotspots.md).
fn table_profiles(seqs: &[Vec<u8>], k: usize, table: &KmerTable) -> Vec<Vec<u16>> {
    if seqs.is_empty() || table.keys.is_empty() {
        return seqs
            .iter()
            .map(|seq| vec![0u16; seq.len().saturating_sub(k - 1)])
            .collect();
    }
    // Collect every valid window as (key, location); per-sequence vectors so
    // the collection itself can run in parallel, then flatten for sorting.
    let key_bytes = k.div_ceil(4);
    let per_seq: Vec<(Vec<u8>, Vec<Loc>)> = seqs
        .par_iter()
        .enumerate()
        .map(|(si, seq)| {
            let mut keys = Vec::with_capacity(seq.len().saturating_sub(k - 1) * key_bytes);
            let mut locs = Vec::with_capacity(seq.len().saturating_sub(k - 1));
            super::canonical_keys(seq, k, |p, key| {
                keys.extend_from_slice(key.to_bytes());
                locs.push(Loc::new(si, p));
            });
            (keys, locs)
        })
        .collect();
    let n_bytes: usize = per_seq.iter().map(|(k, _)| k.len()).sum();
    let n_windows = n_bytes / key_bytes;
    let mut keys: Vec<u8> = Vec::with_capacity(n_bytes);
    let mut locs: Vec<Loc> = Vec::with_capacity(n_windows);
    for (k, l) in per_seq {
        keys.extend_from_slice(&k);
        locs.extend(l);
    }
    radix_sort_bytes_par(&mut keys, key_bytes, &mut locs);
    // Merge sorted windows against the (sorted, deduplicated) table: equal
    // keys receive the table count; everything else stays 0 (N windows were
    // never collected).
    let mut out: Vec<Vec<u16>> = seqs
        .iter()
        .map(|seq| vec![0u16; seq.len().saturating_sub(k - 1)])
        .collect();
    let mut i = 0usize;
    let mut j = 0usize;
    let table_n = table.counts.len();
    while i < n_windows && j < table_n {
        let wk = &keys[i * key_bytes..(i + 1) * key_bytes];
        let tk = &table.keys[j * key_bytes..(j + 1) * key_bytes];
        match wk.cmp(tk) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                let count = table.counts[j].min(PROFILE_CAP as u32) as u16;
                while i < n_windows && keys[i * key_bytes..(i + 1) * key_bytes] == *tk {
                    let (si, pos) = locs[i].split();
                    out[si][pos] = count;
                    i += 1;
                }
                j += 1;
            }
        }
    }
    out
}

/// Location of a collected window: sequence index and window start.
#[derive(Clone, Copy)]
struct Loc {
    seq: u32,
    pos: u64,
}

impl Loc {
    fn new(seq: usize, pos: usize) -> Self {
        Self {
            seq: seq as u32,
            pos: pos as u64,
        }
    }

    fn split(self) -> (usize, usize) {
        (self.seq as usize, self.pos as usize)
    }
}

/// File magic for the `.pkp` per-sequence profile file.
const PKP_MAGIC: &[u8; 4] = b"PKPP";
/// Format version.
const PKP_VERSION: u32 = 1;
/// Bincode-free fixed header size: magic + version + k + n_seqs.
const PKP_HEADER_LEN: usize = 20;

/// Write per-sequence profiles to `path` in the `.pkp` format.
///
/// Layout: header (`magic/version/k/n_seqs`) plus one `u64` length and raw
/// little-endian `u16` values per sequence. The file is self-contained (no
/// external k or contig metadata needed to read it back).
pub fn save_profiles(path: &Path, k: usize, profiles: &[Vec<u16>]) -> anyhow::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    w.write_all(PKP_MAGIC).context("writing pkp header")?;
    w.write_all(&PKP_VERSION.to_le_bytes())
        .context("writing pkp header")?;
    w.write_all(&(k as u32).to_le_bytes())
        .context("writing pkp header")?;
    w.write_all(&(profiles.len() as u64).to_le_bytes())
        .context("writing pkp header")?;
    for p in profiles {
        w.write_all(&(p.len() as u64).to_le_bytes())
            .context("writing pkp length")?;
        for &v in p {
            w.write_all(&v.to_le_bytes())
                .context("writing pkp values")?;
        }
    }
    w.flush().context("flushing pkp file")?;
    Ok(())
}

/// Read a `.pkp` file written by [`save_profiles`], validating magic/version
/// and that the stored `k` matches the requested one.
pub fn load_profiles(path: &Path, k: usize) -> anyhow::Result<Vec<Vec<u16>>> {
    let bytes = std::fs::read(path).context("reading pkp file")?;
    anyhow::ensure!(bytes.len() >= PKP_HEADER_LEN, "truncated pkp header");
    anyhow::ensure!(&bytes[0..4] == PKP_MAGIC, "not a pgr profile (bad magic)");
    let version = u32::from_le_bytes(bytes[4..8].try_into().unwrap());
    anyhow::ensure!(version == PKP_VERSION, "unsupported pkp version {version}");
    let stored_k = u32::from_le_bytes(bytes[8..12].try_into().unwrap()) as usize;
    anyhow::ensure!(
        stored_k == k,
        "pkp k mismatch: file has k={stored_k}, requested {k}"
    );
    let n_seqs = u64::from_le_bytes(bytes[12..20].try_into().unwrap());
    let mut off = PKP_HEADER_LEN;
    let mut profiles = Vec::with_capacity(n_seqs as usize);
    for _ in 0..n_seqs {
        let len_bytes = bytes.get(off..off + 8).context("truncated pkp length")?;
        let len = u64::from_le_bytes(len_bytes.try_into().unwrap());
        off += 8;
        let n_bytes = (len as usize)
            .checked_mul(2)
            .context("pkp length overflow")?;
        let body = bytes
            .get(off..off + n_bytes)
            .context("truncated pkp values")?;
        let mut p = Vec::with_capacity(len as usize);
        for chunk in body.chunks_exact(2) {
            p.push(u16::from_le_bytes(chunk.try_into().unwrap()));
        }
        profiles.push(p);
        off += n_bytes;
    }
    anyhow::ensure!(off == bytes.len(), "trailing bytes in pkp file");
    Ok(profiles)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::kmer::count::build_table;
    use rand::{rngs::StdRng, Rng, SeedableRng};

    /// Deterministic pseudo-random DNA block (same LCG as pgi tests).
    fn random_block(len: usize, seed: u64) -> Vec<u8> {
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

    /// A random block whose k-mer windows are all canonical-unique.
    fn unique_block(k: usize, seed0: u64) -> Vec<u8> {
        (0..100u64)
            .map(|i| random_block(80, seed0 + i))
            .find(|b| {
                build_table(std::slice::from_ref(b), k)
                    .unwrap()
                    .counts
                    .len()
                    == b.len() - k + 1
            })
            .expect("a collision-free block must exist")
    }

    #[test]
    fn self_profile_counts_across_sequences() {
        // Two identical copies as separate sequences: every block k-mer
        // counts 2 in the combined table, so both profiles are all 2.
        let block = unique_block(6, 42);
        let seqs = vec![block.clone(), block.clone()];
        let table = build_table(&seqs, 6).unwrap();
        let profiles = self_profiles(&seqs, 6, &table);
        assert_eq!(profiles.len(), 2);
        assert!(profiles[0].iter().all(|&v| v == 2));
        assert!(profiles[1].iter().all(|&v| v == 2));
    }

    #[test]
    fn profile_zero_at_n_windows() {
        let seq = b"ACGTACGTNNACGTACGT".to_vec();
        let table = build_table(std::slice::from_ref(&seq), 4).unwrap();
        let profiles = self_profiles(std::slice::from_ref(&seq), 4, &table);
        // 15 windows; the 5 touching N (starts 5..9) must be 0, others > 0.
        assert_eq!(profiles[0].len(), 15);
        for (p, &v) in profiles[0].iter().enumerate() {
            assert_eq!(v == 0, (5..10).contains(&p), "window {p} value {v}");
        }
    }

    #[test]
    fn relative_profile_uses_table_count() {
        // Lib holds one copy; genome has two copies. The relative profile
        // must report the table count (1), not the genomic count (2), and 0
        // for windows overlapping the N gap (absent from the table).
        let lib = unique_block(6, 7);
        let k = 6;
        let table = build_table(std::slice::from_ref(&lib), k).unwrap();
        let mut genome = lib.clone();
        genome.extend(std::iter::repeat_n(b'N', k));
        genome.extend_from_slice(&lib);
        let profiles = relative_profiles(&[genome], k, &table);
        let first = lib.len() - k + 1; // windows fully inside the first copy
        let zeros = 2 * k - 1; // windows overlapping the N run
        assert!(profiles[0][..first].iter().all(|&v| v == 1));
        assert!(profiles[0][first..first + zeros].iter().all(|&v| v == 0));
        assert!(profiles[0][first + zeros..].iter().all(|&v| v == 1));

        // A sequence sharing no canonical k-mer with the table profiles as 0.
        let other = (0..100u64)
            .map(|i| random_block(22, 900 + i))
            .find(|b| {
                relative_profiles(std::slice::from_ref(b), k, &table)[0]
                    .iter()
                    .all(|&v| v == 0)
            })
            .expect("an absent block must exist");
        assert_eq!(
            relative_profiles(&[other], k, &table)[0]
                .iter()
                .sum::<u16>(),
            0
        );
    }

    #[test]
    fn sort_merge_matches_binary_search() {
        // Random multi-sequence inputs: the sort+merge implementation must
        // agree with the previous per-window partition_point lookup.
        let mut rng = StdRng::seed_from_u64(20260809);
        for k in [1usize, 4, 8, 17, 40] {
            for _ in 0..10 {
                let seqs: Vec<Vec<u8>> = (0..rng.random_range(1..8))
                    .map(|_| random_block(rng.random_range(20..200), rng.random()))
                    .collect();
                let table = build_table(&seqs, k).unwrap();
                let got = table_profiles(&seqs, k, &table);
                let expected: Vec<Vec<u16>> = seqs
                    .iter()
                    .map(|seq| {
                        let mut out = vec![0u16; seq.len().saturating_sub(k - 1)];
                        crate::libs::kmer::canonical_keys(seq, k, |p, key| {
                            let kb = table.key_bytes();
                            let mut lo = 0usize;
                            let mut hi = table.counts.len();
                            while lo < hi {
                                let mid = (lo + hi) / 2;
                                if &table.keys[mid * kb..(mid + 1) * kb] < key.to_bytes() {
                                    lo = mid + 1;
                                } else {
                                    hi = mid;
                                }
                            }
                            if lo < table.counts.len()
                                && &table.keys[lo * kb..(lo + 1) * kb] == key.to_bytes()
                            {
                                out[p] = table.counts[lo].min(PROFILE_CAP as u32) as u16;
                            }
                        });
                        out
                    })
                    .collect();
                assert_eq!(got, expected, "k={k}");
            }
        }
    }

    #[test]
    fn save_load_profiles_roundtrip() {
        let profiles = vec![
            vec![1, 2, 3, 0, 32767],
            vec![5, 5, 5],
            vec![],
            vec![0u16; 1000],
        ];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.pkp");
        save_profiles(&path, 17, &profiles).unwrap();
        let loaded = load_profiles(&path, 17).unwrap();
        assert_eq!(loaded, profiles);
    }

    #[test]
    fn pkp_layout_and_rejections() {
        let profiles = vec![vec![7u16, 8], vec![9u16]];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.pkp");
        save_profiles(&path, 9, &profiles).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(&bytes[0..4], b"PKPP");
        assert_eq!(&bytes[4..8], &1u32.to_le_bytes());
        assert_eq!(&bytes[8..12], &9u32.to_le_bytes());
        assert_eq!(&bytes[12..20], &2u64.to_le_bytes()); // n_seqs
                                                         // seq 0: len 2 + [7, 8]
        assert_eq!(&bytes[20..28], &2u64.to_le_bytes());
        assert_eq!(&bytes[28..32], &[7, 0, 8, 0]);

        // k mismatch and bad magic must be rejected.
        assert!(load_profiles(&path, 10).is_err());
        let bad = dir.path().join("bad.pkp");
        std::fs::write(&bad, b"XXXX").unwrap();
        assert!(load_profiles(&bad, 9).is_err());
        // Truncated body must be rejected.
        std::fs::write(&path, &bytes[..bytes.len() - 3]).unwrap();
        assert!(load_profiles(&path, 9).is_err());
    }
}
