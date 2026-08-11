//! Build a canonical k-mer count table and persist it as `.pkt`.

use super::KmerTable;
use anyhow::Context;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::path::{Path, PathBuf};

/// File magic for the `.pkt` k-mer table cache.
const PKT_MAGIC: &[u8; 4] = b"PKTT";
/// Format version.
const PKT_VERSION: u32 = 1;

/// Fixed-size `.pkt` header, serialized with bincode (the entry payload is
/// written separately so raw `Vec<u128>` keys are never bincode-encoded).
#[derive(Serialize, Deserialize)]
struct PktHeader {
    magic: [u8; 4],
    version: u32,
    k: u32,
    n_entries: u64,
    key_bytes: u32,
}

/// Bincode byte size of the fixed-size [`PktHeader`] (4+4+4+8+4 LE fields).
const PKT_HEADER_LEN: usize = 24;

/// Build a canonical k-mer count table from `seqs`.
///
/// Every N-free k-mer window contributes its canonical key; counts accumulate
/// across all sequences (FastK `-p` semantics). `k` must be `<= MAX_K`.
pub fn build_table(seqs: &[Vec<u8>], k: usize) -> anyhow::Result<KmerTable> {
    anyhow::ensure!(
        k > 0 && k <= crate::libs::kmer::key::Kmer::MAX_K,
        "k must be in 1..={}, got {k}",
        crate::libs::kmer::key::Kmer::MAX_K
    );
    let per_seq: Vec<Vec<crate::libs::kmer::key::Kmer>> = seqs
        .par_iter()
        .map(|seq| {
            let mut keys = Vec::new();
            super::canonical_keys(seq, k, |_, key| keys.push(key));
            keys
        })
        .collect();
    let n: usize = per_seq.iter().map(Vec::len).sum();
    let mut keys: Vec<u8> = Vec::with_capacity(n * k.div_ceil(4));
    for v in &per_seq {
        for km in v {
            keys.extend_from_slice(km.to_bytes());
        }
    }
    Ok(count_keys(keys, k))
}

/// Sorts a packed raw canonical key list (with duplicates) into a count
/// table.
///
/// The deduplication tail shared by [`build_table`] and the memory-bounded
/// bucket path of `fq norm`. `k` must already be validated.
pub(crate) fn count_keys(mut keys: Vec<u8>, k: usize) -> KmerTable {
    let key_bytes = k.div_ceil(4);
    if keys.is_empty() {
        return KmerTable {
            k,
            keys,
            counts: Vec::new(),
        };
    }
    let n_keys = keys.len() / key_bytes;
    crate::libs::ds::radix_sort::radix_sort_bytes_par(&mut keys, key_bytes, &mut vec![(); n_keys]);
    // Group equal keys (now contiguous after the sort) into counts.
    let mut counts: Vec<u32> = Vec::with_capacity(n_keys);
    let mut i = 0usize;
    let mut w = 0usize;
    while i < n_keys {
        let mut j = i + 1;
        while j < n_keys
            && keys[j * key_bytes..(j + 1) * key_bytes] == keys[i * key_bytes..(i + 1) * key_bytes]
        {
            j += 1;
        }
        if w != i {
            keys.copy_within(i * key_bytes..(i + 1) * key_bytes, w * key_bytes);
        }
        counts.push((j - i).min(u32::MAX as usize) as u32);
        w += 1;
        i = j;
    }
    keys.truncate(w * key_bytes);
    KmerTable { k, keys, counts }
}

/// Write `table` to `path` (`.pkt`) atomically: header (bincode) plus one
/// packed key of `ceil(2k/8)` bytes and a `u32` count per entry.
pub fn save(table: &KmerTable, path: &Path) -> anyhow::Result<()> {
    let key_bytes = table.key_bytes();
    let header = PktHeader {
        magic: *PKT_MAGIC,
        version: PKT_VERSION,
        k: table.k as u32,
        n_entries: table.counts.len() as u64,
        key_bytes: key_bytes as u32,
    };
    let mut buf = bincode::serialize(&header).context("serializing pkt header")?;
    buf.reserve(table.keys.len() + table.counts.len() * 4);
    for (i, &count) in table.counts.iter().enumerate() {
        buf.extend_from_slice(&table.keys[i * key_bytes..(i + 1) * key_bytes]);
        buf.extend_from_slice(&count.to_le_bytes());
    }
    let mut name = path.as_os_str().to_os_string();
    name.push(format!(".tmp.{}", std::process::id()));
    let tmp = PathBuf::from(name);
    {
        let mut w = std::io::BufWriter::new(std::fs::File::create(&tmp)?);
        w.write_all(&buf)?;
        w.flush()?;
    }
    std::fs::rename(&tmp, path)?;
    Ok(())
}

/// K-mer length stored in a `.pkt` file header, without loading the table.
pub fn k_of(path: &Path) -> anyhow::Result<usize> {
    let bytes = std::fs::read(path)?;
    let header_bytes = bytes
        .get(..PKT_HEADER_LEN)
        .context("truncated pkt header")?;
    let header: PktHeader = bincode::deserialize(header_bytes).context("bad pkt header")?;
    if &header.magic != PKT_MAGIC {
        anyhow::bail!("not a pgr k-mer table (bad magic)");
    }
    if header.version != PKT_VERSION {
        anyhow::bail!(
            "unsupported pkt version {} (expected {PKT_VERSION})",
            header.version
        );
    }
    Ok(header.k as usize)
}

/// Read a `.pkt` table written by [`save`], validating magic/version/length
/// and that the stored `k` matches the requested one.
pub fn load(path: &Path, k: usize) -> anyhow::Result<KmerTable> {
    let bytes = std::fs::read(path)?;
    let header_bytes = bytes
        .get(..PKT_HEADER_LEN)
        .context("truncated pkt header")?;
    let header: PktHeader = bincode::deserialize(header_bytes).context("bad pkt header")?;
    if &header.magic != PKT_MAGIC {
        anyhow::bail!("not a pgr k-mer table (bad magic)");
    }
    if header.version != PKT_VERSION {
        anyhow::bail!(
            "unsupported pkt version {} (expected {PKT_VERSION})",
            header.version
        );
    }
    let stored_k = header.k as usize;
    anyhow::ensure!(
        stored_k == k,
        "repeat table k={stored_k} conflicts with -k {k} (rebuild)"
    );
    let key_bytes = header.key_bytes as usize;
    anyhow::ensure!(
        key_bytes == (2 * k).div_ceil(8),
        "repeat table key size {key_bytes} does not match k={k}"
    );
    let n_entries = header.n_entries as usize;
    let entry_len = key_bytes
        .checked_add(4)
        .and_then(|e| e.checked_mul(n_entries))
        .context("pkt entry count overflow")?;
    anyhow::ensure!(
        bytes.len() == PKT_HEADER_LEN + entry_len,
        "truncated pkt table ({} bytes, expected {})",
        bytes.len(),
        PKT_HEADER_LEN + entry_len
    );

    let mut keys = Vec::with_capacity(n_entries * key_bytes);
    let mut counts = Vec::with_capacity(n_entries);
    let mut off = PKT_HEADER_LEN;
    for _ in 0..n_entries {
        keys.extend_from_slice(
            bytes
                .get(off..off + key_bytes)
                .context("truncated pkt entry")?,
        );
        let count_bytes: [u8; 4] = bytes
            .get(off + key_bytes..off + key_bytes + 4)
            .context("truncated pkt count")?
            .try_into()
            .unwrap();
        counts.push(u32::from_le_bytes(count_bytes));
        off += key_bytes + 4;
    }
    Ok(KmerTable { k, keys, counts })
}

#[cfg(test)]
mod tests {
    use super::*;

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

    /// A random block whose k-mer windows are all canonical-unique (needed
    /// when a test asserts exact per-window counts on a duplicated block).
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
    fn build_counts_duplicates() {
        // Each canonical k-mer of the duplicated block appears twice across
        // the two copies; k-mers straddling N gaps are skipped entirely.
        let block = unique_block(6, 42);
        let mut seq = block.clone();
        seq.extend_from_slice(&block);
        let table = build_table(&[seq], 6).unwrap();
        assert_eq!(table.k, 6);
        assert!(!table.keys.is_empty());
        // The first window's canonical key appears once per copy.
        let canonical = crate::libs::kmer::key::Kmer::from_bases(&block[..6], 6)
            .unwrap()
            .canonical();
        let kb = table.key_bytes();
        let mut lo = 0usize;
        let mut hi = table.counts.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if &table.keys[mid * kb..(mid + 1) * kb] < canonical.to_bytes() {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        assert_eq!(&table.keys[lo * kb..(lo + 1) * kb], canonical.to_bytes());
        assert_eq!(table.counts[lo], 2, "duplicated block k-mer must count 2");
        // Total windows = duplicated sequence length - k + 1, all valid.
        let total: u32 = table.counts.iter().sum();
        assert_eq!(total as usize, 2 * block.len() - 6 + 1);
    }

    #[test]
    fn build_counts_hand_checked() {
        // Small sequence with a known canonical count:
        // "AAAA" (k=2) has canonical keys AA (1x), and reverse pairs are
        // merged (e.g. AC == GT canonical).
        let table = build_table(&[b"AAAA".to_vec()], 2).unwrap();
        let mut counts = std::collections::BTreeMap::new();
        for (i, &count) in table.counts.iter().enumerate() {
            counts.insert(table.key_at(i).to_bytes().to_vec(), count);
        }
        let aa = crate::libs::kmer::key::Kmer::from_bases(b"AA", 2).unwrap();
        assert_eq!(
            counts.get(aa.to_bytes()),
            Some(&3),
            "AA appears at 3 positions"
        );
        assert_eq!(table.counts.len(), 1, "one unique canonical k-mer");

        // Case-insensitive: lowercase input merges into the same keys.
        let lower = build_table(&[b"aaaa".to_vec()], 2).unwrap();
        assert_eq!(table.keys, lower.keys);
        assert_eq!(table.counts, lower.counts);
    }

    #[test]
    fn n_runs_split_but_no_keys_inside() {
        // A window touching N contributes nothing; flanks still count.
        let table = build_table(&[b"ACGTACGTNNACGTACGT".to_vec()], 4).unwrap();
        let total: u32 = table.counts.iter().sum();
        // 18 bases, 15 windows; the 5 windows covering N positions 8..9
        // (starts 5..9) are invalid, so 15 - 5 = 10 valid windows.
        assert_eq!(total, 10);
    }

    #[test]
    fn runs_shorter_than_k_emit_nothing() {
        // Regression: an N-free run shorter than k must neither emit nor
        // panic (the old window scan could cross the N gap).
        let table = build_table(&[b"ACGTNNNNNNNNACGTACGTACGTACGT".to_vec()], 8).unwrap();
        let total: u32 = table.counts.iter().sum();
        // Left run "ACGT" (4 < 8): nothing; right run 16 bases: 9 windows.
        assert_eq!(total, 9);
    }

    #[test]
    fn save_load_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("lib.pkt");
        let seqs = vec![
            b"ACGTACGTACGTACGTACGTNNNACGTACGTACGTACGT".to_vec(),
            b"TTTTTGGGGGCCCCCAAAAA".to_vec(),
        ];
        let table = build_table(&seqs, 9).unwrap();
        save(&table, &path).unwrap();
        let loaded = load(&path, 9).unwrap();
        assert_eq!(loaded.k, table.k);
        assert_eq!(loaded.keys, table.keys);
        assert_eq!(loaded.counts, table.counts);
    }

    #[test]
    fn load_rejects_truncated_and_wrong_k() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("lib.pkt");
        let table = build_table(&[b"ACGTACGTACGTACGTACGT".to_vec()], 8).unwrap();
        save(&table, &path).unwrap();

        let full = std::fs::read(&path).unwrap();
        std::fs::write(&path, &full[..full.len() - 5]).unwrap();
        assert!(load(&path, 8).is_err(), "truncated pkt must be rejected");
        std::fs::write(&path, &full).unwrap();
        assert!(
            load(&path, 10).is_err(),
            "k mismatch must be rejected as stale"
        );

        let bad = dir.path().join("bad.pkt");
        std::fs::write(&bad, b"XXXX").unwrap();
        assert!(load(&bad, 8).is_err(), "bad magic must be rejected");
    }
}
