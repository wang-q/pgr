//! Build a `.pgi` index from FASTA or 2bit sequences.

use super::{pack_position, PgiEntry, PgiIndex};
use crate::libs::nt::rc_key;
use crate::libs::syncmer::SyncmerParams;
use anyhow::Context;
use std::collections::{HashSet, VecDeque};

/// Parallel growable buffers for (key, contig, pos, strand) records.
struct RecordBuf {
    keys: Vec<u128>,
    payloads: Vec<u64>,
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
    let factor = crate::libs::syncmer::hash_factor(7);
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
    // Positions already queued. The closed-syncmer rule can select the same
    // position twice (it is the minimal s-mer at the last position of one
    // window and at the first position of the next); dedupe so the index
    // holds exactly one record per k-mer position.
    let mut queued: HashSet<usize> = HashSet::new();

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
            queued.remove(&pos);
            // The rolling key is `seq[i-k+1..=i]`, which equals the k-mer at
            // `pos` only when this iteration completes exactly that window
            // (`pos + k - 1 == i`). Positions selected twice (already
            // emitted) or queued out of order (a window-start minimum is
            // discovered after its k-mer window completed) pop later; their
            // key must be recomputed from the sequence instead of using the
            // already-rolled key.
            let key = if pos + k - 1 == i {
                (kvalid >= k && i + 1 >= k).then_some(kx)
            } else {
                kmer_key_at(seq, pos, k, &codes)
            };
            if let Some(key) = key {
                out.keys.push(key);
                out.payloads.push(pack_position(cid, pos as u32, 0));
                if !no_rev {
                    out.keys.push(rc_key(key, k));
                    out.payloads.push(pack_position(cid, pos as u32, 1));
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
                    if start + k <= n && queued.insert(start) {
                        pending.push_back(start);
                    }
                } else if ch == min_val && j + k <= n && queued.insert(j) {
                    pending.push_back(j);
                }
            }
        }
    }
}

/// 2-bit key of the k-mer starting at `pos`; `None` when it contains an N.
fn kmer_key_at(seq: &[u8], pos: usize, k: usize, codes: &[u64; 256]) -> Option<u128> {
    let mut key = 0u128;
    for &b in &seq[pos..pos + k] {
        let code = codes[b as usize];
        if code == 4 {
            return None;
        }
        key = (key << 2) | code as u128;
    }
    Some(key)
}

/// Build an index from named sequences.
///
/// Each syncmer position seeds a k-mer (both strands unless `no_rev`); the
/// resulting records are sorted by key (in-place MSD radix) and grouped into
/// unique entries, globally ascending by k-mer.
pub fn build_from_seqs(
    mut contigs: Vec<(String, Vec<u8>)>,
    k: usize,
    smer: usize,
    window: usize,
    no_rev: bool,
    mask: bool,
) -> anyhow::Result<PgiIndex> {
    // k must fit `2*k` significant bits in a u128 kmer key; check `k <= 64`
    // directly so an extreme CLI value (e.g. `usize::MAX`) is rejected with a
    // friendly error instead of overflowing `k * 2`.
    anyhow::ensure!(k > 0 && k <= 64, "k must be in 1..=64, got {k}");
    anyhow::ensure!(smer > 0, "smer must be positive");
    anyhow::ensure!(window > 0, "window must be positive");
    if mask {
        // FastGA `-M` semantics: soft-masked (lowercase) bases become N so
        // windows touching them are skipped. The k-mer code is otherwise
        // case-insensitive, so without this a soft-masked copy shares seeds
        // with its uppercase twin but the case-sensitive extension DP fails
        // and the chain falls back to an unscored (all-zero) PSL block.
        harden_soft_mask(&mut contigs);
    }
    // `SeedHit` packs contig ids into `u16`; refuse inputs that would
    // truncate silently.
    anyhow::ensure!(
        contigs.len() <= u16::MAX as usize,
        "too many contigs: {} (max {})",
        contigs.len(),
        u16::MAX
    );
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
    anyhow::ensure!(
        payloads.len() <= u32::MAX as usize,
        "too many k-mer records: {} (max {})",
        payloads.len(),
        u32::MAX
    );

    let mut entries: Vec<PgiEntry> = Vec::with_capacity(keys.len());
    let mut positions: Vec<u64> = Vec::with_capacity(keys.len());
    let mut i = 0usize;
    while i < keys.len() {
        let kmer = keys[i];
        let pos_start = positions.len() as u32;
        let mut j = i;
        // A syncmer position can be selected twice (it is the minimum of two
        // adjacent windows); `collect_one_contig` dedups via `queued` only
        // while the position is still pending, so when the position is
        // flushed before its second selection (small-k parameters) the same
        // (kmer, pos, strand) record is emitted twice. Drop the exact
        // duplicate payloads here so the index holds one record per physical
        // position (a diff frequency would falsely trip the `--freq` filter).
        let mut seen: std::collections::HashSet<u64> =
            std::collections::HashSet::with_capacity((j - i).min(64));
        while j < keys.len() && keys[j] == kmer {
            if seen.insert(payloads[j]) {
                positions.push(payloads[j]);
            }
            j += 1;
        }
        entries.push(PgiEntry {
            kmer,
            pos_start,
            freq: (positions.len() - pos_start as usize) as u32,
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
    let mut reader = crate::libs::fmt::seq::SeqReader::new(path)
        .with_context(|| format!("failed to open FASTA {path}"))?;
    let mut rec = crate::libs::fmt::seq::SeqRecord::new();
    let mut contigs = Vec::new();
    while reader.read_record(&mut rec)? {
        let name = String::from_utf8(rec.name().to_vec()).context("FASTA name utf8")?;
        contigs.push((name, rec.sequence().to_vec()));
    }
    Ok(contigs)
}

/// Read all sequences from a 2bit file; `mask` applies the stored soft-mask
/// blocks (lowercased bases).
pub fn read_2bit(path: &str, mask: bool) -> anyhow::Result<Vec<(String, Vec<u8>)>> {
    let mut tb = crate::libs::fmt::twobit::TwoBitFile::open(path)
        .with_context(|| format!("failed to open 2bit {path}"))?;
    let names = tb.get_sequence_names();
    let mut contigs = Vec::with_capacity(names.len());
    for name in names {
        let seq = tb
            .read_sequence(&name, None, None, !mask)
            .with_context(|| format!("reading {name} from 2bit"))?
            .into_bytes();
        contigs.push((name, seq));
    }
    Ok(contigs)
}

/// Replace soft-masked (lowercase) bases with N so windows touching them are
/// skipped by the syncmer/k-mer scan.
fn harden_soft_mask(contigs: &mut [(String, Vec<u8>)]) {
    for (_, seq) in contigs {
        for b in seq.iter_mut() {
            if b.is_ascii_lowercase() {
                *b = b'N';
            }
        }
    }
}

/// Build an index from a FASTA or 2bit input file (extension decides).
/// `mask` skips soft-masked regions (FASTA lowercase / 2bit mask blocks),
/// matching FastGA `-M` semantics.
pub fn build_from_path(
    path: &str,
    k: usize,
    smer: usize,
    window: usize,
    no_rev: bool,
    mask: bool,
) -> anyhow::Result<PgiIndex> {
    let is_2bit = std::path::Path::new(path)
        .extension()
        .and_then(|e| e.to_str())
        == Some("2bit");
    let mut contigs = if is_2bit {
        read_2bit(path, mask)?
    } else {
        read_fasta(path)?
    };
    if mask {
        harden_soft_mask(&mut contigs);
    }
    build_from_seqs(contigs, k, smer, window, no_rev, false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::fmt::twobit::TwoBitWriter;
    use crate::libs::nt::{rc_key, rolling_kmer_keys};
    use crate::libs::pgi::unpack_position;
    use std::collections::HashSet;

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

    fn write_fasta(dir: &std::path::Path, name: &str, seq: &[u8]) -> String {
        let path = dir.join(name);
        let mut text = String::from(">c1\n");
        for chunk in seq.chunks(60) {
            text.push_str(std::str::from_utf8(chunk).unwrap());
            text.push('\n');
        }
        std::fs::write(&path, text).unwrap();
        path.to_string_lossy().into_owned()
    }

    #[test]
    fn no_duplicate_records_small_k() {
        // Regression: under small-k parameters (k <= window+smer-1) a
        // syncmer position selected twice could be re-emitted after it was
        // already flushed (`collect_one_contig`'s `queued` dedup only holds
        // while the position is pending), corrupting the index with duplicate
        // (kmer, pos, strand) records and falsely inflating the `--freq`
        // filter counts. The index must hold exactly one record per position.
        for (k, smer, window) in [
            (8usize, 8usize, 8usize),
            (10, 8, 5),
            (6, 4, 4),
            (4, 4, 4),
            (40, 8, 5), // defaults: must stay clean too
            (12, 8, 5),
            (20, 16, 6),
        ] {
            for seed in 0..20u64 {
                let seq = pseudo_random_seq(50_000, seed);
                let idx = build_from_seqs(
                    vec![(String::from("c"), seq)],
                    k,
                    smer,
                    window,
                    false,
                    false,
                )
                .unwrap();
                let mut seen: std::collections::HashSet<(u128, u32, u8)> = Default::default();
                for e in &idx.entries {
                    for &rec in
                        &idx.positions[e.pos_start as usize..(e.pos_start + e.freq) as usize]
                    {
                        let (cid, pos, strand) = crate::libs::pgi::unpack_position(rec);
                        assert!(
                            seen.insert((e.kmer, pos, strand)),
                            "duplicate record (k={k}, smer={smer}, window={window}, seed={seed}): \
                             key={:x} pos={pos} strand={strand} cid={cid}",
                            e.kmer
                        );
                    }
                }
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
    fn index_records_match_sequence_positions() {
        // Regression: syncmer positions selected twice (the minimum of two
        // adjacent windows) or queued out of order used to be flushed with
        // the current rolling key, pairing positions with the wrong k-mer
        // and corrupting the index. Every record's key must equal the k-mer
        // starting at its position, and no (key, pos, strand) may repeat.
        let seq = pseudo_random_seq(100_000, 42);
        let idx = build_from_seqs(
            vec![(String::from("c1"), seq.clone())],
            40,
            8,
            5,
            false,
            false,
        )
        .unwrap();
        let codes: [u64; 256] = {
            let mut c = [4u64; 256];
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
        let mut seen: HashSet<(u128, u32, u8)> = HashSet::new();
        for e in &idx.entries {
            for &rec in &idx.positions[e.pos_start as usize..(e.pos_start + e.freq) as usize] {
                let (cid, pos, strand) = unpack_position(rec);
                assert_eq!(cid, 0);
                assert!(
                    pos as usize + 40 <= seq.len(),
                    "position out of range: {pos}"
                );
                let kmer = kmer_key_at(&seq, pos as usize, 40, &codes)
                    .expect("indexed position must be N-free");
                let expected = if strand == 0 {
                    e.kmer
                } else {
                    rc_key(e.kmer, 40)
                };
                assert_eq!(
                    kmer, expected,
                    "key mismatch at pos {pos} (strand {strand})"
                );
                assert!(
                    seen.insert((e.kmer, pos, strand)),
                    "duplicate record at pos {pos}"
                );
            }
        }
        assert!(idx.n_positions() > 0);
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
                .map(|(key, rec)| {
                    let (cid, pos, strand) = unpack_position(rec);
                    (key, cid, pos, strand)
                })
                .collect();
            single_recs.sort_unstable();
            // The closed-syncmer rule can select one position twice (it is
            // the minimum of two adjacent windows); dedupe both sides.
            single_recs.dedup();

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
            ref_recs.dedup();
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
            .map(|(key, rec)| {
                let (cid, pos, strand) = unpack_position(rec);
                (key, cid, pos, strand)
            })
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
        single_recs.dedup();
        ref_recs.dedup();
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
            .map(|(key, rec)| {
                let (cid, pos, strand) = unpack_position(rec);
                (key, cid, pos, strand)
            })
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
        single_recs.dedup();
        ref_recs.dedup();
        assert_eq!(single_recs, ref_recs, "single-pass mismatch with N");
    }

    #[test]
    fn randomized_single_pass_matches_reference() {
        // Property test: on random sequences with randomly placed N runs
        // (including at the very start/end) and different sampling
        // parameters, the single-pass collection must emit exactly the
        // reference (syncmer_dna + rolling keys) records.
        let (k, smer, window) = (40usize, 8usize, 5usize);
        let params = SyncmerParams {
            smer,
            window,
            seed: 7,
        };
        use rand::SeedableRng;
        let mut rng = rand::rngs::StdRng::seed_from_u64(2024);
        use rand::Rng;
        for trial in 0..30 {
            let n = 300 + rng.random_range(0..3000);
            let mut seq: Vec<u8> = (0..n).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();
            // 0-3 N runs of 1-80 bp at random positions (may overlap).
            for _ in 0..rng.random_range(0..4) {
                let start = rng.random_range(0..n);
                let len = 1 + rng.random_range(0..80);
                for b in &mut seq[start..n.min(start + len)] {
                    *b = b'N';
                }
            }

            let mut buf = RecordBuf {
                keys: Vec::new(),
                payloads: Vec::new(),
            };
            collect_one_contig(&seq, 0, k, &params, false, &mut buf);
            let mut single_recs: Vec<(u128, u32, u32, u8)> = buf
                .keys
                .into_iter()
                .zip(buf.payloads)
                .map(|(key, rec)| {
                    let (cid, pos, strand) = unpack_position(rec);
                    (key, cid, pos, strand)
                })
                .collect();
            single_recs.sort_unstable();
            single_recs.dedup();

            let sm = crate::libs::syncmer::syncmer_dna(&seq, &params).unwrap();
            let keys = crate::libs::nt::rolling_kmer_keys(&seq, k);
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
            ref_recs.dedup();
            assert_eq!(
                single_recs, ref_recs,
                "single-pass mismatch at trial {trial} (n={n})"
            );
        }
    }

    #[test]
    fn mask_skips_lowercase_fasta_kmers() {
        // An uppercase zone plus a lowercase (soft-masked) zone with a
        // different sequence: masked zones must not contribute k-mers.
        let upper = pseudo_random_seq(200, 1);
        let lower: Vec<u8> = pseudo_random_seq(200, 2)
            .iter()
            .map(|b| b.to_ascii_lowercase())
            .collect();
        let mut seq = upper.clone();
        seq.extend_from_slice(&lower);
        let dir = tempfile::TempDir::new().unwrap();
        let path = write_fasta(dir.path(), "g.fa", &seq);
        let masked = build_from_path(&path, 10, 4, 2, false, true).unwrap();
        let plain = build_from_path(&path, 10, 4, 2, false, false).unwrap();
        assert!(
            masked.n_unique() < plain.n_unique(),
            "masking must drop k-mers ({} vs {})",
            masked.n_unique(),
            plain.n_unique()
        );
        let plain_keys: HashSet<u128> = plain.entries.iter().map(|e| e.kmer).collect();
        assert!(
            masked.entries.iter().all(|e| plain_keys.contains(&e.kmer)),
            "masked k-mers must be a subset of the unmasked ones"
        );
    }

    #[test]
    fn mask_skips_2bit_mask_blocks() {
        let upper = pseudo_random_seq(200, 1);
        let lower: Vec<u8> = pseudo_random_seq(200, 2)
            .iter()
            .map(|b| b.to_ascii_lowercase())
            .collect();
        let mut dna = String::from_utf8(upper).unwrap();
        dna.push_str(&String::from_utf8(lower).unwrap());
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("g.2bit");
        let mut file = std::fs::File::create(&path).unwrap();
        TwoBitWriter::new(&mut file)
            .write(&[("c1", &dna)], true)
            .unwrap();
        drop(file);

        let path = path.to_string_lossy().into_owned();
        // Without mask the 2bit reads unmasked (all uppercase); with mask the
        // stored mask blocks come back as lowercase and are skipped by build.
        let unmasked = read_2bit(&path, false).unwrap();
        let masked = read_2bit(&path, true).unwrap();
        assert!(
            unmasked[0].1.iter().any(|b| b.is_ascii_uppercase()),
            "unmasked 2bit read is uppercase"
        );
        assert!(
            masked[0].1.iter().any(|b| b.is_ascii_lowercase()),
            "masked 2bit read applies soft mask"
        );
        let idx_masked = build_from_path(&path, 10, 4, 2, false, true).unwrap();
        let idx_plain = build_from_path(&path, 10, 4, 2, false, false).unwrap();
        assert!(idx_masked.n_unique() < idx_plain.n_unique());
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
        let idx = build_from_seqs(vec![(String::from("c"), seq)], 40, 8, 5, false, false).unwrap();
        assert!(idx.entries.len() > 10_000, "too few unique k-mers");
        for w in idx.entries.windows(2) {
            assert!(w[0].kmer <= w[1].kmer, "entries not sorted by k-mer");
        }
    }

    #[test]
    fn no_rev_halves_strands() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let both = build_from_seqs(
            vec![(String::from("c1"), seq.clone())],
            10,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let fwd = build_from_seqs(vec![(String::from("c1"), seq)], 10, 4, 2, true, false).unwrap();
        assert!(both.n_positions() >= fwd.n_positions());
    }
}
