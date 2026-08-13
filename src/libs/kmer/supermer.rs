//! FastK-style super-mer / minimizer two-stage k-mer counting (stage B).
//!
//! Stage 1 partitions each sequence into super-mers: maximal runs of k-mer
//! windows whose canonical m-mer values stay above the run's minimum (a new
//! strictly smaller m-mer cuts the run, as does a `k - m + 1` window gap to
//! the minimum position, `MAX_SUPER`). Each span is packed into a fixed-size
//! record (span bytes, canonical-m-mer orientation, plus the window count)
//! and the records are sorted, collapsing identical spans. Stage 2 expands
//! every unique span into canonical k-mers weighted by the span multiplicity,
//! sorts those, and accumulates counts. The output is byte-identical to
//! [`super::count::count_keys`].
//!
//! On high-coverage low-error data the stage-1 collapse shrinks the sorting
//! volume by roughly the coverage factor (FastK's design target, see
//! `notes/references/fastk.md` §3); on fully unique data the extra stage-1
//! pass costs only a few percent. The minimizer is a fixed-length canonical
//! m-mer (FastK adapts `PAD_LEN` per input; a fixed default matches its
//! typical trained range and keeps the implementation simple).

use super::KmerTable;
use rayon::prelude::*;

/// Default minimizer length (FastK's adaptive `PAD_LEN` typically lands in
/// 10..=13 after training).
pub const DEFAULT_M: usize = 12;

/// Minimizer length for a given `k`. `m <= k - 1` keeps at least two m-mers
/// in every k-mer window; the two-stage path is gated on `k >= 3`.
pub fn minimizer_len(k: usize) -> usize {
    DEFAULT_M.min(k.saturating_sub(1))
}

/// Two-stage super-mer count table with the default minimizer, byte-identical
/// to [`super::count::count_keys`].
pub fn build_table(seqs: &[Vec<u8>], k: usize) -> anyhow::Result<KmerTable> {
    build_table_with_m(seqs, k, minimizer_len(k))
}

/// [`build_table`] with an explicit minimizer length.
pub fn build_table_with_m(seqs: &[Vec<u8>], k: usize, m: usize) -> anyhow::Result<KmerTable> {
    anyhow::ensure!(
        k > 0 && k <= super::key::Kmer::MAX_K,
        "k must be in 1..={}, got {k}",
        super::key::Kmer::MAX_K
    );
    anyhow::ensure!(
        m >= 2 && m < k,
        "minimizer length m must be in 2..=k-1, got m={m} for k={k}"
    );
    let key_bytes = k.div_ceil(4);
    // A span is at most 2k - m bases (m_pos - s <= k - m and the run closes
    // once q - m_pos >= k - m + 1); the +1 leaves slack for the re-scan.
    let max_span = 2 * k - m + 1;
    let span_bytes = max_span.div_ceil(4);
    let rec_bytes = span_bytes + 2;

    let per_seq: Vec<Vec<u8>> = seqs
        .par_iter()
        .map(|seq| pack_sequence(seq, k, m, span_bytes, rec_bytes))
        .collect();
    let n: usize = per_seq.iter().map(Vec::len).sum();
    let mut records: Vec<u8> = Vec::with_capacity(n);
    for mut rec in per_seq {
        records.append(&mut rec);
    }
    let n_records = records.len() / rec_bytes;
    if n_records == 0 {
        return Ok(KmerTable {
            k,
            keys: Vec::new(),
            counts: Vec::new(),
        });
    }
    // Stage 1: sort super-mer records (span + window count), grouping
    // identical spans.
    crate::libs::ds::radix_sort::radix_sort_bytes_par(
        &mut records,
        rec_bytes,
        &mut vec![(); n_records],
    );

    // Stage 2: expand each unique span into weighted canonical k-mers.
    let mut keys: Vec<u8> = Vec::new();
    let mut weights: Vec<u32> = Vec::new();
    let mut i = 0usize;
    while i < n_records {
        let mut j = i + 1;
        while j < n_records
            && records[j * rec_bytes..(j + 1) * rec_bytes]
                == records[i * rec_bytes..(i + 1) * rec_bytes]
        {
            j += 1;
        }
        expand_span(
            &records[i * rec_bytes..(i + 1) * rec_bytes],
            k,
            key_bytes,
            span_bytes,
            (j - i) as u32,
            &mut keys,
            &mut weights,
        );
        i = j;
    }
    let n_entries = weights.len();
    if n_entries == 0 {
        return Ok(KmerTable {
            k,
            keys: Vec::new(),
            counts: Vec::new(),
        });
    }
    crate::libs::ds::radix_sort::radix_sort_bytes_par(&mut keys, key_bytes, &mut weights);
    // Sum weights of equal canonical keys (deduplicating the key buffer, same
    // compaction as `count::count_keys`).
    let mut counts: Vec<u32> = Vec::with_capacity(n_entries);
    let mut w = 0usize;
    let mut idx = 0usize;
    while idx < n_entries {
        let mut j = idx + 1;
        while j < n_entries
            && keys[j * key_bytes..(j + 1) * key_bytes]
                == keys[idx * key_bytes..(idx + 1) * key_bytes]
        {
            j += 1;
        }
        let sum: u64 = weights[idx..j].iter().map(|&c| c as u64).sum();
        if w != idx {
            keys.copy_within(idx * key_bytes..(idx + 1) * key_bytes, w * key_bytes);
        }
        counts.push(sum.min(u32::MAX as u64) as u32);
        w += 1;
        idx = j;
    }
    keys.truncate(w * key_bytes);
    Ok(KmerTable { k, keys, counts })
}

/// Fixed parameters shared by the per-sequence packers.
struct PackCtx<'a> {
    k: usize,
    m: usize,
    span_bytes: usize,
    rec_bytes: usize,
    codes: &'a [u64; 256],
}

/// Pack all super-mers of one sequence into fixed-size records (N-free runs
/// are partitioned independently, matching `canonical_keys`).
fn pack_sequence(seq: &[u8], k: usize, m: usize, span_bytes: usize, rec_bytes: usize) -> Vec<u8> {
    let codes = super::base_codes();
    let ctx = PackCtx {
        k,
        m,
        span_bytes,
        rec_bytes,
        codes: &codes,
    };
    let mut records = Vec::with_capacity(seq.len().saturating_sub(k) / (k - m).max(1) * rec_bytes);
    let n = seq.len();
    let mut start = 0usize;
    while start < n {
        while start < n && codes[seq[start] as usize] == 4 {
            start += 1;
        }
        if start >= n {
            break;
        }
        let mut end = start;
        while end < n && codes[seq[end] as usize] != 4 {
            end += 1;
        }
        if end - start >= k {
            pack_run(seq, start, end, &ctx, &mut records);
        }
        start = end;
    }
    records
}

/// Partition one N-free run into super-mers and append their records.
fn pack_run(seq: &[u8], run_start: usize, run_end: usize, ctx: &PackCtx, records: &mut Vec<u8>) {
    let PackCtx { k, m, codes, .. } = *ctx;
    let l = run_end - run_start;
    let n_windows = l - k + 1;
    let win_m = k - m;
    // Rolling canonical m-mer values: `fwd` packs the m-mer starting at the
    // current position, `rc` its reverse complement, so the minimum of the
    // two is strand-invariant (FastK's `flp` comparison).
    let mask = if m >= 16 {
        u32::MAX
    } else {
        (1u32 << (2 * m)) - 1
    };
    let mut fwd = 0u32;
    let mut rc = 0u32;
    for j in 0..m {
        fwd = (fwd << 2) | codes[seq[run_start + j] as usize] as u32;
    }
    for j in 0..m {
        rc = (rc << 2) | (3 - codes[seq[run_start + m - 1 - j] as usize]) as u32;
    }
    let mut mval = Vec::with_capacity(l - m + 1);
    let mut flp = Vec::with_capacity(l - m + 1);
    for q in 0..=l - m {
        mval.push(fwd.min(rc));
        flp.push(rc < fwd);
        if q + m < l {
            let leave = codes[seq[run_start + q] as usize] as u32;
            let enter = codes[seq[run_start + q + m] as usize] as u32;
            fwd = ((fwd << 2) | enter) & mask;
            rc = (rc >> 2) | ((3 - leave) << (2 * m - 2));
        }
    }
    // The first k-mer window contains m-mers at 0..=win_m.
    let mut mc = u32::MAX;
    let mut m_pos = 0usize;
    for (q, &v) in mval.iter().enumerate().take(win_m + 1) {
        if v < mc {
            mc = v;
            m_pos = q;
        }
    }
    let max_super = k - m + 1;
    let mut s = 0usize; // run start window
    for i in 1..n_windows {
        let q = i + win_m;
        let mp = mval[q];
        if mp < mc || q - m_pos >= max_super {
            // Close the run on windows [s, i-1]: span bases [s, i-1+k).
            emit_span(
                seq,
                run_start + s,
                (i - 1) + k - s,
                flp[m_pos],
                ctx,
                records,
            );
            if mp < mc {
                mc = mp;
                m_pos = q;
            } else {
                // Force cut: the new minimum is the min of the remaining
                // m-mers (when both conditions hold, mp < mc makes the
                // re-scan land on q, so the branches agree).
                mc = u32::MAX;
                for (x, &v) in mval.iter().enumerate().skip(m_pos + 1).take(q - m_pos) {
                    if v < mc {
                        mc = v;
                        m_pos = x;
                    }
                }
            }
            s = i;
        }
    }
    // Final run: windows [s, n_windows-1], span [s, l).
    emit_span(seq, run_start + s, l - s, flp[m_pos], ctx, records);
}

/// Append one super-mer record: packed span in the defining m-mer's canonical
/// orientation, zero-padded to `span_bytes`, then the window count (u16 BE).
fn emit_span(
    seq: &[u8],
    base_start: usize,
    span_len: usize,
    flip: bool,
    ctx: &PackCtx,
    records: &mut Vec<u8>,
) {
    let PackCtx {
        k,
        span_bytes,
        rec_bytes,
        ..
    } = *ctx;
    let base = records.len();
    records.resize(base + rec_bytes, 0);
    pack_span(
        &mut records[base..base + span_bytes],
        seq,
        base_start,
        span_len,
        flip,
    );
    let sln = span_len - k + 1;
    records[base + span_bytes] = (sln >> 8) as u8;
    records[base + span_bytes + 1] = sln as u8;
}

/// Pack `len` bases of `seq` from `start` into `out`, 2 bits per base with
/// the 5'-most base in the high bits (FastK byte layout), either forward or
/// reverse-complemented.
fn pack_span(out: &mut [u8], seq: &[u8], start: usize, len: usize, flip: bool) {
    for i in 0..len {
        let b = seq[start + if flip { len - 1 - i } else { i }];
        let c = match b {
            b'A' | b'a' => 0u8,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            _ => 3, // T/t: the caller only reaches N-free runs
        };
        let c = if flip { 3 - c } else { c };
        out[i / 4] |= c << (2 * (3 - i % 4));
    }
}

/// Emit every canonical k-mer of a unique span, weighted by its multiplicity.
fn expand_span(
    rec: &[u8],
    k: usize,
    key_bytes: usize,
    span_bytes: usize,
    ct: u32,
    keys: &mut Vec<u8>,
    weights: &mut Vec<u32>,
) {
    let sln = ((rec[span_bytes] as usize) << 8) | rec[span_bytes + 1] as usize;
    // The packed span continues past the first k-mer, so the low pad bits of
    // the copied key bytes hold the span's next bases; clear them to the
    // FastK zero-pad layout before rolling (the direct path keeps pads zero
    // by construction).
    let mut key_buf = [0u8; super::key::Kmer::MAX_K / 4];
    key_buf[..key_bytes].copy_from_slice(&rec[..key_bytes]);
    let pad = 8 * key_bytes - 2 * k;
    if pad < 8 {
        key_buf[key_bytes - 1] &= 0xFFu8 << pad;
    }
    let mut win = super::key::Kmer::from_bytes(k, &key_buf[..key_bytes]);
    let mut win_rc = win.rc();
    emit_canonical(&win, &win_rc, ct, keys, weights);
    for o in 1..sln {
        let idx = k - 1 + o;
        let b = (rec[idx / 4] >> (2 * (3 - idx % 4))) & 3;
        win.push_right(b);
        win_rc.push_left(3 - b);
        emit_canonical(&win, &win_rc, ct, keys, weights);
    }
}

/// Canonical emit matching `canonical_keys` (first-half-byte comparison).
fn emit_canonical(
    win: &super::key::Kmer,
    win_rc: &super::key::Kmer,
    ct: u32,
    keys: &mut Vec<u8>,
    weights: &mut Vec<u32>,
) {
    let half = win.key_bytes().div_ceil(2);
    if win.to_bytes()[..half] <= win_rc.to_bytes()[..half] {
        keys.extend_from_slice(win.to_bytes());
    } else {
        keys.extend_from_slice(win_rc.to_bytes());
    }
    weights.push(ct);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::kmer::count;

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

    /// Random DNA with Ns sprinkled at fixed offsets.
    fn noisy_block(len: usize, seed: u64) -> Vec<u8> {
        let mut b = random_block(len, seed);
        for i in (7..len).step_by(31) {
            b[i] = b'N';
        }
        b
    }

    fn assert_same_table(direct: &KmerTable, supermer: &KmerTable) {
        assert_eq!(direct.k, supermer.k);
        assert_eq!(
            direct.keys, supermer.keys,
            "packed keys must match (k={})",
            direct.k
        );
        assert_eq!(
            direct.counts, supermer.counts,
            "counts must match (k={})",
            direct.k
        );
    }

    #[test]
    fn matches_direct_on_random_data() {
        for k in [5usize, 8, 17, 31, 64, 100] {
            let seqs = vec![
                random_block(500, 1),
                random_block(333, 2),
                noisy_block(700, 3),
                random_block(50, 4),
            ];
            let direct = count::build_table(&seqs, k).unwrap();
            let supermer = build_table(&seqs, k).unwrap();
            assert_same_table(&direct, &supermer);
        }
    }

    #[test]
    fn matches_direct_with_duplicates() {
        let block = random_block(80, 42);
        let mut seqs = Vec::new();
        for i in 0..50u64 {
            seqs.push(block.clone());
            seqs.push(random_block(60, i + 7));
        }
        for k in [5usize, 17, 31, 100] {
            let direct = count::build_table(&seqs, k).unwrap();
            let supermer = build_table(&seqs, k).unwrap();
            assert_same_table(&direct, &supermer);
        }
    }

    #[test]
    fn multiplicity_above_u16() {
        // 70,000 identical reads collapse to one super-mer whose k-mers each
        // carry a weight above u16::MAX (the direct path must agree).
        let seqs = vec![b"ACGTACGTACGT".to_vec(); 70_000];
        let k = 5usize;
        let direct = count::build_table(&seqs, k).unwrap();
        let supermer = build_table(&seqs, k).unwrap();
        assert_same_table(&direct, &supermer);
        assert_eq!(
            supermer.counts.iter().map(|&c| c as usize).sum::<usize>(),
            70_000 * 8
        );
    }

    #[test]
    fn matches_direct_on_minimizer_edges() {
        // Smallest k (m = 2) and m = k-1 (two m-mers per window).
        let seqs = vec![random_block(300, 9), noisy_block(250, 10)];
        for (k, m) in [(3usize, 2usize), (7, 6), (13, 12)] {
            let direct = count::build_table(&seqs, k).unwrap();
            let supermer = build_table_with_m(&seqs, k, m).unwrap();
            assert_same_table(&direct, &supermer);
        }
    }

    #[test]
    fn matches_direct_across_k_sweep() {
        // Sweep k across key-byte boundaries (4/5, 8/9, 12/13, ...) with a
        // mix of clean and noisy reads.
        let seqs = vec![
            random_block(600, 21),
            noisy_block(450, 22),
            random_block(90, 23),
            b"ACGTACGTACGTNNACGTACGTACGTACGTACGTACGTNNACGTACGTACGT".to_vec(),
        ];
        for k in 3..=40usize {
            let direct = count::build_table(&seqs, k).unwrap();
            let supermer = build_table(&seqs, k).unwrap();
            assert_same_table(&direct, &supermer);
        }
    }

    #[test]
    fn kmer_shared_across_spans_sums_weights() {
        // A canonical k-mer appearing inside two different super-mers (with
        // different multiplicities) must accumulate both weights.
        let seqs = vec![
            b"ACGTACGTACGTACGTACGT".to_vec(), // one span, k-mer at offset 4
            b"ACGTACGTACGTACGTACGT".to_vec(), // duplicate -> weight 2
            b"TTACGTACGTACGTACGTTT".to_vec(), // same k-mer in another span
        ];
        let k = 8usize;
        let direct = count::build_table(&seqs, k).unwrap();
        let supermer = build_table(&seqs, k).unwrap();
        assert_same_table(&direct, &supermer);
    }

    #[test]
    fn reverse_complement_reads_merge() {
        // The same region seen from both strands must count identically (the
        // span orientation flip makes the opposite-strand spans merge).
        let fwd = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let rev = fwd
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect::<Vec<_>>();
        for k in [7usize, 17, 31] {
            let direct = count::build_table(&[fwd.clone(), rev.clone()], k).unwrap();
            let supermer = build_table(&[fwd.clone(), rev.clone()], k).unwrap();
            assert_same_table(&direct, &supermer);
            // Two reads -> every window counted twice (canonical merges).
            assert_eq!(
                supermer.counts.iter().map(|&c| c as usize).sum::<usize>(),
                2 * (fwd.len() - k + 1)
            );
        }
    }

    #[test]
    fn empty_and_short_inputs() {
        let direct = count::build_table(&[], 17).unwrap();
        let supermer = build_table(&[], 17).unwrap();
        assert_same_table(&direct, &supermer);

        let short = vec![
            Vec::new(),
            b"ACG".to_vec(),
            b"NNNN".to_vec(),
            b"ACGTACGTACGT".to_vec(),
        ];
        for k in [5usize, 17] {
            let direct = count::build_table(&short, k).unwrap();
            let supermer = build_table(&short, k).unwrap();
            assert_same_table(&direct, &supermer);
        }
    }

    #[test]
    fn case_insensitive() {
        let lower = vec![b"acgtacgtacgtnnacgtacgtacgt".to_vec()];
        let upper = vec![b"ACGTACGTACGTNNACGTACGTACGT".to_vec()];
        for k in [5usize, 17] {
            let a = build_table(&lower, k).unwrap();
            let b = build_table(&upper, k).unwrap();
            assert_same_table(&a, &b);
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(build_table_with_m(&[b"ACGT".to_vec()], 3, 1).is_err());
        assert!(build_table_with_m(&[b"ACGT".to_vec()], 3, 3).is_err());
        assert!(build_table_with_m(&[b"ACGT".to_vec()], 129, 12).is_err());
    }
}
