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
//! pass costs only a few percent. The minimizer is a canonical m-mer whose
//! length follows the anchr-measured heuristic `min(12, max(5, k/4))`
//! (`anchr notes/design/asm-assemble.md` §12.3: k=31 best at m=8, k=63/99
//! at m=12), matching FastK's per-input training in the common range.

use super::KmerTable;
use rayon::prelude::*;

/// Minimizer length cap (FastK's adaptive `PAD_LEN` typically lands in
/// 10..=13 after training; the heuristic starts at 5 for small k).
pub const DEFAULT_M: usize = 12;

/// Minimizer length for a given `k`: `min(12, max(5, ceil(k/4)))`, bounded
/// by `m <= k - 1` so every k-mer window keeps at least two m-mers.
pub fn minimizer_len(k: usize) -> usize {
    DEFAULT_M.min(5.max(k.div_ceil(4))).min(k.saturating_sub(1))
}

/// Two-stage super-mer count table with the default minimizer, byte-identical
/// to [`super::count::count_keys`].
pub fn build_table(seqs: &[Vec<u8>], k: usize) -> anyhow::Result<KmerTable> {
    build_table_with_m(seqs, k, minimizer_len(k))
}

/// [`build_table`] with an explicit minimizer length.
pub fn build_table_with_m(seqs: &[Vec<u8>], k: usize, m: usize) -> anyhow::Result<KmerTable> {
    build_impl(seqs, k, m)
}

/// [`build_table`] over borrowed sequence slices, so callers with a
/// streaming path (e.g. anchr's `TadpoleTable`) can feed records without
/// materializing `Vec<Vec<u8>>`.
pub fn build_table_slices(seqs: &[&[u8]], k: usize) -> anyhow::Result<KmerTable> {
    build_impl(seqs, k, minimizer_len(k))
}

/// [`build_table_slices`] with an explicit minimizer length.
pub fn build_table_slices_with_m(seqs: &[&[u8]], k: usize, m: usize) -> anyhow::Result<KmerTable> {
    build_impl(seqs, k, m)
}

/// [`build_table_slices`] with a sliding-window quality gate (anchr
/// `min_prob` semantics): a window counts only when the product of its base
/// correctness probabilities reaches `min_prob`; `min_prob <= 0.0` or empty
/// qualities disable the gate.
pub fn build_table_slices_qual(
    seqs: &[&[u8]],
    quals: &[&[u8]],
    k: usize,
    min_prob: f32,
) -> anyhow::Result<KmerTable> {
    build_impl_qual(seqs, quals, k, minimizer_len(k), min_prob)
}

/// [`build_table_slices_qual`] with an explicit minimizer length.
pub fn build_table_slices_qual_with_m(
    seqs: &[&[u8]],
    quals: &[&[u8]],
    k: usize,
    m: usize,
    min_prob: f32,
) -> anyhow::Result<KmerTable> {
    build_impl_qual(seqs, quals, k, m, min_prob)
}

/// Shared two-stage implementation over any slice of byte sequences.
fn build_impl<S: AsRef<[u8]> + Sync>(seqs: &[S], k: usize, m: usize) -> anyhow::Result<KmerTable> {
    anyhow::ensure!(
        k > 0 && k <= super::key::Kmer::MAX_K,
        "k must be in 1..={}, got {k}",
        super::key::Kmer::MAX_K
    );
    anyhow::ensure!(
        (2..=16).contains(&m) && m < k,
        "minimizer length m must be in 2..=min(16, k-1), got m={m} for k={k}"
    );
    // A span is at most 2k - m bases (m_pos - s <= k - m and the run closes
    // once q - m_pos >= k - m + 1); the +1 leaves slack for the re-scan.
    let max_span = 2 * k - m + 1;
    let span_bytes = max_span.div_ceil(4);
    let rec_bytes = span_bytes + 2;
    let codes = super::base_codes();
    let ctx = PackCtx {
        k,
        m,
        span_bytes,
        rec_bytes,
        codes: &codes,
    };
    // Pack in coarse chunks into one contiguous buffer per chunk (FastK
    // packs per IO block): far fewer allocations than one `Vec` per read.
    const PACK_CHUNK: usize = 4096;
    let t0 = std::time::Instant::now();
    let per_chunk: Vec<Vec<u8>> = seqs
        .par_chunks(PACK_CHUNK)
        .map(|chunk| {
            let est = chunk.iter().map(|s| s.as_ref().len()).sum::<usize>();
            let mut records =
                Vec::with_capacity(est.saturating_sub(k) / (k - m).max(1) * rec_bytes + rec_bytes);
            for seq in chunk {
                pack_sequence_into(seq.as_ref(), &ctx, &mut records);
            }
            records
        })
        .collect();
    let pack = t0.elapsed();
    let t1 = std::time::Instant::now();
    let table = finish_records(per_chunk, &ctx)?;
    let finish = t1.elapsed();
    if std::env::var_os("PGR_SUPERMER_TIMING").is_some() {
        eprintln!(
            "supermer k={k} m={m}: pack={:.3}s finish={:.3}s total={:.3}s",
            pack.as_secs_f64(),
            finish.as_secs_f64(),
            (pack + finish).as_secs_f64()
        );
    }
    Ok(table)
}

/// Sort packed super-mer records, group identical spans, and expand them
/// into weighted canonical k-mers (shared by plain and quality-gated paths).
fn finish_records(per_chunk: Vec<Vec<u8>>, ctx: &PackCtx) -> anyhow::Result<KmerTable> {
    let k = ctx.k;
    let key_bytes = k.div_ceil(4);
    let span_bytes = ctx.span_bytes;
    let rec_bytes = ctx.rec_bytes;
    let n: usize = per_chunk.iter().map(Vec::len).sum();
    let t0 = std::time::Instant::now();
    let mut records: Vec<u8> = Vec::with_capacity(n);
    for mut rec in per_chunk {
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
    let sort1 = t0.elapsed();
    // Stage 2: group identical spans (adjacent after the sort), then expand
    // each unique span into weighted canonical k-mers. The boundary scan is
    // a parallel filter over adjacent records; the expansion is parallelized
    // over coarse blocks of groups (each block appends into its own buffers,
    // so the writes never overlap).
    let t1 = std::time::Instant::now();
    let boundaries: Vec<usize> = (1..n_records)
        .into_par_iter()
        .filter(|&i| {
            records[i * rec_bytes..(i + 1) * rec_bytes]
                != records[(i - 1) * rec_bytes..i * rec_bytes]
        })
        .collect();
    // Build the group table from the boundary indices without shifting the
    // array (an `insert(0, ..)` here would move millions of elements).
    let mut groups: Vec<(usize, u32, usize)> = Vec::with_capacity(boundaries.len() + 1); // (rec, ct, sln)
    let mut prev = 0usize;
    for &start in std::iter::once(&0)
        .chain(boundaries.iter())
        .chain(std::iter::once(&n_records))
    {
        if start != prev {
            let rec = &records[prev * rec_bytes..(prev + 1) * rec_bytes];
            let sln = ((rec[span_bytes] as usize) << 8) | rec[span_bytes + 1] as usize;
            groups.push((prev, (start - prev) as u32, sln));
            prev = start;
        }
    }
    if groups.is_empty() {
        return Ok(KmerTable {
            k,
            keys: Vec::new(),
            counts: Vec::new(),
        });
    }
    const EXPAND_CHUNK: usize = 1 << 13; // groups per parallel block
    let n_blocks = groups.len().div_ceil(EXPAND_CHUNK);
    let per_block: Vec<(Vec<u8>, Vec<u32>)> = (0..n_blocks)
        .into_par_iter()
        .map(|b| {
            let start = b * EXPAND_CHUNK;
            let end = (start + EXPAND_CHUNK).min(groups.len());
            let est: usize = groups[start..end].iter().map(|&(_, _, sln)| sln).sum();
            let mut keys = Vec::with_capacity(est * key_bytes);
            let mut weights = Vec::with_capacity(est);
            for &(ri, ct, _) in &groups[start..end] {
                expand_span(
                    &records[ri * rec_bytes..(ri + 1) * rec_bytes],
                    k,
                    key_bytes,
                    span_bytes,
                    ct,
                    &mut keys,
                    &mut weights,
                );
            }
            (keys, weights)
        })
        .collect();
    let expand = t1.elapsed();
    let n_entries: usize = per_block.iter().map(|(k, _)| k.len() / key_bytes).sum();
    let t2 = std::time::Instant::now();
    let mut keys: Vec<u8> = Vec::with_capacity(n_entries * key_bytes);
    let mut weights: Vec<u32> = Vec::with_capacity(n_entries);
    for (mut kb, mut wb) in per_block {
        keys.append(&mut kb);
        weights.append(&mut wb);
    }
    debug_assert_eq!(keys.len(), weights.len() * key_bytes);
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
    let sort2 = t2.elapsed();
    if std::env::var_os("PGR_SUPERMER_TIMING").is_some() {
        eprintln!(
            "  sort1={:.3}s expand={:.3}s sort2={:.3}s spans={n_records} emitted={n_entries}",
            sort1.as_secs_f64(),
            expand.as_secs_f64(),
            sort2.as_secs_f64()
        );
    }
    Ok(KmerTable { k, keys, counts })
}

/// Quality-gated two-stage count; `min_prob <= 0.0` or empty qualities fall
/// back to the ungated path.
fn build_impl_qual(
    seqs: &[&[u8]],
    quals: &[&[u8]],
    k: usize,
    m: usize,
    min_prob: f32,
) -> anyhow::Result<KmerTable> {
    anyhow::ensure!(
        k > 0 && k <= super::key::Kmer::MAX_K,
        "k must be in 1..={}, got {k}",
        super::key::Kmer::MAX_K
    );
    anyhow::ensure!(
        (2..=16).contains(&m) && m < k,
        "minimizer length m must be in 2..=min(16, k-1), got m={m} for k={k}"
    );
    if min_prob <= 0.0 || quals.is_empty() {
        return build_impl(seqs, k, m);
    }
    anyhow::ensure!(
        quals.len() == seqs.len(),
        "{} sequences but {} quality strings",
        seqs.len(),
        quals.len()
    );
    for (seq, qual) in seqs.iter().zip(quals) {
        anyhow::ensure!(
            qual.is_empty() || qual.len() == seq.len(),
            "sequence length {} does not match quality length {}",
            seq.len(),
            qual.len()
        );
    }
    let max_span = 2 * k - m + 1;
    let span_bytes = max_span.div_ceil(4);
    let rec_bytes = span_bytes + 2;
    let codes = super::base_codes();
    let ctx = PackCtx {
        k,
        m,
        span_bytes,
        rec_bytes,
        codes: &codes,
    };
    let (prob_correct, prob_correct_inv) = prob_tables();
    // Pack in coarse chunks into one contiguous buffer per chunk (FastK
    // packs per IO block): far fewer allocations than one `Vec` per read.
    const PACK_CHUNK: usize = 4096;
    let per_chunk: Vec<Vec<u8>> = seqs
        .par_chunks(PACK_CHUNK)
        .zip(quals.par_chunks(PACK_CHUNK))
        .map(|(seq_chunk, qual_chunk)| {
            let est = seq_chunk.iter().map(|s| s.len()).sum::<usize>();
            let mut records =
                Vec::with_capacity(est.saturating_sub(k) / (k - m).max(1) * rec_bytes + rec_bytes);
            for (seq, qual) in seq_chunk.iter().zip(qual_chunk) {
                pack_sequence_into_qual(
                    seq,
                    qual,
                    &ctx,
                    min_prob,
                    &prob_correct,
                    &prob_correct_inv,
                    &mut records,
                );
            }
            records
        })
        .collect();
    finish_records(per_chunk, &ctx)
}

/// Fixed parameters shared by the per-sequence packers.
struct PackCtx<'a> {
    k: usize,
    m: usize,
    span_bytes: usize,
    rec_bytes: usize,
    codes: &'a [u64; 256],
}

/// Append the super-mer records of one sequence to `records` (N-free runs
/// are partitioned independently, matching `canonical_keys`).
fn pack_sequence_into(seq: &[u8], ctx: &PackCtx, records: &mut Vec<u8>) {
    let codes = ctx.codes;
    let k = ctx.k;
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
            pack_run(seq, start, end, ctx, records);
        }
        start = end;
    }
}

/// BBTools `QualityTools.PROB_ERROR` (q=0 -> 0.75, q=1 -> 0.7, else 10^-0.1q).
fn prob_error() -> [f32; 128] {
    let mut r = [0f32; 128];
    for (i, v) in r.iter_mut().enumerate() {
        *v = (10f64.powf(-0.1 * i as f64)) as f32;
    }
    r[0] = 0.75;
    r[1] = 0.7;
    r
}

/// Base correctness probability tables (BBTools `PROB_CORRECT[_INVERSE]`).
fn prob_tables() -> ([f32; 128], [f32; 128]) {
    let err = prob_error();
    let mut correct = [0f32; 128];
    let mut inverse = [0f32; 128];
    for i in 0..128 {
        let c = 1.0 - err[i];
        correct[i] = c;
        inverse[i] = 1.0 / c;
    }
    (correct, inverse)
}

/// Like [`pack_sequence_into`], but skips windows whose base correctness
/// probability falls below `min_prob` (anchr `minprob` semantics); a skipped
/// window cuts the super-mer run like an N boundary.
fn pack_sequence_into_qual(
    seq: &[u8],
    qual: &[u8],
    ctx: &PackCtx,
    min_prob: f32,
    prob_correct: &[f32; 128],
    prob_correct_inv: &[f32; 128],
    records: &mut Vec<u8>,
) {
    if qual.is_empty() {
        // Empty quality disables the gate for this read (anchr per-read rule).
        pack_sequence_into(seq, ctx, records);
        return;
    }
    let codes = ctx.codes;
    let k = ctx.k;
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
            // Sliding window probability in the same f32 order as anchr
            // `count_read_kmers_packed`; invalid windows split the stretch
            // into contiguous valid segments.
            let mut seg_start = start;
            let mut prob = 1.0f32;
            let mut len = 0usize;
            let mut i = start;
            while i < end {
                let q = (qual[i] as usize).min(127);
                prob *= prob_correct[q];
                if len >= k {
                    prob *= prob_correct_inv[(qual[i - k] as usize).min(127)];
                }
                len += 1;
                i += 1;
                if len >= k && prob < min_prob {
                    let win_start = i - k;
                    if win_start > seg_start {
                        // Segment windows [seg_start, win_start), covering
                        // bases [seg_start, win_start + k - 1); the end must
                        // include the last window's final k-1 bases.
                        pack_run(seq, seg_start, win_start + k - 1, ctx, records);
                    }
                    seg_start = win_start + 1;
                }
            }
            if end > seg_start && end - seg_start >= k {
                pack_run(seq, seg_start, end, ctx, records);
            }
        }
        start = end;
    }
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
    fn long_k_supermer_sorts_large_records() {
        // Regression: k > 128 makes super-mer records exceed 64 bytes
        // (rec_bytes = ceil((2k-m+1)/4) + 2; k=160, m=12 -> 80 bytes), which
        // overflowed the fixed [0u8; 64] scratch buffers in radix_sort_bytes.
        let seqs = vec![random_block(600, 160), noisy_block(450, 161)];
        for k in [141usize, 160, 256] {
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
        assert!(build_table_with_m(&[b"ACGT".to_vec()], 257, 12).is_err());
    }

    #[test]
    fn slices_api_matches_vec_api() {
        let seqs = vec![random_block(500, 31), noisy_block(400, 32)];
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        for k in [17usize, 31, 100] {
            let a = build_table(&seqs, k).unwrap();
            let b = build_table_slices(&refs, k).unwrap();
            let c = build_table_slices_with_m(&refs, k, 8).unwrap();
            assert_same_table(&a, &b);
            let d = build_table_with_m(&seqs, k, 8).unwrap();
            assert_same_table(&c, &d);
        }
    }

    /// The minimizer length only affects stage-1 collapse, never the output:
    /// every valid window is expanded exactly once regardless of the run
    /// partition, so any legal `m` yields a byte-identical table. This is
    /// the basis for reusing a fixed-m minimizer extraction across k.
    #[test]
    fn output_independent_of_minimizer() {
        let seqs = [random_block(500, 91), noisy_block(400, 92)];
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        for k in [13usize, 21, 31, 61] {
            let base = build_table_slices(&refs, k).unwrap();
            // m <= 16 keeps the packed m-mer inside u32 (2m bits); larger m
            // is rejected by the API's validation.
            for m in 2..=k.min(16).min(k - 1) {
                let t = build_table_slices_with_m(&refs, k, m).unwrap();
                assert_same_table(&base, &t);
            }
        }
    }

    /// Deterministic quality bytes in the Phred+33 range (30..=69).
    fn random_quals(len: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..len)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                30 + ((x >> 33) as u8) % 40
            })
            .collect()
    }

    /// Reference: quality-gated emission + `count_keys`, mirroring anchr
    /// `count_read_kmers_packed` (same sliding f32 order) but independent
    /// of the super-mer record pipeline.
    fn reference_qual_table(seqs: &[&[u8]], quals: &[&[u8]], k: usize, min_prob: f32) -> KmerTable {
        let (pc, pci) = prob_tables();
        let mut keys: Vec<u8> = Vec::new();
        for (&seq, &qual) in seqs.iter().zip(quals) {
            if qual.is_empty() {
                crate::libs::kmer::canonical_keys(seq, k, |_, km| {
                    keys.extend_from_slice(km.to_bytes());
                });
                continue;
            }
            let mut prob = 1.0f32;
            let mut len = 0usize;
            let mut i = 0usize;
            while i < seq.len() {
                if crate::libs::kmer::base_codes()[seq[i] as usize] == 4 {
                    len = 0;
                    prob = 1.0;
                    i += 1;
                    continue;
                }
                let q = (qual[i] as usize).min(127);
                prob *= pc[q];
                if len >= k {
                    prob *= pci[(qual[i - k] as usize).min(127)];
                }
                len += 1;
                i += 1;
                if len >= k && prob >= min_prob {
                    let start = i - k;
                    let km = crate::libs::kmer::key::Kmer::from_bases(&seq[start..start + k], k)
                        .expect("N-free window");
                    keys.extend_from_slice(km.canonical().to_bytes());
                }
            }
        }
        count::count_keys(keys, k)
    }

    #[test]
    fn qual_gated_matches_reference() {
        let seqs: Vec<Vec<u8>> = (0..8u64)
            .map(|s| random_block(250 + s as usize * 37, s + 11))
            .collect();
        let quals: Vec<Vec<u8>> = seqs
            .iter()
            .enumerate()
            .map(|(i, s)| random_quals(s.len(), 100 + i as u64))
            .collect();
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        let qrefs: Vec<&[u8]> = quals.iter().map(Vec::as_slice).collect();
        for min_prob in [0.5f32, 0.9, 0.99, 1.0] {
            for k in [7usize, 13, 31] {
                let got = build_table_slices_qual(&refs, &qrefs, k, min_prob).unwrap();
                let expected = reference_qual_table(&refs, &qrefs, k, min_prob);
                assert_same_table(&expected, &got);
            }
        }
    }

    #[test]
    fn qual_gated_zero_min_prob_matches_ungated() {
        let seqs = [random_block(500, 21), noisy_block(400, 22)];
        let quals: Vec<Vec<u8>> = seqs.iter().map(|s| random_quals(s.len(), 7)).collect();
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        let qrefs: Vec<&[u8]> = quals.iter().map(Vec::as_slice).collect();
        for k in [13usize, 31] {
            let a = build_table_slices(&refs, k).unwrap();
            let b = build_table_slices_qual(&refs, &qrefs, k, 0.0).unwrap();
            assert_same_table(&a, &b);
        }
    }

    #[test]
    fn qual_gated_empty_quals_matches_ungated() {
        let seqs = [random_block(500, 31), noisy_block(400, 32)];
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        for k in [13usize, 31] {
            let a = build_table_slices(&refs, k).unwrap();
            let b = build_table_slices_qual(&refs, &[], k, 0.9).unwrap();
            assert_same_table(&a, &b);
        }
    }

    #[test]
    fn qual_gated_length_mismatch_errors() {
        let seqs = [random_block(100, 1), random_block(100, 2)];
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        let quals = [random_quals(100, 3)];
        let qrefs: Vec<&[u8]> = quals.iter().map(Vec::as_slice).collect();
        assert!(build_table_slices_qual(&refs, &qrefs, 13, 0.9).is_err());

        let one = [random_block(100, 4)];
        let refs1: Vec<&[u8]> = one.iter().map(Vec::as_slice).collect();
        let short = [random_quals(99, 5)];
        let qshort: Vec<&[u8]> = short.iter().map(Vec::as_slice).collect();
        assert!(build_table_slices_qual(&refs1, &qshort, 13, 0.9).is_err());
    }

    #[test]
    fn qual_gated_with_m_matches_default() {
        let seqs = [random_block(500, 41), noisy_block(400, 42)];
        let quals: Vec<Vec<u8>> = seqs.iter().map(|s| random_quals(s.len(), 9)).collect();
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        let qrefs: Vec<&[u8]> = quals.iter().map(Vec::as_slice).collect();
        for k in [13usize, 31] {
            let a = build_table_slices_qual(&refs, &qrefs, k, 0.9).unwrap();
            let b =
                build_table_slices_qual_with_m(&refs, &qrefs, k, minimizer_len(k), 0.9).unwrap();
            assert_same_table(&a, &b);
        }
    }

    #[test]
    fn qual_gated_all_low_is_empty() {
        let seq = random_block(300, 51);
        let refs = [seq.as_slice()];
        let quals = [vec![0u8; 300]];
        let qrefs = [quals[0].as_slice()];
        // q=0 -> correct 0.25; 0.25^13 < 0.99 gates out every window.
        let table = build_table_slices_qual(&refs, &qrefs, 13, 0.99).unwrap();
        assert!(table.keys.is_empty());
        assert!(table.counts.is_empty());
    }

    #[test]
    fn qual_gated_mixed_empty_quals() {
        let seqs = [random_block(300, 61), random_block(400, 62)];
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        // One read without qualities: its windows must all count.
        let quals: Vec<Vec<u8>> = vec![Vec::new(), random_quals(400, 63)];
        let qrefs: Vec<&[u8]> = quals.iter().map(Vec::as_slice).collect();
        for min_prob in [0.9f32, 0.99] {
            let got = build_table_slices_qual(&refs, &qrefs, 13, min_prob).unwrap();
            let expected = reference_qual_table(&refs, &qrefs, 13, min_prob);
            assert_same_table(&expected, &got);
            // The ungated read alone contributes all 300-13+1 windows.
            let total: u32 = got.counts.iter().sum();
            assert!(total >= (300 - 13 + 1) as u32);
        }
    }

    /// Alternating valid/invalid windows: a q=0 base drags the k windows
    /// containing it below `min_prob`, producing adjacent valid -> invalid
    /// -> valid transitions and single-window segments. The ungated
    /// high-quality windows are every `k+1`-th window, all others gated out.
    #[test]
    fn qual_gated_alternating_validity() {
        let seq = random_block(400, 71);
        let k = 7usize;
        let mut qual = vec![40u8; seq.len()];
        for i in (k..seq.len()).step_by(k + 1) {
            qual[i] = 0;
        }
        let refs = [seq.as_slice()];
        let qrefs = [qual.as_slice()];
        for min_prob in [0.9f32, 0.95, 0.99] {
            let got = build_table_slices_qual(&refs, &qrefs, k, min_prob).unwrap();
            let expected = reference_qual_table(&refs, &qrefs, k, min_prob);
            assert_same_table(&expected, &got);
            // Valid windows are the multiples of k+1 (0, 8, 16, ... <= len-k).
            let total: u32 = got.counts.iter().sum();
            assert_eq!(total as usize, (seq.len() - k) / (k + 1) + 1);
        }
    }
}
