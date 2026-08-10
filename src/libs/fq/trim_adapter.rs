//! BBTools 39.38 `bbduk`-compatible adapter trimming and k-mer filtering.
//!
//! Ports `jgi.BBDuk` for the anchr trim pipeline parameter set: k-mer
//! right-trimming (`ktrim=r`, `mink`, `hdist`), overlap trimming (`tbo`),
//! even pair trimming (`tpe`), quality trimming (`qtrim=r`), `maxns`,
//! `minlen`, `ftm`, and `tossbrokenreads`. Output is byte-identical to
//! `bbduk.sh ... ordered=t` (see tests/bbtools/Lambda/README.md).

use crate::libs::fmt::fq::write_fq;
use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use anyhow::{Context, Result};
use std::collections::HashSet;
use std::io::Write;
use std::sync::Arc;

/// Bits per base (DNA).
const BPB: u32 = 2;

/// Options for a BBDuk-style trimming pass.
#[derive(Debug, Clone)]
pub struct AdapterTrimOptions {
    /// Full k-mer size (bbduk `k`).
    pub k: usize,
    /// Minimum short k-mer size (bbduk `mink`).
    pub mink: usize,
    /// Reference hamming distance (bbduk `hdist`).
    pub hdist: usize,
    /// Right-trim at the first matching k-mer (`ktrim=r`).
    pub ktrim_right: bool,
    /// Trim implied adapters from mate overlap (`tbo`).
    pub tbo: bool,
    /// Trim both mates to equal length (`tpe`).
    pub tpe: bool,
    /// Right quality trim (`qtrim=r`).
    pub qtrim_right: bool,
    /// Quality threshold for `qtrim=r` (`trimq`).
    pub trimq: u8,
    /// Discard reads shorter than this (`minlen`).
    pub minlen: usize,
    /// Discard reads with more than this many N bases (`maxns`).
    pub maxns: i64,
    /// Right-trim lengths to a multiple of this (`ftm`; 0 disables).
    pub ftm: usize,
    /// Discard a pair when either mate fails (`tossbrokenreads`).
    pub toss_broken_reads: bool,
    /// Reference FASTA of adapters/contaminants (`ref`).
    pub ref_file: String,
    /// Input quality ASCII offset (33 or 64).
    pub quality_base: u8,
    /// Discard reads with more than this many matching k-mers (filter mode;
    /// used when `ktrim_right` is false).
    pub max_bad_kmers: usize,
}

/// Precomputed bit masks for a k-mer size.
struct Masks {
    mask: i64,
    shift2: u32,
    length_mask: Vec<i64>,
    clear_mask: Vec<i64>,
    set_mask: [[i64; 4]; 64],
    right_mask: Vec<i64>,
    kmask: i64,
}

impl Masks {
    fn new(k: usize) -> Self {
        let mut length_mask = vec![0i64; k + 1];
        let mut clear_mask = vec![0i64; k];
        let mut right_mask = vec![0i64; k + 1];
        let mut set_mask = [[0i64; 4]; 64];
        for i in 0..k {
            let shift = (BPB * i as u32) as usize;
            clear_mask[i] = !(3i64 << shift);
            right_mask[i] = !((-1i64) << shift);
            length_mask[i] = 1i64 << shift;
            for (j, v) in set_mask[i].iter_mut().enumerate() {
                *v = (j as i64) << shift;
            }
        }
        right_mask[k] = !((-1i64) << (BPB * k as u32));
        length_mask[k] = 1i64 << (BPB * k as u32);
        let mask = if 2 * k >= 64 {
            -1
        } else {
            !((-1i64) << (BPB * k as u32))
        };
        let kmask = length_mask[k];
        Self {
            mask,
            shift2: BPB * k as u32 - 2,
            length_mask,
            clear_mask,
            set_mask,
            right_mask,
            kmask,
        }
    }
}

/// 2-bit code tables (dna.AminoAcid `symbolToNumber0` / complement).
fn symbol_to_number0(b: u8) -> i64 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 0,
    }
}

fn symbol_to_complement_number0(b: u8) -> i64 {
    match b {
        b'A' | b'a' => 3,
        b'C' | b'c' => 2,
        b'G' | b'g' => 1,
        b'T' | b't' | b'U' | b'u' => 0,
        _ => 0,
    }
}

fn symbol_to_number(b: u8) -> i64 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => -1,
    }
}

fn is_fully_defined(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'U' | b'u'
    )
}

/// 8-bit reverse-complement lookup table (AminoAcid.rcompBinaryTable).
fn rcomp_table() -> [u8; 256] {
    let mut r = [0u8; 256];
    for (i, v) in r.iter_mut().enumerate() {
        let mut x = (!(i as i64)) & 0xFF;
        let mut out = 0u8;
        for _ in 0..4 {
            out = (out << 2) | ((x & 3) as u8);
            x >>= 2;
        }
        *v = out;
    }
    r
}

/// Reverse-complement a 2-bit k-mer (AminoAcid.reverseComplementBinaryFast).
fn rcomp(mut kmer: i64, len: usize, table: &[u8; 256]) -> i64 {
    let mut out = 0i64;
    let extra = len & 3;
    let mut k = len;
    for _ in 0..extra {
        out = (out << 2) | ((!kmer) & 3);
        kmer >>= 2;
    }
    k -= extra;
    while k > 0 {
        out = (out << 8) | (table[(kmer & 0xFF) as usize] as i64);
        kmer >>= 8;
        k -= 4;
    }
    out
}

/// Canonical table key: max(kmer, rkmer) with the length bit set.
fn to_value(kmer: i64, rkmer: i64, length_mask: i64) -> i64 {
    (kmer.max(rkmer)) | length_mask
}

/// Builds the reference k-mer set (full k-mers + short end k-mers + hdist
/// single-substitution variants), matching BBDuk's LoadThread.
fn build_table(ref_file: &str, opts: &AdapterTrimOptions) -> Result<HashSet<i64>> {
    let mut reader = SeqReader::new(ref_file)
        .with_context(|| format!("Failed to open reader for {}", ref_file))?;
    let masks = Masks::new(opts.k);
    let rc = rcomp_table();
    let mut table = HashSet::new();
    let mut rec = SeqRecord::new();
    while reader.read_record(&mut rec)? {
        add_ref_sequence(rec.sequence(), opts, &masks, &rc, &mut table);
    }
    Ok(table)
}

/// Stores all k-mers (and variants) of one reference sequence.
#[allow(clippy::too_many_arguments)]
fn add_ref_sequence(
    seq: &[u8],
    opts: &AdapterTrimOptions,
    masks: &Masks,
    rc: &[u8; 256],
    table: &mut HashSet<i64>,
) {
    if seq.len() < opts.k {
        return;
    }
    let k = opts.k;
    let mut kmer = 0i64;
    let mut rkmer = 0i64;
    let mut len = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        let x = symbol_to_number0(b);
        let x2 = symbol_to_complement_number0(b);
        kmer = ((kmer << BPB) | x) & masks.mask;
        rkmer = ((rkmer >> BPB) | (x2 << masks.shift2)) & masks.mask;
        if is_fully_defined(b) {
            len += 1;
        } else {
            len = 0;
            rkmer = 0;
        }
        if len >= k {
            let extra_base = if i + 1 >= seq.len() {
                -1
            } else {
                symbol_to_number(seq[i + 1])
            };
            add_kmer(kmer, rkmer, k, extra_base, opts.hdist, masks, rc, table);
            if opts.mink > 0 && opts.mink < k {
                if i == k - 1 {
                    add_right_shift(kmer, rkmer, opts, masks, rc, table);
                }
                if i == seq.len() - 1 {
                    add_left_shift(kmer, rkmer, extra_base, opts, masks, rc, table);
                }
            }
        }
    }
}

/// Adds a k-mer (with `dist` hamming mutations) to the table.
#[allow(clippy::only_used_in_recursion, clippy::too_many_arguments)]
fn add_kmer(
    kmer: i64,
    rkmer: i64,
    len: usize,
    extra_base: i64,
    dist: usize,
    masks: &Masks,
    rc: &[u8; 256],
    table: &mut HashSet<i64>,
) {
    let key = to_value(kmer, rkmer, masks.length_mask[len]);
    table.insert(key);
    if dist > 0 {
        for j in 0..4 {
            for i in 0..len {
                let temp = (kmer & masks.clear_mask[i]) | masks.set_mask[i][j];
                if temp != kmer {
                    let rtemp = rcomp(temp, len, rc);
                    add_kmer(temp, rtemp, len, extra_base, dist - 1, masks, rc, table);
                }
            }
        }
        // editDistance=0 in this pipeline: no insertion/deletion variants.
    }
}

/// Short k-mers (mink..k) of the reference prefix (addToMapRightShift).
fn add_right_shift(
    mut kmer: i64,
    mut rkmer: i64,
    opts: &AdapterTrimOptions,
    masks: &Masks,
    rc: &[u8; 256],
    table: &mut HashSet<i64>,
) {
    for i in (opts.mink..opts.k).rev() {
        let extra_base = kmer & 3;
        kmer >>= BPB;
        rkmer &= masks.right_mask[i];
        add_kmer(kmer, rkmer, i, extra_base, opts.hdist, masks, rc, table);
    }
}

/// Short k-mers (mink..k) of the reference suffix (addToMapLeftShift).
fn add_left_shift(
    mut kmer: i64,
    mut rkmer: i64,
    extra_base: i64,
    opts: &AdapterTrimOptions,
    masks: &Masks,
    rc: &[u8; 256],
    table: &mut HashSet<i64>,
) {
    for i in (opts.mink..opts.k).rev() {
        kmer &= masks.right_mask[i];
        rkmer >>= BPB;
        add_kmer(kmer, rkmer, i, extra_base, opts.hdist, masks, rc, table);
    }
}

/// Exact table lookup (bbduk query hamming distance is 0 in this pipeline).
fn get_value(kmer: i64, rkmer: i64, length_mask: i64, table: &HashSet<i64>) -> bool {
    let key = to_value(kmer, rkmer, length_mask);
    table.contains(&key)
}

/// Counts matching k-mers of a read (BBDuk countSetKmers); stops early once
/// the count exceeds `max_bad_kmers`.
fn count_set_kmers(seq: &[u8], opts: &AdapterTrimOptions, table: &HashSet<i64>) -> usize {
    if seq.len() < opts.k || table.is_empty() {
        return 0;
    }
    let masks = Masks::new(opts.k);
    let mut kmer = 0i64;
    let mut rkmer = 0i64;
    let mut len = 0usize;
    let mut found = 0usize;
    for &b in seq {
        let x = symbol_to_number0(b);
        let x2 = symbol_to_complement_number0(b);
        kmer = ((kmer << BPB) | x) & masks.mask;
        rkmer = ((rkmer >> BPB) | (x2 << masks.shift2)) & masks.mask;
        len += 1; // forbidNs=false with hdist>0
        if len >= opts.k && get_value(kmer, rkmer, masks.kmask, table) {
            found += 1;
            if found > opts.max_bad_kmers {
                return found;
            }
        }
    }
    found
}

/// A read with its buffers, mimicking `stream.Read` trimming.
struct ReadBuf {
    seq: Vec<u8>,
    qual: Vec<u8>,
    discarded: bool,
}

impl ReadBuf {
    fn len(&self) -> usize {
        self.seq.len()
    }

    fn count_undefined(&self) -> usize {
        self.seq
            .iter()
            .filter(|&&b| symbol_to_number(b) < 0)
            .count()
    }

    /// Trim to [left, right] inclusive (TrimRead.trimToPosition).
    fn trim_to_position(&mut self, left: usize, right: usize, min_result: usize) {
        let len = self.len();
        let left_trim = left.min(len);
        let mut right_trim = len.saturating_sub(right.saturating_add(1));
        let min_result = len.min(min_result);
        if left_trim + right_trim + min_result > len {
            right_trim = 1.max(len.saturating_sub(min_result));
        }
        if left_trim + right_trim > 0 {
            let keep = len - left_trim - right_trim;
            self.seq.drain(keep..);
            self.seq.drain(0..left_trim);
            if self.qual.len() >= left_trim + right_trim {
                self.qual.drain(keep..);
                self.qual.drain(0..left_trim);
            } else {
                self.qual.clear();
            }
        }
    }
}

/// k-mer right trimming (bbduk `ktrim=r`). Returns bases trimmed.
fn ktrim(read: &mut ReadBuf, opts: &AdapterTrimOptions, table: &HashSet<i64>) -> usize {
    let min_len = 1.max(opts.k.min(opts.mink));
    if read.len() < min_len || table.is_empty() {
        return 0;
    }
    let masks = Masks::new(opts.k);
    let bases = read.seq.clone();
    let mut kmer = 0i64;
    let mut rkmer = 0i64;
    let mut len = 0usize;
    let mut min_loc: i64 = i64::MAX;
    let mut found = 0usize;
    for (i, &b) in bases.iter().enumerate() {
        let x = symbol_to_number0(b);
        let x2 = symbol_to_complement_number0(b);
        kmer = ((kmer << BPB) | x) & masks.mask;
        rkmer = ((rkmer >> BPB) | (x2 << masks.shift2)) & masks.mask;
        len += 1; // forbidNs=false with hdist>0
        if len >= opts.k && get_value(kmer, rkmer, masks.kmask, table) {
            min_loc = min_loc.min(i as i64 - opts.k as i64 + 1);
            found += 1;
        }
    }
    if found == 0 && opts.mink > 0 && opts.mink < opts.k {
        // Short k-mers at the read's right end.
        kmer = 0;
        rkmer = 0;
        len = 0;
        let stop = bases.len();
        let lim = stop as i64 - opts.k as i64;
        let mut i = stop as i64 - 1;
        while i > lim {
            let b = bases[i as usize];
            let x = symbol_to_number0(b);
            let x2 = symbol_to_complement_number0(b);
            kmer |= x << (BPB * len as u32);
            rkmer = ((rkmer << BPB) | x2) & masks.mask;
            len += 1;
            if len >= opts.mink && get_value(kmer, rkmer, masks.length_mask[len], table) {
                min_loc = i;
                found += 1;
            }
            i -= 1;
        }
    }
    if found == 0 {
        return 0;
    }
    let before = read.len();
    if opts.ktrim_right {
        // Keep [0, min_loc-1].
        let right = (min_loc - 1).max(0) as usize;
        read.trim_to_position(0, right, 1);
    }
    before - read.len()
}

/// BBMergeOverlapper.mateByOverlapRatioJava (no-quality path, strict mode).
///
/// `r2` must already be reverse-complemented. Returns the best insert size or
/// -1; `ambig` reports an ambiguous overlap.
#[allow(clippy::too_many_arguments)]
fn mate_by_overlap_ratio(
    a: &[u8],
    b: &[u8],
    min_overlap0: usize,
    min_overlap: usize,
    min_insert0: usize,
    min_insert: usize,
    max_ratio: f32,
    min_second_ratio: f32,
    margin: f32,
    offset: f32,
    g_incr: f32,
    b_incr: f32,
    ambig: &mut bool,
) -> i32 {
    let alen = a.len() as i32;
    let blen = b.len() as i32;
    let min_length = alen.min(blen);
    let n = b'N';
    let mut max_ratio = max_ratio;
    let x = find_best_ratio(
        a,
        b,
        min_overlap0,
        min_overlap,
        min_insert,
        max_ratio,
        offset,
        g_incr,
        b_incr,
    );
    if x > max_ratio {
        return -1;
    }
    max_ratio = max_ratio.min(x);
    let margin2 = (margin + offset) / min_length as f32;
    let mut best_insert = -1i32;
    let mut best_ratio = 1f32;
    let mut best_ambig = false;
    let mut second_best_ratio = 1f32;
    let extra_mult = 1.2f32;
    let extra_badlimit = 20f32;
    let largest = alen + blen - min_overlap0 as i32;
    let smallest = min_insert0 as i32;
    let mut insert = largest;
    while insert >= smallest {
        let istart = if insert <= blen { 0 } else { insert - blen };
        let jstart = if insert >= blen { 0 } else { blen - insert };
        let overlap_length = (alen - istart).min(blen - jstart).min(insert);
        let badlimit = extra_mult * (best_ratio.min(max_ratio) * margin * overlap_length as f32)
            + 1.0
            + extra_badlimit;
        let mut good = 0f32;
        let mut bad = 0f32;
        let imax = istart + overlap_length;
        let mut i = istart;
        let mut j = jstart;
        while i < imax && bad <= badlimit {
            let ca = a[i as usize];
            let cb = b[j as usize];
            if ca == cb {
                if ca != n {
                    good += g_incr;
                }
            } else {
                bad += b_incr;
            }
            i += 1;
            j += 1;
        }
        if bad <= badlimit {
            if bad == 0.0 && good > min_overlap0 as f32 && good < min_overlap as f32 {
                *ambig = true;
                return -1;
            }
            let ratio = (bad + offset) / overlap_length as f32;
            if ratio < best_ratio * margin {
                let this_ambig = ratio * margin >= best_ratio || good < min_overlap as f32;
                if ratio < best_ratio {
                    second_best_ratio = best_ratio;
                    best_insert = insert;
                    best_ratio = ratio;
                } else if ratio < second_best_ratio {
                    second_best_ratio = ratio;
                }
                best_ambig = this_ambig;
                if (best_ambig && best_ratio < margin2) || second_best_ratio < min_second_ratio {
                    *ambig = true;
                    return -1;
                }
            }
        }
        insert -= 1;
    }
    if second_best_ratio < min_second_ratio {
        best_ambig = true;
    }
    if !best_ambig && best_ratio > max_ratio {
        best_insert = -1;
    }
    *ambig = best_ambig;
    if best_insert < 0 {
        -1
    } else {
        best_insert
    }
}

/// findBestRatio: fast pre-screen for overlap detection.
#[allow(clippy::too_many_arguments)]
fn find_best_ratio(
    a: &[u8],
    b: &[u8],
    min_overlap0: usize,
    min_overlap: usize,
    min_insert: usize,
    max_ratio: f32,
    offset: f32,
    g_incr: f32,
    b_incr: f32,
) -> f32 {
    let alen = a.len() as i32;
    let blen = b.len() as i32;
    let n = b'N';
    let mut best_ratio = max_ratio + 0.0001;
    let halfmax = max_ratio * 0.5;
    let largest = alen + blen - min_overlap as i32;
    let smallest = min_insert as i32;
    let mut insert = largest;
    while insert >= smallest {
        let istart = if insert <= blen { 0 } else { insert - blen };
        let jstart = if insert >= blen { 0 } else { blen - insert };
        let overlap_length = (alen - istart).min(blen - jstart).min(insert);
        let badlimit = best_ratio * overlap_length as f32 + 20.0;
        let mut good = 0f32;
        let mut bad = 0f32;
        let imax = istart + overlap_length;
        let mut i = istart;
        let mut j = jstart;
        while i < imax && bad <= badlimit {
            let ca = a[i as usize];
            let cb = b[j as usize];
            if ca == cb {
                if ca != n {
                    good += g_incr;
                }
            } else {
                bad += b_incr;
            }
            i += 1;
            j += 1;
        }
        if bad <= badlimit {
            if bad == 0.0 && good > min_overlap0 as f32 && good < min_overlap as f32 {
                return 100.0;
            }
            let ratio = (bad + offset) / overlap_length as f32;
            if ratio < best_ratio {
                best_ratio = ratio;
                if good >= min_overlap as f32 && ratio < halfmax {
                    return best_ratio;
                }
            }
        }
        insert -= 1;
    }
    best_ratio
}

/// Expected errors of a read (Read.expectedErrors, countUndefined=false).
fn expected_errors(seq: &[u8], qual: &[u8], prob: &[f32; 128]) -> f32 {
    let mut sum = 0f32;
    for (i, &b) in seq.iter().enumerate() {
        if is_fully_defined(b) {
            let q = qual.get(i).copied().unwrap_or(0) as usize;
            sum += prob[q.min(127)];
        }
    }
    sum
}

fn prob_error() -> [f32; 128] {
    let mut r = [0f32; 128];
    for (i, v) in r.iter_mut().enumerate() {
        *v = (10f64.powf(-0.1 * i as f64)) as f32;
    }
    r[0] = 0.75;
    r[1] = 0.7;
    r
}

/// Quality right trim (TrimRead.trimFast with optimalMode).
fn qtrim_right(read: &mut ReadBuf, trimq: u8, prob: &[f32; 128]) -> usize {
    if read.len() == 0 {
        return 0;
    }
    let bases = read.seq.clone();
    let qual = read.qual.clone();
    if qual.is_empty() {
        return 0;
    }
    // avgErrorRate = phredToProbError(trimq), like Parser.trimE().
    let avg_error_rate = 10f64.powf(-0.1 * trimq as f64) as f32;
    let nprob = (avg_error_rate * 1.1).clamp(0.75, 1.0);
    let mut max_score = 0f32;
    let mut score = 0f32;
    let mut max_loc = -1i64;
    let mut max_count = -1i64;
    let mut count = 0i64;
    for (i, &b) in bases.iter().enumerate() {
        let q = qual[i];
        let prob_error = if b == b'N' || q < 1 {
            nprob
        } else {
            prob[q as usize]
        };
        let delta = avg_error_rate - prob_error;
        score += delta;
        if score > 0.0 {
            count += 1;
            if score > max_score || (score == max_score && count > max_count) {
                max_score = score;
                max_count = count;
                max_loc = i as i64;
            }
        } else {
            score = 0.0;
            count = 0;
        }
    }
    let right = if max_score > 0.0 {
        read.len() - max_loc as usize - 1
    } else {
        read.len()
    };
    let before = read.len();
    if right > 0 {
        read.trim_to_position(0, read.len().saturating_sub(right + 1), 1);
    }
    before - read.len()
}

/// Process one interleaved pair (or single); returns surviving pair.
fn process_pair(
    r1: &mut ReadBuf,
    mut r2: Option<&mut ReadBuf>,
    opts: &AdapterTrimOptions,
    table: &HashSet<i64>,
    prob: &[f32; 128],
) {
    let minlen = opts.minlen;
    // forceTrimModulo (ftm) before k-mer trimming.
    if opts.ftm > 0 {
        for r in [Some(&mut *r1), r2.as_deref_mut()].into_iter().flatten() {
            if r.len() > 0 {
                let b0 = r.len() - 1 - r.len() % opts.ftm;
                r.trim_to_position(0, b0, 1);
            }
        }
    }
    if r1.len() < minlen {
        r1.discarded = true;
    }
    if let Some(r2) = r2.as_deref_mut() {
        if r2.len() < minlen {
            r2.discarded = true;
        }
    }

    let mut xsum = 0usize;
    if opts.ktrim_right {
        if !r1.discarded {
            xsum += ktrim(r1, opts, table);
            if r1.len() < minlen {
                r1.discarded = true;
            }
        }
        if let Some(r2) = r2.as_deref_mut() {
            if !r2.discarded {
                xsum += ktrim(r2, opts, table);
                if r2.len() < minlen {
                    r2.discarded = true;
                }
            }
        }
        // tpe: trim the longer mate to the shorter length.
        if opts.tpe && xsum > 0 {
            if let Some(r2) = r2.as_deref_mut() {
                if r1.len() != r2.len() {
                    if r1.len() > r2.len() {
                        r1.trim_to_position(0, r2.len() - 1, 1);
                    } else {
                        r2.trim_to_position(0, r1.len() - 1, 1);
                    }
                }
            }
        }
    } else if !table.is_empty() {
        // Filter mode: discard reads with more than max_bad_kmers matches.
        if !r1.discarded && count_set_kmers(&r1.seq, opts, table) > opts.max_bad_kmers {
            r1.discarded = true;
        }
        if let Some(r2) = r2.as_deref_mut() {
            if !r2.discarded && count_set_kmers(&r2.seq, opts, table) > opts.max_bad_kmers {
                r2.discarded = true;
            }
        }
    }
    let remove = if opts.toss_broken_reads {
        r1.discarded || r2.as_ref().is_some_and(|r| r.discarded)
    } else {
        r1.discarded && r2.as_ref().is_none_or(|r| r.discarded)
    };

    // tbo: overlap trimming (only when both reads survive).
    if !remove && opts.tbo {
        if let Some(r2) = r2.as_deref_mut() {
            let e1 = expected_errors(&r1.seq, &r1.qual, prob);
            let e2 = expected_errors(&r2.seq, &r2.qual, prob);
            if e1.max(e2) < 15.0 {
                // mateByOverlapRatio expects r2 reverse-complemented.
                r2.seq.reverse();
                for (i, &b) in r2.seq.clone().iter().enumerate() {
                    r2.seq[i] = match b {
                        b'A' | b'a' => b'T',
                        b'C' | b'c' => b'G',
                        b'G' | b'g' => b'C',
                        b'T' | b't' | b'U' | b'u' => b'A',
                        _ => b'N',
                    };
                }
                let mut ambig = false;
                let best_insert = mate_by_overlap_ratio(
                    &r1.seq, &r2.seq, 7, 14, 16, 40, 0.05, 0.12, 9.0, 0.5, 0.95, 0.95, &mut ambig,
                );
                r2.seq.reverse();
                for (i, &b) in r2.seq.clone().iter().enumerate() {
                    r2.seq[i] = match b {
                        b'A' | b'a' => b'T',
                        b'C' | b'c' => b'G',
                        b'G' | b'g' => b'C',
                        b'T' | b't' | b'U' | b'u' => b'A',
                        _ => b'N',
                    };
                }
                if best_insert > 0 && !ambig {
                    if (best_insert as usize) < r1.len() {
                        r1.trim_to_position(0, best_insert as usize - 1, 1);
                    }
                    if (best_insert as usize) < r2.len() {
                        r2.trim_to_position(0, best_insert as usize - 1, 1);
                    }
                }
            }
        }
    }

    // Quality trim + final filters.
    if opts.qtrim_right {
        if !r1.discarded {
            qtrim_right(r1, opts.trimq, prob);
            if r1.len() < minlen {
                r1.discarded = true;
            }
        }
        if let Some(r2) = r2.as_deref_mut() {
            if !r2.discarded {
                qtrim_right(r2, opts.trimq, prob);
                if r2.len() < minlen {
                    r2.discarded = true;
                }
            }
        }
    }
    if !r1.discarded && opts.maxns >= 0 && r1.count_undefined() > opts.maxns as usize {
        r1.discarded = true;
    }
    if let Some(r2) = r2.as_deref_mut() {
        if !r2.discarded && opts.maxns >= 0 && r2.count_undefined() > opts.maxns as usize {
            r2.discarded = true;
        }
    }
    if opts.toss_broken_reads {
        if r1.discarded {
            if let Some(r2) = r2 {
                r2.discarded = true;
            }
        } else if let Some(r2) = r2 {
            if r2.discarded {
                r1.discarded = true;
            }
        }
    }
}

/// Runs adapter trimming over interleaved pairs (1 file) or R1/R2 (2 files).
pub fn trim_adapter<W: Write>(
    infiles: &[String],
    out: &mut W,
    opts: &AdapterTrimOptions,
    parallel: usize,
) -> Result<()> {
    let table = Arc::new(build_table(&opts.ref_file, opts)?);
    let prob = prob_error();
    let opts = Arc::new(opts.clone());
    let quality_base = opts.quality_base;
    let pairs = crate::libs::fq::pairs::PairReader::new(infiles)?;
    crate::libs::par::ordered_map(
        pairs,
        parallel,
        move |pair| {
            let pair = pair?;
            Ok(process_one(&pair, &opts, &table, &prob))
        },
        |t| {
            let Some(t) = t else {
                return Ok(());
            };
            write_record(out, &t.r1, &t.r1_seq, &t.r1_qual, quality_base)?;
            if let (Some(r2), Some(s2), Some(q2)) = (&t.r2, &t.r2_seq, &t.r2_qual) {
                write_record(out, r2, s2, q2, quality_base)?;
            }
            Ok(())
        },
    )
}

/// A surviving pair with trimmed sequence/quality buffers.
#[derive(Clone)]
struct TrimmedPair {
    /// Original R1 record (header preserved).
    pub r1: SeqRecord,
    /// Trimmed R1 sequence.
    pub r1_seq: Vec<u8>,
    /// Trimmed R1 quality (phred).
    pub r1_qual: Vec<u8>,
    /// Original R2 record.
    pub r2: Option<SeqRecord>,
    /// Trimmed R2 sequence.
    pub r2_seq: Option<Vec<u8>>,
    /// Trimmed R2 quality (phred).
    pub r2_qual: Option<Vec<u8>>,
}

/// Processes one pair through the trim/filter pipeline.
fn process_one(
    pair: &(SeqRecord, Option<SeqRecord>),
    opts: &AdapterTrimOptions,
    table: &HashSet<i64>,
    prob: &[f32; 128],
) -> Option<TrimmedPair> {
    let (rec1, rec2) = pair;
    let mut r1 = make_read_buf(rec1, opts.quality_base);
    let mut r2 = rec2.as_ref().map(|r| make_read_buf(r, opts.quality_base));
    process_pair(&mut r1, r2.as_mut(), opts, table, prob);
    if r1.discarded {
        return None;
    }
    Some(TrimmedPair {
        r1: rec1.clone(),
        r1_seq: r1.seq,
        r1_qual: r1.qual,
        r2: rec2.clone(),
        r2_seq: r2.as_ref().and_then(|r| {
            if r.discarded {
                None
            } else {
                Some(r.seq.clone())
            }
        }),
        r2_qual: r2.as_ref().and_then(|r| {
            if r.discarded {
                None
            } else {
                Some(r.qual.clone())
            }
        }),
    })
}

/// Writes a trimmed FASTQ record with the original header.
fn write_record<W: Write>(
    w: &mut W,
    rec: &SeqRecord,
    seq: &[u8],
    qual: &[u8],
    quality_base: u8,
) -> anyhow::Result<()> {
    let comment = rec.comment();
    let header = if comment.is_empty() {
        rec.name().to_string()
    } else {
        format!("{} {}", rec.name(), comment)
    };
    let mut out_qual = qual.to_vec();
    for q in out_qual.iter_mut() {
        *q = q.saturating_add(quality_base);
    }
    write_fq(w, &header, seq, &out_qual)?;
    Ok(())
}

/// Builds a working read buffer with phred (offset-subtracted) quality.
fn make_read_buf(rec: &SeqRecord, quality_base: u8) -> ReadBuf {
    ReadBuf {
        seq: rec.sequence().to_vec(),
        qual: rec
            .quality_scores()
            .iter()
            .map(|&q| q.saturating_sub(quality_base))
            .collect(),
        discarded: false,
    }
}
