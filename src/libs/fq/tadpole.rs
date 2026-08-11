//! Tadpole-compatible k-mer error correction, read extension, and tossing.
//!
//! Ports the BBTools `tadpole.sh` correct/extend/discard modes: a canonical
//! k-mer count table (quality-filtered by `minprob`), per-read error
//! correction by local reassembly through the k-mer graph, conservative read
//! extension that stops at branches, and junk/low-depth read discarding.

use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use crate::libs::fq::qual::{from_phred, to_phred};
use crate::libs::nt::rev_comp;
use anyhow::Result;
use std::collections::HashMap;
use std::io::Write;

/// Default k-mer length (tadpole.sh `k`).
pub const DEFAULT_K: usize = 31;

/// Options mirroring the tadpole.sh defaults used by the anchr merge flow.
#[derive(Debug, Clone)]
pub struct TadpoleOptions {
    /// K-mer length.
    pub k: usize,
    /// Ignore k-mers whose probability of being error-free is below this.
    pub min_prob: f32,
    /// Minimum k-mer depth to seed an extension.
    pub min_count_seed: usize,
    /// Minimum k-mer depth to continue an extension.
    pub min_count_extend: usize,
    /// Branch ratio at high depth (branchmult1).
    pub branch_mult1: f32,
    /// Branch ratio at low depth (branchmult2).
    pub branch_mult2: f32,
    /// Second-highest depth considered "low" (branchlower).
    pub branch_lower_const: usize,
    /// Error ratio multiplier (errormult1).
    pub error_mult1: f32,
    /// Alternative error ratio multiplier (errormult2).
    pub error_mult2: f32,
    /// Quality factor for the error multiplier.
    pub error_mult_q_factor: f32,
    /// Max second-highest depth for the low-depth error rule (errorlowerconst).
    pub error_lower_const: usize,
    /// Minimum depth of a k-mer to be considered correct (mincountcorrect).
    pub min_count_correct: usize,
    /// Absolute path-similarity tolerance (pathsimilarityconstant).
    pub path_similarity_constant: usize,
    /// Fractional path-similarity tolerance (pathsimilarityfraction).
    pub path_similarity_fraction: f32,
    /// K-mers to verify after an error in reassembly (errorextensionreassemble).
    pub error_extension_reassemble: usize,
    /// K-mers to verify in pincer mode (errorextensionpincer).
    pub error_extension_pincer: usize,
    /// K-mers to verify in tail mode (errorextensiontail).
    pub error_extension_tail: usize,
    /// Do not correct bases within this distance of read ends (deadzone).
    pub dead_zone: usize,
    /// Sliding-window length for reassembly quality filtering (window).
    pub window_len: usize,
    /// Max corrections in a window (windowcount).
    pub window_count: usize,
    /// Max quality sum in a window (qualsum).
    pub window_qual_sum: usize,
    /// Undo corrections that lower k-mer coverage (eccrollback).
    pub ecc_rollback: bool,
    /// Run k-mer reassembly error correction (tadpole `ecc`; off in extend
    /// and discard-only modes, matching Java's per-mode default).
    pub ecc: bool,
    /// Require both directions to agree in the read middle (requirebidirectional).
    pub ecc_require_bidirectional: bool,
    /// Extend to the right by at most this many bases.
    pub extend_right: usize,
    /// Extend to the left by at most this many bases.
    pub extend_left: usize,
    /// Trim random trailing bases of partial extensions (extendrollback).
    pub extension_rollback: usize,
    /// Discard reads that cannot be used for assembly (tossjunk).
    pub toss_junk: bool,
    /// Discard reads containing k-mers at or below this depth (tossdepth).
    pub toss_depth: i64,
    /// Discard reads with uncorrectable errors (tossuncorrectable).
    pub toss_uncorrectable: bool,
    /// Minimum fraction of low-depth k-mers to discard a read (lowdepthfraction).
    pub low_depth_discard_fraction: f32,
    /// Only discard a pair if both reads fail (requirebothbad).
    pub require_both_bad: bool,
}

impl Default for TadpoleOptions {
    fn default() -> Self {
        Self {
            k: DEFAULT_K,
            min_prob: 0.5,
            min_count_seed: 3,
            min_count_extend: 2,
            branch_mult1: 20.0,
            branch_mult2: 3.0,
            branch_lower_const: 3,
            error_mult1: 16.0,
            error_mult2: 2.6,
            error_mult_q_factor: 0.002,
            error_lower_const: 4,
            min_count_correct: 3,
            path_similarity_constant: 3,
            path_similarity_fraction: 0.45,
            error_extension_reassemble: 5,
            error_extension_pincer: 5,
            error_extension_tail: 9,
            dead_zone: 0,
            window_len: 12,
            window_count: 6,
            window_qual_sum: 80,
            ecc_rollback: true,
            ecc: false,
            ecc_require_bidirectional: true,
            extend_right: 0,
            extend_left: 0,
            extension_rollback: 3,
            toss_junk: false,
            toss_depth: -1,
            toss_uncorrectable: false,
            low_depth_discard_fraction: 0.0,
            require_both_bad: false,
        }
    }
}

/// Canonical k-mer count table (2-bit encoding, forward vs reverse-complement).
#[derive(Debug, Clone, Default)]
pub struct TadpoleTable {
    map: HashMap<Kmer, u32>,
}

impl TadpoleTable {
    /// Builds the table from reads (bases, phred qualities) with `minprob`
    /// quality filtering, mirroring `KmerTableSetU.addKmersToTable`.
    pub fn build(reads: &[(Vec<u8>, Vec<u8>)], k: usize, min_prob: f32) -> Self {
        let mut map: HashMap<Kmer, u32> = HashMap::new();
        let (prob_correct, prob_correct_inv) = prob_tables();
        for (bases, quals) in reads {
            count_read_kmers(
                &mut map,
                bases,
                quals,
                k,
                min_prob,
                &prob_correct,
                &prob_correct_inv,
            );
        }
        Self { map }
    }

    /// Count of the canonical form of `kmer` (0 when absent).
    pub(crate) fn get_count(&self, kmer: &Kmer) -> u32 {
        self.map.get(&kmer.canonical()).copied().unwrap_or(0)
    }

    /// Iterates the distinct canonical k-mers and their counts.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&Kmer, u32)> {
        self.map.iter().map(|(k, &c)| (k, c))
    }

    /// Counts of the four right-extensions of `kmer`.
    pub(crate) fn fill_right_counts(&self, kmer: &Kmer) -> [u32; 4] {
        let mut out = [0u32; 4];
        for (i, c) in out.iter_mut().enumerate() {
            let mut x = kmer.clone();
            x.push_right(i as u8);
            *c = self.get_count(&x);
        }
        out
    }

    /// Counts of the four left-extensions of `kmer`.
    pub(crate) fn fill_left_counts(&self, kmer: &Kmer) -> [u32; 4] {
        let mut out = [0u32; 4];
        for (i, c) in out.iter_mut().enumerate() {
            let mut x = kmer.clone();
            x.push_left(i as u8);
            *c = self.get_count(&x);
        }
        out
    }
}

/// Multi-word k-mer key (2 bits per base, base 0 in the low bits, up to
/// `2k` bits), mirroring BBTools `Kmer` long arrays for k > 64.
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct Kmer {
    k: usize,
    words: Vec<u64>,
}

impl Kmer {
    pub(crate) fn new(k: usize) -> Self {
        let n = (2 * k).div_ceil(64);
        Self {
            k,
            words: vec![0; n],
        }
    }

    pub(crate) fn reset(&mut self) {
        for w in &mut self.words {
            *w = 0;
        }
    }

    /// Bit index of base `i` (0-based from the 5' end).
    fn bit_pos(i: usize) -> (usize, u32) {
        let bit = 2 * i;
        (bit / 64, (bit % 64) as u32)
    }

    pub(crate) fn base_at(&self, i: usize) -> u8 {
        debug_assert!(i < self.k);
        let (w, b) = Self::bit_pos(i);
        ((self.words[w] >> b) & 3) as u8
    }

    fn set_base(&mut self, i: usize, x: u8) {
        debug_assert!(i < self.k);
        let (w, b) = Self::bit_pos(i);
        self.words[w] = (self.words[w] & !(3u64 << b)) | ((x as u64) << b);
    }

    /// Shift left by one base (2 bits) and append `x` at the 3' end.
    pub(crate) fn push_right(&mut self, x: u8) {
        // Shift the whole 2k-bit window left by 2; only the highest word is
        // masked so bits beyond the window are dropped.
        let n = self.words.len();
        // Number of valid bits in the highest word: 2k mod 64 (64 when the
        // window exactly fills whole words).
        let top_bits = (2 * self.k) % 64;
        let top_mask = if top_bits == 0 {
            u64::MAX
        } else {
            (1u64 << top_bits) - 1
        };
        let mut carry = 0u64;
        for (i, w) in self.words.iter_mut().enumerate() {
            let old = *w;
            let shifted = (old << 2) | carry;
            *w = if i == n - 1 {
                shifted & top_mask
            } else {
                shifted
            };
            carry = old >> 62;
        }
        self.words[0] |= (x as u64) & 3;
    }

    /// Shift right by one base (2 bits) and prepend `x` at the 5' end.
    pub(crate) fn push_left(&mut self, x: u8) {
        // Shift the whole 2k-bit window right by 2: the low 2 bits of word
        // i+1 become the top 2 bits of word i, and the low 2 bits of word 0
        // are discarded. Then the new base is placed at the window top.
        let n = self.words.len();
        for i in 0..n - 1 {
            self.words[i] = (self.words[i] >> 2) | (self.words[i + 1] << 62);
        }
        self.words[n - 1] >>= 2;
        let (w, b) = Self::bit_pos(self.k - 1);
        self.words[w] = (self.words[w] & !(3u64 << b)) | (((x as u64) & 3) << b);
    }

    /// Reverse complement, mirroring `rc_kmer` over the base sequence.
    pub(crate) fn rc(&self) -> Self {
        let mut r = Self::new(self.k);
        for i in 0..self.k {
            let x = self.base_at(i);
            r.set_base(self.k - 1 - i, 3 - x);
        }
        r
    }

    /// Lexicographic comparison of the base sequences (5' to 3').
    pub(crate) fn cmp_bases(&self, other: &Self) -> std::cmp::Ordering {
        for i in (0..self.k).rev() {
            match self.base_at(i).cmp(&other.base_at(i)) {
                std::cmp::Ordering::Equal => {}
                o => return o,
            }
        }
        std::cmp::Ordering::Equal
    }

    /// Canonical key (lexicographically smaller of forward / reverse-complement).
    pub(crate) fn canonical(&self) -> Self {
        let r = self.rc();
        if r.cmp_bases(self).is_lt() {
            r
        } else {
            self.clone()
        }
    }
}

/// 2-bit base code, mirroring `AminoAcid.baseToNumber`: A=0, C=1, G=2,
/// T/U=3, and -1 for everything else (N and ambiguity reset the k-mer window).
pub(crate) fn base_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 0,
    }
}

/// Reverse-complement code, mirroring `AminoAcid.baseToComplementNumber`:
/// A=3, C=2, G=1, T/U=0, and -1 for everything else (including N).
fn base_comp_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 3,
        b'C' | b'c' => 2,
        b'G' | b'g' => 1,
        b'T' | b't' | b'U' | b'u' => 0,
        _ => 0,
    }
}

/// `AminoAcid.baseToNumber >= 0`: A/C/G/T/U count as defined (baseToNumber
/// is filled with -1 and only ACGTU are overwritten).
pub(crate) fn base_defined(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'U' | b'u'
    )
}

/// `QualityTools.PROB_ERROR`: phred quality to error probability.
fn prob_error(q: u8) -> f32 {
    match q {
        0 => 0.75,
        1 => 0.7,
        _ => (10f64.powf(-0.1 * q as f64)) as f32,
    }
}

/// `QualityTools.PROB_CORRECT` and `PROB_CORRECT_INVERSE`, precomputed like
/// the Java arrays so the sliding `minprob` product uses exactly the same
/// float operations (multiply by a precomputed inverse, never divide).
fn prob_tables() -> (Vec<f32>, Vec<f32>) {
    let mut correct = Vec::with_capacity(128);
    let mut inverse = Vec::with_capacity(128);
    for q in 0..128u16 {
        let c = 1.0 - prob_error(q as u8);
        correct.push(c);
        inverse.push(1.0 / c);
    }
    (correct, inverse)
}

/// Counts the k-mers of one read, mirroring `KmerTableSetU.addKmersToTable`
/// (canonical keys, sliding `minprob` quality gate, N resets the window).
fn count_read_kmers(
    map: &mut HashMap<Kmer, u32>,
    bases: &[u8],
    quals: &[u8],
    k: usize,
    min_prob: f32,
    prob_correct: &[f32],
    prob_correct_inv: &[f32],
) {
    if bases.len() < k {
        return;
    }
    let min_prob2 = if min_prob > 0.0 && !quals.is_empty() {
        min_prob
    } else {
        0.0
    };
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    let mut prob = 1f32;
    for (i, &b) in bases.iter().enumerate() {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            if min_prob2 > 0.0 {
                prob *= prob_correct[quals[i] as usize];
                if len >= k {
                    let oldq = quals[i - k];
                    prob *= prob_correct_inv[oldq as usize];
                }
            }
            len += 1;
        } else {
            len = 0;
            kmer.reset();
            prob = 1.0;
        }
        if len >= k && prob >= min_prob2 {
            *map.entry(kmer.canonical()).or_insert(0) += 1;
        }
    }
}

/// The forward k-mers of a read (position-wise, `None` for invalid windows),
/// mirroring `KmerTableSet.fillKmers`.
fn fill_kmers(bases: &[u8], k: usize) -> Vec<Option<Kmer>> {
    let mut out = Vec::with_capacity(bases.len().saturating_sub(k - 1));
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    let min = k - 1;
    for (i, &b) in bases.iter().enumerate() {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if i >= min {
            if len >= k {
                out.push(Some(kmer.clone()));
            } else {
                out.push(None);
            }
        }
    }
    out
}

/// Fill the per-window canonical k-mer counts of a read.
fn fill_counts(kmers: &[Option<Kmer>], table: &TadpoleTable) -> Vec<i64> {
    kmers
        .iter()
        .map(|k| {
            if let Some(kmer) = k {
                raw_count(kmer, table)
            } else {
                0
            }
        })
        .collect()
}

/// Raw table count, mirroring Java `getCount`: -1 when the k-mer is absent.
fn raw_count(kmer: &Kmer, table: &TadpoleTable) -> i64 {
    let c = table.get_count(kmer);
    if c == 0 {
        -1
    } else {
        c as i64
    }
}

/// `KmerTableSet.regenerateCounts`: recompute window counts starting at `ca`
/// after a base change, resetting at undefined bases (count 0 for invalid
/// windows, raw -1 for absent k-mers).
fn regenerate_counts(bases: &[u8], counts: &mut [i64], table: &TadpoleTable, k: usize, ca: usize) {
    let b = ca + k - 1;
    let lim = bases.len().min(b + k + 1);
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    for (j, &base) in bases[ca..lim].iter().enumerate() {
        let i = ca + j;
        if base_defined(base) {
            kmer.push_right(base_code(base));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if i >= b {
            let idx = i + 1 - k;
            if len >= k {
                counts[idx] = raw_count(&kmer, table);
            } else {
                counts[idx] = 0;
            }
        }
    }
}

/// Error-correction statistics (thread-local equivalents in BBTools).
#[derive(Debug, Default, Clone)]
pub struct ErrorTracker {
    pub suspected: usize,
    pub detected_reassemble: usize,
    pub corrected_reassemble_inner: usize,
    pub corrected_reassemble_outer: usize,
    pub rollback: bool,
}

impl ErrorTracker {
    fn corrected(&self) -> usize {
        self.corrected_reassemble_inner + self.corrected_reassemble_outer
    }
}

/// `isSimilar`: two k-mer depths are similar within absolute/fractional tolerances.
fn is_similar(a: i64, b: i64, opts: &TadpoleOptions) -> bool {
    let min = a.min(b);
    let max = a.max(b);
    let dif = max - min;
    (dif as f32) < opts.path_similarity_constant as f32
        || (dif as f32) < (max as f32) * opts.path_similarity_fraction
}

/// `isError(high, low)` (errorPath=1).
fn is_error2(high: i64, low: i64, opts: &TadpoleOptions) -> bool {
    let em1 = opts.error_mult1;
    (low as f32) * em1 < high as f32
        || (low <= opts.error_lower_const as i64
            && (high as f32) >= (opts.min_count_correct as f32).max(low as f32 * opts.error_mult2))
}

/// `isError(high, low, q)` (errorPath=1, quality-weighted).
fn is_error3(high: i64, low: i64, q: u8, opts: &TadpoleOptions) -> bool {
    let em1 = opts.error_mult1 * (1.0 + q as f32 * opts.error_mult_q_factor);
    (low as f32) * em1 < high as f32
        || (low <= opts.error_lower_const as i64
            && (high as f32) >= (opts.min_count_correct as f32).max(low as f32 * opts.error_mult2))
}

/// `isErrorBidirectional`.
fn is_error_bidirectional(a: i64, b: i64, qa: u8, qb: u8, opts: &TadpoleOptions) -> bool {
    if a >= b {
        is_error3(a, b, qb, opts)
    } else {
        is_error3(b, a, qa, opts)
    }
}

/// `isSubstitution`: isolated 1bp substitution candidate.
fn is_substitution(
    ca: usize,
    error_extension: usize,
    qb: u8,
    counts: &[i64],
    k: usize,
    opts: &TadpoleOptions,
) -> bool {
    let cb = ca + 1;
    let a_count = counts[ca];
    let b_count = counts[cb];
    if is_error3(a_count, b_count, qb, opts)
        && similar_range(
            a_count,
            ca as isize - error_extension as isize,
            ca as isize - 1,
            counts,
            opts,
        )
        && error_range(a_count, ca + 2, ca + k, counts, opts)
    {
        let cc = ca + k;
        let cd = cc + 1;
        if cd < counts.len() {
            let c_count = counts[cc];
            let d_count = counts[cd];
            is_error2(a_count, d_count, opts) || is_error3(d_count, c_count, qb, opts)
        } else {
            true
        }
    } else {
        false
    }
}

fn similar_range(a: i64, loc1: isize, loc2: isize, counts: &[i64], opts: &TadpoleOptions) -> bool {
    if loc2 < 0 {
        // Java clamps loc2 to -1 and the loop body never runs (empty range).
        return true;
    }
    let lo = loc1.max(0) as usize;
    let hi = (loc2 as usize).min(counts.len() - 1);
    if lo > hi {
        return true;
    }
    counts[lo..=hi].iter().all(|&c| is_similar(a, c, opts))
}

fn error_range(a: i64, loc1: usize, loc2: usize, counts: &[i64], opts: &TadpoleOptions) -> bool {
    let hi = loc2.min(counts.len() - 1);
    if loc1 > hi {
        return true;
    }
    counts[loc1..=hi].iter().all(|&c| is_error2(a, c, opts))
}

/// `countErrors`: count error positions, skipping `k` after each hit.
fn count_errors(counts: &[i64], quals: Option<&[u8]>, k: usize, opts: &TadpoleOptions) -> usize {
    let mut possible = 0usize;
    let mut i = 1usize;
    while i < counts.len() {
        let (a, b) = (counts[i - 1], counts[i]);
        let error = match quals {
            Some(q) => is_error_bidirectional(a, b, q[i - 1], q[i + k - 1], opts),
            None => is_error_bidirectional(a, b, 20, 20, opts),
        };
        if error {
            possible += 1;
            i += k;
        } else {
            i += 1;
        }
    }
    possible
}

/// `hasErrorsFast`: sampled k-mer depth screen for likely errors.
fn has_errors_fast(kmers: &[Option<Kmer>], table: &TadpoleTable, opts: &TadpoleOptions) -> bool {
    if kmers.is_empty() {
        return false;
    }
    let incr = (opts.k / 2).clamp(1, 9);
    let mcc = opts.min_count_correct as i64;
    let mut prev = -1i64;
    let mut i = 0usize;
    while i < kmers.len() {
        let count = match &kmers[i] {
            Some(kmer) => raw_count(kmer, table),
            None => return true,
        };
        let min = count.min(prev);
        let max = count.max(prev);
        if count < mcc || (i > 0 && is_error2(max + 1, min - 1, opts)) {
            return true;
        }
        prev = count;
        i += incr;
    }
    if let Some(kmer) = kmers.last() {
        let count = match kmer {
            Some(kmer) => raw_count(kmer, table),
            None => return true,
        };
        let min = count.min(prev);
        let max = count.max(prev);
        return count < mcc || is_error2(max + 1, min - 1, opts);
    }
    false
}

/// `isJunction(max, second)` with branch-resolution thresholds.
pub(crate) fn is_junction(max: u32, second: u32, opts: &TadpoleOptions) -> bool {
    if second < 1
        || (second as f32) * opts.branch_mult1 < max as f32
        || (second <= opts.branch_lower_const as u32
            && (max as f32)
                >= (opts.min_count_extend as f32).max(second as f32 * opts.branch_mult2))
    {
        return false;
    }
    true
}

/// Extends a sequence to the right by at most `distance` bases, mirroring
/// `Tadpole1.extendToRight2`. Returns the number of bases added.
#[allow(clippy::too_many_arguments)]
fn extend_to_right2(
    bases: &mut Vec<u8>,
    table: &TadpoleTable,
    opts: &TadpoleOptions,
    distance: usize,
    include_junction_base: bool,
    use_left: bool,
) -> usize {
    let k = opts.k;
    let initial = bases.len();
    if initial < k {
        return 0;
    }
    // Build the rightmost k-mer.
    let mut kmer = Kmer::new(k);
    let mut rkmer = Kmer::new(k);
    let mut len = 0usize;
    for &b in &bases[initial - k..initial] {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            rkmer.push_left(base_comp_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
            rkmer.reset();
        }
    }
    if len < k {
        return 0;
    }
    let count = table.get_count(&kmer);
    if count < opts.min_count_seed as u32 {
        return 0;
    }

    let mut left_max_pos = 0usize;
    let mut left_max = opts.min_count_extend as u32;
    let mut left_second = 0u32;
    if use_left {
        let lc = table.fill_left_counts(&kmer);
        left_max_pos = argmax2(&lc, &mut left_max);
        left_second = lc[second_highest_position(&lc)];
    }

    let rc = table.fill_right_counts(&kmer);
    let mut right_max = 0u32;
    let mut right_max_pos = argmax2(&rc, &mut right_max);
    let mut right_second_pos = second_highest_position(&rc);
    let mut right_second = rc[right_second_pos];

    if right_max < opts.min_count_extend as u32 {
        return 0;
    }
    if is_junction(right_max, right_second, opts)
        || (use_left && is_junction(left_max, left_second, opts))
    {
        return 0;
    }

    let max_len = initial + distance;
    let mut added = 0usize;
    // Tadpole1 (k<=31) appends the junction base when the forward k-mer is
    // the canonical maximum; Tadpole2 (k>31) canonicalizes to the minimum,
    // so the condition flips.
    let canonical_is_rc = k > 31;
    while bases.len() < max_len {
        let b = right_max_pos as u8;
        let x = right_max_pos as u8;
        let x2 = 3 - x;
        let evicted = kmer.base_at(k - 1);
        kmer.push_right(x);
        rkmer.push_left(x2);

        if use_left {
            let lc = table.fill_left_counts(&kmer);
            left_max_pos = argmax2(&lc, &mut left_max);
            left_second = lc[second_highest_position(&lc)];
        }
        let rc = table.fill_right_counts(&kmer);
        right_max_pos = argmax2(&rc, &mut right_max);
        right_second_pos = second_highest_position(&rc);
        right_second = rc[right_second_pos];

        let junc_r = is_junction(right_max, right_second, opts);
        let junc_l = use_left && is_junction(left_max, left_second, opts);
        // Tadpole2 (k>31) appends the junction base when the k-mer's
        // canonical orientation is the forward one (`key()==array1` in
        // BBTools; the reverse-complement key is the other branch).
        let kmer_is_rc = kmer.cmp_bases(&rkmer).is_lt();
        if junc_r || junc_l {
            if include_junction_base
                && if canonical_is_rc {
                    kmer_is_rc
                } else {
                    !kmer_is_rc
                }
            {
                bases.push(number_to_base(b));
                added += 1;
            }
            break;
        }
        if use_left && left_max_pos != evicted as usize {
            if include_junction_base
                && if canonical_is_rc {
                    kmer_is_rc
                } else {
                    !kmer_is_rc
                }
            {
                bases.push(number_to_base(b));
                added += 1;
            }
            break;
        }
        bases.push(number_to_base(b));
        added += 1;
        if right_max < opts.min_count_extend as u32 {
            break;
        }
    }
    added
}

pub(crate) fn number_to_base(n: u8) -> u8 {
    match n {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

pub(crate) fn argmax2(a: &[u32; 4], max: &mut u32) -> usize {
    let mut pos = 0usize;
    *max = a[0];
    for (i, &x) in a.iter().enumerate().skip(1) {
        if x > *max {
            *max = x;
            pos = i;
        }
    }
    pos
}

pub(crate) fn second_highest_position(a: &[u32; 4]) -> usize {
    let (mut p, mut p2) = if a[0] >= a[1] { (0, 1) } else { (1, 0) };
    for i in 2..a.len() {
        let x = a[i];
        if x > a[p2] {
            if x >= a[p] {
                p2 = p;
                p = i;
            } else {
                p2 = i;
            }
        }
    }
    p2
}

/// `isJunk`: read cannot be used for assembly.
pub fn is_junk(bases: &[u8], table: &TadpoleTable, opts: &TadpoleOptions, paired: bool) -> bool {
    let k = opts.k;
    let blen = bases.len();
    if blen < k {
        return true;
    }
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    for &b in &bases[..k] {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
    }
    if len >= k {
        let lc = table.fill_left_counts(&kmer);
        let max_pos = argmax2(&lc, &mut 0);
        if lc[max_pos] > 0 {
            return false;
        }
    }
    let mut max_depth = 0u32;
    for &b in &bases[k..] {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if len < k {
            continue;
        }
        {
            let depth = table.get_count(&kmer);
            if depth > max_depth {
                max_depth = depth;
                if max_depth > 1 && (!paired || max_depth > 2) {
                    return false;
                }
            }
        }
    }
    if len >= k && !paired {
        let rc = table.fill_right_counts(&kmer);
        let max_pos = argmax2(&rc, &mut 0);
        if rc[max_pos] > 0 {
            return false;
        }
    }
    true
}

/// `hasKmersAtOrBelow`: does the read have enough low-depth k-mers to toss?
pub fn has_kmers_at_or_below(
    bases: &[u8],
    table: &TadpoleTable,
    opts: &TadpoleOptions,
    too_low: u32,
    fraction: f32,
) -> bool {
    let k = opts.k;
    let blen = bases.len();
    if blen < k {
        return true;
    }
    let mut kmer = Kmer::new(k);
    let limit = ((blen - k + 1) as f32 * fraction).round().max(1.0) as usize;
    let mut valid = 0usize;
    let mut invalid = 0usize;
    let mut len = 0usize;
    for &b in bases.iter() {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if len >= k {
            let depth = table.get_count(&kmer);
            if depth > too_low {
                valid += 1;
            } else {
                invalid += 1;
                if invalid >= limit {
                    return true;
                }
            }
        }
    }
    let limit2 = ((valid + invalid) as f32 * fraction).round().max(1.0) as usize;
    valid < 1 || invalid >= limit2
}

/// `Read.expectedErrors` (phred qualities, countUndefined=true).
pub fn expected_errors(quals: &[u8]) -> f32 {
    quals.iter().map(|&q| prob_error(q)).sum()
}

/// Error-corrects one read in place (reassemble-only path), mirroring
/// `Tadpole1.errorCorrect` + `Tadpole.reassemble`. Returns corrections applied.
pub fn error_correct(
    bases: &mut Vec<u8>,
    quals: &mut Vec<u8>,
    table: &TadpoleTable,
    opts: &TadpoleOptions,
    tracker: &mut ErrorTracker,
) -> usize {
    tracker.suspected = 0;
    tracker.detected_reassemble = 0;
    tracker.corrected_reassemble_inner = 0;
    tracker.corrected_reassemble_outer = 0;
    tracker.rollback = false;

    let kmers = fill_kmers(bases, opts.k);
    let valid = kmers.len();
    if valid < 2 {
        return 0;
    }
    let has_undefined = bases.iter().any(|&b| !base_defined(b));
    if !has_undefined && !has_errors_fast(&kmers, table, opts) {
        return 0;
    }
    let mut counts = fill_counts(&kmers, table);
    // FASTA input has no qualities; BBTools substitutes a fixed quality 20.
    let qs = if quals.is_empty() {
        None
    } else {
        Some(quals.as_slice())
    };
    let possible_errors = count_errors(&counts, qs, opts.k, opts);
    tracker.suspected = possible_errors;
    let expected = expected_errors(quals);
    let counts0 = counts.clone();
    let bases0 = bases.clone();
    let quals0 = quals.clone();

    let corrected = reassemble(
        bases,
        quals,
        table,
        opts,
        &mut counts,
        tracker,
        opts.error_extension_reassemble,
    );
    debug_assert_eq!(corrected, tracker.corrected());

    if opts.ecc_rollback && (tracker.corrected() > 0 || tracker.rollback) {
        if !tracker.rollback && tracker.corrected() > 3 {
            let mult = (0.5f32 * (0.5 + 0.01 * bases.len() as f32)).max(1.0);
            let ce = count_errors(&counts, Some(quals), opts.k, opts);
            let c1 = ce > 0 && tracker.corrected() as f32 > mult + expected;
            let c2 = tracker.corrected() as f32 > 2.5 * mult + expected;
            if c1 || c2 {
                tracker.rollback = true;
            }
        }
        if !tracker.rollback {
            for i in 0..counts.len() {
                // Java clamps both sides to 0 before the rollback comparison.
                let a = counts0[i].max(0);
                let b = counts[i].max(0);
                if b < a - 1 && !is_similar(a, b, opts) {
                    tracker.rollback = true;
                }
            }
        }
        if tracker.rollback {
            *bases = bases0;
            *quals = quals0;
            tracker.corrected_reassemble_inner = 0;
            tracker.corrected_reassemble_outer = 0;
            return 0;
        }
    }
    tracker.corrected()
}

/// `reassemble`: multi-pass local reassembly error correction.
#[allow(clippy::too_many_arguments)]
fn reassemble(
    bases: &mut [u8],
    quals: &mut [u8],
    table: &TadpoleTable,
    opts: &TadpoleOptions,
    counts: &mut Vec<i64>,
    tracker: &mut ErrorTracker,
    error_extension: usize,
) -> usize {
    if bases.len() < opts.k + 1 + opts.dead_zone {
        return 0;
    }
    let mut corrected = 0usize;
    let mut corrected_incr;
    let mut detected_incr;
    let mut uncorrected;
    let detected0 = tracker.detected_reassemble;
    corrected_incr = reassemble_pass(bases, quals, table, opts, counts, tracker, error_extension);
    corrected += corrected_incr;
    detected_incr = tracker.detected_reassemble - detected0;
    uncorrected = detected_incr.saturating_sub(corrected_incr);
    let mut passes = 1usize;
    while passes < 6 && corrected_incr > 0 && uncorrected > 0 {
        tracker.detected_reassemble -= uncorrected;
        let detected0 = tracker.detected_reassemble;
        corrected_incr =
            reassemble_pass(bases, quals, table, opts, counts, tracker, error_extension);
        corrected += corrected_incr;
        detected_incr = tracker.detected_reassemble - detected0;
        uncorrected = detected_incr.saturating_sub(corrected_incr);
        passes += 1;
    }
    corrected
}

/// `reassemble_pass`: forward + reverse passes, window filtering, consensus.
#[allow(clippy::too_many_arguments)]
fn reassemble_pass(
    bases: &mut [u8],
    quals: &mut [u8],
    table: &TadpoleTable,
    opts: &TadpoleOptions,
    counts: &mut Vec<i64>,
    tracker: &mut ErrorTracker,
    error_extension: usize,
) -> usize {
    if bases.len() < opts.k + 1 + opts.dead_zone {
        return 0;
    }
    let mut from_left = bases.to_vec();
    let mut from_right = bases.to_vec();
    let mut counts2 = counts.clone();
    reassemble_inner(
        &mut from_left,
        quals,
        table,
        opts,
        &mut counts2,
        error_extension,
    );

    from_right = rev_comp(&from_right).collect();
    let qr: Vec<u8> = quals.iter().rev().copied().collect();
    counts2 = counts.clone();
    counts2.reverse();
    reassemble_inner(
        &mut from_right,
        &qr,
        table,
        opts,
        &mut counts2,
        error_extension,
    );
    from_right = rev_comp(&from_right).collect();

    let mut corrected_inner = 0usize;
    let mut corrected_outer = 0usize;
    let mut detected_inner = 0usize;
    let mut detected_outer = 0usize;
    let mut rollback = false;
    for i in 0..bases.len() {
        let a = bases[i];
        let b = from_left[i];
        let c = from_right[i];
        if a != b || a != c {
            if b == c {
                detected_inner += 1;
            } else {
                detected_outer += 1;
                if a != b && a != c {
                    rollback = true;
                }
            }
        }
        if b == a {
            from_left[i] = 0;
        }
        if c == a {
            from_right[i] = 0;
        }
    }
    let detected = detected_inner + detected_outer;
    tracker.detected_reassemble += detected;
    if rollback || detected == 0 {
        return 0;
    }

    clear_window2(&mut from_left, quals, opts);
    // Java clears fromRight while it is in reversed orientation with the
    // reversed qualities; clearing the forward-oriented copy with reversed
    // qualities would mis-weight the two read ends.
    from_right.reverse();
    clear_window2(&mut from_right, &qr, opts);
    from_right.reverse();

    for i in 0..bases.len() {
        let a = bases[i];
        let b = from_left[i];
        let c = from_right[i];
        let mut d = a;
        if b == 0 && c == 0 {
            // nothing
        } else if b == c {
            d = b;
        } else if b == 0 {
            d = c;
        } else if c == 0 {
            d = b;
        } else if b != c {
            // keep a
        }
        if opts.ecc_require_bidirectional && b != c && i >= opts.k && i < bases.len() - opts.k {
            d = a;
        }
        if d != a {
            let mut q = if quals.is_empty() { 30 } else { quals[i] };
            if b == c {
                corrected_inner += 1;
                q = q.saturating_add(8).clamp(24, 32);
            } else {
                corrected_outer += 1;
                q = q.saturating_add(4).clamp(20, 28);
            }
            if !rollback {
                bases[i] = d;
                if !quals.is_empty() {
                    quals[i] = q;
                }
            }
        }
    }
    if rollback && corrected_inner + corrected_outer > 0 {
        tracker.rollback = true;
        return 0;
    }
    tracker.corrected_reassemble_inner += corrected_inner;
    tracker.corrected_reassemble_outer += corrected_outer;
    let corrected = corrected_inner + corrected_outer;
    if corrected > 0 {
        // Regenerate counts for all windows.
        let kmers = fill_kmers(bases, opts.k);
        *counts = fill_counts(&kmers, table);
    }
    corrected
}

/// `reassemble_inner`: per-position substitution detection and correction.
#[allow(clippy::too_many_arguments)]
fn reassemble_inner(
    bases: &mut [u8],
    quals: &[u8],
    table: &TadpoleTable,
    opts: &TadpoleOptions,
    counts: &mut [i64],
    error_extension: usize,
) -> usize {
    let k = opts.k;
    let length = bases.len();
    if length < k + 1 + opts.dead_zone {
        return 0;
    }
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    let mut corrected = 0usize;
    let lim = length - opts.dead_zone - 1;
    for a in 0..lim {
        if base_defined(bases[a]) {
            kmer.push_right(base_code(bases[a]));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if len >= k {
            let b = a + 1;
            // len>=k implies a+1>=k, so this cannot underflow.
            let ca = a + 1 - k;
            let a_count = counts[ca];
            let qb = if quals.is_empty() { 20 } else { quals[b] };
            if is_substitution(ca, error_extension, qb, counts, k, opts) {
                let rc = table.fill_right_counts(&kmer);
                let right_max_pos = argmax2(&rc, &mut 0);
                let right_max = rc[right_max_pos];
                let right_second_pos = second_highest_position(&rc);
                let right_second = rc[right_second_pos];
                let base = bases[b];
                // Java `baseToNumber` is -1 for N, so an N never matches the
                // preferred extension and always goes through the correction.
                let num = if base_defined(base) {
                    base_code(base) as i64
                } else {
                    -1
                };
                if right_max >= opts.min_count_extend as u32 {
                    // BBTools compares the base code to the *count* here
                    // (`if(num==rightMax)`), not to the position index; the
                    // base is treated as already-correct when they coincide.
                    if num == right_max as i64 {
                    } else if (is_error3(right_max as i64, right_second as i64, qb, opts)
                        || !is_junction(right_max, right_second, opts))
                        && is_similar(a_count, right_max as i64, opts)
                    {
                        bases[b] = number_to_base(right_max_pos as u8);
                        corrected += 1;
                        // Regenerate counts for windows ca+1..=ca+k (those
                        // containing the changed base at ca+k).
                        regenerate_counts(bases, counts, table, k, ca);
                    }
                }
            }
        }
    }
    corrected
}

/// `clearWindow2`: sliding-window quality filter over correction candidates.
fn clear_window2(bb: &mut [u8], quals: &[u8], opts: &TadpoleOptions) -> usize {
    let len = bb.len();
    let window = opts.window_len as isize;
    let mut cleared = 0usize;
    let mut count = 0usize;
    let mut qsum = 0usize;
    for (i, prev) in (0..len as isize).zip((-window)..) {
        let b = bb[i as usize];
        if b != 0 && (quals.is_empty() || quals[i as usize] > 0) {
            count += 1;
            if !quals.is_empty() {
                qsum += quals[i as usize] as usize;
            }
            if count > opts.window_count || qsum > opts.window_qual_sum {
                let start = (i - window).max(0) as usize;
                for b in &mut bb[start..] {
                    if *b != 0 {
                        *b = 0;
                        cleared += 1;
                    }
                }
                return cleared;
            }
        }
        if prev >= 0 && bb[prev as usize] > 0 && (quals.is_empty() || quals[prev as usize] > 0) {
            count -= 1;
            if !quals.is_empty() {
                qsum -= quals[prev as usize] as usize;
            }
        }
    }
    cleared
}

/// Extends one read in place (both ends), mirroring `processRead` extension.
pub fn extend_read(
    bases: &mut Vec<u8>,
    quals: &mut Vec<u8>,
    table: &TadpoleTable,
    opts: &TadpoleOptions,
    numeric_id: u64,
) -> usize {
    let mut extension_right = 0usize;
    let mut extension_left = 0usize;
    if opts.extend_right > 0 {
        extension_right = extend_read_one_side(bases, quals, table, opts, opts.extend_right);
    }
    if opts.extend_left > 0 {
        // Reverse-complement, extend, reverse back.
        let mut rc: Vec<u8> = rev_comp(bases).collect();
        let mut rq: Vec<u8> = quals.iter().rev().copied().collect();
        extension_left = extend_read_one_side(&mut rc, &mut rq, table, opts, opts.extend_left);
        *bases = rev_comp(&rc).collect();
        *quals = rq.iter().rev().copied().collect();
    }
    let mut extension = extension_right + extension_left;
    if opts.extension_rollback > 0 {
        let mut left_mod = 0usize;
        let mut right_mod = 0usize;
        if extension_left > 0 && extension_left < opts.extend_left {
            left_mod =
                extension_left.min((numeric_id % (opts.extension_rollback as u64 + 1)) as usize);
            extension_left -= left_mod;
        }
        if extension_right > 0 && extension_right < opts.extend_right {
            right_mod =
                extension_right.min((numeric_id % (opts.extension_rollback as u64 + 1)) as usize);
            extension_right -= right_mod;
        }
        if left_mod > 0 || right_mod > 0 {
            // Trim left_mod bases from the 5' end and right_mod from 3'.
            let keep_from = left_mod.min(bases.len());
            let keep_to = bases.len().saturating_sub(right_mod);
            if keep_from < keep_to {
                *bases = bases[keep_from..keep_to].to_vec();
                if !quals.is_empty() {
                    *quals = quals[keep_from..keep_to].to_vec();
                }
            } else {
                bases.clear();
                quals.clear();
            }
        }
        extension = extension_left + extension_right;
    }
    extension
}

/// Extends one read's 3' end in place, mirroring `Tadpole.extendToRight2`
/// as called by BBMerge (`extendAndMerge` / `extendRead`): right junction
/// only, no left-branch check.
pub fn extend_read_right(
    bases: &mut Vec<u8>,
    quals: &mut Vec<u8>,
    table: &TadpoleTable,
    opts: &TadpoleOptions,
    distance: usize,
) -> usize {
    // BBMerge `extendRead` calls `extendToRight2(..., false)` (no junction
    // base); the standalone `fq extend` path keeps the junction base.
    let initial = bases.len();
    if initial < opts.k {
        return 0;
    }
    let added = extend_to_right2(bases, table, opts, distance, false, false);
    if added > 0 && !quals.is_empty() {
        quals.resize(bases.len(), 30);
    }
    added
}

/// Extends one read end (3' after any RC flip) by up to `distance` bases.
fn extend_read_one_side(
    bases: &mut Vec<u8>,
    quals: &mut Vec<u8>,
    table: &TadpoleTable,
    opts: &TadpoleOptions,
    distance: usize,
) -> usize {
    let initial = bases.len();
    if initial < opts.k {
        return 0;
    }
    // BBTools never initializes the left-counts buffer for read extension
    // (`ExtendThread.leftCounts` stays null), so only the right junction is
    // considered; the left-branch check is disabled.
    let added = extend_to_right2(bases, table, opts, distance, true, false);
    if added > 0 && !quals.is_empty() {
        quals.resize(bases.len(), 30);
    }
    added
}

/// Per-read processing outcome counters (subset of tadpole.sh stats).
#[derive(Debug, Default, Clone)]
pub struct TadpoleStats {
    pub reads_in: u64,
    pub bases_in: u64,
    pub bases_extended: u64,
    pub reads_extended: u64,
    pub reads_corrected: u64,
    pub bases_corrected: u64,
    pub reads_detected: u64,
    pub bases_detected: u64,
    pub reads_fully_corrected: u64,
    pub reads_discarded: u64,
    pub bases_discarded: u64,
    pub rollbacks: u64,
}

/// Main per-read processing pipeline, mirroring `ExtendThread.processRead`.
#[allow(clippy::too_many_arguments)]
pub fn process_read(
    bases: &mut Vec<u8>,
    quals: &mut Vec<u8>,
    table: &TadpoleTable,
    opts: &TadpoleOptions,
    stats: &mut TadpoleStats,
    numeric_id: u64,
    mate: Option<usize>,
    discard_mate: &mut bool,
) -> bool {
    let initial_len = bases.len();
    let mut tracker = ErrorTracker::default();
    if opts.ecc && (opts.toss_uncorrectable || opts.ecc_rollback) {
        let corrected = error_correct(bases, quals, table, opts, &mut tracker);
        if tracker.rollback {
            stats.rollbacks += 1;
        }
        let detected = tracker.detected_reassemble;
        if detected > 0 {
            stats.reads_detected += 1;
            stats.bases_detected += detected as u64;
            if corrected > 0 {
                stats.reads_corrected += 1;
                stats.bases_corrected += corrected as u64;
            }
            if corrected == detected
                || (corrected > 0 && count_errors_from(bases, quals, table, opts) == 0)
            {
                stats.reads_fully_corrected += 1;
            } else if opts.toss_uncorrectable {
                if mate.is_some() && !opts.require_both_bad {
                    *discard_mate = true;
                }
                return true; // discard this read
            }
        }
    }

    if opts.toss_junk && is_junk(bases, table, opts, mate.is_some_and(|m| m >= opts.k)) {
        return true;
    }

    if opts.toss_depth >= 0
        && has_kmers_at_or_below(
            bases,
            table,
            opts,
            opts.toss_depth as u32,
            opts.low_depth_discard_fraction,
        )
    {
        if mate.is_some() && !opts.require_both_bad {
            *discard_mate = true;
        }
        return true;
    }

    if opts.extend_right > 0 || opts.extend_left > 0 {
        let ext = extend_read(bases, quals, table, opts, numeric_id);
        if ext > 0 {
            stats.bases_extended += ext as u64;
            stats.reads_extended += 1;
        }
    }
    stats.bases_in += initial_len as u64;
    false
}

fn count_errors_from(
    bases: &[u8],
    quals: &[u8],
    table: &TadpoleTable,
    opts: &TadpoleOptions,
) -> usize {
    let kmers = fill_kmers(bases, opts.k);
    let counts = fill_counts(&kmers, table);
    let qs = if quals.is_empty() { None } else { Some(quals) };
    count_errors(&counts, qs, opts.k, opts)
}

/// Applies the BBTools phred round-trip to a record's quality scores.
fn canonicalize_quality(rec: &mut SeqRecord) {
    if rec.quality_scores().is_empty() {
        return;
    }
    let seq = rec.sequence().to_vec();
    let raw = rec.quality_scores().to_vec();
    let phred = to_phred(&seq, &raw);
    rec.set_quality(from_phred(&phred));
}

/// Runs the tadpole correct/extend/discard pipeline over FASTQ input
/// (1 interleaved file or 2 paired files), mirroring `tadpole.sh`.
pub fn run<W: Write>(
    infiles: &[String],
    out: &mut W,
    opts: &TadpoleOptions,
) -> Result<TadpoleStats> {
    anyhow::ensure!(
        opts.k >= 1,
        "k-mer length must be at least 1, got {}",
        opts.k
    );
    // Pass 1: read all records into memory, canonicalizing qualities.
    let mut records: Vec<SeqRecord> = Vec::new();
    let mut reader1 = SeqReader::new(&infiles[0])?;
    let mut reader2 = if infiles.len() > 1 {
        Some(SeqReader::new(&infiles[1])?)
    } else {
        None
    };
    let mut rec = SeqRecord::new();
    loop {
        if !reader1.read_record(&mut rec)? {
            break;
        }
        canonicalize_quality(&mut rec);
        records.push(rec.clone());
        if let Some(r) = reader2.as_mut() {
            if !r.read_record(&mut rec)? {
                anyhow::bail!("unpaired trailing read in {}", infiles[0]);
            }
            canonicalize_quality(&mut rec);
            records.push(rec.clone());
        } else if !reader1.read_record(&mut rec)? {
            anyhow::bail!("unpaired trailing read in {}", infiles[0]);
        } else {
            canonicalize_quality(&mut rec);
            records.push(rec.clone());
        }
    }

    // Pass 2: count k-mers from the canonicalized (phred) qualities.
    let reads: Vec<(Vec<u8>, Vec<u8>)> = records
        .iter()
        .map(|r| {
            (
                r.sequence().to_vec(),
                to_phred(r.sequence(), r.quality_scores()),
            )
        })
        .collect();
    let table = TadpoleTable::build(&reads, opts.k, opts.min_prob);

    // Pass 3: process pairs and write surviving reads.
    let mut stats = TadpoleStats {
        reads_in: records.len() as u64,
        ..Default::default()
    };
    let mut i = 0usize;
    while i < records.len() {
        let r1 = records[i].clone();
        let r2 = if i + 1 < records.len() {
            Some(records[i + 1].clone())
        } else {
            None
        };
        // BBTools assigns one numeric ID per pair (both mates share it).
        let id = (i / 2) as u64;
        let mut bases1 = r1.sequence().to_vec();
        let mut quals1 = to_phred(&bases1, r1.quality_scores());
        let mut bases2 = r2
            .as_ref()
            .map(|r| r.sequence().to_vec())
            .unwrap_or_default();
        let mut quals2 = r2
            .as_ref()
            .map(|r| to_phred(r.sequence(), r.quality_scores()))
            .unwrap_or_default();
        let mut discard_mate = false;
        let d1 = process_read(
            &mut bases1,
            &mut quals1,
            &table,
            opts,
            &mut stats,
            id,
            r2.as_ref().map(|r| r.sequence().len()),
            &mut discard_mate,
        );
        let d2 = if discard_mate {
            true
        } else if r2.is_some() {
            let mate_len = bases1.len(); // r1 length after its own processing
            process_read(
                &mut bases2,
                &mut quals2,
                &table,
                opts,
                &mut stats,
                id,
                Some(mate_len),
                &mut discard_mate,
            )
        } else {
            true
        };
        // Either read's processing may discard the other as its mate
        // (tossdepth / tossuncorrectable without requireBothBad).
        let d1 = d1 || discard_mate;
        // A pair is dropped only when both reads are discarded; otherwise
        // both are written (discarded mates keep their processed state).
        if d1 && d2 {
            stats.reads_discarded += 1 + r2.is_some() as u64;
            stats.bases_discarded += bases1.len() as u64 + bases2.len() as u64;
        } else {
            write_record(out, &r1, &bases1, &from_phred(&quals1))?;
            if let Some(r) = r2.as_ref() {
                write_record(out, r, &bases2, &from_phred(&quals2))?;
            }
        }
        i += 2;
    }
    Ok(stats)
}

fn write_record<W: Write>(w: &mut W, rec: &SeqRecord, seq: &[u8], qual: &[u8]) -> Result<()> {
    let header = if rec.comment().is_empty() {
        rec.name().to_string()
    } else {
        format!("{} {}", rec.name(), rec.comment())
    };
    if qual.is_empty() {
        crate::libs::fmt::fq::write_fa(w, &header, seq)?;
    } else {
        crate::libs::fmt::fq::write_fq(w, &header, seq, qual)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_rc_is_identity() {
        // ACGT (0,1,2,3) RC is ACGT itself.
        let mut k = Kmer::new(4);
        for (i, &b) in b"ACGT".iter().enumerate() {
            k.set_base(i, base_code(b));
        }
        assert_eq!(k.rc().cmp_bases(&k), std::cmp::Ordering::Equal);
        // The canonical key is the lexicographically smaller orientation.
        let mut r = Kmer::new(4);
        for (i, &b) in b"GTCA".iter().enumerate() {
            r.set_base(i, base_code(b));
        }
        // RC(GTCA) = TGAC > GTCA, so canonical(r) must equal r itself.
        assert_eq!(r.canonical().cmp_bases(&r), std::cmp::Ordering::Equal);
        assert_eq!(r.canonical().cmp_bases(&r.rc()), std::cmp::Ordering::Less);
    }

    #[test]
    fn rolling_kmers_match_set_base_layout() {
        for k in [4usize, 31, 62, 81] {
            // Build "ACGT" repeated (truncated to k) by rolling push_right
            // and by set_base in the same orientation (base i of the window
            // occupies the low 2 bits after the window is full).
            let seq: Vec<u8> = (0..k).map(|i| b"ACGT"[i % 4]).collect();
            let mut direct = Kmer::new(k);
            for (i, &b) in seq.iter().rev().enumerate() {
                direct.set_base(i, base_code(b));
            }
            let mut rolled = Kmer::new(k);
            for &b in &seq {
                rolled.push_right(base_code(b));
            }
            eprintln!("k={k} rolled={:?} direct={:?}", rolled.words, direct.words);
            assert_eq!(
                rolled.cmp_bases(&direct),
                std::cmp::Ordering::Equal,
                "k={k}"
            );

            // Rolling rc (push complements forward) must equal rc() of the
            // rolled kmer: this is the invariant extend_to_right2 relies on.
            let mut rolled_rc = Kmer::new(k);
            for &b in &seq {
                rolled_rc.push_left(base_comp_code(b));
            }
            assert_eq!(
                rolled_rc.cmp_bases(&rolled.rc()),
                std::cmp::Ordering::Equal,
                "rc-invariant k={k}"
            );

            // Canonical key must equal the lexicographically smaller of the
            // forward and rc orientations (u128 reference for k <= 62).
            if k <= 62 {
                let f = kmer_to_u128(&rolled, k);
                let r = kmer_to_u128(&rolled.rc(), k);
                assert_eq!(
                    kmer_to_u128(&rolled.canonical(), k),
                    f.min(r),
                    "canonical k={k}"
                );
            }

            // One rolling push_left (window full) must drop the oldest base
            // and prepend the new one.
            let mut rl = rolled_rc.clone();
            rl.push_left(1);
            assert_eq!(rl.base_at(k - 1), 1, "push_left top k={k}");
            for i in 0..k - 1 {
                assert_eq!(
                    rl.base_at(i),
                    rolled_rc.base_at(i + 1),
                    "push_left shift k={k} i={i}"
                );
            }

            let mut rr = rolled.clone();
            rr.push_right(2);
            assert_eq!(rr.base_at(0), 2, "push_right bottom k={k}");
            for i in 1..k {
                assert_eq!(
                    rr.base_at(i),
                    rolled.base_at(i - 1),
                    "push_right shift k={k} i={i}"
                );
            }
        }
    }
}

#[cfg(test)]
fn kmer_to_u128(k: &Kmer, kmer_len: usize) -> u128 {
    let mut x = 0u128;
    for i in (0..kmer_len).rev() {
        x = (x << 2) | k.base_at(i) as u128;
    }
    x
}

#[test]
fn counting_matches_simple_expected() {
    let reads = vec![(b"ACGTACGT".to_vec(), vec![40; 8])];
    let table = TadpoleTable::build(&reads, 4, 0.5);
    // "ACGT" kmer and its RC (ACGT) are the same key; counts = 5 windows.
    assert_eq!(table.map.values().sum::<u32>(), 5);
}

#[test]
fn junk_detects_short_read() {
    let opts = TadpoleOptions::default();
    let table = TadpoleTable::build(&[], opts.k, opts.min_prob);
    assert!(is_junk(b"ACGT".as_ref(), &table, &opts, false));
}

#[test]
fn read36_left_extension_matches_golden() {
    // Reproduce the golden `fq extend` run on the committed subset and
    // check that read 36's left extension matches BBTools.
    let infile = "tests/bbtools/Lambda/golden/ecco_sub.fq.gz";
    let mut reader = SeqReader::new(infile).unwrap();
    let mut rec = SeqRecord::new();
    let mut reads: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    while reader.read_record(&mut rec).unwrap() {
        let seq = rec.sequence().to_vec();
        let quals = to_phred(&seq, rec.quality_scores());
        reads.push((seq, quals));
    }
    let k = 62usize;
    let table = TadpoleTable::build(&reads, k, 0.5);
    let opts = TadpoleOptions {
        k,
        extend_left: 20,
        extend_right: 20,
        ..TadpoleOptions::default()
    };
    let r36 = &reads[35].0;
    let mut bases = r36.clone();
    let mut quals = vec![30; bases.len()];
    // The `run` pipeline assigns one numeric ID per pair; read 36 is in
    // pair 17 (0-based), which drives the extension-rollback trim.
    extend_read(&mut bases, &mut quals, &table, &opts, 17);
    let seq = String::from_utf8_lossy(&bases).into_owned();
    // Golden: input + GTGGAA on the left + GAAGGCATTAACGCCTCTGC right.
    let golden = format!(
        "GTGGAA{}{}",
        String::from_utf8_lossy(r36),
        "GAAGGCATTAACGCCTCTGC"
    );
    assert_eq!(seq, golden, "read36 extension mismatch");
}

#[test]
fn k81_table_counts_match_bruteforce() {
    // The merge phase-4 extension uses k=81 (Tadpole2 long-k path);
    // verify the k-mer table counts against a brute-force scan.
    let infile = "tests/bbtools/Lambda/golden/ext_sub.fq.gz";
    let mut reader = SeqReader::new(infile).unwrap();
    let mut rec = SeqRecord::new();
    let mut reads: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    while reader.read_record(&mut rec).unwrap() {
        let seq = rec.sequence().to_vec();
        let quals = to_phred(&seq, rec.quality_scores());
        reads.push((seq, quals));
    }
    let k = 81usize;
    let table = TadpoleTable::build(&reads, k, 0.5);
    // Probe: the last 81 bases of the first read.
    let seq = &reads[0].0;
    let mut probe = Kmer::new(k);
    for &b in &seq[seq.len() - k..] {
        probe.push_right(base_code(b));
    }
    let probe_canon = probe.canonical();
    let (pc, pci) = prob_tables();
    let mut expected = 0u64;
    for (r, q) in &reads {
        let mut kk = Kmer::new(k);
        let mut len = 0usize;
        let mut prob = 1f32;
        for (i, &bb) in r.iter().enumerate() {
            if base_defined(bb) {
                kk.push_right(base_code(bb));
                prob *= pc[q[i] as usize];
                if len >= k {
                    prob *= pci[q[i - k] as usize];
                }
                len += 1;
            } else {
                len = 0;
                kk.reset();
                prob = 1.0;
            }
            if len >= k
                && prob >= 0.5
                && kk.canonical().cmp_bases(&probe_canon) == std::cmp::Ordering::Equal
            {
                expected += 1;
            }
        }
    }
    assert_eq!(
        table.get_count(&probe),
        expected as u32,
        "k=81 table count mismatch"
    );
}

#[test]
fn seed_kmer_count_symmetric() {
    // First Lambda read from ecco_sub.fq.gz (108 bp).
    let seq = b"AGAGATTCTTGGCGGAGAAACCATAATTGCATCTACTCGTCGCGAACCGCTTTCATCCGGCACAGTATCAAGGTATTTTATGCGCGCACGAAAAGCATC".to_vec();
    let quals = vec![40; seq.len()];
    let k = 62usize;
    let table = TadpoleTable::build(&[(seq.clone(), quals.clone())], k, 0.5);
    let rc: Vec<u8> = rev_comp(&seq).collect();
    for (label, s) in [("forward", &seq), ("rc", &rc)] {
        let mut kmer = Kmer::new(k);
        for &b in &s[s.len() - k..] {
            kmer.push_right(base_code(b));
        }
        eprintln!(
            "{label} tail kmer count={} words={:?} canonical={:?}",
            table.get_count(&kmer),
            kmer.words,
            kmer.canonical().words
        );
    }
    // Directly compare the two canonical forms.
    let mut f = Kmer::new(k);
    for &b in &seq[seq.len() - k..] {
        f.push_right(base_code(b));
    }
    let mut r = Kmer::new(k);
    for &b in &rc[rc.len() - k..] {
        r.push_right(base_code(b));
    }
    eprintln!(
        "f.canonical={:?} r.canonical={:?} rc_of_f={:?}",
        f.canonical().words,
        r.canonical().words,
        f.rc().words
    );
    // Canonical must be orientation-invariant: canonical(f) == canonical(rc(f)).
    let f_rc = f.rc();
    assert_eq!(
        f.canonical().cmp_bases(&f_rc.canonical()),
        std::cmp::Ordering::Equal,
        "canonical orientation-invariance broken"
    );

    // String-level check: rc() of a kmer must equal the reverse
    // complement of the sequence it encodes.
    for (label, s) in [("forward", &seq), ("rc", &rc)] {
        let mut kmer = Kmer::new(k);
        for &b in &s[s.len() - k..] {
            kmer.push_right(base_code(b));
        }
        let kmer_seq: Vec<u8> = (0..k).map(|i| kmer.base_at(i)).collect();
        let rc_seq: Vec<u8> = (0..k).map(|i| kmer.rc().base_at(i)).collect();
        eprintln!("{label} kmer_seq={kmer_seq:?} rc_seq={rc_seq:?}");
        // The rc of the kmer's base sequence (base_at order is 5'->3' as
        // built; verify the reversal+complement relationship explicitly).
        for i in 0..k {
            assert_eq!(
                rc_seq[i],
                3 - kmer_seq[k - 1 - i],
                "{label} rc mismatch at {i}"
            );
        }
    }
}
