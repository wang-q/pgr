/// CIGAR (Compact Idiosyncratic Gapped Alignment Report) operations.
///
/// Implements bit-packed CIGAR storage, coordinate projection, and
/// statistical summaries.
///
/// # Bit-packing
///
/// Each `CigarOp` packs an op code and length into a single `u32`:
/// - bits[31:29] = op code (0: '=', 1: 'X', 2: 'I', 3: 'D', 4: 'M')
/// - bits[28:0]  = length (max 512 Mbp)
///
/// This is alignment-friendly, memory-efficient (4 bytes per op), and
/// enables branch-free coordinate projection via `target_delta`/`query_delta`.
use std::fmt;

// ── Op code constants ────────────────────────────────────────────
const OP_EQ: u32 = 0; // '='
const OP_X: u32 = 1; // 'X'
const OP_I: u32 = 2; // 'I'
const OP_D: u32 = 3; // 'D'
const OP_M: u32 = 4; // 'M'

/// Bit-packed CIGAR operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp(pub(crate) u32);

impl CigarOp {
    /// Create a new `CigarOp` from length and op character.
    ///
    /// Returns an error if the op character is invalid.
    pub fn try_new(len: u32, op: char) -> anyhow::Result<Self> {
        let code = match op {
            '=' => OP_EQ,
            'X' => OP_X,
            'I' => OP_I,
            'D' => OP_D,
            'M' => OP_M,
            _ => anyhow::bail!("invalid CIGAR op: '{op}'"),
        };
        Ok(Self((code << 29) | (len & 0x1FFF_FFFF)))
    }

    /// Create a new `CigarOp` from length and op character (unchecked).
    ///
    /// # Panics
    /// Panics if `op` is not one of '=', 'X', 'I', 'D', 'M'.
    pub(crate) fn new(len: u32, op: char) -> Self {
        let code = match op {
            '=' => OP_EQ,
            'X' => OP_X,
            'I' => OP_I,
            'D' => OP_D,
            'M' => OP_M,
            _ => panic!("invalid CIGAR op: '{op}'"),
        };
        Self((code << 29) | (len & 0x1FFF_FFFF))
    }

    /// Reconstruct from a raw bit-packed u32 (deserialization).
    pub fn from_raw(val: u32) -> Self {
        Self(val)
    }

    /// Decode the op character.
    ///
    /// Returns '?' for invalid op codes (from corrupted `from_raw` inputs).
    pub fn op(self) -> char {
        match self.0 >> 29 {
            OP_EQ => '=',
            OP_X => 'X',
            OP_I => 'I',
            OP_D => 'D',
            OP_M => 'M',
            _ => '?', // invalid op code from corrupted raw value
        }
    }

    /// Decode the length.
    ///
    /// Note: `CigarOp` represents a single op (not a collection), so there is
    /// no meaningful `is_empty` — length is always >= 1 by construction.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(self) -> u32 {
        self.0 & 0x1FFF_FFFF
    }

    /// Advance on the target axis.
    ///
    /// 'I' contributes 0 (insertion in query = gap in target),
    /// all other ops contribute their length.
    pub fn target_delta(self) -> u32 {
        match self.op() {
            'I' => 0,
            _ => self.len(),
        }
    }

    /// Advance on the query axis.
    ///
    /// 'D' contributes 0 (deletion in query = gap in query),
    /// all other ops contribute their length.
    pub fn query_delta(self) -> u32 {
        match self.op() {
            'D' => 0,
            _ => self.len(),
        }
    }
}

impl fmt::Display for CigarOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.len(), self.op())
    }
}

// ── String ↔ Vec<CigarOp> ────────────────────────────────────────

/// Parse a CIGAR string into a vector of `CigarOp`.
///
/// Returns an error if the string contains an invalid op character.
pub fn parse_cigar(s: &str) -> anyhow::Result<Vec<CigarOp>> {
    let mut ops = Vec::new();
    let mut len: u32 = 0;

    for c in s.chars() {
        if c.is_ascii_digit() {
            len = len
                .saturating_mul(10)
                .saturating_add((c as u8 - b'0') as u32);
        } else {
            if !matches!(c, '=' | 'X' | 'I' | 'D' | 'M') {
                anyhow::bail!("invalid CIGAR op: '{c}'");
            }
            ops.push(CigarOp::new(len, c));
            len = 0;
        }
    }

    if len > 0 {
        anyhow::bail!("trailing digits without CIGAR op: {s}");
    }

    Ok(ops)
}

/// Format a slice of `CigarOp` into a CIGAR string.
pub fn format_cigar(ops: &[CigarOp]) -> String {
    let mut s = String::new();
    for op in ops {
        use fmt::Write;
        // `fmt::Write` for `String` is infallible (capacity grows as needed).
        let _ = write!(&mut s, "{op}");
    }
    s
}

/// Extract and parse the `cg:Z:` tag from a PAF tag list. Empty if absent.
pub fn extract_cigar(tags: &[String]) -> anyhow::Result<Vec<CigarOp>> {
    for tag in tags {
        if let Some(s) = tag.strip_prefix("cg:Z:") {
            return parse_cigar(s);
        }
    }
    Ok(Vec::new())
}

// ── Reversal (for bidirectional index) ───────────────────────────

/// Reverse a CIGAR operation list, swapping `I` and `D`.
///
/// When an alignment is viewed from the query's perspective (instead of the
/// target's), the CIGAR must be read backwards and insertions/deletions
/// swapped: an insertion in the original query becomes a deletion in the
/// mirrored record (and vice versa). `=`/`X`/`M` ops are unchanged.
pub fn reverse_cigar(ops: &[CigarOp]) -> Vec<CigarOp> {
    ops.iter()
        .rev()
        .map(|&op| {
            let new_op = match op.op() {
                'I' => 'D',
                'D' => 'I',
                c => c,
            };
            CigarOp::new(op.len(), new_op)
        })
        .collect()
}

/// Extract the sub-CIGAR corresponding to a target sub-interval `[ts, te)`.
///
/// The input CIGAR is assumed to start at `target_start` on the target axis.
/// Insertions (`I`) are included when their target anchor lies inside the
/// interval; deletions and match/mismatch ops are clipped to the interval.
/// Consecutive identical ops are merged.
pub fn slice_cigar_by_target(
    cigar: &[CigarOp],
    target_start: i32,
    ts: i32,
    te: i32,
) -> Vec<CigarOp> {
    let mut out = Vec::new();
    let mut ct = target_start;
    for op in cigar {
        let td = op.target_delta() as i32;
        match op.op() {
            '=' | 'X' | 'M' | 'D' => {
                let os = ts.max(ct);
                let oe = te.min(ct + td);
                if os < oe {
                    push_or_merge(&mut out, op.op(), (oe - os) as u32);
                }
            }
            'I' if ct >= ts && ct < te => {
                push_or_merge(&mut out, 'I', op.len());
            }
            _ => {}
        }
        ct += td;
    }
    out
}

fn push_or_merge(ops: &mut Vec<CigarOp>, op: char, len: u32) {
    if len == 0 {
        return;
    }
    if let Some(last) = ops.last_mut() {
        if last.op() == op {
            *last = CigarOp::new(last.len() + len, op);
            return;
        }
    }
    ops.push(CigarOp::new(len, op));
}

// ── Statistics ───────────────────────────────────────────────────

/// Summary statistics computed from a CIGAR operation list.
///
/// Provides both per‑event and per‑base counts for insertions and
/// deletions, matching the two identity metrics (gi / bi).
#[derive(Debug, Clone, Default)]
pub struct CigarStats {
    /// Matching bases (`M` and `=`).
    pub matches: u32,
    /// Mismatching bases (`X`).
    pub mismatches: u32,
    /// Insertion events (one per `I` op).
    pub ins_events: u32,
    /// Insertion bases (sum of `I` op lengths).
    pub ins_bp: u32,
    /// Deletion events (one per `D` op).
    pub del_events: u32,
    /// Deletion bases (sum of `D` op lengths).
    pub del_bp: u32,
}

/// Compute `CigarStats` from a slice of `CigarOp`.
pub fn cigar_stats(ops: &[CigarOp]) -> CigarStats {
    let mut s = CigarStats::default();
    for op in ops {
        let len = op.len();
        match op.op() {
            'M' | '=' => s.matches += len,
            'X' => s.mismatches += len,
            'I' => {
                s.ins_events += 1;
                s.ins_bp += len;
            }
            'D' => {
                s.del_events += 1;
                s.del_bp += len;
            }
            _ => {}
        }
    }
    s
}

/// Total alignment block length (all bases including indels).
pub fn block_length(stats: &CigarStats) -> u32 {
    stats.matches + stats.mismatches + stats.ins_bp + stats.del_bp
}

// ── Identity ──────────────────────────────────────────────────────

/// Gap-compressed identity.
///
/// `gi = matches / (matches + mismatches + #indel_events)`
///
/// Each indel counts as **one event** regardless of its length,
/// making this metric lenient toward long indels (evaluates homology).
pub fn gap_compressed_identity(ops: &[CigarOp]) -> f64 {
    let s = cigar_stats(ops);
    let total = s.matches + s.mismatches + s.ins_events + s.del_events;
    if total == 0 {
        0.0
    } else {
        s.matches as f64 / total as f64
    }
}

/// Block identity.
///
/// `bi = matches / (matches + mismatches + indel_bp_total)`
///
/// Each indel base counts as a difference, making this metric strict
/// (evaluates sequence identity).
pub fn block_identity(ops: &[CigarOp]) -> f64 {
    let s = cigar_stats(ops);
    let total = s.matches + s.mismatches + s.ins_bp + s.del_bp;
    if total == 0 {
        0.0
    } else {
        s.matches as f64 / total as f64
    }
}

// ── MAF alignment → CIGAR (pgr‑specific) ─────────────────────────

/// Per-column class masks for a pairwise alignment: four `u32` words per
/// 32 columns (insertion / deletion / match / mismatch). Columns where both
/// sequences gap (`-` vs `-`) are unset in every mask and skipped by scans.
pub struct AlignmentMask {
    pub ins: Vec<u32>,
    pub del: Vec<u32>,
    pub m: Vec<u32>,
    pub x: Vec<u32>,
}

/// Classifies every column of an alignment into I/D/=/X masks, so callers
/// that need both CIGAR and `cs:Z` (e.g. `maf_block_to_paf`) scan once.
pub fn classify_alignment(r#ref: &[u8], qry: &[u8]) -> anyhow::Result<AlignmentMask> {
    if r#ref.len() != qry.len() {
        anyhow::bail!("alignment vectors must have equal length");
    }
    let n_words = r#ref.len().div_ceil(32);
    let mut mask = AlignmentMask {
        ins: vec![0; n_words],
        del: vec![0; n_words],
        m: vec![0; n_words],
        x: vec![0; n_words],
    };
    let (ref_chunks, ref_rem) = r#ref.as_chunks::<32>();
    let (qry_chunks, _) = qry.as_chunks::<32>();
    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx2") {
        // SAFETY: gated on runtime AVX2 support.
        unsafe { avx2::classify_chunks_avx2(&mut mask, ref_chunks, qry_chunks) };
        classify_remainder(
            &mut mask,
            ref_rem,
            &qry[ref_chunks.len() * 32..],
            ref_chunks.len(),
        );
        return Ok(mask);
    }
    for (wi, (rc, qc)) in ref_chunks.iter().zip(qry_chunks).enumerate() {
        for k in 0..32 {
            set_class(&mut mask, wi, k, class_of(rc[k], qc[k]));
        }
    }
    classify_remainder(
        &mut mask,
        ref_rem,
        &qry[ref_chunks.len() * 32..],
        ref_chunks.len(),
    );
    Ok(mask)
}

/// Classifies the final partial chunk (scalar for all paths).
fn classify_remainder(mask: &mut AlignmentMask, rem: &[u8], qry_rem: &[u8], word: usize) {
    for (k, (&rc, &qc)) in rem.iter().zip(qry_rem).enumerate() {
        set_class(mask, word, k, class_of(rc, qc));
    }
}

/// Column class as a bit per mask: `0` skip, `1` I, `2` D, `3` =, `4` X.
#[inline]
fn class_of(rc: u8, qc: u8) -> u8 {
    match (rc, qc) {
        (b'-', b'-') => 0,
        (b'-', _) => 1,
        (_, b'-') => 2,
        _ if rc.eq_ignore_ascii_case(&qc) => 3,
        _ => 4,
    }
}

#[inline]
fn set_class(mask: &mut AlignmentMask, wi: usize, k: usize, class: u8) {
    let bit = 1 << k;
    match class {
        1 => mask.ins[wi] |= bit,
        2 => mask.del[wi] |= bit,
        3 => mask.m[wi] |= bit,
        4 => mask.x[wi] |= bit,
        _ => {}
    }
}

#[cfg(target_arch = "x86_64")]
mod avx2 {
    use super::*;
    use std::arch::x86_64::*;

    #[inline]
    unsafe fn load(ptr: *const u8) -> __m256i {
        _mm256_loadu_si256(ptr as *const __m256i)
    }

    #[inline]
    unsafe fn set1(b: u8) -> __m256i {
        _mm256_set1_epi8(b as i8)
    }

    #[inline]
    unsafe fn is_letter(v: __m256i) -> __m256i {
        // Lowercased ASCII letter check: 'a' <= v <= 'z'.
        let ge_a = _mm256_cmpeq_epi8(_mm256_max_epu8(v, set1(0x61)), v);
        let le_z = _mm256_cmpeq_epi8(_mm256_min_epu8(v, set1(0x7A)), v);
        _mm256_and_si256(ge_a, le_z)
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn classify_chunks_avx2(
        mask: &mut AlignmentMask,
        ref_chunks: &[[u8; 32]],
        qry_chunks: &[[u8; 32]],
    ) {
        let dash = set1(b'-');
        let or20 = set1(0x20);
        for (wi, (rc, qc)) in ref_chunks.iter().zip(qry_chunks).enumerate() {
            let r = load(rc.as_ptr());
            let q = load(qc.as_ptr());
            let gap_r = _mm256_cmpeq_epi8(r, dash);
            let gap_q = _mm256_cmpeq_epi8(q, dash);
            let ins = _mm256_andnot_si256(gap_q, gap_r);
            let del = _mm256_andnot_si256(gap_r, gap_q);
            let skip = _mm256_and_si256(gap_r, gap_q);
            let rl = _mm256_or_si256(r, or20);
            let ql = _mm256_or_si256(q, or20);
            let eq_ci = _mm256_and_si256(
                _mm256_and_si256(_mm256_cmpeq_epi8(rl, ql), is_letter(rl)),
                is_letter(ql),
            );
            let eq_case = _mm256_cmpeq_epi8(r, q);
            let m = _mm256_andnot_si256(skip, _mm256_or_si256(eq_case, eq_ci));
            let x = _mm256_andnot_si256(
                _mm256_or_si256(_mm256_or_si256(skip, ins), _mm256_or_si256(del, m)),
                _mm256_set1_epi8(-1),
            );
            mask.ins[wi] = _mm256_movemask_epi8(ins) as u32;
            mask.del[wi] = _mm256_movemask_epi8(del) as u32;
            mask.m[wi] = _mm256_movemask_epi8(m) as u32;
            mask.x[wi] = _mm256_movemask_epi8(x) as u32;
        }
    }
}

/// Merges consecutive equal-op columns into `CigarOp`s (mask-driven).
pub fn scan_cigar_ops(mask: &AlignmentMask) -> Vec<CigarOp> {
    let mut ops: Vec<CigarOp> = Vec::new();
    for wi in 0..mask.ins.len() {
        let mut rem = mask.m[wi];
        let mut col = 0usize;
        while rem != 0 {
            let tz = rem.trailing_zeros() as usize;
            for k in col..col + tz {
                push_op(&mut ops, mask, wi, k);
            }
            let ones = (rem >> tz).trailing_ones() as usize;
            merge_op(&mut ops, '=', ones as u32);
            col = col + tz + ones;
            rem = if tz + ones >= 32 {
                0
            } else {
                rem >> (tz + ones)
            };
        }
        for k in col..32 {
            push_op(&mut ops, mask, wi, k);
        }
    }
    ops
}

#[inline]
fn push_op(ops: &mut Vec<CigarOp>, mask: &AlignmentMask, wi: usize, k: usize) {
    let bit = 1 << k;
    let op_char = if mask.ins[wi] & bit != 0 {
        'I'
    } else if mask.del[wi] & bit != 0 {
        'D'
    } else if mask.x[wi] & bit != 0 {
        'X'
    } else {
        return; // gap-gap column
    };
    merge_op(ops, op_char, 1);
}

#[inline]
fn merge_op(ops: &mut Vec<CigarOp>, op_char: char, len: u32) {
    match ops.last_mut() {
        Some(last) if last.op() == op_char => {
            let new_len = last.len() + len;
            *last = CigarOp::new(new_len, op_char);
        }
        _ => ops.push(CigarOp::new(len, op_char)),
    }
}

/// `cs:Z` compact string from masks plus the original alignment (I/D/X
/// columns emit their bases).
pub fn scan_cs(mask: &AlignmentMask, r#ref: &[u8], qry: &[u8]) -> String {
    let mut cs = String::new();
    let mut run = 0usize;
    let n = r#ref.len();
    let flush = |run: &mut usize, cs: &mut String| {
        if *run > 0 {
            cs.push(':');
            cs.push_str(&run.to_string());
            *run = 0;
        }
    };
    for (wi, _) in mask.ins.iter().enumerate() {
        let base = wi * 32;
        let n_in_word = 32.min(n - base);
        let word_bits = if n_in_word == 32 {
            u32::MAX
        } else {
            (1u32 << n_in_word) - 1
        };
        let mut rem = mask.m[wi] & word_bits;
        let mut col = 0usize;
        while rem != 0 {
            let tz = rem.trailing_zeros() as usize;
            for k in col..col + tz {
                cs_push_col(&mut cs, &mut run, &flush, mask, r#ref, qry, base + k);
            }
            let ones = (rem >> tz).trailing_ones() as usize;
            run += ones;
            col = col + tz + ones;
            rem = if tz + ones >= 32 {
                0
            } else {
                rem >> (tz + ones)
            };
        }
        for k in col..n_in_word {
            cs_push_col(&mut cs, &mut run, &flush, mask, r#ref, qry, base + k);
        }
    }
    flush(&mut run, &mut cs);
    cs
}

#[inline]
fn cs_push_col(
    cs: &mut String,
    run: &mut usize,
    flush: &impl Fn(&mut usize, &mut String),
    mask: &AlignmentMask,
    r#ref: &[u8],
    qry: &[u8],
    i: usize,
) {
    let bit = 1 << (i % 32);
    let rc = r#ref[i];
    let qc = qry[i];
    let wi = i / 32;
    if mask.ins[wi] & bit != 0 {
        flush(run, cs);
        cs.push('+');
        cs.push(qc.to_ascii_uppercase() as char);
    } else if mask.del[wi] & bit != 0 {
        flush(run, cs);
        cs.push('-');
        cs.push(rc.to_ascii_uppercase() as char);
    } else if mask.x[wi] & bit != 0 {
        flush(run, cs);
        cs.push('*');
        cs.push(rc.to_ascii_uppercase() as char);
        cs.push(qc.to_ascii_uppercase() as char);
    }
}

/// Build CIGAR from two MAF `s`-line alignment strings (byte slices).
///
/// Each position is compared (case-insensitive, so soft-masked bases count):
/// - `ref[i] == '-' && qry[i] != '-'` → `I` (insertion in query)
/// - `ref[i] != '-' && qry[i] == '-'` → `D` (deletion in query)
/// - `ref[i] eq_ignore_ascii_case qry[i]` → `=` (match)
/// - otherwise → `X` (mismatch)
///
/// Consecutive identical ops are merged.
pub fn cigar_from_alignment(r#ref: &[u8], qry: &[u8]) -> anyhow::Result<Vec<CigarOp>> {
    let mask = classify_alignment(r#ref, qry)?;
    Ok(scan_cigar_ops(&mask))
}

/// Compact reversible CIGAR (`cs:Z`) from an aligned pair, FastGA `-pafs`
/// style: `:N` match runs, `*<ref><qry>` mismatches, `+<qry>` insertions,
/// `-<ref>` deletions.
pub fn cs_from_alignment(r#ref: &[u8], qry: &[u8]) -> anyhow::Result<String> {
    let mask = classify_alignment(r#ref, qry)?;
    Ok(scan_cs(&mask, r#ref, qry))
}

// ── Tests ─────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::StdRng, Rng, SeedableRng};

    fn random_aln(rng: &mut StdRng, len: usize) -> Vec<u8> {
        let alphabet = *b"ACGTacgtNMRWSYKVHDX-*_ 0123456789";
        (0..len)
            .map(|_| alphabet[rng.random_range(0..alphabet.len())])
            .collect()
    }

    #[test]
    fn masks_match_scalar_classification() {
        // Random alignments: the SIMD/scalar classify masks must agree with
        // the per-column reference at every position and stay disjoint.
        let mut rng = StdRng::seed_from_u64(20260809);
        for len in [0usize, 1, 31, 32, 33, 63, 64, 65, 500, 1000] {
            for _ in 0..20 {
                let r = random_aln(&mut rng, len);
                let q = random_aln(&mut rng, len);
                let mask = classify_alignment(&r, &q).unwrap();
                for i in 0..len {
                    let wi = i / 32;
                    let bit = 1 << (i % 32);
                    let set = [
                        (mask.ins[wi] & bit != 0, 1usize),
                        (mask.del[wi] & bit != 0, 2),
                        (mask.m[wi] & bit != 0, 3),
                        (mask.x[wi] & bit != 0, 4),
                    ];
                    let n_set = set.iter().filter(|(s, _)| *s).count();
                    assert_eq!(n_set, usize::from(class_of(r[i], q[i]) != 0), "col {i}");
                    let class = set.iter().find(|(s, _)| *s).map(|(_, c)| *c).unwrap_or(0);
                    assert_eq!(class, class_of(r[i], q[i]) as usize, "col {i} len {len}");
                }
            }
        }
    }

    // ── CigarOp bit-packing ───────────────────────────────────

    #[test]
    fn test_cigar_op_roundtrip() {
        for (len, op_char) in [(10, '='), (5, 'I'), (3, 'D'), (1, 'X'), (0, 'M')] {
            let op = CigarOp::new(len, op_char);
            assert_eq!(op.op(), op_char, "op mismatch");
            assert_eq!(op.len(), len, "len mismatch");
        }
    }

    #[test]
    fn cs_from_alignment_mixed_ops() {
        // ref:  ACGTACGTACGT
        // qry:  ACGTAAGTA-GT
        // 5 matches, mismatch C/A, 3 matches, deletion C, 2 matches
        let r = b"ACGTACGTACGT";
        let q = b"ACGTAAGTA-GT";
        assert_eq!(
            cs_from_alignment(r, q).unwrap(),
            ":5*CA:3-C:2",
            "FastGA -pafs style cs string"
        );
    }

    #[test]
    fn cs_from_alignment_indels_and_length_check() {
        // Insertions and deletions carry their bases; the string reproduces
        // the alignment columns.
        let r = b"AC-GTACGT";
        let q = b"ACGGTACGT";
        assert_eq!(cs_from_alignment(r, q).unwrap(), ":2+G:6");
        assert!(cs_from_alignment(r, &q[..5]).is_err(), "length mismatch");
    }

    #[test]
    fn test_target_delta() {
        assert_eq!(CigarOp::new(10, '=').target_delta(), 10);
        assert_eq!(CigarOp::new(5, 'I').target_delta(), 0);
        assert_eq!(CigarOp::new(3, 'D').target_delta(), 3);
        assert_eq!(CigarOp::new(7, 'M').target_delta(), 7);
    }

    #[test]
    fn test_query_delta() {
        assert_eq!(CigarOp::new(10, '=').query_delta(), 10);
        assert_eq!(CigarOp::new(5, 'I').query_delta(), 5);
        assert_eq!(CigarOp::new(3, 'D').query_delta(), 0);
        assert_eq!(CigarOp::new(7, 'M').query_delta(), 7);
    }

    #[test]
    fn test_zero_len_op() {
        let op = CigarOp::new(0, 'I');
        assert_eq!(op.target_delta(), 0);
        assert_eq!(op.query_delta(), 0);
    }

    #[test]
    #[should_panic(expected = "invalid CIGAR op")]
    fn test_invalid_op_panics() {
        CigarOp::new(10, 'Q');
    }

    // ── String ↔ Vec<CigarOp> ─────────────────────────────────

    #[test]
    fn test_parse_cigar_basic() {
        let ops = parse_cigar("10=5I5D").unwrap();
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], CigarOp::new(10, '='));
        assert_eq!(ops[1], CigarOp::new(5, 'I'));
        assert_eq!(ops[2], CigarOp::new(5, 'D'));
    }

    #[test]
    fn test_parse_cigar_empty() {
        let ops = parse_cigar("").unwrap();
        assert!(ops.is_empty());
    }

    #[test]
    fn test_parse_cigar_digits_only() {
        assert!(parse_cigar("10").is_err());
        assert!(parse_cigar("10=5").is_err());
    }

    #[test]
    fn test_parse_cigar_zero_len() {
        let ops = parse_cigar("0=5I").unwrap();
        assert_eq!(ops.len(), 2);
        assert_eq!(ops[0], CigarOp::new(0, '='));
        assert_eq!(ops[1], CigarOp::new(5, 'I'));
    }

    #[test]
    fn test_format_cigar_roundtrip() {
        let cases = ["10=5I5D", "3M1I2D", "", "100="];
        for case in cases {
            let ops = parse_cigar(case).unwrap();
            let formatted = format_cigar(&ops);
            assert_eq!(formatted, case, "roundtrip failed for '{case}'");
        }
    }

    // ── reverse_cigar ────────────────────────────────────────

    #[test]
    fn test_reverse_cigar_basic() {
        // 10M5I3D → reversed: 3I5D10M
        let ops = parse_cigar("10M5I3D").unwrap();
        let rev = reverse_cigar(&ops);
        assert_eq!(format_cigar(&rev), "3I5D10M");
    }

    #[test]
    fn test_reverse_cigar_no_indels() {
        // 10=2X8= → reversed: 8=2X10= (no I/D swap, just order reversed)
        let ops = parse_cigar("10=2X8=").unwrap();
        let rev = reverse_cigar(&ops);
        assert_eq!(format_cigar(&rev), "8=2X10=");
    }

    #[test]
    fn test_reverse_cigar_empty() {
        let rev = reverse_cigar(&[]);
        assert!(rev.is_empty());
    }

    #[test]
    fn test_reverse_cigar_double_reversal() {
        // reverse(reverse(x)) == x (I↔D swapped twice = identity)
        let ops = parse_cigar("5M3I2D7=").unwrap();
        let rev2 = reverse_cigar(&reverse_cigar(&ops));
        assert_eq!(format_cigar(&rev2), format_cigar(&ops));
    }

    #[test]
    fn test_reverse_cigar_only_indels() {
        // 5I3D → reversed: 3I5D
        let ops = parse_cigar("5I3D").unwrap();
        let rev = reverse_cigar(&ops);
        assert_eq!(format_cigar(&rev), "3I5D");
    }

    #[test]
    fn test_reverse_cigar_preserves_lengths() {
        let ops = parse_cigar("100M1I99M1D200=").unwrap();
        let rev = reverse_cigar(&ops);
        // Total length consumed should be preserved per-axis
        let orig_query: u32 = ops.iter().map(|o| o.query_delta()).sum();
        let rev_query: u32 = rev.iter().map(|o| o.query_delta()).sum();
        assert_eq!(orig_query, rev_query, "query-axis length changed");
        let orig_target: u32 = ops.iter().map(|o| o.target_delta()).sum();
        let rev_target: u32 = rev.iter().map(|o| o.target_delta()).sum();
        assert_eq!(orig_target, rev_target, "target-axis length changed");
    }

    // ── Statistics ────────────────────────────────────────────

    #[test]
    fn test_cigar_stats_basic() {
        let ops = parse_cigar("10=5I3D").unwrap();
        let s = cigar_stats(&ops);
        assert_eq!(s.matches, 10);
        assert_eq!(s.mismatches, 0);
        assert_eq!(s.ins_events, 1);
        assert_eq!(s.ins_bp, 5);
        assert_eq!(s.del_events, 1);
        assert_eq!(s.del_bp, 3);
    }

    #[test]
    fn test_cigar_stats_with_mismatch() {
        let ops = parse_cigar("5=2X3I").unwrap();
        let s = cigar_stats(&ops);
        assert_eq!(s.matches, 5);
        assert_eq!(s.mismatches, 2);
        assert_eq!(s.ins_events, 1);
        assert_eq!(s.ins_bp, 3);
    }

    #[test]
    fn test_cigar_stats_multiple_events() {
        let ops = parse_cigar("3I5=2D4=1I").unwrap();
        let s = cigar_stats(&ops);
        assert_eq!(s.matches, 9);
        assert_eq!(s.ins_events, 2);
        assert_eq!(s.ins_bp, 4);
        assert_eq!(s.del_events, 1);
        assert_eq!(s.del_bp, 2);
    }

    #[test]
    fn test_block_length() {
        let ops = parse_cigar("10=5I3D").unwrap();
        let s = cigar_stats(&ops);
        assert_eq!(block_length(&s), 18); // 10 + 0 + 5 + 3
    }

    // ── Identity ──────────────────────────────────────────────

    #[test]
    fn test_gi_pure_match() {
        let ops = parse_cigar("10=").unwrap();
        assert!((gap_compressed_identity(&ops) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_gi_with_insertion() {
        let ops = parse_cigar("10=5I").unwrap();
        let gi = gap_compressed_identity(&ops);
        let expected = 10.0 / (10.0 + 0.0 + 1.0);
        assert!((gi - expected).abs() < 1e-6);
    }

    #[test]
    fn test_gi_with_deletion() {
        let ops = parse_cigar("10=5D").unwrap();
        let gi = gap_compressed_identity(&ops);
        let expected = 10.0 / (10.0 + 0.0 + 1.0);
        assert!((gi - expected).abs() < 1e-6);
    }

    #[test]
    fn test_gi_mixed() {
        let ops = parse_cigar("10=2X3I4D").unwrap();
        let gi = gap_compressed_identity(&ops);
        let expected = 10.0 / (10.0 + 2.0 + 2.0);
        assert!((gi - expected).abs() < 1e-6);
    }

    #[test]
    fn test_gi_empty() {
        assert!((gap_compressed_identity(&[]) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_bi_with_insertion() {
        let ops = parse_cigar("10=5I").unwrap();
        let bi = block_identity(&ops);
        let expected = 10.0 / (10.0 + 0.0 + 5.0);
        assert!((bi - expected).abs() < 1e-6);
    }

    #[test]
    fn test_bi_empty() {
        assert!((block_identity(&[]) - 0.0).abs() < 1e-9);
    }

    // ── MAF alignment → CIGAR (pgr‑specific) ─────────────────

    #[test]
    fn test_cigar_from_alignment_all_match() {
        let ops = cigar_from_alignment(b"ACGT", b"ACGT").unwrap();
        assert_eq!(ops, vec![CigarOp::new(4, '=')]);
    }

    #[test]
    fn test_cigar_from_alignment_mismatches() {
        // ACGT vs AGGT → = X = =
        let ops = cigar_from_alignment(b"ACGT", b"AGGT").unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::new(1, '='),
                CigarOp::new(1, 'X'),
                CigarOp::new(2, '='),
            ]
        );
    }

    #[test]
    fn test_cigar_from_alignment_case_insensitive() {
        // Soft-masked bases (lowercase) count as match
        let ops = cigar_from_alignment(b"acgt", b"ACGT").unwrap();
        assert_eq!(ops, vec![CigarOp::new(4, '=')]);
    }

    #[test]
    fn test_cigar_from_alignment_ref_gap() {
        let ops = cigar_from_alignment(b"ACG-", b"ACGT").unwrap();
        assert_eq!(ops, vec![CigarOp::new(3, '='), CigarOp::new(1, 'I')]);
    }

    #[test]
    fn test_cigar_from_alignment_qry_gap() {
        let ops = cigar_from_alignment(b"ACGT", b"ACG-").unwrap();
        assert_eq!(ops, vec![CigarOp::new(3, '='), CigarOp::new(1, 'D')]);
    }

    #[test]
    fn test_cigar_from_alignment_interleaved() {
        // AC-TG vs ACGT- → = = I = D
        let ops = cigar_from_alignment(b"AC-TG", b"ACGT-").unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::new(2, '='),
                CigarOp::new(1, 'I'),
                CigarOp::new(1, '='),
                CigarOp::new(1, 'D'),
            ]
        );
    }

    #[test]
    fn test_cigar_from_alignment_terminal_gaps() {
        // -ACGT- vs TACGTA → I = = = = I
        let ops = cigar_from_alignment(b"-ACGT-", b"TACGTA").unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::new(1, 'I'),
                CigarOp::new(4, '='),
                CigarOp::new(1, 'I'),
            ]
        );
    }

    #[test]
    fn test_cigar_from_alignment_all_gaps() {
        let ops = cigar_from_alignment(b"---", b"---").unwrap();
        assert!(ops.is_empty());
    }

    #[test]
    fn test_cigar_from_alignment_merge_consecutive() {
        let ops = cigar_from_alignment(b"ACG--T", b"ACGTT-").unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::new(3, '='),
                CigarOp::new(2, 'I'),
                CigarOp::new(1, 'D'),
            ]
        );
    }

    #[test]
    fn test_format_cigar_only() {
        // Direct format without parse dependency
        let ops = vec![CigarOp::new(10, 'M'), CigarOp::new(1, 'I')];
        assert_eq!(format_cigar(&ops), "10M1I");
    }

    #[test]
    fn test_cigar_stats_all_ops() {
        // Cover all five CIGAR op types
        let ops = parse_cigar("5M3=2X4I1D").unwrap();
        let s = cigar_stats(&ops);
        assert_eq!(s.matches, 8); // 5M + 3=
        assert_eq!(s.mismatches, 2);
        assert_eq!(s.ins_events, 1);
        assert_eq!(s.ins_bp, 4);
        assert_eq!(s.del_events, 1);
        assert_eq!(s.del_bp, 1);
    }

    #[test]
    fn test_cigar_from_alignment_mixed_gaps() {
        // ref: A-CG--T, qry: A-CGTT-, col 2 both-gap skipped
        let ops = cigar_from_alignment(b"A-CG--T", b"A-CGTT-").unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::new(3, '='),
                CigarOp::new(2, 'I'),
                CigarOp::new(1, 'D'),
            ]
        );
    }

    #[test]
    fn test_cigar_from_alignment_length_mismatch() {
        assert!(cigar_from_alignment(b"ACG", b"ACGT").is_err());
    }
}
