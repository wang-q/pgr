/// CIGAR (Compact Idiosyncratic Gapped Alignment Report) operations.
///
/// Implements bit-packed CIGAR storage, coordinate projection, and
/// statistical summaries.  Modeled after impg's `CigarOp`
/// (`impg-0.4.1/src/impg.rs:73-138`) with additional statistics
/// inspired by wgatools' `Cigar` struct.
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
    /// Panics on invalid op (internal invariant — callers must validate).
    pub fn new(len: u32, op: char) -> Self {
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
    pub fn op(self) -> char {
        match self.0 >> 29 {
            OP_EQ => '=',
            OP_X => 'X',
            OP_I => 'I',
            OP_D => 'D',
            OP_M => 'M',
            _ => unreachable!(),
        }
    }

    /// Decode the length.
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
/// Mirrors impg's `parse_cigar_to_delta` (`impg.rs:2884-2898`).
///
/// # Panics
///
/// Panics if the string contains an invalid op character.
pub fn parse_cigar(s: &str) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut len: u32 = 0;

    for c in s.chars() {
        if c.is_ascii_digit() {
            len = len
                .saturating_mul(10)
                .saturating_add((c as u8 - b'0') as u32);
        } else {
            ops.push(CigarOp::new(len, c));
            len = 0;
        }
    }

    ops
}

/// Format a slice of `CigarOp` into a CIGAR string.
pub fn format_cigar(ops: &[CigarOp]) -> String {
    let mut s = String::new();
    for op in ops {
        use fmt::Write;
        write!(&mut s, "{op}").unwrap();
    }
    s
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

/// Build CIGAR from two MAF `s`-line alignment strings (byte slices).
///
/// Each position is compared:
/// - `ref[i] == '-' && qry[i] != '-'` → `I` (insertion in query)
/// - `ref[i] != '-' && qry[i] == '-'` → `D` (deletion in query)
/// - otherwise → `M` (match/mismatch, not distinguished in v1)
///
/// Consecutive identical ops are merged.
pub fn cigar_from_alignment(r#ref: &[u8], qry: &[u8]) -> Vec<CigarOp> {
    assert_eq!(
        r#ref.len(),
        qry.len(),
        "alignment vectors must have equal length"
    );

    let mut ops: Vec<CigarOp> = Vec::new();

    for (&rc, &qc) in r#ref.iter().zip(qry.iter()) {
        let op_char = match (rc, qc) {
            (b'-', b'-') => continue, // both gaps — degenerate, skip
            (b'-', _) => 'I',
            (_, b'-') => 'D',
            _ => 'M',
        };

        match ops.last_mut() {
            Some(last) if last.op() == op_char => {
                let new_len = last.len() + 1;
                *last = CigarOp::new(new_len, op_char);
            }
            _ => ops.push(CigarOp::new(1, op_char)),
        }
    }

    ops
}

// ── Tests ─────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

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
        // mirrors impg test_parse_cigar_to_delta_basic
        let ops = parse_cigar("10=5I5D");
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], CigarOp::new(10, '='));
        assert_eq!(ops[1], CigarOp::new(5, 'I'));
        assert_eq!(ops[2], CigarOp::new(5, 'D'));
    }

    #[test]
    fn test_parse_cigar_empty() {
        let ops = parse_cigar("");
        assert!(ops.is_empty());
    }

    #[test]
    fn test_parse_cigar_digits_only() {
        let ops = parse_cigar("10");
        assert!(ops.is_empty());
    }

    #[test]
    fn test_parse_cigar_zero_len() {
        let ops = parse_cigar("0=5I");
        assert_eq!(ops.len(), 2);
        assert_eq!(ops[0], CigarOp::new(0, '='));
        assert_eq!(ops[1], CigarOp::new(5, 'I'));
    }

    #[test]
    fn test_format_cigar_roundtrip() {
        let cases = ["10=5I5D", "3M1I2D", "", "100="];
        for case in cases {
            let ops = parse_cigar(case);
            let formatted = format_cigar(&ops);
            assert_eq!(formatted, case, "roundtrip failed for '{case}'");
        }
    }

    // ── Statistics ────────────────────────────────────────────

    #[test]
    fn test_cigar_stats_basic() {
        let ops = parse_cigar("10=5I3D");
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
        let ops = parse_cigar("5=2X3I");
        let s = cigar_stats(&ops);
        assert_eq!(s.matches, 5);
        assert_eq!(s.mismatches, 2);
        assert_eq!(s.ins_events, 1);
        assert_eq!(s.ins_bp, 3);
    }

    #[test]
    fn test_cigar_stats_multiple_events() {
        let ops = parse_cigar("3I5=2D4=1I");
        let s = cigar_stats(&ops);
        assert_eq!(s.matches, 9);
        assert_eq!(s.ins_events, 2);
        assert_eq!(s.ins_bp, 4);
        assert_eq!(s.del_events, 1);
        assert_eq!(s.del_bp, 2);
    }

    #[test]
    fn test_block_length() {
        let ops = parse_cigar("10=5I3D");
        let s = cigar_stats(&ops);
        assert_eq!(block_length(&s), 18); // 10 + 0 + 5 + 3
    }

    // ── Identity ──────────────────────────────────────────────

    #[test]
    fn test_gi_pure_match() {
        let ops = parse_cigar("10=");
        assert!((gap_compressed_identity(&ops) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_gi_with_insertion() {
        let ops = parse_cigar("10=5I");
        let gi = gap_compressed_identity(&ops);
        let expected = 10.0 / (10.0 + 0.0 + 1.0);
        assert!((gi - expected).abs() < 1e-6);
    }

    #[test]
    fn test_gi_with_deletion() {
        let ops = parse_cigar("10=5D");
        let gi = gap_compressed_identity(&ops);
        let expected = 10.0 / (10.0 + 0.0 + 1.0);
        assert!((gi - expected).abs() < 1e-6);
    }

    #[test]
    fn test_gi_mixed() {
        let ops = parse_cigar("10=2X3I4D");
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
        let ops = parse_cigar("10=5I");
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
        let ops = cigar_from_alignment(b"ACGT", b"ACGT");
        assert_eq!(ops, vec![CigarOp::new(4, 'M')]);
    }

    #[test]
    fn test_cigar_from_alignment_ref_gap() {
        let ops = cigar_from_alignment(b"ACG-", b"ACGT");
        assert_eq!(ops, vec![CigarOp::new(3, 'M'), CigarOp::new(1, 'I')]);
    }

    #[test]
    fn test_cigar_from_alignment_qry_gap() {
        let ops = cigar_from_alignment(b"ACGT", b"ACG-");
        assert_eq!(ops, vec![CigarOp::new(3, 'M'), CigarOp::new(1, 'D')]);
    }

    #[test]
    fn test_cigar_from_alignment_interleaved() {
        // AC-TG vs ACGT- →  M M I M D
        let ops = cigar_from_alignment(b"AC-TG", b"ACGT-");
        assert_eq!(
            ops,
            vec![
                CigarOp::new(2, 'M'),
                CigarOp::new(1, 'I'),
                CigarOp::new(1, 'M'),
                CigarOp::new(1, 'D'),
            ]
        );
    }

    #[test]
    fn test_cigar_from_alignment_terminal_gaps() {
        // -ACGT- vs TACGTA →  I M M M M I
        let ops = cigar_from_alignment(b"-ACGT-", b"TACGTA");
        assert_eq!(
            ops,
            vec![
                CigarOp::new(1, 'I'),
                CigarOp::new(4, 'M'),
                CigarOp::new(1, 'I'),
            ]
        );
    }

    #[test]
    fn test_cigar_from_alignment_all_gaps() {
        let ops = cigar_from_alignment(b"---", b"---");
        assert!(ops.is_empty());
    }

    #[test]
    fn test_cigar_from_alignment_merge_consecutive() {
        let ops = cigar_from_alignment(b"ACG--T", b"ACGTT-");
        assert_eq!(
            ops,
            vec![
                CigarOp::new(3, 'M'),
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
        let ops = parse_cigar("5M3=2X4I1D");
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
        let ops = cigar_from_alignment(b"A-CG--T", b"A-CGTT-");
        assert_eq!(
            ops,
            vec![
                CigarOp::new(3, 'M'),
                CigarOp::new(2, 'I'),
                CigarOp::new(1, 'D'),
            ]
        );
    }

    #[test]
    #[should_panic(expected = "alignment vectors must have equal length")]
    fn test_cigar_from_alignment_length_mismatch() {
        cigar_from_alignment(b"ACG", b"ACGT");
    }
}
