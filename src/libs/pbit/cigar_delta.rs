//! CIGAR-driven delta codec: bit-packed CIGAR ops + X/I base stream.
//!
//! Pairs with [`super::lz_diff`] to provide an alternative delta encoding for
//! PAF-covered sample segments. The packed layout (before flate2 compression)
//! is:
//!
//! ```text
//! u32 cigar_op_count
//! [u32; cigar_op_count]   // raw bit-packed CigarOp values
//! u32 xi_base_count
//! [u8; xi_base_count]     // X/I bases (ASCII, in CIGAR forward-traversal order)
//! ```

use anyhow::{bail, Result};

use crate::libs::paf::cigar::CigarOp;

/// Pack CIGAR ops + X/I bases into a flate2-compressed byte buffer.
pub fn pack_cigar(ops: &[CigarOp], xi_bases: &[u8]) -> Result<Vec<u8>> {
    use std::io::Write;
    let mut raw = Vec::with_capacity(8 + ops.len() * 4 + xi_bases.len());
    raw.extend_from_slice(&(ops.len() as u32).to_le_bytes());
    for op in ops {
        // CigarOp.0 is pub(crate), accessible from pbit (same crate).
        raw.extend_from_slice(&op.0.to_le_bytes());
    }
    raw.extend_from_slice(&(xi_bases.len() as u32).to_le_bytes());
    raw.extend_from_slice(xi_bases);
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(&raw)?;
    Ok(encoder.finish()?)
}

/// Unpack a flate2-compressed buffer into (CIGAR ops, X/I bases).
pub fn unpack_cigar(packed: &[u8]) -> Result<(Vec<CigarOp>, Vec<u8>)> {
    use std::io::Read;
    let mut decoder = flate2::read::GzDecoder::new(packed);
    let mut raw = Vec::new();
    decoder.read_to_end(&mut raw)?;
    let mut cursor = std::io::Cursor::new(raw);
    let mut buf4 = [0u8; 4];
    cursor.read_exact(&mut buf4)?;
    let op_count = u32::from_le_bytes(buf4) as usize;
    let mut ops = Vec::with_capacity(op_count);
    for _ in 0..op_count {
        cursor.read_exact(&mut buf4)?;
        ops.push(CigarOp::from_raw(u32::from_le_bytes(buf4)));
    }
    cursor.read_exact(&mut buf4)?;
    let xi_count = u32::from_le_bytes(buf4) as usize;
    let mut xi_bases = vec![0u8; xi_count];
    cursor.read_exact(&mut xi_bases)?;
    Ok((ops, xi_bases))
}

/// Apply CIGAR to a reference slice, consuming X/I bases, producing the
/// sample sequence (no '-' gap insertion, no coordinate trimming).
///
/// Simplified variant of `build_maf_block`: produces raw sample sequence
/// directly. `M` ops are not expected here (compressor splits M into `=/X`);
/// encountering `M` indicates corrupt or unsupported packed CIGAR data.
pub fn apply_cigar(ref_seq: &[u8], ops: &[CigarOp], xi_bases: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    let mut rt: usize = 0; // ref cursor
    let mut xi: usize = 0; // X/I base stream cursor
    for op in ops {
        let len = op.len() as usize;
        match op.op() {
            '=' => {
                if rt + len > ref_seq.len() {
                    bail!("CIGAR '=' exceeds reference length");
                }
                out.extend_from_slice(&ref_seq[rt..rt + len]);
                rt += len;
            }
            'M' => {
                bail!("unexpected CIGAR op 'M' in pbit delta (should have been split to =/X)");
            }
            'X' => {
                if xi + len > xi_bases.len() {
                    bail!("CIGAR 'X' exceeds X/I base stream");
                }
                if rt + len > ref_seq.len() {
                    bail!("CIGAR 'X' exceeds reference length");
                }
                out.extend_from_slice(&xi_bases[xi..xi + len]);
                xi += len;
                rt += len; // X consumes reference too
            }
            'I' => {
                if xi + len > xi_bases.len() {
                    bail!("CIGAR 'I' exceeds X/I base stream");
                }
                out.extend_from_slice(&xi_bases[xi..xi + len]);
                xi += len;
                // I does not advance reference
            }
            'D' => {
                if rt + len > ref_seq.len() {
                    bail!("CIGAR 'D' exceeds reference length");
                }
                rt += len; // D advances reference only
            }
            other => bail!("invalid CIGAR op: '{}'", other),
        }
    }
    if xi != xi_bases.len() {
        bail!(
            "CIGAR consumed {} X/I bases but {} were packed",
            xi,
            xi_bases.len()
        );
    }
    if rt != ref_seq.len() {
        bail!(
            "CIGAR consumed {} reference bases but {} available",
            rt,
            ref_seq.len()
        );
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::paf::cigar::CigarOp;

    fn ops(specs: &[(u32, char)]) -> Vec<CigarOp> {
        specs
            .iter()
            .map(|(len, op)| CigarOp::try_new(*len, *op).unwrap())
            .collect()
    }

    #[test]
    fn test_pack_unpack_roundtrip_empty() {
        let packed = pack_cigar(&[], &[]).unwrap();
        let (ops, xi) = unpack_cigar(&packed).unwrap();
        assert!(ops.is_empty());
        assert!(xi.is_empty());
    }

    #[test]
    fn test_pack_unpack_roundtrip_match_only() {
        let ops_in = ops(&[(100, '='), (200, '=')]);
        let packed = pack_cigar(&ops_in, &[]).unwrap();
        let (ops_out, xi) = unpack_cigar(&packed).unwrap();
        assert_eq!(ops_out.len(), 2);
        assert_eq!(ops_out[0].op(), '=');
        assert_eq!(ops_out[0].len(), 100);
        assert_eq!(ops_out[1].op(), '=');
        assert_eq!(ops_out[1].len(), 200);
        assert!(xi.is_empty());
    }

    #[test]
    fn test_pack_unpack_roundtrip_with_xi() {
        let ops_in = ops(&[(10, '='), (3, 'X'), (5, 'I'), (2, 'D'), (4, '=')]);
        // 3 X bases + 5 I bases = 8 bases total
        let xi_correct: Vec<u8> = b"ACGTTTGC".to_vec();
        let packed = pack_cigar(&ops_in, &xi_correct).unwrap();
        let (ops_out, xi_out) = unpack_cigar(&packed).unwrap();
        assert_eq!(ops_out.len(), 5);
        assert_eq!(xi_out, xi_correct);
    }

    #[test]
    fn test_apply_cigar_pure_match() {
        let ref_seq = b"ACGTACGTACGT";
        let ops_in = ops(&[(12, '=')]);
        let out = apply_cigar(ref_seq, &ops_in, &[]).unwrap();
        assert_eq!(out, ref_seq);
    }

    #[test]
    fn test_apply_cigar_with_snps() {
        let ref_seq = b"ACGTACGTACGT";
        // =4 X1 =4 X1 =2: SNP at ref[4] (A→G) and ref[9] (C→N)
        let ops_in = ops(&[(4, '='), (1, 'X'), (4, '='), (1, 'X'), (2, '=')]);
        let xi = b"GN";
        let out = apply_cigar(ref_seq, &ops_in, xi).unwrap();
        // ref: A(0) C(1) G(2) T(3) A(4) C(5) G(6) T(7) A(8) C(9) G(10) T(11)
        // =4: ACGT (ref[0..4])
        // X1: G (xi), consumes ref[4..5]="A"
        // =4: CGTA (ref[5..9])
        // X1: N (xi), consumes ref[9..10]="C"
        // =2: GT (ref[10..12])
        // result: ACGT + G + CGTA + N + GT = ACGTGCGTANGT
        assert_eq!(out, b"ACGTGCGTANGT");
    }

    #[test]
    fn test_apply_cigar_with_insertion() {
        let ref_seq = b"ACGTACGT";
        // =4 I3 =4: insert "TTT" between pos 4 and 5
        let ops_in = ops(&[(4, '='), (3, 'I'), (4, '=')]);
        let xi = b"TTT";
        let out = apply_cigar(ref_seq, &ops_in, xi).unwrap();
        assert_eq!(out, b"ACGTTTTACGT");
    }

    #[test]
    fn test_apply_cigar_with_deletion() {
        let ref_seq = b"ACGTACGTACGT";
        // =4 D2 =6: delete ref[4..6]="AC"
        let ops_in = ops(&[(4, '='), (2, 'D'), (6, '=')]);
        let out = apply_cigar(ref_seq, &ops_in, &[]).unwrap();
        // ref: A(0) C(1) G(2) T(3) A(4) C(5) G(6) T(7) A(8) C(9) G(10) T(11)
        // =4: ACGT (ref[0..4])
        // D2: skip ref[4..6]="AC"
        // =6: GTACGT (ref[6..12])
        // result: ACGT + GTACGT = ACGTGTACGT
        assert_eq!(out, b"ACGTGTACGT");
    }

    #[test]
    fn test_apply_cigar_combined() {
        let ref_seq = b"ACGTACGTACGT";
        // =2 X1 =2 I2 D2 =5: X at ref[2], insert "GG" after ref[4], delete ref[5..7]
        let ops_in = ops(&[(2, '='), (1, 'X'), (2, '='), (2, 'I'), (2, 'D'), (5, '=')]);
        let xi = b"TGG"; // 1 X + 2 I
        let out = apply_cigar(ref_seq, &ops_in, xi).unwrap();
        // ref: A(0) C(1) G(2) T(3) A(4) C(5) G(6) T(7) A(8) C(9) G(10) T(11)
        // =2 (ref[0..2]):  AC
        // X1 (ref[2..3]):  T (xi)
        // =2 (ref[3..5]):  TA
        // I2:              GG (xi)
        // D2 (ref[5..7]):  skip CG
        // =5 (ref[7..12]): TACGT
        // result: AC + T + TA + GG + TACGT = ACTTAGGTACGT
        assert_eq!(out, b"ACTTAGGTACGT");
    }

    #[test]
    fn test_apply_cigar_with_n_bases() {
        let ref_seq = b"ACGTACGT";
        // =3 X2 =3: replace ref[3..5] with NN
        let ops_in = ops(&[(3, '='), (2, 'X'), (3, '=')]);
        let xi = b"NN";
        let out = apply_cigar(ref_seq, &ops_in, xi).unwrap();
        // =3: ACG (ref[0..3])
        // X2: NN (xi), consumes ref[3..5]="TA"
        // =3: ref[5..8]="CGT"
        // result: ACG + NN + CGT = ACGNNCGT
        assert_eq!(out, b"ACGNNCGT");
    }

    #[test]
    fn test_apply_cigar_m_rejected() {
        let ref_seq = b"ACGTACGT";
        let ops_in = ops(&[(8, 'M')]);
        let err = apply_cigar(ref_seq, &ops_in, &[]).unwrap_err();
        assert!(err.to_string().contains("unexpected CIGAR op 'M'"));
    }

    #[test]
    fn test_apply_cigar_error_ref_overflow() {
        let ref_seq = b"ACGT";
        let ops_in = ops(&[(5, '=')]);
        let err = apply_cigar(ref_seq, &ops_in, &[]).unwrap_err();
        assert!(err.to_string().contains("exceeds reference length"));
    }

    #[test]
    fn test_apply_cigar_error_xi_overflow() {
        let ref_seq = b"ACGTACGT";
        let ops_in = ops(&[(2, 'X')]);
        let err = apply_cigar(ref_seq, &ops_in, b"A").unwrap_err();
        assert!(err.to_string().contains("exceeds X/I base stream"));
    }

    #[test]
    fn test_apply_cigar_error_xi_excess() {
        // xi_bases has one more base than the CIGAR consumes; must be rejected
        // rather than silently discarded (corrupted-archive guard).
        let ref_seq = b"ACGTACGT";
        let ops_in = ops(&[(2, 'X')]);
        let err = apply_cigar(ref_seq, &ops_in, b"ABC").unwrap_err();
        assert!(err
            .to_string()
            .contains("consumed 2 X/I bases but 3 were packed"));
    }

    #[test]
    fn test_pack_unpack_roundtrip_with_n() {
        let ops_in = ops(&[(3, '='), (2, 'X'), (3, 'I')]);
        let xi_in = b"NNACG";
        let packed = pack_cigar(&ops_in, xi_in).unwrap();
        let (ops_out, xi_out) = unpack_cigar(&packed).unwrap();
        assert_eq!(ops_out.len(), 3);
        assert_eq!(xi_out, xi_in);
    }

    #[test]
    fn test_full_roundtrip_via_apply() {
        // Compose ref + sample, derive CIGAR+XI, pack, unpack, apply, compare.
        let ref_seq = b"ACGTACGTACGTACGT";
        // ops: =3 X1 =6 I2 D2 =4 (1 X + 2 I = 3 XI bases)
        let ops_in = ops(&[(3, '='), (1, 'X'), (6, '='), (2, 'I'), (2, 'D'), (4, '=')]);
        let xi_in = b"NTT";
        let packed = pack_cigar(&ops_in, xi_in).unwrap();
        let (ops_out, xi_out) = unpack_cigar(&packed).unwrap();
        let reconstructed = apply_cigar(ref_seq, &ops_out, &xi_out).unwrap();
        // ref: A(0) C(1) G(2) T(3) A(4) C(5) G(6) T(7) A(8) C(9) G(10) T(11) A(12) C(13) G(14) T(15)
        // =3 (ref[0..3]):  ACG
        // X1 (ref[3..4]):  N (xi)
        // =6 (ref[4..10]): ACGTAC
        // I2:              TT (xi)
        // D2 (ref[10..12]): skip GT
        // =4 (ref[12..16]):ACGT
        // result: ACG + N + ACGTAC + TT + ACGT = ACGNACGTACTTACGT
        assert_eq!(reconstructed, b"ACGNACGTACTTACGT");
    }

    #[test]
    fn test_apply_cigar_error_x_ref_overflow() {
        // X consumes reference too; an X op that overruns ref_seq must be
        // rejected (previously only the X/I base stream was bounds-checked).
        let ref_seq = b"ACGT";
        let ops_in = ops(&[(3, '='), (2, 'X')]);
        let err = apply_cigar(ref_seq, &ops_in, b"GG").unwrap_err();
        assert!(err.to_string().contains("exceeds reference length"));
    }

    #[test]
    fn test_apply_cigar_error_ref_underconsumed() {
        // CIGAR covers only part of the reference; the leftover ref bases must
        // be rejected rather than silently producing a truncated sample.
        let ref_seq = b"ACGTACGT";
        let ops_in = ops(&[(4, '=')]);
        let err = apply_cigar(ref_seq, &ops_in, &[]).unwrap_err();
        assert!(err
            .to_string()
            .contains("consumed 4 reference bases but 8 available"));
    }
}
