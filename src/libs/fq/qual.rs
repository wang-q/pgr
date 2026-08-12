//! FASTQ quality conversion and Phred offset detection (base layer, kept in
//! pgr when the fq command group migrates to anchr).

use crate::libs::fmt::seq::SeqRecord;

/// Phred offset for Sanger / Illumina 1.8+ quality encoding.
pub const PHRED33: u8 = 33;

/// Phred offset for Illumina 1.3-1.7 quality encoding (Solexa approximated).
pub const PHRED64: u8 = 64;

/// Auto-detects the Phred offset (+33 or +64) from a record sample using the
/// BBDuk flip-flop heuristic (BBTools-40.01 `stream/FASTQ.java` `testQuality`):
/// a quality implying Q>54 under +33 (or an N base carrying `@`/`B`) flips to
/// +64, a quality implying Q<-5 under +64 flips back, reads >= 200 bp force
/// +33, and two flips fall back to +33.
pub fn detect_quality_base(sample: &[SeqRecord]) -> u8 {
    const QUAL_THRESH: i32 = 54;
    const FORCE_PHRED33_LEN: usize = 200;

    let mut offset = PHRED33 as i32;
    let mut flips = 0;
    let mut detect = true;
    let mut junk_chars = 0;

    for rec in sample {
        let seq = rec.sequence();
        let qual = rec.quality_scores();
        if detect && seq.len() >= FORCE_PHRED33_LEN {
            if offset != PHRED33 as i32 {
                offset = PHRED33 as i32;
            }
            detect = false;
        }
        for (i, &qc) in qual.iter().enumerate() {
            if qc < PHRED33 {
                junk_chars += 1;
            }
            if !detect {
                continue;
            }
            let q = qc as i32 - offset;
            let flip_to_64 = offset == PHRED33 as i32
                && (q > QUAL_THRESH || (seq[i] == b'N' && (q == 31 || q == 33)));
            let flip_to_33 = offset == PHRED64 as i32 && q < -5;
            if flip_to_64 {
                offset = PHRED64 as i32;
                flips += 1;
            } else if flip_to_33 {
                offset = PHRED33 as i32;
                flips += 1;
            }
            if flips == 2 {
                detect = false;
                offset = PHRED33 as i32;
            }
        }
    }
    if junk_chars > 0 {
        PHRED33
    } else {
        offset as u8
    }
}

/// A/C/G/T/U to a 2-bit code; `None` for anything else.
pub fn base_to_number(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// ASCII-33 FASTQ quality to phred (BBTools `Vector.applyQualOffset` with
/// delta -33): 0 for non-ACGT bases, otherwise at least 2.
pub fn to_phred(bases: &[u8], quals: &[u8]) -> Vec<u8> {
    if quals.is_empty() {
        return Vec::new();
    }
    quals
        .iter()
        .enumerate()
        .map(|(i, &q)| {
            let q = q as i16 - 33;
            if base_to_number(bases[i]).is_some() {
                q.max(2) as u8
            } else {
                0
            }
        })
        .collect()
}

/// Phred to ASCII-33 FASTQ quality.
pub fn from_phred(quals: &[u8]) -> Vec<u8> {
    quals.iter().map(|&q| q.saturating_add(33)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::fmt::seq::SeqReader;
    use std::io::Cursor;

    fn records_from(fastq: &str) -> Vec<SeqRecord> {
        let mut reader = SeqReader::from_reader(Box::new(Cursor::new(fastq.as_bytes())));
        let mut rec = SeqRecord::new();
        let mut out = Vec::new();
        while reader.read_record(&mut rec).unwrap() {
            out.push(rec.clone());
        }
        out
    }

    #[test]
    fn detect_phred33_from_high_quality_sanger() {
        // Chars up to 'J' (74) never exceed QUAL_THRESH (87).
        let recs = records_from("@r\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
        assert_eq!(detect_quality_base(&recs), 33);
    }

    #[test]
    fn detect_phred64_from_old_illumina() {
        // Chars above 87 (e.g. 'Y' = 89) imply Q>54 under +33 -> flip to 64.
        let recs = records_from("@r\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\nYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY\n");
        assert_eq!(detect_quality_base(&recs), 64);
    }

    #[test]
    fn detect_forces_phred33_for_long_reads() {
        // 200+ bp with chars that would otherwise flip to +64.
        let seq = "A".repeat(200);
        let qual = "Y".repeat(200);
        let recs = records_from(&format!("@r\n{seq}\n+\n{qual}\n"));
        assert_eq!(detect_quality_base(&recs), 33);
    }

    #[test]
    fn detect_n_with_at_or_b_quality_flips_to_64() {
        // N with '@' (64) or 'B' (66) is a +64 signal.
        let recs = records_from("@r\nANAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\nA@AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
        assert_eq!(detect_quality_base(&recs), 64);
    }
}
