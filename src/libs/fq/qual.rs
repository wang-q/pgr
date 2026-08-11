//! FASTQ quality conversion (BBTools `applyQualOffset` semantics).

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
