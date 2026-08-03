use crate::libs::nt::{self, Nt, NT_VAL};

/// block -> row -> column
pub static AA_TAB: &[[[char; 4]; 4]; 4] = &[
    [
        ['K', 'N', 'K', 'N'], // AAA, AAC, AAG, AAU/AAT
        ['T', 'T', 'T', 'T'], // ACA, ACC, ACG, ACU/ACT
        ['R', 'S', 'R', 'S'], // AGA, AGC, AGG, AGU/AGT
        ['I', 'I', 'M', 'I'], // AUA/ATA, AUC/ATC, AUG/ATG, AUU/ATT
    ],
    [
        ['Q', 'H', 'Q', 'H'], // CAA, CAC, CAG, CAU/CAT
        ['P', 'P', 'P', 'P'], // CCA, CCC, CCG, CCU/CCT
        ['R', 'R', 'R', 'R'], // CGA, CGC, CGG, CGU/CGT
        ['L', 'L', 'L', 'L'], // CUA/CTA, CUC/CTC, CUG/CTG, CUU/CTT
    ],
    [
        ['E', 'D', 'E', 'D'], // GAA, GAC, GAG, GAU/GAT
        ['A', 'A', 'A', 'A'], // GCA, GCC, GCG, GCU/GCT
        ['G', 'G', 'G', 'G'], // GGA, GGC, GGG, GGU/GGT
        ['V', 'V', 'V', 'V'], // GUA/GTA, GUC/GTC, GUG/GTG, GUU/GTT
    ],
    [
        ['*', 'Y', '*', 'Y'], // UAA/TAA, UAC/TAC, UAG/TAG, UAU/TAT
        ['S', 'S', 'S', 'S'], // UCA/TCA, UCC/TCC, UCG/TCG, UCU/TCT
        ['*', 'C', 'W', 'C'], // UGA/TGA, UGC/TGC, UGG/TGG, UGU/TGT
        ['L', 'F', 'L', 'F'], // UUA/TTA, UUC/TTC, UUG/TTG, UUU/TTT
    ],
];

/// ```
/// let dna = b"GCTAGTCGTATCGTAGCTAGTC";
/// assert_eq!(&pgr::libs::translate::translate(dna), "ASRIVAS");
///
/// let rna = b"GCUAGUCGUAUCGUAGCUAGUC";
/// assert_eq!(&pgr::libs::translate::translate(rna), "ASRIVAS");
///
/// // To shift the reading frame, pass in a slice
/// assert_eq!(&pgr::libs::translate::translate(&dna[1..]), "LVVS*LV");
/// assert_eq!(&pgr::libs::translate::translate(&dna[2..]), "*SYRS*");
/// ```
// https://github.com/dweb0/protein-translate/blob/master/src/lib.rs
pub fn translate(seq: &[u8]) -> String {
    let mut peptide = String::with_capacity(seq.len() / 3);

    for triplet in seq.chunks_exact(3) {
        let c1 = NT_VAL[triplet[0] as usize];
        let c2 = NT_VAL[triplet[1] as usize];
        let c3 = NT_VAL[triplet[2] as usize];

        if c1 >= Nt::N as usize || c2 >= Nt::N as usize || c3 >= Nt::N as usize {
            peptide.push('X');
        } else {
            peptide.push(AA_TAB[c1][c2][c3]);
        }
    }
    peptide
}

/// Detect ORFs in a translated protein sequence
///
/// # Examples
///
/// ```
/// let protein = "MGGMGG*AGG";
/// let orfs = pgr::libs::translate::find_orfs(protein);
/// assert_eq!(orfs, vec![
///     ("MGGMGG*".to_string(), 0, 7),
///     ("AGG".to_string(), 7, 10)
/// ]);
/// ```
pub fn find_orfs(protein: &str) -> Vec<(String, usize, usize)> {
    let mut orfs = Vec::new();
    let mut start = 0;

    while start < protein.len() {
        // Find the start aa (not just M)
        if let Some(orf_start) = protein[start..].find(|c: char| c != 'X' && c != 'x') {
            let orf_start_pos = start + orf_start;
            // Find the stop codon '*'
            if let Some(orf_end) = protein[orf_start_pos..].find('*') {
                let orf_end_pos = orf_start_pos + orf_end;
                let orf_seq = &protein[orf_start_pos..=orf_end_pos];
                // Include the stop codon
                orfs.push((orf_seq.to_string(), orf_start_pos, orf_end_pos + 1));
                start = orf_end_pos + 1; // Continue searching for the next ORF
            } else {
                // If no stop codon is found, treat the remaining sequence as an ORF
                let orf_seq = &protein[orf_start_pos..];
                orfs.push((orf_seq.to_string(), orf_start_pos, protein.len()));
                break;
            }
        } else {
            break; // No start codon found, end the search
        }
    }

    orfs
}

/// Translate DNA in all six frames.
///
/// Returns `Vec<(protein, frame, is_reverse)>` where `frame` is the offset
/// (0, 1, or 2) and `is_reverse` indicates whether the protein came from the
/// reverse-complement strand. Frames 0-2 are forward, 3-5 are reverse.
pub fn six_frame_translation(dna: &[u8]) -> Vec<(String, usize, bool)> {
    let mut translations = Vec::new();

    // Forward frames.
    for frame in 0..3 {
        let frame_dna = &dna[frame..];
        let protein = translate(frame_dna);
        translations.push((protein, frame, false));
    }

    // Reverse-complement frames.
    let dna_rc = nt::rev_comp(dna).collect::<Vec<_>>();
    for frame in 0..3 {
        let frame_dna = &dna_rc[frame..];
        let protein = translate(frame_dna);
        translations.push((protein, frame, true));
    }

    translations
}

/// Filter ORFs and convert protein coordinates to 1-based DNA coordinates.
/// Returns (orf_start, orf_end, orf_seq) tuples in 1-based DNA coords.
pub fn filter_and_convert_orfs(
    orfs: &[(String, usize, usize)],
    dna_len: usize,
    frame: usize,
    is_reverse: bool,
    min_len: usize,
    require_start_m: bool,
    require_end_star: bool,
) -> Vec<(usize, usize, String)> {
    let dna_start = if is_reverse { dna_len - frame } else { frame };

    let mut result = Vec::new();
    for (orf_seq, start, end) in orfs {
        if orf_seq.len() < min_len {
            continue;
        }
        if require_start_m && !orf_seq.starts_with('M') {
            continue;
        }
        if require_end_star && !orf_seq.ends_with('*') {
            continue;
        }

        let (orf_start, orf_end) = if is_reverse {
            (dna_start - end * 3 + 1, dna_start - start * 3)
        } else {
            (dna_start + start * 3 + 1, dna_start + end * 3)
        };

        result.push((orf_start, orf_end, orf_seq.clone()));
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_basic() {
        let dna = b"GCTAGTCGTATCGTAGCTAGTC";
        assert_eq!(translate(dna), "ASRIVAS");
    }

    #[test]
    fn test_translate_with_n() {
        // N in first position
        assert_eq!(translate(b"NCT"), "X");
        // N in second position
        assert_eq!(translate(b"CNT"), "X");
        // N in third position
        assert_eq!(translate(b"CTN"), "X");
        // All N
        assert_eq!(translate(b"NNN"), "X");
    }

    #[test]
    fn test_translate_invalid_chars() {
        // Invalid char 'Z' (not in NT_VAL map for A/C/G/T/N/etc, usually Invalid)
        // 'Z' is ascii 90. NT_VAL[90] is 255 (Invalid).
        assert_eq!(translate(b"ZCT"), "X");
        assert_eq!(translate(b"CZT"), "X");
        assert_eq!(translate(b"CTZ"), "X");
    }

    #[test]
    fn test_translate_lower_case() {
        assert_eq!(translate(b"gct"), "A");
    }

    #[test]
    fn test_translate_mixed_case() {
        assert_eq!(translate(b"Gct"), "A");
    }

    #[test]
    fn test_translate_iupac() {
        // R (A or G) -> should be X in translation?
        // NT_VAL[R] = 4 (N). So it should be X.
        assert_eq!(translate(b"RCT"), "X");
    }

    #[test]
    fn test_translate_non_ascii() {
        // 255 is non-ascii
        assert_eq!(translate(&[255, b'A', b'T']), "X");
    }
}
