//! Elementary SD decomposition from cluster FASTA (BISER decompose, k-mer based).

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

/// k-mer length used for elementary SD seeding.
pub const KMER: usize = 10;
/// Merge shared k-mer positions within this gap into one elementary SD.
pub const MAX_GAP: usize = 50;
/// Minimum elementary SD length in bp.
pub const MIN_LEN: usize = 100;

/// One elementary SD fragment within a cluster sequence.
#[derive(Debug, Clone)]
pub struct ElemSd {
    pub species: String,
    pub chrom: String,
    /// 0-based half-open coordinates within the cluster sequence.
    pub begin: usize,
    pub end: usize,
    pub set_id: u32,
    pub length: usize,
    /// Number of shared k-mer hits inside the fragment.
    pub score: u32,
    pub strand: char,
}

/// Parsed cluster FASTA record: (species, chrom, strand, sequence).
type SeqRecord = (String, String, char, Vec<u8>);

/// Parse a cluster FASTA header `{species}#{chrom}{strand}#{start}#{end}`.
fn parse_header(header: &str) -> Option<(String, String, char)> {
    let mut parts = header.split('#');
    let species = parts.next()?.to_string();
    let chrom_strand = parts.next()?;
    let strand = chrom_strand.chars().last()?;
    if strand != '+' && strand != '-' {
        return None;
    }
    let chrom = chrom_strand[..chrom_strand.len() - 1].to_string();
    Some((species, chrom, strand))
}

/// 2-bit rolling k-mer hash (A=0, C=1, G=2, T=3; non-ACGT breaks the window).
fn kmer_hashes(seq: &[u8]) -> Vec<u64> {
    let mask = (1u64 << (2 * KMER)) - 1;
    let mut out = Vec::with_capacity(seq.len().saturating_sub(KMER) + 1);
    let mut h: u64 = 0;
    let mut valid = 0usize;
    for &b in seq {
        let code = match b {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        };
        if code == 4 {
            valid = 0;
            h = 0;
            continue;
        }
        h = ((h << 2) | code) & mask;
        valid += 1;
        if valid >= KMER {
            out.push(h);
        }
    }
    out
}

/// Decompose a cluster FASTA into elementary SD fragments and write BED rows.
///
/// Shared k-mer seeds (a k-mer present in >= 2 distinct sequences) are merged
/// into fragments with a gap tolerance of [`MAX_GAP`]; fragments shorter than
/// [`MIN_LEN`] are dropped. `set_id` is assigned per fragment in output order
/// (simplified vs BISER's cross-copy set grouping).
pub fn decompose_fasta<R: BufRead, W: Write>(reader: R, writer: &mut W) -> anyhow::Result<()> {
    let seqs = parse_fasta(reader)?;
    if seqs.is_empty() {
        return Ok(());
    }

    // k-mer index: hash -> positions (seq_idx, pos).
    let mut index: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
    for (si, (_, _, _, s)) in seqs.iter().enumerate() {
        let hs = kmer_hashes(s);
        for (pos, h) in hs.iter().enumerate() {
            index.entry(*h).or_default().push((si, pos));
        }
    }

    // Mark shared k-mer positions (hash seen in >= 2 distinct sequences).
    let mut shared: Vec<Vec<bool>> = seqs
        .iter()
        .map(|(_, _, _, s)| vec![false; s.len()])
        .collect();
    for positions in index.values() {
        let distinct: HashSet<usize> = positions.iter().map(|&(si, _)| si).collect();
        if distinct.len() < 2 {
            continue;
        }
        for &(si, pos) in positions {
            shared[si][pos] = true;
        }
    }

    // Merge shared runs into fragments (gap tolerance MAX_GAP).
    let mut set_id = 0u32;
    for (si, (species, chrom, strand, s)) in seqs.iter().enumerate() {
        let mut i = 0usize;
        while i < s.len() {
            if !shared[si][i] {
                i += 1;
                continue;
            }
            let begin = i;
            let mut score = 0u32;
            let mut last = i;
            while i < s.len() && (shared[si][i] || i - last <= MAX_GAP) {
                if shared[si][i] {
                    score += 1;
                    last = i;
                }
                i += 1;
            }
            let end = last + 1;
            if end - begin >= MIN_LEN {
                set_id += 1;
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    species,
                    chrom,
                    begin,
                    end,
                    set_id,
                    end - begin,
                    score,
                    strand
                )?;
            }
        }
    }
    Ok(())
}

/// Read a cluster FASTA into `SeqRecord`s.
fn parse_fasta<R: BufRead>(reader: R) -> anyhow::Result<Vec<SeqRecord>> {
    let mut seqs: Vec<SeqRecord> = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if !name.is_empty() {
                if let Some((sp, chr, strand)) = parse_header(&name) {
                    seqs.push((sp, chr, strand, std::mem::take(&mut seq)));
                }
            }
            name = rest.trim().to_string();
            seq.clear();
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !name.is_empty() {
        if let Some((sp, chr, strand)) = parse_header(&name) {
            seqs.push((sp, chr, strand, seq));
        }
    }
    Ok(seqs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_cluster_header() {
        let (sp, chr, st) = parse_header("mg1655#NC_000913+#100#200").unwrap();
        assert_eq!(sp, "mg1655");
        assert_eq!(chr, "NC_000913");
        assert_eq!(st, '+');
        let (_, _, st) = parse_header("mg1655#NC_000913-#100#200").unwrap();
        assert_eq!(st, '-');
    }

    #[test]
    fn kmer_hash_len() {
        let seq = b"ACGTACGTACGTACGTACGT";
        assert_eq!(kmer_hashes(seq).len(), seq.len() - KMER + 1);
    }

    #[test]
    fn decompose_detects_shared_fragment() {
        // Two sequences sharing a 150 bp identical fragment (shared k-mers),
        // plus distinct flanks so only the shared part survives MIN_LEN.
        let shared = "ACGT".repeat(38); // 152 bp
        let a = format!("TTTT{}AAAA", shared);
        let b = format!("CCCC{}GGGG", shared);
        let fa = format!(
            ">sp#chr+#0#{}\n{}\n>sp#chr+#0#{}\n{}\n",
            a.len(),
            a,
            b.len(),
            b
        );
        let mut out = Vec::new();
        decompose_fasta(std::io::Cursor::new(fa.as_bytes()), &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let rows: Vec<&str> = text.lines().collect();
        assert_eq!(
            rows.len(),
            2,
            "expected one fragment per sequence, got {text}"
        );
        for row in &rows {
            let fields: Vec<&str> = row.split('\t').collect();
            let len: usize = fields[5].parse().unwrap();
            assert!(len >= 100, "fragment should cover the shared region: {row}");
        }
    }
}
