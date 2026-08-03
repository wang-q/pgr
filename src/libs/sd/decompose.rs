//! Elementary SD decomposition from cluster FASTA (BISER decompose, k-mer based).

use crate::libs::ds::Dsu;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

/// k-mer length used for elementary SD seeding.
pub const KMER: usize = 10;
/// Merge shared k-mer positions within this gap into one elementary SD.
pub const MAX_GAP: usize = 50;
/// Minimum elementary SD length in bp.
pub const MIN_LEN: usize = 100;
/// Minimum number of distinct shared k-mers for two fragments to be grouped
/// into the same elementary SD set (avoids over-grouping via conserved kmers).
pub const MIN_SHARED_KMERS: u32 = 5;

/// One elementary SD fragment, projected to genome coordinates.
#[derive(Debug, Clone)]
pub struct ElemSd {
    pub species: String,
    pub chrom: String,
    /// 0-based half-open genome coordinates (strand-corrected).
    pub begin: usize,
    pub end: usize,
    pub set_id: u32,
    pub length: usize,
    /// Number of shared k-mer hits inside the fragment.
    pub score: u32,
    pub strand: char,
}

/// Parsed cluster FASTA record: (species, chrom, strand, gstart, gend, sequence).
type SeqRecord = (String, String, char, usize, usize, Vec<u8>);

/// Parse a cluster FASTA header `{species}#{chrom}{strand}#{start}#{end}`.
fn parse_header(header: &str) -> Option<(String, String, char, usize, usize)> {
    let mut parts = header.split('#');
    let species = parts.next()?.to_string();
    let chrom_strand = parts.next()?;
    let strand = chrom_strand.chars().last()?;
    if strand != '+' && strand != '-' {
        return None;
    }
    let chrom = chrom_strand[..chrom_strand.len() - 1].to_string();
    let start: usize = parts.next()?.parse().ok()?;
    let end: usize = parts.next()?.parse().ok()?;
    if end < start {
        return None; // reversed/empty interval
    }
    Some((species, chrom, strand, start, end))
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
/// [`MIN_LEN`] are dropped. Fragments that share a k-mer are unioned into the
/// same elementary SD set (BISER `set_id` semantics); coordinates are
/// projected back to the genome using the cluster header's start/strand.
pub fn decompose_fasta<R: BufRead, W: Write>(reader: R, writer: &mut W) -> anyhow::Result<()> {
    let seqs = parse_fasta(reader)?;
    if seqs.is_empty() {
        return Ok(());
    }

    // k-mer index: hash -> positions (seq_idx, pos).
    let mut index: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
    for (si, (_, _, _, _, _, s)) in seqs.iter().enumerate() {
        let hs = kmer_hashes(s);
        for (pos, h) in hs.iter().enumerate() {
            index.entry(*h).or_default().push((si, pos));
        }
    }

    // Mark shared k-mer positions (hash seen in >= 2 distinct sequences).
    let mut shared: Vec<Vec<bool>> = seqs
        .iter()
        .map(|(_, _, _, _, _, s)| vec![false; s.len()])
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
    let mut frags: Vec<(usize, usize, usize, u32)> = Vec::new();
    let mut pos_to_frag: Vec<Vec<Option<usize>>> = seqs
        .iter()
        .map(|(_, _, _, _, _, s)| vec![None; s.len()])
        .collect();
    for (si, (_, _, _, _, _, s)) in seqs.iter().enumerate() {
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
                let frag_id = frags.len();
                for p in begin..end {
                    if shared[si][p] {
                        pos_to_frag[si][p] = Some(frag_id);
                    }
                }
                frags.push((si, begin, end, score));
            }
        }
    }

    // Fragments sharing >= MIN_SHARED_KMERS distinct k-mers belong to the same
    // elementary SD set. A single conserved k-mer is not enough (it would
    // over-group every fragment in the cluster).
    let mut dsu = Dsu::new(frags.len());
    let mut kmer_frags: HashMap<u64, Vec<usize>> = HashMap::new();
    for (h, positions) in &index {
        let mut fs: Vec<usize> = positions
            .iter()
            .filter_map(|&(si, pos)| pos_to_frag[si][pos])
            .collect();
        fs.sort_unstable();
        fs.dedup();
        if fs.len() >= 2 {
            kmer_frags.insert(*h, fs);
        }
    }
    let mut pair_count: HashMap<(usize, usize), u32> = HashMap::new();
    for fs in kmer_frags.values() {
        for i in 0..fs.len() {
            for j in (i + 1)..fs.len() {
                *pair_count.entry((fs[i], fs[j])).or_default() += 1;
            }
        }
    }
    for ((a, b), c) in pair_count {
        if c >= MIN_SHARED_KMERS {
            dsu.union(a, b);
        }
    }

    // Assign set_ids by connected component (order of first occurrence).
    let mut set_of: Vec<u32> = vec![0; frags.len()];
    let mut next_set = 0u32;
    let mut seen_set: HashMap<usize, u32> = HashMap::new();
    for (i, sid_out) in set_of.iter_mut().enumerate() {
        let root = dsu.find(i);
        *sid_out = *seen_set.entry(root).or_insert_with(|| {
            next_set += 1;
            next_set
        });
    }

    // Project each fragment to genome coordinates and emit BED rows.
    for (frag_id, &(si, begin, end, score)) in frags.iter().enumerate() {
        let (species, chrom, strand, gstart, gend, _s) = &seqs[si];
        let (gb, ge) = if *strand == '-' {
            // Saturating projection: a malformed header (span smaller than
            // the sequence, or coordinates beyond the contig) must not
            // underflow into a panic (Zero Panic).
            (gend.saturating_sub(end), gend.saturating_sub(begin))
        } else {
            (gstart.saturating_add(begin), gstart.saturating_add(end))
        };
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            species,
            chrom,
            gb,
            ge,
            set_of[frag_id],
            ge - gb,
            score,
            strand
        )?;
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
                if let Some((sp, chr, strand, gs, ge)) = parse_header(&name) {
                    seqs.push((sp, chr, strand, gs, ge, std::mem::take(&mut seq)));
                } else {
                    log::warn!("skipping cluster record with unparseable header: {name}");
                }
            }
            name = rest.trim().to_string();
            seq.clear();
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !name.is_empty() {
        if let Some((sp, chr, strand, gs, ge)) = parse_header(&name) {
            seqs.push((sp, chr, strand, gs, ge, seq));
        } else {
            log::warn!("skipping cluster record with unparseable header: {name}");
        }
    }
    Ok(seqs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_cluster_header() {
        let (sp, chr, st, gs, ge) = parse_header("mg1655#NC_000913+#100#200").unwrap();
        assert_eq!(sp, "mg1655");
        assert_eq!(chr, "NC_000913");
        assert_eq!(st, '+');
        assert_eq!((gs, ge), (100, 200));
        let (_, _, st, _, _) = parse_header("mg1655#NC_000913-#100#200").unwrap();
        assert_eq!(st, '-');
    }

    #[test]
    fn kmer_hash_len() {
        let seq = b"ACGTACGTACGTACGTACGT";
        assert_eq!(kmer_hashes(seq).len(), seq.len() - KMER + 1);
    }

    #[test]
    fn decompose_detects_shared_fragment() {
        // Two sequences sharing a 152 bp non-periodic fragment (deterministic
        // LCG); cluster spans genome [100, 100+len).
        let mut x = 12345u64;
        let shared: String = (0..152)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4] as char
            })
            .collect();
        let a = format!("TTTT{}AAAA", shared);
        let b = format!("CCCC{}GGGG", shared);
        let fa = format!(
            ">sp#chr+#100#{}\n{}\n>sp#chr+#100#{}\n{}\n",
            100 + a.len(),
            a,
            100 + b.len(),
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
        // Both fragments belong to the same elementary SD set.
        let set_ids: HashSet<&str> = rows.iter().map(|r| r.split('\t').nth(4).unwrap()).collect();
        assert_eq!(set_ids.len(), 1, "expected shared set_id, got {text}");
        for row in &rows {
            let fields: Vec<&str> = row.split('\t').collect();
            let begin: usize = fields[2].parse().unwrap();
            let end: usize = fields[3].parse().unwrap();
            assert!(begin >= 100 && end - begin >= 100, "unexpected row: {row}");
        }
    }

    #[test]
    fn malformed_header_does_not_panic() {
        // Regression: a minus-strand header with start > end (or a header
        // span smaller than the sequence) used to underflow the projection
        // (`gend - end`) into an arithmetic panic. It must be rejected or
        // clamped, never panic.
        let shared: String = (0..120)
            .map(|i| b"ACGT"[(i % 4) as usize] as char)
            .collect();
        // start=100 > end=50 on the minus strand.
        let reversed = format!(
            ">sp#chr-#100#50\nAAAA{}TTTT\n>sp#chr+#0#200\nCCCC{}GGGG\n",
            shared, shared
        );
        let mut out = Vec::new();
        decompose_fasta(std::io::Cursor::new(reversed.as_bytes()), &mut out).unwrap();
        // The reversed header is rejected at parse time.
        assert!(out.is_empty(), "reversed header must be skipped: {out:?}");

        // Header span (10 bp) smaller than the 200 bp sequence: the shared
        // fragment [40, 160) projects beyond the declared span; the row must
        // stay in range (clamped), not panic.
        let short_span = format!(
            ">sp#chr+#0#10\nAAAA{}TTTT\n>sp#chr+#0#200\nCCCC{}GGGG\n",
            shared, shared
        );
        let mut out = Vec::new();
        decompose_fasta(std::io::Cursor::new(short_span.as_bytes()), &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        for row in text.lines() {
            let fields: Vec<&str> = row.split('\t').collect();
            let begin: usize = fields[2].parse().unwrap();
            let end: usize = fields[3].parse().unwrap();
            assert!(end >= begin, "clamped row must stay ordered: {row}");
        }
    }
}
