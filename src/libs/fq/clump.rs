//! BBTools-compatible read clumping (clumpify sort order).
//!
//! Ports BBTools 39.38 `clump.Clumpify` / `KmerSort1` with default options
//! (k=31, seed=1, hashes=4, border=1, no dedupe): interleaved pairs are
//! sorted by the pivot k-mer of R1, which is byte-identical to
//! `clumpify.sh` output (see tests/bbtools/Lambda/README.md).

use crate::libs::fmt::fq::write_fq;
use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use anyhow::Result;
use std::cmp::Ordering;
use std::io::Write;

/// 62-bit k-mer mask for k=31 (2*k bits).
const KMER_MASK: i64 = (1i64 << 62) - 1;
/// Default BBTools k-mer size.
pub const DEFAULT_K: usize = 31;
/// Default BBTools comparator seed.
pub const DEFAULT_SEED: u64 = 1;
/// Default BBTools hash count.
const DEFAULT_HASHES: usize = 4;
/// Default BBTools border size.
const DEFAULT_BORDER: usize = 1;

/// Sorts paired reads by R1 pivot k-mer and writes them in the resulting
/// order, reproducing `clumpify.sh` default output. `infiles` is one
/// interleaved FASTQ or two files (R1, R2).
pub fn clump<W: Write>(infiles: &[String], out: &mut W, k: usize, seed: u64) -> Result<()> {
    let pairs = read_pairs(infiles)?;
    let sorted = sort_pairs(pairs, k, seed);
    for (r1, r2) in &sorted {
        write_record(out, r1)?;
        if let Some(r2) = r2 {
            write_record(out, r2)?;
        }
    }
    Ok(())
}

/// Reads interleaved pairs (1 file) or R1/R2 pairs (2 files).
pub(crate) fn read_pairs(infiles: &[String]) -> Result<Vec<(SeqRecord, Option<SeqRecord>)>> {
    let mut pairs: Vec<(SeqRecord, Option<SeqRecord>)> = Vec::new();
    if infiles.len() == 2 {
        let mut reader1 = SeqReader::new(&infiles[0])?;
        let mut reader2 = SeqReader::new(&infiles[1])?;
        let mut rec1 = SeqRecord::new();
        let mut rec2 = SeqRecord::new();
        loop {
            if !reader1.read_record(&mut rec1)? {
                break;
            }
            let has2 = reader2.read_record(&mut rec2)?;
            pairs.push((rec1.clone(), has2.then(|| rec2.clone())));
            if !has2 {
                break;
            }
        }
    } else {
        let mut reader = SeqReader::new(&infiles[0])?;
        let mut rec1 = SeqRecord::new();
        let mut rec2 = SeqRecord::new();
        loop {
            if !reader.read_record(&mut rec1)? {
                break;
            }
            let has2 = reader.read_record(&mut rec2)?;
            pairs.push((rec1.clone(), has2.then(|| rec2.clone())));
            if !has2 {
                break;
            }
        }
    }
    Ok(pairs)
}

/// Sorts pairs by R1 pivot k-mer (clumpify-compatible order).
fn sort_pairs(
    pairs: Vec<(SeqRecord, Option<SeqRecord>)>,
    k: usize,
    seed: u64,
) -> Vec<(SeqRecord, Option<SeqRecord>)> {
    let codes = make_codes(seed);
    let prob = prob_error();
    let mut keyed: Vec<(ReadKey, usize)> = pairs
        .iter()
        .enumerate()
        .map(|(i, (r1, _))| (fill_max(r1.sequence(), k, &codes), i))
        .collect();
    keyed.sort_by(|(ka, ia), (kb, ib)| compare(&pairs, *ia, *ib, ka, kb, &prob));
    keyed.into_iter().map(|(_, i)| pairs[i].clone()).collect()
}

/// Pivot k-mer of a read: the highest-hash canonical k-mer window.
#[derive(Debug, Clone, Copy)]
struct ReadKey {
    kmer: i64,
    position: usize,
    minus: bool,
}

/// Scans a read for its maximum-hash canonical k-mer (KmerComparator.fillMax).
fn fill_max(seq: &[u8], k: usize, codes: &[[i64; 256]]) -> ReadKey {
    if seq.len() < k {
        return fill_short(seq, k);
    }
    let border = if seq.len() > k + 4 * DEFAULT_BORDER {
        DEFAULT_BORDER
    } else {
        0
    };
    let shift2 = (2 * k - 2) as u32;
    let mut kmer = 0i64;
    let mut rkmer = 0i64;
    let mut len = 0usize;
    let mut top_code = i64::MIN;
    let mut top = None;
    let max = seq.len() - border;
    for (i, &b) in seq.iter().enumerate().take(max).skip(border) {
        let x = base_to_number(b);
        let x2 = base_to_complement(b);
        kmer = ((kmer << 2) | x) & KMER_MASK;
        rkmer = ((rkmer >> 2) | (x2 << shift2)) & KMER_MASK;
        if x < 0 {
            len = 0;
        } else {
            len += 1;
        }
        if len >= k {
            let kmax = kmer.max(rkmer);
            let code = hash(kmax, codes);
            if code > top_code {
                top_code = code;
                top = Some((kmax, i, kmax != kmer));
            }
        }
    }
    match top {
        Some((kmer, position, minus)) => ReadKey {
            kmer,
            position,
            minus,
        },
        None => fill_short(seq, k),
    }
}

/// Key for reads shorter than k (KmerComparator.fillShort).
fn fill_short(seq: &[u8], k: usize) -> ReadKey {
    let max = seq.len().min(k);
    let shift2 = (2 * k - 2) as u32;
    let mut kmer = 0i64;
    let mut rkmer = 0i64;
    for &b in &seq[..max] {
        let x = base_to_number0(b);
        let x2 = base_to_complement0(b);
        kmer = ((kmer << 2) | x) & KMER_MASK;
        rkmer = ((rkmer >> 2) | (x2 << shift2)) & KMER_MASK;
    }
    let kmax = kmer.max(rkmer);
    ReadKey {
        kmer: kmax,
        position: max.saturating_sub(1),
        minus: kmax != kmer,
    }
}

/// BBTools k-mer hash: XOR of `hashes` code bytes (KmerComparator.hash).
fn hash(kmer: i64, codes: &[[i64; 256]]) -> i64 {
    let mut code = kmer;
    let mut k = kmer;
    for row in codes.iter().take(DEFAULT_HASHES) {
        let x = (k & 0xFF) as usize;
        k >>= 8;
        code ^= row[x];
    }
    code & i64::MAX
}

/// Full comparator: key, sequence, expected errors, then header id.
fn compare(
    pairs: &[(SeqRecord, Option<SeqRecord>)],
    ia: usize,
    ib: usize,
    ka: &ReadKey,
    kb: &ReadKey,
    prob: &[f32; 128],
) -> Ordering {
    let order = compare_key(ka, kb);
    if order != Ordering::Equal {
        return order;
    }
    let a = &pairs[ia].0;
    let b = &pairs[ib].0;
    let a2 = pairs[ia].1.as_ref();
    let b2 = pairs[ib].1.as_ref();
    let order = compare_sequence(a, a2, b, b2);
    if order != Ordering::Equal {
        return order;
    }
    let ea = expected_errors(a, a2, prob);
    let eb = expected_errors(b, b2, prob);
    if ea != eb {
        return if ea > eb {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }
    id(a).cmp(&id(b))
}

/// ReadKey ordering: bigger k-mer first, plus strand first, bigger position first.
fn compare_key(a: &ReadKey, b: &ReadKey) -> Ordering {
    if a.kmer != b.kmer {
        return b.kmer.cmp(&a.kmer);
    }
    if a.minus != b.minus {
        return if a.minus {
            Ordering::Greater
        } else {
            Ordering::Less
        };
    }
    b.position.cmp(&a.position)
}

/// Sequence comparison: longer first, then bytes, then mate, then quality sum.
fn compare_sequence(
    a1: &SeqRecord,
    a2: Option<&SeqRecord>,
    b1: &SeqRecord,
    b2: Option<&SeqRecord>,
) -> Ordering {
    compare_bytes(a1.sequence(), b1.sequence())
        .then_with(|| match (a2, b2) {
            (Some(a), Some(b)) => compare_bytes(a.sequence(), b.sequence()),
            (None, Some(_)) => Ordering::Greater,
            (Some(_), None) => Ordering::Less,
            (None, None) => Ordering::Equal,
        })
        .then_with(|| {
            // BBTools sums the entry (R1) quality only, not the mate.
            let qa = quality_sum(a1);
            let qb = quality_sum(b1);
            qb.cmp(&qa)
        })
}

/// Longer sequence first, then byte-wise ascending.
fn compare_bytes(a: &[u8], b: &[u8]) -> Ordering {
    match a.len().cmp(&b.len()) {
        Ordering::Less => Ordering::Greater,
        Ordering::Greater => Ordering::Less,
        Ordering::Equal => a.cmp(b),
    }
}

/// Sum of raw quality bytes.
fn quality_sum(rec: &SeqRecord) -> u32 {
    rec.quality_scores().iter().map(|&q| q as u32).sum()
}

/// Expected error count of a pair (float accumulation like BBTools).
fn expected_errors(a: &SeqRecord, a2: Option<&SeqRecord>, prob: &[f32; 128]) -> f32 {
    let mut sum = 0f32;
    let mut add = |rec: &SeqRecord| {
        for (i, &b) in rec.sequence().iter().enumerate() {
            let q = if b == b'A' || b == b'C' || b == b'G' || b == b'T' {
                rec.quality_scores().get(i).copied().unwrap_or(0) as usize
            } else {
                0
            };
            sum += prob[q.min(127)];
        }
    };
    add(a);
    if let Some(a2) = a2 {
        add(a2);
    }
    sum
}

/// Header line without the leading `@`, used as the final tie-breaker.
fn id(rec: &SeqRecord) -> Vec<u8> {
    let comment = rec.comment();
    if comment.is_empty() {
        rec.name().to_vec()
    } else {
        let mut v = rec.name().to_vec();
        v.push(b' ');
        v.extend_from_slice(comment);
        v
    }
}

/// Writes a FASTQ record, preserving the `name comment` header layout.
fn write_record<W: Write>(w: &mut W, rec: &SeqRecord) -> anyhow::Result<()> {
    let comment = rec.comment();
    let header = if comment.is_empty() {
        rec.name().to_string()
    } else {
        format!("{} {}", rec.name(), comment)
    };
    write_fq(w, &header, rec.sequence(), rec.quality_scores())?;
    Ok(())
}

/// 2-bit code tables (dna.AminoAcid).
fn base_to_number(b: u8) -> i64 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => -1,
    }
}

fn base_to_complement(b: u8) -> i64 {
    match b {
        b'A' | b'a' => 3,
        b'C' | b'c' => 2,
        b'G' | b'g' => 1,
        b'T' | b't' | b'U' | b'u' => 0,
        _ => -1,
    }
}

fn base_to_number0(b: u8) -> i64 {
    match base_to_number(b) {
        -1 => 0,
        x => x,
    }
}

fn base_to_complement0(b: u8) -> i64 {
    match base_to_complement(b) {
        -1 => 0,
        x => x,
    }
}

/// `align2.QualityTools.PROB_ERROR` (f32, 128 entries).
fn prob_error() -> [f32; 128] {
    let mut r = [0f32; 128];
    for (i, v) in r.iter_mut().enumerate() {
        *v = (10f64.powf(-0.1 * i as f64)) as f32;
    }
    r[0] = 0.75;
    r[1] = 0.7;
    r
}

/// java.util.Random, byte-compatible with the JVM implementation.
struct JavaRandom {
    seed: u64,
}

impl JavaRandom {
    fn new(seed: u64) -> Self {
        Self {
            seed: (seed ^ 0x5DEECE66D) & ((1u64 << 48) - 1),
        }
    }

    fn next(&mut self, bits: u32) -> u32 {
        self.seed = self.seed.wrapping_mul(0x5DEECE66D).wrapping_add(0xB) & ((1u64 << 48) - 1);
        (self.seed >> (48 - bits)) as u32
    }

    fn next_long(&mut self) -> i64 {
        // Java: ((long)next(32) << 32) + next(32); the low half is added as a
        // sign-extended int, so a set high bit subtracts 1 from the top half.
        let hi = self.next(32) as i32 as i64;
        let lo = self.next(32) as i32 as i64;
        (hi << 32).wrapping_add(lo)
    }

    fn next_int(&mut self, bound: u32) -> u32 {
        if bound.is_power_of_two() {
            ((bound as u64 * self.next(31) as u64) >> 31) as u32
        } else {
            loop {
                let bits = self.next(31) as i64;
                let val = bits % bound as i64;
                if bits - val + (bound as i64 - 1) >= 0 {
                    return val as u32;
                }
            }
        }
    }
}

/// `SketchObject.makeCodes(8, 256, seed, true)` with antialiasing.
fn make_codes(seed: u64) -> [[i64; 256]; 8] {
    let mut rng = JavaRandom::new(seed);
    let mut r = [[0i64; 256]; 8];
    for row in r.iter_mut() {
        for v in row.iter_mut() {
            *v = rng.next_long() & i64::MAX;
        }
    }
    for _ in 0..3 {
        antialias(&mut r, &mut rng);
    }
    r
}

fn antialias(r: &mut [[i64; 256]; 8], rng: &mut JavaRandom) {
    for row in r.iter_mut() {
        for bit in 0..64 {
            antialias_numbers(row, rng);
            antialias_bit(row, rng, bit);
        }
    }
}

fn antialias_numbers(row: &mut [i64; 256], rng: &mut JavaRandom) {
    for v in row.iter_mut() {
        while v.count_ones() < 31 {
            *v |= 1i64 << rng.next_int(64);
        }
        while v.count_ones() > 33 {
            *v &= !(1i64 << rng.next_int(64));
        }
    }
}

fn antialias_bit(row: &mut [i64; 256], rng: &mut JavaRandom, bit: u32) {
    let or_mask = 1i64 << bit;
    let and_mask = !or_mask;
    let mut ones = row.iter().map(|&v| (v >> bit) & 1).sum::<i64>();
    while ones < 127 {
        let mut loc = rng.next_int(256) as usize;
        while row[loc] & or_mask != 0 {
            loc = rng.next_int(256) as usize;
        }
        row[loc] |= or_mask;
        ones += 1;
    }
    while ones > 129 {
        let mut loc = rng.next_int(256) as usize;
        while row[loc] & or_mask == 0 {
            loc = rng.next_int(256) as usize;
        }
        row[loc] &= and_mask;
        ones -= 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn java_random_matches_jvm_sequence() {
        let mut rng = JavaRandom::new(1);
        assert_eq!(rng.next_long(), 0xbb1ad57319b89cd8u64 as i64);
        assert_eq!(rng.next_long(), 0x68fb0e6f684df992u64 as i64);
        assert_eq!(rng.next_long(), 0x352cccfc0946b8f0u64 as i64);
        // The low half has its high bit set, exercising the sign-extended
        // addition in Java's nextLong (top half is E4, not E5).
        assert_eq!(rng.next_long(), 0x552cf1e4a8ab85ddu64 as i64);
        assert_eq!(rng.next_int(256), 247);
        assert_eq!(rng.next_int(64), 45);
    }

    #[test]
    fn prob_error_table_matches_bbtools() {
        let p = prob_error();
        assert_eq!(p[0], 0.75);
        assert_eq!(p[1], 0.7);
        assert!((p[2] - 0.63095736).abs() < 1e-6);
        assert!((p[30] - 0.001).abs() < 1e-7);
    }
}
