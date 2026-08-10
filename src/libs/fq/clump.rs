//! BBTools-compatible read clumping (clumpify sort order).
//!
//! Ports BBTools 39.38 `clump.Clumpify` / `KmerSort1` with default options
//! (k=31, seed=1, hashes=4, border=1, no dedupe): interleaved pairs are
//! sorted by the pivot k-mer of R1, which is byte-identical to
//! `clumpify.sh` output (see tests/bbtools/Lambda/README.md).

use crate::libs::fmt::fq::write_fq;
use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use anyhow::{Context, Result};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

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
/// Input quality ASCII offset (BBTools stores phred internally).
const QUALITY_BASE: u8 = 33;
/// Maximum non-matching reads scanned past during dedupe (`scanlimit`).
const SCAN_LIMIT: usize = 5;
/// Options for the clumpify-compatible sort.
#[derive(Debug, Clone)]
pub struct ClumpOptions {
    /// K-mer size (clumpify `k`).
    pub k: usize,
    /// Comparator seed (clumpify `seed`).
    pub seed: u64,
    /// Remove duplicate read pairs.
    pub dedupe: bool,
    /// Maximum substitutions allowed in a duplicate (clumpify `dupesubs`).
    pub dupesubs: usize,
    /// User memory cap in bytes (`--mem`); default 2 GiB.
    pub mem: Option<u64>,
    /// Override the external-path bucket count.
    pub buckets: Option<usize>,
    /// Force the sorting path; `Auto` decides from the memory budget.
    pub mode: SortMode,
    /// Worker threads for the parallel sort/buckets.
    pub parallel: usize,
}

/// Sorting path selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortMode {
    /// Pick the in-memory path when the estimate fits the budget.
    Auto,
    /// Always sort globally in memory.
    Global,
    /// Always use the external hash-bucket path.
    Bucket,
}

/// Sorts paired reads by R1 pivot k-mer and writes them in the resulting
/// order, reproducing `clumpify.sh` default output. `infiles` is one
/// interleaved FASTQ or two files (R1, R2). With `dedupe`, whole-pair
/// duplicates (R1 and R2 both matching within `dupesubs` substitutions) are
/// removed, keeping the higher-quality copy.
pub fn clump<W: Write + Send>(infiles: &[String], out: &mut W, opts: &ClumpOptions) -> Result<()> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.parallel.max(1))
        .build()
        .context("failed to build rayon pool")?;
    pool.install(|| {
        let estimate = estimate_fastq_bytes(infiles);
        let cap = crate::libs::sys::mem_cap(opts.mem);
        match opts.mode {
            SortMode::Auto => {
                if estimate <= cap {
                    clump_in_memory(infiles, out, opts)?;
                } else {
                    clump_buckets(infiles, out, opts, cap, estimate)?;
                }
            }
            SortMode::Global => clump_in_memory(infiles, out, opts)?,
            SortMode::Bucket => clump_buckets(infiles, out, opts, cap, estimate)?,
        }
        Ok(())
    })
}

/// In-memory path: load all pairs, sort, optional dedupe.
fn clump_in_memory<W: Write>(infiles: &[String], out: &mut W, opts: &ClumpOptions) -> Result<()> {
    let pairs = read_pairs(infiles)?;
    let keyed = sort_keyed(pairs, opts.k, opts.seed);
    let kept = if opts.dedupe {
        dedupe_pairs(keyed, opts.dupesubs)
    } else {
        keyed.into_iter().map(|k| k.pair).collect()
    };
    for (r1, r2) in &kept {
        write_record(out, r1)?;
        if let Some(r2) = r2 {
            write_record(out, r2)?;
        }
    }
    Ok(())
}

/// External path: hash pairs into buckets by pivot k-mer, sort each bucket in
/// memory, and emit buckets in order. Deterministic for fixed options; the
/// order is bucket-concatenated (documented divergence from the in-memory
/// global order, matching BBTools' own large-data behavior).
fn clump_buckets<W: Write>(
    infiles: &[String],
    out: &mut W,
    opts: &ClumpOptions,
    cap: u64,
    estimate: u64,
) -> Result<()> {
    let buckets = opts.buckets.unwrap_or_else(|| {
        let per = ((cap as f64 * 0.8).max(1.0)) as u64;
        ((estimate.div_ceil(per)).max(2) as usize).min(4096)
    });
    let tmp = temp_dir_for();
    std::fs::create_dir_all(&tmp)
        .with_context(|| format!("failed to create temp dir {}", tmp.display()))?;

    let mut writers: Vec<Option<BufWriter<File>>> = Vec::with_capacity(buckets);
    writers.resize_with(buckets, || None);
    let codes = make_codes(opts.seed);
    let write_result = for_each_pair(infiles, |r1, r2| {
        let key = fill_max(r1.sequence(), opts.k, &codes);
        let bucket = (key.kmer as u64 % buckets as u64) as usize;
        if writers[bucket].is_none() {
            let path = tmp.join(format!("bucket_{bucket:05}.fq"));
            let f = File::create(&path)
                .with_context(|| format!("failed to create {}", path.display()))?;
            writers[bucket] = Some(BufWriter::new(f));
        }
        let w = writers[bucket].as_mut().unwrap();
        write_record(w, r1)?;
        if let Some(r2) = r2 {
            write_record(w, r2)?;
        }
        Ok(())
    });
    for w in writers.iter_mut().flatten() {
        w.flush()?;
    }
    drop(writers);
    write_result?;

    // Process buckets in parallel, but in memory-bounded waves so concurrent
    // bucket sorts cannot exceed the memory budget.
    let bucket_estimate = (estimate / buckets as u64).max(1);
    let wave = ((cap as f64 * 0.8) / bucket_estimate as f64)
        .ceil()
        .max(1.0) as usize;
    let process_result = (|| -> Result<()> {
        let mut b = 0usize;
        while b < buckets {
            let end = (b + wave).min(buckets);
            let chunk: Vec<usize> = (b..end).collect();
            let results: Vec<(usize, Result<Vec<u8>>)> = chunk
                .par_iter()
                .map(|&bi| {
                    let path = tmp.join(format!("bucket_{bi:05}.fq"));
                    let bytes = (|| -> Result<Vec<u8>> {
                        if !path.exists() {
                            return Ok(Vec::new());
                        }
                        let pairs = read_pairs(&[path.to_string_lossy().into_owned()])?;
                        let keyed = sort_keyed(pairs, opts.k, opts.seed);
                        let kept = if opts.dedupe {
                            dedupe_pairs(keyed, opts.dupesubs)
                        } else {
                            keyed.into_iter().map(|k| k.pair).collect()
                        };
                        let mut buf = Vec::new();
                        for (r1, r2) in &kept {
                            write_record(&mut buf, r1)?;
                            if let Some(r2) = r2 {
                                write_record(&mut buf, r2)?;
                            }
                        }
                        Ok(buf)
                    })();
                    (bi, bytes)
                })
                .collect();
            for (_, bytes) in results {
                out.write_all(&bytes?)?;
            }
            b = end;
        }
        Ok(())
    })();
    let _ = std::fs::remove_dir_all(&tmp);
    process_result
}

/// Creates a unique temporary directory for bucket files.
fn temp_dir_for() -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!("pgr-clump-{}-{nanos}", std::process::id()))
}

/// Conservative in-memory footprint estimate for FASTQ inputs: gzipped inputs
/// are expanded ~4x and records carry ~2x overhead.
fn estimate_fastq_bytes(infiles: &[String]) -> u64 {
    infiles
        .iter()
        .map(|f| {
            let bytes = std::fs::metadata(f).map(|m| m.len()).unwrap_or(0);
            if f.ends_with(".gz") {
                bytes * 8
            } else {
                bytes * 2
            }
        })
        .sum()
}

/// Streams interleaved pairs (1 file) or R1/R2 pairs (2 files).
fn for_each_pair(
    infiles: &[String],
    mut f: impl FnMut(&SeqRecord, Option<&SeqRecord>) -> Result<()>,
) -> Result<()> {
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
            f(&rec1, has2.then_some(&rec2))?;
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
            f(&rec1, has2.then_some(&rec2))?;
            if !has2 {
                break;
            }
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

/// A pair with its pivot key and expected-error estimate.
struct KeyedPair {
    pair: (SeqRecord, Option<SeqRecord>),
    key: ReadKey,
    errors: f32,
}

/// Sorts pairs by R1 pivot k-mer (clumpify-compatible order).
fn sort_keyed(pairs: Vec<(SeqRecord, Option<SeqRecord>)>, k: usize, seed: u64) -> Vec<KeyedPair> {
    let codes = make_codes(seed);
    let prob = prob_error();
    let mut keyed: Vec<KeyedPair> = pairs
        .into_iter()
        .map(|pair| {
            let key = fill_max(pair.0.sequence(), k, &codes);
            let errors = expected_errors(&pair.0, pair.1.as_ref(), &prob);
            KeyedPair { pair, key, errors }
        })
        .collect();
    keyed.par_sort_by(compare);
    keyed
}

/// Removes whole-pair duplicates within clumps (Clump.removeDuplicates with
/// `dupesubs` exact matching semantics).
fn dedupe_pairs(keyed: Vec<KeyedPair>, dupesubs: usize) -> Vec<(SeqRecord, Option<SeqRecord>)> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < keyed.len() {
        let kmer = keyed[i].key.kmer;
        let mut j = i + 1;
        while j < keyed.len() && keyed[j].key.kmer == kmer {
            j += 1;
        }
        let mut discarded = vec![false; j - i];
        dedupe_clump(&keyed[i..j], &mut discarded, dupesubs);
        for (a, d) in discarded.iter().enumerate() {
            if !d {
                out.push(keyed[i + a].pair.clone());
            }
        }
        i = j;
    }
    out
}

/// Removes duplicates within one clump (Clump.removeDuplicates with the
/// BBTools scan parameters: scan=0 for exact dupesubs=0, otherwise scan=5
/// with a wider retry when more than `maxDiscarded` reads are removed).
fn dedupe_clump(clump: &[KeyedPair], discarded: &mut [bool], dupesubs: usize) {
    let mut scan = if dupesubs < 1 { 0 } else { SCAN_LIMIT };
    let mut max_discarded = scan + 10;
    loop {
        let removed = dedupe_pass(clump, discarded, dupesubs, scan, max_discarded);
        if !(dupesubs > 0 && removed > max_discarded) {
            break;
        }
        scan += 10;
        max_discarded = max_discarded * 2 + 20;
    }
}

/// One dedupe scan pass over a clump (Clump.removeDuplicates_inner).
fn dedupe_pass(
    clump: &[KeyedPair],
    discarded: &mut [bool],
    dupesubs: usize,
    scan_limit: usize,
    max_discarded: usize,
) -> usize {
    let mut removed = 0usize;
    for a in 0..clump.len() {
        if discarded[a] {
            continue;
        }
        let mut unequals = 0usize;
        let mut discarded_seen = 0usize;
        let mut b = a + 1;
        while b < clump.len() && unequals <= scan_limit && discarded_seen <= max_discarded {
            if discarded[b] {
                discarded_seen += 1;
            } else {
                if !key_equal(&clump[a].key, &clump[b].key) {
                    break;
                }
                if pair_equal(&clump[a].pair, &clump[b].pair, dupesubs) {
                    if clump[b].errors >= clump[a].errors {
                        discarded[b] = true;
                        removed += 1;
                        unequals = 0;
                    } else {
                        discarded[a] = true;
                        removed += 1;
                        break;
                    }
                } else {
                    unequals += 1;
                }
            }
            b += 1;
        }
    }
    removed
}

/// ReadKey equality: kmer, strand, and window position must all match.
fn key_equal(a: &ReadKey, b: &ReadKey) -> bool {
    a.kmer == b.kmer && a.minus == b.minus && a.position == b.position
}

/// Whole-pair equality within `max_subs` substitutions (N is a wildcard).
fn pair_equal(
    a: &(SeqRecord, Option<SeqRecord>),
    b: &(SeqRecord, Option<SeqRecord>),
    max_subs: usize,
) -> bool {
    if !seq_equal(a.0.sequence(), b.0.sequence(), max_subs) {
        return false;
    }
    match (&a.1, &b.1) {
        (Some(a2), Some(b2)) => seq_equal(a2.sequence(), b2.sequence(), max_subs),
        _ => true,
    }
}

/// Byte equality within `max_subs` substitutions; N matches anything.
fn seq_equal(a: &[u8], b: &[u8], max_subs: usize) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut subs = 0usize;
    for (x, y) in a.iter().zip(b) {
        if x != y && *x != b'N' && *y != b'N' {
            subs += 1;
            if subs > max_subs {
                return false;
            }
        }
    }
    true
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
fn compare(a: &KeyedPair, b: &KeyedPair) -> Ordering {
    let order = compare_key(&a.key, &b.key);
    if order != Ordering::Equal {
        return order;
    }
    // BBTools KmerComparator.compareSequence defaults to true; clumpify
    // always orders identical keys by sequence, then expected errors.
    let order =
        compare_sequence_records(&a.pair.0, a.pair.1.as_ref(), &b.pair.0, b.pair.1.as_ref());
    if order != Ordering::Equal {
        return order;
    }
    if a.errors != b.errors {
        return if a.errors > b.errors {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }
    id(&a.pair.0).cmp(&id(&b.pair.0))
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
fn compare_sequence_records(
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

/// Sum of phred quality bytes (order-equivalent to BBTools quality sums).
fn quality_sum(rec: &SeqRecord) -> u32 {
    rec.quality_scores()
        .iter()
        .map(|&q| q.saturating_sub(QUALITY_BASE) as u32)
        .sum()
}

/// Expected error count of a pair (float accumulation like BBTools).
fn expected_errors(a: &SeqRecord, a2: Option<&SeqRecord>, prob: &[f32; 128]) -> f32 {
    let mut sum = 0f32;
    let mut add = |rec: &SeqRecord| {
        for (i, &b) in rec.sequence().iter().enumerate() {
            let q = if b == b'A' || b == b'C' || b == b'G' || b == b'T' {
                rec.quality_scores()
                    .get(i)
                    .copied()
                    .unwrap_or(0)
                    .saturating_sub(QUALITY_BASE) as usize
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
