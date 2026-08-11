//! Perfect-match read mapping (bbmap `perfectmode` replacement).
//!
//! Builds a canonical k-mer position index over a reference and maps every
//! read by verifying full-length exact matches (no mismatches, no gaps) at
//! all candidate positions. Design: `notes/design/asm-map.md`.

use crate::libs::ds::radix_sort::radix_sort_bytes;
use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use crate::libs::kmer::canonical_keys;
use crate::libs::kmer::key::Kmer;
use crate::libs::nt::rev_comp;
use anyhow::Result;
use rayon::prelude::*;
use std::io::Write;

/// One reference contig (name + sequence).
pub struct RefRecord {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Canonical k-mer position index over the reference: sorted keys with
/// packed `(contig_id, pos)` payloads, binary-searched per seed.
pub struct MapIndex {
    /// Packed canonical k-mer keys (FastK bytes), `key_bytes` each.
    keys: Vec<u8>,
    payloads: Vec<u64>,
}

/// Options for [`map_files`].
pub struct MapOptions {
    pub k: usize,
    pub outm: Option<String>,
    pub outu: Option<String>,
    /// Interleave the two read files as R1/R2 pairs and emit paired SAM
    /// (FLAG 0x1/0x2/0x40/0x80, RNEXT/PNEXT, TLEN) for mapped pairs.
    pub paired: bool,
    /// Stop after processing this many read records (pairs count as two).
    pub max_reads: Option<u64>,
}

/// The two SAM outputs (mapped / unmapped).
struct SamOutputs {
    outm: Option<Box<dyn Write + Send>>,
    outu: Option<Box<dyn Write + Send>>,
}

/// Mapping statistics.
#[derive(Debug, Default, Clone)]
pub struct MapStats {
    pub reads_in: u64,
    pub mapped: u64,
    pub unmapped: u64,
    pub hits: u64,
}

/// Reads FASTA/FASTQ files (plain or gzipped) into reference records.
pub fn read_fasta(paths: &[String]) -> Result<Vec<RefRecord>> {
    let mut refs = Vec::new();
    for path in paths {
        let mut reader = SeqReader::new(path)?;
        let mut rec = SeqRecord::new();
        while reader.read_record(&mut rec)? {
            refs.push(RefRecord {
                name: String::from_utf8_lossy(rec.name()).into_owned(),
                seq: rec.sequence().to_vec(),
            });
        }
    }
    Ok(refs)
}

/// Builds the canonical k-mer index (`k` up to `Kmer::MAX_K`).
pub fn build_index(refs: &[RefRecord], k: usize) -> Result<MapIndex> {
    anyhow::ensure!(
        k > 0 && k <= Kmer::MAX_K,
        "k-mer length must be in 1..={}, got {k}",
        Kmer::MAX_K
    );
    let key_bytes = k.div_ceil(4);
    let mut keys = Vec::new();
    let mut payloads = Vec::new();
    for (cid, r) in refs.iter().enumerate() {
        let cid = cid as u64;
        canonical_keys(&r.seq, k, |pos, key| {
            keys.extend_from_slice(key.to_bytes());
            payloads.push((cid << 32) | pos as u64);
        });
    }
    radix_sort_bytes(&mut keys, key_bytes, &mut payloads);
    Ok(MapIndex { keys, payloads })
}

/// First index whose packed key is `>= key`.
fn lower_bound(keys: &[u8], key_bytes: usize, key: &Kmer) -> usize {
    let mut lo = 0usize;
    let mut hi = keys.len() / key_bytes;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if &keys[mid * key_bytes..(mid + 1) * key_bytes] < key.to_bytes() {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

/// First index whose packed key is `> key`.
fn upper_bound(keys: &[u8], key_bytes: usize, key: &Kmer) -> usize {
    let mut lo = 0usize;
    let mut hi = keys.len() / key_bytes;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if &keys[mid * key_bytes..(mid + 1) * key_bytes] <= key.to_bytes() {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

/// Maps all reads from `read_paths` and writes the SAM outputs.
///
/// Reads are processed in parallel (rayon) and written in input order, so
/// the output is deterministic across runs.
pub fn map_files(refs: &[RefRecord], read_paths: &[String], opts: &MapOptions) -> Result<MapStats> {
    let index = build_index(refs, opts.k)?;

    let mut stats = MapStats::default();
    let mut outputs = SamOutputs {
        outm: opts.outm.as_ref().map(|p| pgr_writer(p)).transpose()?,
        outu: opts.outu.as_ref().map(|p| pgr_writer(p)).transpose()?,
    };
    // Both SAM outputs carry the same header (bbmap writes outm/outu with
    // the shared header), so line-count ratios stay symmetric.
    for w in [outputs.outm.as_mut(), outputs.outu.as_mut()]
        .into_iter()
        .flatten()
    {
        write_sam_header(w, refs)?;
    }

    let block_size = 100_000usize;
    if opts.paired {
        map_paired(
            refs,
            read_paths,
            opts,
            &index,
            &mut outputs,
            &mut stats,
            block_size,
        )?;
    } else {
        map_single(
            refs,
            read_paths,
            opts,
            &index,
            &mut outputs,
            &mut stats,
            block_size,
        )?;
    }
    drop(outputs);
    Ok(stats)
}

/// Single-end path: stream reads in blocks (bounded memory), map each block
/// in parallel, and write results in input order.
fn map_single(
    refs: &[RefRecord],
    read_paths: &[String],
    opts: &MapOptions,
    index: &MapIndex,
    outputs: &mut SamOutputs,
    stats: &mut MapStats,
    block_size: usize,
) -> Result<()> {
    let mut block: Vec<SeqRecord> = Vec::with_capacity(block_size);
    let mut processed = 0u64;
    for path in read_paths {
        let mut reader = SeqReader::new(path)?;
        let mut rec = SeqRecord::new();
        while reader.read_record(&mut rec)? {
            if opts.max_reads.is_some_and(|m| processed >= m) {
                break;
            }
            processed += 1;
            block.push(rec.clone());
            if block.len() >= block_size {
                write_block(&block, opts.k, refs, index, outputs, stats)?;
                block.clear();
            }
        }
    }
    if !block.is_empty() {
        write_block(&block, opts.k, refs, index, outputs, stats)?;
    }
    Ok(())
}

/// Paired path: interleave the two read files (same length) into R1/R2
/// pairs, map both ends, and write paired SAM records in input order.
fn map_paired(
    refs: &[RefRecord],
    read_paths: &[String],
    opts: &MapOptions,
    index: &MapIndex,
    outputs: &mut SamOutputs,
    stats: &mut MapStats,
    block_size: usize,
) -> Result<()> {
    anyhow::ensure!(
        read_paths.len() == 2,
        "paired mode requires exactly 2 read files (R1, R2), got {}",
        read_paths.len()
    );
    let mut r1_reader = SeqReader::new(&read_paths[0])?;
    let mut r2_reader = SeqReader::new(&read_paths[1])?;
    let mut block: Vec<(SeqRecord, SeqRecord)> = Vec::with_capacity(block_size / 2);
    let mut processed = 0u64;
    loop {
        if opts.max_reads.is_some_and(|m| processed >= m) {
            break;
        }
        let mut a = SeqRecord::new();
        let mut b = SeqRecord::new();
        let a_ok = r1_reader.read_record(&mut a)?;
        let b_ok = r2_reader.read_record(&mut b)?;
        if !a_ok && !b_ok {
            break;
        }
        anyhow::ensure!(
            a_ok == b_ok,
            "paired files have different read counts: {} vs {}",
            read_paths[0],
            read_paths[1]
        );
        processed += 2;
        block.push((a, b));
        if block.len() >= block_size / 2 {
            write_pair_block(&block, opts.k, refs, index, outputs, stats)?;
            block.clear();
        }
    }
    if !block.is_empty() {
        write_pair_block(&block, opts.k, refs, index, outputs, stats)?;
    }
    Ok(())
}

/// Maps one block of reads (parallel) and writes the SAM records in input
/// order.
fn write_block(
    block: &[SeqRecord],
    k: usize,
    refs: &[RefRecord],
    index: &MapIndex,
    outputs: &mut SamOutputs,
    stats: &mut MapStats,
) -> Result<()> {
    let results: Vec<ReadResult> = block
        .par_iter()
        .map(|rec| map_one(rec, k, index, refs))
        .collect();
    stats.reads_in += results.len() as u64;
    for r in &results {
        if r.hits.is_empty() {
            stats.unmapped += 1;
            if let Some(w) = outputs.outu.as_mut() {
                write_unmapped(w, r)?;
            }
        } else {
            stats.mapped += 1;
            stats.hits += r.hits.len() as u64;
            if let Some(w) = outputs.outm.as_mut() {
                for &(cid, pos, rc) in &r.hits {
                    write_mapped(w, r, &refs[cid as usize].name, pos, rc)?;
                }
            }
        }
    }
    Ok(())
}

/// Maps one block of R1/R2 pairs (parallel) and writes paired SAM records in
/// input order.
fn write_pair_block(
    block: &[(SeqRecord, SeqRecord)],
    k: usize,
    refs: &[RefRecord],
    index: &MapIndex,
    outputs: &mut SamOutputs,
    stats: &mut MapStats,
) -> Result<()> {
    let results: Vec<(ReadResult, ReadResult)> = block
        .par_iter()
        .map(|(a, b)| (map_one(a, k, index, refs), map_one(b, k, index, refs)))
        .collect();
    stats.reads_in += (results.len() * 2) as u64;
    for (r1, r2) in &results {
        match (r1.hits.first(), r2.hits.first()) {
            (Some(&h1), Some(&h2)) => {
                stats.mapped += 2;
                stats.hits += 2;
                if let Some(w) = outputs.outm.as_mut() {
                    write_pair_mapped(w, r1, h1, r2, h2, refs)?;
                }
            }
            _ => {
                stats.unmapped += 2;
                if let Some(w) = outputs.outu.as_mut() {
                    write_pair_unmapped(w, r1, r2)?;
                }
            }
        }
    }
    Ok(())
}

/// Maps one read into a [`ReadResult`] (name, sequence, qualities, hits).
fn map_one(rec: &SeqRecord, k: usize, index: &MapIndex, refs: &[RefRecord]) -> ReadResult {
    let read = rec.sequence();
    let hits = map_read(read, k, index, refs);
    ReadResult {
        name: String::from_utf8_lossy(rec.name()).into_owned(),
        seq: read.to_vec(),
        qual: rec.quality_scores().to_vec(),
        hits,
    }
}

/// One read's mapping outcome (all exact-match positions).
struct ReadResult {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    hits: Vec<(u32, u32, bool)>,
}

/// Exact-match mapping of one read: seed on its first k-mer, then verify
/// the full length (forward or reverse) at every candidate position.
fn map_read(read: &[u8], k: usize, index: &MapIndex, refs: &[RefRecord]) -> Vec<(u32, u32, bool)> {
    if read.len() < k {
        return Vec::new();
    }
    let Some(seed) = Kmer::from_bases(&read[..k], k) else {
        return Vec::new();
    };
    let seed = seed.canonical();
    let key_bytes = k.div_ceil(4);
    let lo = lower_bound(&index.keys, key_bytes, &seed);
    let hi = upper_bound(&index.keys, key_bytes, &seed);
    if lo == hi {
        return Vec::new();
    }
    let rc: Vec<u8> = rev_comp(read).collect();
    let mut hits = Vec::new();
    for &p in &index.payloads[lo..hi] {
        let cid = (p >> 32) as usize;
        let q = (p & 0xffff_ffff) as usize;
        let seq = &refs[cid].seq;
        let l = read.len();
        // Forward: the seed k-mer sits at the read start.
        if q + l <= seq.len() && &seq[q..q + l] == read {
            hits.push((cid as u32, q as u32, false));
        }
        // Reverse: the seed k-mer is the rc of the reference window ending
        // at q, i.e. the read starts at q - (l - k).
        else if q >= l - k && q + k <= seq.len() && seq[q - (l - k)..q + k] == rc[..] {
            hits.push((cid as u32, (q - (l - k)) as u32, true));
        }
    }
    // Deterministic output order: by reference position (radix sort does not
    // promise a stable order among equal keys).
    hits.sort_unstable();
    hits
}

fn write_sam_header<W: Write>(w: &mut W, refs: &[RefRecord]) -> Result<()> {
    writeln!(w, "@HD\tVN:1.6\tSO:unknown")?;
    for r in refs {
        writeln!(w, "@SQ\tSN:{}\tLN:{}", r.name, r.seq.len())?;
    }
    Ok(())
}

fn qual_string(q: &[u8]) -> String {
    if q.is_empty() {
        "*".to_string()
    } else {
        String::from_utf8_lossy(q).into_owned()
    }
}

fn write_mapped<W: Write>(
    w: &mut W,
    r: &ReadResult,
    rname: &str,
    pos: u32,
    rc: bool,
) -> Result<()> {
    let flag = if rc { 16 } else { 0 };
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t255\t{}M\t*\t0\t0\t{}\t{}",
        r.name,
        flag,
        rname,
        pos + 1,
        r.seq.len(),
        String::from_utf8_lossy(&r.seq),
        qual_string(&r.qual),
    )?;
    Ok(())
}

fn write_unmapped<W: Write>(w: &mut W, r: &ReadResult) -> Result<()> {
    writeln!(
        w,
        "{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}",
        r.name,
        String::from_utf8_lossy(&r.seq),
        qual_string(&r.qual),
    )?;
    Ok(())
}

/// A mapped R1/R2 pair: FLAG 0x1 + 0x40/0x80 (+ 0x2 when the pair is a
/// proper FR pair), RNEXT/PNEXT mate coordinates, and signed TLEN.
fn write_pair_mapped<W: Write>(
    w: &mut W,
    r1: &ReadResult,
    h1: (u32, u32, bool),
    r2: &ReadResult,
    h2: (u32, u32, bool),
    refs: &[RefRecord],
) -> Result<()> {
    let (c1, p1, rc1) = h1;
    let (c2, p2, rc2) = h2;
    let (flag1, flag2, tlen1, tlen2) = match proper_pair_insert(h1, h2, r1.seq.len(), r2.seq.len())
    {
        Some((insert, left_is_r1)) => {
            let t1 = if left_is_r1 {
                insert as i64
            } else {
                -(insert as i64)
            };
            (0x2, 0x2, t1, -t1)
        }
        None => (0, 0, 0, 0),
    };
    let flag1 = 0x1 | 0x40 | flag1 | (rc1 as u16 * 0x10) | (rc2 as u16 * 0x20);
    let flag2 = 0x1 | 0x80 | flag2 | (rc2 as u16 * 0x10) | (rc1 as u16 * 0x20);
    let same = c1 == c2;
    let rnext1: &str = if same { "=" } else { &refs[c2 as usize].name };
    let rnext2: &str = if same { "=" } else { &refs[c1 as usize].name };
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t255\t{}M\t{}\t{}\t{}\t{}\t{}",
        r1.name,
        flag1,
        refs[c1 as usize].name,
        p1 + 1,
        r1.seq.len(),
        rnext1,
        p2 + 1,
        tlen1,
        String::from_utf8_lossy(&r1.seq),
        qual_string(&r1.qual),
    )?;
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t255\t{}M\t{}\t{}\t{}\t{}\t{}",
        r2.name,
        flag2,
        refs[c2 as usize].name,
        p2 + 1,
        r2.seq.len(),
        rnext2,
        p1 + 1,
        tlen2,
        String::from_utf8_lossy(&r2.seq),
        qual_string(&r2.qual),
    )?;
    Ok(())
}

/// An unmapped pair: both records FLAG 0x4 (segmented, mate unmapped).
fn write_pair_unmapped<W: Write>(w: &mut W, r1: &ReadResult, r2: &ReadResult) -> Result<()> {
    for (flag, r) in [(0x1 | 0x4 | 0x8 | 0x40, r1), (0x1 | 0x4 | 0x8 | 0x80, r2)] {
        writeln!(
            w,
            "{}\t{}\t*\t0\t0\t*\t*\t0\t0\t{}\t{}",
            r.name,
            flag,
            String::from_utf8_lossy(&r.seq),
            qual_string(&r.qual),
        )?;
    }
    Ok(())
}

/// Insert size and leftmost read of a proper FR pair (`None` otherwise):
/// same contig, opposite strands, reads pointing inward.
fn proper_pair_insert(
    h1: (u32, u32, bool),
    h2: (u32, u32, bool),
    len1: usize,
    len2: usize,
) -> Option<(usize, bool)> {
    let (c1, p1, rc1) = h1;
    let (c2, p2, rc2) = h2;
    if c1 != c2 || rc1 == rc2 {
        return None;
    }
    let (left, right, left_is_r1, right_len) = if !rc1 && rc2 && p1 < p2 {
        (p1, p2, true, len2)
    } else if rc1 && !rc2 && p2 < p1 {
        (p2, p1, false, len1)
    } else {
        return None; // outward orientation
    };
    Some((right as usize + right_len - left as usize, left_is_r1))
}

fn pgr_writer(path: &str) -> Result<Box<dyn Write + Send>> {
    crate::libs::io::writer(path).map(|w| Box::new(w) as Box<dyn Write + Send>)
}
