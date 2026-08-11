//! Perfect-match read mapping (bbmap `perfectmode` replacement).
//!
//! Builds a canonical k-mer position index over a reference and maps every
//! read by verifying full-length exact matches (no mismatches, no gaps) at
//! all candidate positions. Design: `notes/design/asm-map.md`.

use crate::libs::ds::radix_sort::radix_sort_u128;
use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use crate::libs::kmer::canonical_keys;
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
    keys: Vec<u128>,
    payloads: Vec<u64>,
}

/// Options for [`map_files`].
pub struct MapOptions {
    pub k: usize,
    pub outm: Option<String>,
    pub outu: Option<String>,
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

/// Builds the canonical k-mer index (`k` in 1..=64, u128 key).
pub fn build_index(refs: &[RefRecord], k: usize) -> Result<MapIndex> {
    anyhow::ensure!(
        (1..=64).contains(&k),
        "k-mer length must be in 1..=64, got {k}"
    );
    let mut keys = Vec::new();
    let mut payloads = Vec::new();
    for (cid, r) in refs.iter().enumerate() {
        let cid = cid as u64;
        canonical_keys(&r.seq, k, |pos, key| {
            keys.push(key);
            payloads.push((cid << 32) | pos as u64);
        });
    }
    radix_sort_u128(&mut keys, &mut payloads, 2 * k as u32);
    Ok(MapIndex { keys, payloads })
}

/// Maps all reads from `read_paths` and writes the SAM outputs.
///
/// Reads are processed in parallel (rayon) and written in input order, so
/// the output is deterministic across runs.
pub fn map_files(refs: &[RefRecord], read_paths: &[String], opts: &MapOptions) -> Result<MapStats> {
    let index = build_index(refs, opts.k)?;

    let mut stats = MapStats::default();
    let mut outm = opts.outm.as_ref().map(|p| pgr_writer(p)).transpose()?;
    let mut outu = opts.outu.as_ref().map(|p| pgr_writer(p)).transpose()?;
    // Both SAM outputs carry the same header (bbmap writes outm/outu with
    // the shared header), so line-count ratios stay symmetric.
    for w in [outm.as_mut(), outu.as_mut()].into_iter().flatten() {
        write_sam_header(w, refs)?;
    }

    // Stream reads in blocks (bounded memory) and process each block in
    // parallel; results are written in input order.
    let block_size = 100_000usize;
    let mut block: Vec<SeqRecord> = Vec::with_capacity(block_size);
    for path in read_paths {
        let mut reader = SeqReader::new(path)?;
        let mut rec = SeqRecord::new();
        while reader.read_record(&mut rec)? {
            block.push(rec.clone());
            if block.len() >= block_size {
                write_block(
                    &block, opts.k, refs, &index, &mut outm, &mut outu, &mut stats,
                )?;
                block.clear();
            }
        }
    }
    if !block.is_empty() {
        write_block(
            &block, opts.k, refs, &index, &mut outm, &mut outu, &mut stats,
        )?;
    }
    drop(outm);
    drop(outu);
    Ok(stats)
}

/// Maps one block of reads (parallel) and writes the SAM records in input
/// order.
fn write_block(
    block: &[SeqRecord],
    k: usize,
    refs: &[RefRecord],
    index: &MapIndex,
    outm: &mut Option<Box<dyn Write + Send>>,
    outu: &mut Option<Box<dyn Write + Send>>,
    stats: &mut MapStats,
) -> Result<()> {
    let results: Vec<ReadResult> = block
        .par_iter()
        .map(|rec| {
            let read = rec.sequence();
            let hits = map_read(read, k, index, refs);
            ReadResult {
                name: String::from_utf8_lossy(rec.name()).into_owned(),
                seq: read.to_vec(),
                qual: rec.quality_scores().to_vec(),
                hits,
            }
        })
        .collect();
    stats.reads_in += results.len() as u64;
    for r in &results {
        if r.hits.is_empty() {
            stats.unmapped += 1;
            if let Some(w) = outu.as_mut() {
                write_unmapped(w, r)?;
            }
        } else {
            stats.mapped += 1;
            stats.hits += r.hits.len() as u64;
            if let Some(w) = outm.as_mut() {
                for &(cid, pos, rc) in &r.hits {
                    write_mapped(w, r, &refs[cid as usize].name, pos, rc)?;
                }
            }
        }
    }
    Ok(())
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
    let Some(seed) = kmer_canonical(&read[..k], k) else {
        return Vec::new();
    };
    let lo = index.keys.partition_point(|&x| x < seed);
    let hi = index.keys.partition_point(|&x| x <= seed);
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

/// Canonical (2-bit, strand-minimized) key of a k-mer window; `None` when
/// the window contains an undefined base (N or ambiguity).
fn kmer_canonical(seq: &[u8], k: usize) -> Option<u128> {
    let codes = crate::libs::kmer::base_codes();
    let kmask = (1u128 << (2 * k)) - 1;
    let rc_top = (2 * k - 2) as u32;
    let mut fwd = 0u128;
    let mut rev = 0u128;
    for &b in seq {
        let c = codes[b as usize] as u128;
        if c == 4 {
            return None;
        }
        fwd = ((fwd << 2) | c) & kmask;
        rev = (rev >> 2) | ((3 - c) << rc_top);
    }
    Some(fwd.min(rev))
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

fn pgr_writer(path: &str) -> Result<Box<dyn Write + Send>> {
    crate::libs::io::writer(path).map(|w| Box::new(w) as Box<dyn Write + Send>)
}
