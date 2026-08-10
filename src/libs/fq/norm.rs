//! BBNorm-style low-depth read filtering.
//!
//! Discards reads whose k-mer coverage is below a minimum depth, following
//! BBTools 39.38 `bbnorm.sh passes=1 bits=16 min=<n> target=9999999` read
//! decision logic but with an exact k-mer count table instead of the
//! approximate `bits=16` hash table.

use crate::libs::fmt::fq::write_fq;
use crate::libs::fmt::seq::SeqRecord;
use crate::libs::fq::clump::temp_dir_for;
use crate::libs::fq::pairs::PairReader;
use crate::libs::kmer::{count, KmerTable};
use anyhow::{Context, Result};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

/// Peak bytes per k-mer while counting one bucket (key + radix scratch +
/// count-vector capacity).
const COUNT_BYTES_PER_KMER: f64 = 21.0;
/// Serialized table bytes per k-mer (key + count).
const TABLE_BYTES_PER_KMER: f64 = 20.0;
/// Peak bytes per k-mer inside a scoring chunk (bucketed key + coverage slot
/// + wrapped record overhead).
const CHUNK_BYTES_PER_KMER: f64 = 24.0;
/// Fraction of the memory cap reserved for one bucket's count working set.
const BUCKET_FRACTION: f64 = 0.35;
/// Fraction of the memory cap reserved for one scoring chunk.
const CHUNK_FRACTION: f64 = 0.4;
/// Maximum bucket count (mirrors the clump external path).
const MAX_BUCKETS: usize = 4096;
/// Estimated peak memory per input base for the in-memory path (canonical
/// keys + sort scratch + count table + record overhead).
const MEM_BYTES_PER_BASE: u64 = 48;
/// Assumed uncompressed-to-gzip ratio when estimating bases from .gz files.
const GZ_EXPANSION: u64 = 8;

/// Options for the k-mer normalization cutoff.
#[derive(Debug, Clone)]
pub struct NormOptions {
    /// K-mer size (bbnorm `k`).
    pub k: usize,
    /// Minimum k-mer depth (bbnorm `min`).
    pub min_depth: usize,
    /// User memory cap in bytes (`--mem`); default 2 GiB.
    pub mem: Option<u64>,
}

/// Filtering results for one read.
struct ReadStats {
    true_depth: i64,
    depth_al: i64,
}

/// Filters reads by k-mer depth; writes survivors in input order.
///
/// Uses the exact count table in memory when the estimated footprint fits
/// `--mem`; larger inputs are counted through external hash buckets and
/// scored in bounded-memory chunks (same output, byte for byte).
pub fn norm<W: Write + Send>(
    infiles: &[String],
    out: &mut W,
    opts: &NormOptions,
    parallel: usize,
) -> Result<()> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(parallel.max(1))
        .build()
        .context("failed to build rayon pool")?;
    pool.install(|| {
        let cap = crate::libs::sys::mem_cap(opts.mem);
        let (mem_est, kmer_est) = estimate_norm(infiles);
        if mem_est <= cap {
            norm_in_memory(infiles, out, opts, parallel)
        } else {
            norm_buckets(infiles, out, opts, cap, kmer_est)
        }
    })
}

/// In-memory path: exact canonical counts, then parallel per-pair scoring.
fn norm_in_memory<W: Write>(
    infiles: &[String],
    out: &mut W,
    opts: &NormOptions,
    parallel: usize,
) -> Result<()> {
    let table = {
        let mut seqs: Vec<Vec<u8>> = Vec::new();
        for pair in PairReader::new(infiles)? {
            let (r1, r2) = pair?;
            seqs.push(r1.sequence().to_vec());
            if let Some(r2) = r2 {
                seqs.push(r2.sequence().to_vec());
            }
        }
        Arc::new(count::build_table(&seqs, opts.k)?)
    };
    let opts = Arc::new(opts.clone());
    // Pass 2: score every read against the table in parallel, writing
    // survivors in input order.
    let pairs = PairReader::new(infiles)?;
    crate::libs::par::ordered_map(
        pairs,
        parallel,
        move |pair| {
            let (r1, r2) = pair?;
            let s1 = read_stats(&r1, &table, &opts);
            let s2 = r2.as_ref().map(|r| read_stats(r, &table, &opts));
            Ok((!pair_tossed(&s1, s2.as_ref(), opts.min_depth)).then_some((r1, r2)))
        },
        |pair| {
            if let Some((r1, r2)) = pair {
                write_record(out, &r1)?;
                if let Some(r2) = r2 {
                    write_record(out, &r2)?;
                }
            }
            Ok(())
        },
    )
}

/// Whether a pair is tossed given both mates' coverage stats.
fn pair_tossed(s1: &ReadStats, s2: Option<&ReadStats>, min_depth: usize) -> bool {
    let min_al = match (s1, s2) {
        (s1, Some(s2)) => {
            if s1.depth_al >= 0 && s2.depth_al >= 0 {
                s1.depth_al.min(s2.depth_al)
            } else if s1.depth_al >= 0 {
                s1.depth_al
            } else {
                s2.depth_al
            }
        }
        (s1, None) => s1.depth_al,
    };
    let max_true = match (s1, s2) {
        (s1, Some(s2)) => s1.true_depth.max(s2.true_depth),
        (s1, None) => s1.true_depth,
    };
    min_al < 0 || max_true < min_depth as i64
}

/// Per-read coverage quantiles (KmerNormalize truedepth/depthAL).
fn read_stats(rec: &SeqRecord, table: &KmerTable, opts: &NormOptions) -> ReadStats {
    let mut cov: Vec<u32> = Vec::new();
    crate::libs::kmer::canonical_keys(rec.sequence(), opts.k, |_, key| {
        let idx = table.keys.partition_point(|&x| x < key);
        let c = if idx < table.keys.len() && table.keys[idx] == key {
            table.counts[idx]
        } else {
            0
        };
        cov.push(c);
    });
    score_cov(&cov, opts)
}

/// Coverage quantiles from a read's k-mer count list (order-independent).
fn score_cov(cov: &[u32], opts: &NormOptions) -> ReadStats {
    if cov.is_empty() {
        return ReadStats {
            true_depth: -1,
            depth_al: -1,
        };
    }
    let mut cov = cov.to_vec();
    cov.sort_unstable();
    let covlast = cov.len() - 1;
    let high = cov[((covlast as f64) * 0.10) as usize];
    let true_depth = cov[((covlast as f64) * 0.46) as usize] as i64;
    let mindepth = opts.min_depth.max((high / 125) as usize);
    let mut above_limit = covlast as i64;
    while above_limit >= 0 && cov[above_limit as usize] < mindepth as u32 {
        above_limit -= 1;
    }
    let mut depth_al = -1;
    let min_kmers = 15usize;
    if above_limit + 1 >= min_kmers as i64 || (above_limit >= 0 && min_kmers > cov.len()) {
        depth_al = cov[((above_limit as f64) * 0.46) as usize] as i64;
    }
    ReadStats {
        true_depth,
        depth_al,
    }
}

/// External path: bucket-count canonical k-mers on disk, then score reads in
/// memory-bounded chunks. Decisions and output order match the in-memory path.
fn norm_buckets<W: Write>(
    infiles: &[String],
    out: &mut W,
    opts: &NormOptions,
    cap: u64,
    est_kmers: u64,
) -> Result<()> {
    let buckets = bucket_count(cap, est_kmers);
    let tmp = temp_dir_for();
    std::fs::create_dir_all(&tmp)
        .with_context(|| format!("failed to create temp dir {}", tmp.display()))?;
    let result = (|| -> Result<()> {
        count_buckets(infiles, opts, buckets, cap, est_kmers, &tmp)?;
        score_in_chunks(infiles, out, opts, buckets, cap, &tmp)
    })();
    let _ = std::fs::remove_dir_all(&tmp);
    result
}

/// Bucket count from the memory cap and the input estimate: one bucket's
/// counting working set stays within `BUCKET_FRACTION` of the cap.
fn bucket_count(cap: u64, est_kmers: u64) -> usize {
    let per_bucket = ((cap as f64 * BUCKET_FRACTION) / COUNT_BYTES_PER_KMER).max(1.0) as u64;
    ((est_kmers.div_ceil(per_bucket)).max(2) as usize).min(MAX_BUCKETS)
}

/// Peak-memory and total-k-mer estimates for the norm workload.
///
/// The k-mer table (not the FASTQ bytes) dominates: ~1 key per base at
/// `MEM_BYTES_PER_BASE` bytes each. Gzip inputs assume `GZ_EXPANSION`x
/// expansion; plain files carry ~2.3 bytes per base.
fn estimate_norm(infiles: &[String]) -> (u64, u64) {
    let mut bases = 0u64;
    for f in infiles {
        let bytes = std::fs::metadata(f).map(|m| m.len()).unwrap_or(0);
        bases += if f.ends_with(".gz") {
            bytes * GZ_EXPANSION
        } else {
            bytes * 3 / 7
        };
    }
    (bases.saturating_mul(MEM_BYTES_PER_BASE), bases)
}

/// Passes A+B: stream canonical keys into per-bucket files, then count each
/// bucket in memory-bounded parallel waves and write sorted (key, count)
/// tables to disk.
fn count_buckets(
    infiles: &[String],
    opts: &NormOptions,
    buckets: usize,
    cap: u64,
    est_kmers: u64,
    tmp: &Path,
) -> Result<()> {
    let mut writers: Vec<Option<BufWriter<File>>> = (0..buckets).map(|_| None).collect();
    for pair in PairReader::new(infiles)? {
        let (r1, r2) = pair?;
        write_bucket_keys(r1.sequence(), opts.k, buckets, tmp, &mut writers)?;
        if let Some(r2) = r2 {
            write_bucket_keys(r2.sequence(), opts.k, buckets, tmp, &mut writers)?;
        }
    }
    for w in writers.iter_mut().flatten() {
        w.flush()?;
    }
    drop(writers);

    // A wave's concurrent sort working sets plus their serialized tables stay
    // within the memory cap; results are collected per wave and written after.
    let bucket_est = (est_kmers.div_ceil(buckets as u64)).max(1);
    let wave = ((cap as f64 * 0.8)
        / (bucket_est as f64 * (COUNT_BYTES_PER_KMER + TABLE_BYTES_PER_KMER)))
        .floor()
        .max(1.0) as usize;
    let mut b = 0usize;
    while b < buckets {
        let end = (b + wave).min(buckets);
        let results: Vec<(usize, Result<Vec<u8>>)> = (b..end)
            .into_par_iter()
            .map(|bi| {
                let path = tmp.join(format!("bucket_{bi:05}.kmer"));
                let bytes = (|| -> Result<Vec<u8>> {
                    if !path.exists() {
                        return Ok(Vec::new());
                    }
                    let keys = read_keys(&path)?;
                    let table = count::count_keys(keys, opts.k);
                    serialize_table(&table)
                })();
                (bi, bytes)
            })
            .collect();
        for (bi, bytes) in results {
            std::fs::write(tmp.join(format!("table_{bi:05}.tbl")), bytes?)?;
        }
        b = end;
    }
    Ok(())
}

/// Pass C: score reads in chunks whose bucketed keys + coverage slots fit the
/// budget, loading one count table at a time, and emit survivors in order.
fn score_in_chunks<W: Write>(
    infiles: &[String],
    out: &mut W,
    opts: &NormOptions,
    buckets: usize,
    cap: u64,
    tmp: &Path,
) -> Result<()> {
    let chunk_budget = ((cap as f64 * CHUNK_FRACTION) / CHUNK_BYTES_PER_KMER).max(1.0) as usize;
    let mut reader = PairReader::new(infiles)?;
    loop {
        let mut chunk: Vec<(SeqRecord, Option<SeqRecord>)> = Vec::new();
        let mut kmers = 0usize;
        for pair in reader.by_ref() {
            let (r1, r2) = pair?;
            kmers += n_kmers(r1.sequence(), opts.k);
            if let Some(r2) = r2.as_ref() {
                kmers += n_kmers(r2.sequence(), opts.k);
            }
            chunk.push((r1, r2));
            if kmers >= chunk_budget {
                break;
            }
        }
        if chunk.is_empty() {
            break;
        }
        score_chunk(&chunk, opts, buckets, tmp, out)?;
    }
    Ok(())
}

/// Distributes a chunk's k-mers to per-bucket lists, fills each read's
/// coverage slots from the counted tables, and writes surviving pairs.
fn score_chunk<W: Write>(
    chunk: &[(SeqRecord, Option<SeqRecord>)],
    opts: &NormOptions,
    buckets: usize,
    tmp: &Path,
    out: &mut W,
) -> Result<()> {
    let mut bucketed: Vec<Vec<(u32, u64)>> = vec![Vec::new(); buckets];
    for (i, (r1, r2)) in chunk.iter().enumerate() {
        let base = 2 * i as u32;
        collect_bucket_keys(r1.sequence(), opts.k, buckets, base, &mut bucketed);
        if let Some(r2) = r2 {
            collect_bucket_keys(r2.sequence(), opts.k, buckets, base + 1, &mut bucketed);
        }
    }
    let mut covs: Vec<Vec<u32>> = vec![Vec::new(); 2 * chunk.len()];
    for (b, keys) in bucketed.iter().enumerate() {
        if keys.is_empty() {
            continue;
        }
        let table = load_table(&tmp.join(format!("table_{b:05}.tbl")))?;
        let hits: Vec<(u32, u32)> = keys
            .par_iter()
            .map(|&(ri, key)| {
                let key = key as u128;
                let idx = table.keys.partition_point(|&x| x < key);
                let c = if idx < table.keys.len() && table.keys[idx] == key {
                    table.counts[idx]
                } else {
                    0
                };
                (ri, c)
            })
            .collect();
        for (ri, c) in hits {
            covs[ri as usize].push(c);
        }
    }
    for (i, (r1, r2)) in chunk.iter().enumerate() {
        let s1 = score_cov(&covs[2 * i], opts);
        let s2 = r2.as_ref().map(|_| score_cov(&covs[2 * i + 1], opts));
        if !pair_tossed(&s1, s2.as_ref(), opts.min_depth) {
            write_record(out, r1)?;
            if let Some(r2) = r2 {
                write_record(out, r2)?;
            }
        }
    }
    Ok(())
}

/// Writes every canonical k-mer of `seq` to its bucket file.
fn write_bucket_keys(
    seq: &[u8],
    k: usize,
    buckets: usize,
    tmp: &Path,
    writers: &mut [Option<BufWriter<File>>],
) -> Result<()> {
    let mut keys = Vec::new();
    crate::libs::kmer::canonical_keys(seq, k, |_, key| keys.push(key));
    for key in keys {
        let b = bucket_of(key, buckets);
        if writers[b].is_none() {
            let path = tmp.join(format!("bucket_{b:05}.kmer"));
            let f = File::create(&path)
                .with_context(|| format!("failed to create {}", path.display()))?;
            writers[b] = Some(BufWriter::new(f));
        }
        // k <= 31 in norm, so 2*k <= 62 bits: the canonical key fits u64.
        writers[b]
            .as_mut()
            .unwrap()
            .write_all(&(key as u64).to_le_bytes())?;
    }
    Ok(())
}

/// Adds `(read index, key)` for every canonical k-mer of `seq`.
fn collect_bucket_keys(
    seq: &[u8],
    k: usize,
    buckets: usize,
    read_idx: u32,
    bucketed: &mut [Vec<(u32, u64)>],
) {
    crate::libs::kmer::canonical_keys(seq, k, |_, key| {
        // k <= 31 in norm, so 2*k <= 62 bits: the canonical key fits u64.
        bucketed[bucket_of(key, buckets)].push((read_idx, key as u64));
    });
}

/// Deterministic bucket index for a canonical k-mer.
fn bucket_of(key: u128, buckets: usize) -> usize {
    (key % buckets as u128) as usize
}

/// Number of N-free k-mer windows in `seq` (what `canonical_keys` emits).
fn n_kmers(seq: &[u8], k: usize) -> usize {
    let mut total = 0usize;
    let mut run = 0usize;
    for &b in seq {
        if matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            run += 1;
        } else {
            if run >= k {
                total += run - k + 1;
            }
            run = 0;
        }
    }
    if run >= k {
        total += run - k + 1;
    }
    total
}

/// Reads 8-byte little-endian keys written by `write_bucket_keys`.
fn read_keys(path: &Path) -> Result<Vec<u128>> {
    let bytes = std::fs::read(path)?;
    anyhow::ensure!(
        bytes.len() % 8 == 0,
        "corrupt bucket file {} ({} bytes)",
        path.display(),
        bytes.len()
    );
    let mut keys = Vec::with_capacity(bytes.len() / 8);
    for chunk in bytes.chunks_exact(8) {
        keys.push(u64::from_le_bytes(chunk.try_into().unwrap()) as u128);
    }
    Ok(keys)
}

/// Serializes a count table as interleaved (key u64 LE, count u32 LE).
fn serialize_table(table: &KmerTable) -> Result<Vec<u8>> {
    let mut buf = Vec::with_capacity(table.keys.len() * 12);
    for (&key, &c) in table.keys.iter().zip(&table.counts) {
        buf.extend_from_slice(&(key as u64).to_le_bytes());
        buf.extend_from_slice(&c.to_le_bytes());
    }
    Ok(buf)
}

/// Loads a table written by `serialize_table` (empty for an empty file).
fn load_table(path: &Path) -> Result<KmerTable> {
    let bytes = std::fs::read(path)?;
    anyhow::ensure!(
        bytes.len() % 12 == 0,
        "corrupt table file {} ({} bytes)",
        path.display(),
        bytes.len()
    );
    let mut keys = Vec::with_capacity(bytes.len() / 12);
    let mut counts = Vec::with_capacity(bytes.len() / 12);
    for chunk in bytes.chunks_exact(12) {
        keys.push(u64::from_le_bytes(chunk[..8].try_into().unwrap()) as u128);
        counts.push(u32::from_le_bytes(chunk[8..].try_into().unwrap()));
    }
    Ok(KmerTable { k: 0, keys, counts })
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

#[cfg(test)]
mod tests {
    use super::*;

    fn mixed_input() -> String {
        let mut input = String::new();
        for i in 0..100 {
            input.push_str(&format!(
                "@hi{i}\n{}\n+\n{}\n",
                "ACGT".repeat(20),
                "I".repeat(80)
            ));
        }
        for i in 0..20 {
            input.push_str(&format!(
                "@mid{i}\n{}\n+\n{}\n",
                "GATTACA".repeat(12),
                "I".repeat(84)
            ));
        }
        input.push_str(&format!(
            "@lo1\n{}\n+\n{}\n",
            "GATCCTAGACGTTCGATCGGTACCTAGCATGCAGTTACGTACGATCGTAGCTAGCGGATCGATC",
            "I".repeat(64)
        ));
        input.push_str(&format!(
            "@lo2\n{}\n+\n{}\n",
            "TTTTTCCCCGGGAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA",
            "I".repeat(64)
        ));
        input
    }

    #[test]
    fn external_path_matches_in_memory_path() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("in.fq");
        std::fs::write(&path, mixed_input()).unwrap();
        let infiles = vec![path.to_str().unwrap().to_string()];

        let mut mem_out = Vec::new();
        norm(
            &infiles,
            &mut mem_out,
            &NormOptions {
                k: 31,
                min_depth: 3,
                mem: None,
            },
            4,
        )
        .unwrap();
        // A 1 KiB cap forces the external bucket path on this input.
        let mut ext_out = Vec::new();
        norm(
            &infiles,
            &mut ext_out,
            &NormOptions {
                k: 31,
                min_depth: 3,
                mem: Some(1 << 10),
            },
            4,
        )
        .unwrap();
        assert_eq!(ext_out, mem_out);
        assert!(mem_out.windows(4).any(|w| w == b"@hi9"));
        assert!(!mem_out.windows(4).any(|w| w == b"@lo1"));
    }
}
