//! BBNorm-style low-depth read filtering.
//!
//! Discards reads whose k-mer coverage is below a minimum depth, following
//! BBTools 39.38 `bbnorm.sh passes=1 bits=16 min=<n> target=9999999` read
//! decision logic but with an exact k-mer count table instead of the
//! approximate `bits=16` hash table.

use crate::libs::fmt::fq::write_fq;
use crate::libs::fmt::seq::SeqRecord;
use crate::libs::kmer::{count, KmerTable};
use anyhow::Result;
use std::collections::HashMap;
use std::io::Write;

use super::clump::read_pairs;

/// Options for the k-mer normalization cutoff.
#[derive(Debug, Clone)]
pub struct NormOptions {
    /// K-mer size (bbnorm `k`).
    pub k: usize,
    /// Minimum k-mer depth (bbnorm `min`).
    pub min_depth: usize,
}

/// Filtering results for one read.
struct ReadStats {
    true_depth: i64,
    depth_al: i64,
}

/// Filters reads by k-mer depth; writes survivors in input order.
pub fn norm<W: Write>(infiles: &[String], out: &mut W, opts: &NormOptions) -> Result<()> {
    let pairs = read_pairs(infiles)?;
    // Pass 1: exact canonical k-mer counts.
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    for (r1, r2) in &pairs {
        seqs.push(r1.sequence().to_vec());
        if let Some(r2) = r2 {
            seqs.push(r2.sequence().to_vec());
        }
    }
    let table = count::build_table(&seqs, opts.k)?;
    let counts: HashMap<u128, u32> = table
        .keys
        .iter()
        .zip(&table.counts)
        .map(|(&k, &c)| (k, c))
        .collect();

    for (r1, r2) in &pairs {
        let s1 = read_stats(r1, &table, &counts, opts);
        let s2 = r2.as_ref().map(|r| read_stats(r, &table, &counts, opts));
        let min_al = match (&s1, &s2) {
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
        let max_true = match (&s1, &s2) {
            (s1, Some(s2)) => s1.true_depth.max(s2.true_depth),
            (s1, None) => s1.true_depth,
        };
        let toss = min_al < 0 || max_true < opts.min_depth as i64;
        if !toss {
            write_record(out, r1)?;
            if let Some(r2) = r2 {
                write_record(out, r2)?;
            }
        }
    }
    Ok(())
}

/// Per-read coverage quantiles (KmerNormalize truedepth/depthAL).
fn read_stats(
    rec: &SeqRecord,
    table: &KmerTable,
    counts: &HashMap<u128, u32>,
    opts: &NormOptions,
) -> ReadStats {
    let seq = rec.sequence();
    if seq.len() < opts.k {
        return ReadStats {
            true_depth: -1,
            depth_al: -1,
        };
    }
    let mut cov: Vec<u32> = Vec::new();
    crate::libs::kmer::canonical_keys(seq, opts.k, |_, key| {
        cov.push(counts.get(&key).copied().unwrap_or(0));
    });
    if cov.is_empty() {
        return ReadStats {
            true_depth: -1,
            depth_al: -1,
        };
    }
    let _ = table;
    cov.sort_unstable();
    let covlast = cov.len() - 1;
    let high = cov[((covlast as f64) * 0.10) as usize];
    let low = cov[((covlast as f64) * 0.75) as usize];
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
    let _ = low;
    ReadStats {
        true_depth,
        depth_al,
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
