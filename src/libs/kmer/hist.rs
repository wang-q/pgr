//! FASTK-compatible k-mer frequency histogram (`.hist`).

use super::khist;
use super::KmerTable;
use anyhow::Context;
use std::io::Write;
use std::path::Path;

/// Highest histogram bin; counts above it are folded into this bin.
///
/// Matches FastK's fixed `low=1` / `high=32767` histogram range so the
/// written file is byte-compatible with `Histex` / `KatGC` / GenomeScope.
pub const HIST_HIGH: u32 = 0x7fff;

/// Fixed `.hist` file size in bytes: 28-byte header + `32767 * 8` payload.
pub const HIST_FILE_LEN: usize = 28 + HIST_HIGH as usize * 8;

/// K-mer frequency histogram with FastK-compatible fixed bins `1..=32767`.
#[derive(Debug, Clone)]
pub struct Histogram {
    /// K-mer length (bp).
    pub k: usize,
    /// `hist[i-1]` = number of distinct k-mers with count `i` (`i in 1..=32767`).
    pub hist: Vec<u64>,
    /// Instance-mode count folded into the top bin (FastK `max_inst`).
    pub max_inst: u64,
}

/// Coverage statistics and a genome-size estimate from a k-mer table.
#[derive(Debug, Clone, Copy)]
pub struct GenomeEstimate {
    /// Coverage at the histogram peak (main k-mer frequency mode).
    pub peak_cov: u64,
    /// Total number of distinct k-mers.
    pub total_distinct: u64,
    /// Total k-mer instances (sum of all counts).
    pub total_kmers: u128,
    /// Estimated genome size = total k-mer instances / peak coverage.
    pub genome_size: f64,
}

/// Estimate coverage and genome size from a count table.
///
/// `peak_cov` is the main k-mer frequency peak (BBTools CallPeaks, skipping
/// the error k-mers that dominate real reads); `genome_size =
/// total_kmers / peak_cov` approximates the haploid length (cheap
/// pre-GenomeScope estimate; model fitting is a separate step).
pub fn estimate(table: &KmerTable) -> GenomeEstimate {
    let total_distinct = table.counts.len() as u64;
    let total_kmers: u128 = table.counts.iter().map(|&c| c as u128).sum();
    let peak_cov = khist::call_peaks(&khist::histogram(table, khist::HIST_MAX))
        .into_iter()
        .max_by_key(|p| p.volume)
        .map(|p| p.center as u64)
        .unwrap_or(0);
    let genome_size = if peak_cov > 0 {
        total_kmers as f64 / peak_cov as f64
    } else {
        0.0
    };
    GenomeEstimate {
        peak_cov,
        total_distinct,
        total_kmers,
        genome_size,
    }
}

/// Aggregate `table.counts` into the fixed `1..=32767` frequency bins.
///
/// Counts above 32767 are folded into the top bin and their instances are
/// accumulated into `max_inst`, mirroring FastK's counting semantics.
pub fn from_table(table: &KmerTable) -> Histogram {
    let mut hist = vec![0u64; HIST_HIGH as usize];
    let mut max_inst = 0u64;
    for &c in &table.counts {
        let c = c as usize;
        if c >= HIST_HIGH as usize {
            hist[HIST_HIGH as usize - 1] += 1;
            max_inst += c as u64;
        } else if c > 0 {
            hist[c - 1] += 1;
        }
    }
    Histogram {
        k: table.k,
        hist,
        max_inst,
    }
}

/// Write `hist` in the FastK `.hist` binary layout (all fields little-endian).
///
/// Header is `int32 k, int32 low=1, int32 high=32767, int64 ilowcnt
/// (= bin 1 count), int64 max_inst`, followed by the `32767` `int64` bins.
pub fn write(path: &Path, hist: &Histogram) -> anyhow::Result<()> {
    anyhow::ensure!(
        hist.hist.len() == HIST_HIGH as usize,
        "histogram must have {HIST_HIGH} bins"
    );
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    w.write_all(&(hist.k as u32).to_le_bytes())
        .context("writing hist header")?;
    w.write_all(&1u32.to_le_bytes())
        .context("writing hist header")?;
    w.write_all(&HIST_HIGH.to_le_bytes())
        .context("writing hist header")?;
    w.write_all(&hist.hist[0].to_le_bytes())
        .context("writing hist header")?;
    w.write_all(&hist.max_inst.to_le_bytes())
        .context("writing hist header")?;
    for &v in &hist.hist {
        w.write_all(&v.to_le_bytes()).context("writing hist bins")?;
    }
    w.flush().context("flushing hist file")?;
    Ok(())
}

/// Read a FastK-compatible `.hist` file.
pub fn load(path: &Path) -> anyhow::Result<Histogram> {
    let bytes = std::fs::read(path).context("reading hist file")?;
    anyhow::ensure!(
        bytes.len() == HIST_FILE_LEN,
        "not a fixed-layout .hist file ({} bytes, expected {HIST_FILE_LEN})",
        bytes.len()
    );
    let k = u32::from_le_bytes(bytes[0..4].try_into().unwrap()) as usize;
    let low = u32::from_le_bytes(bytes[4..8].try_into().unwrap());
    let high = u32::from_le_bytes(bytes[8..12].try_into().unwrap());
    anyhow::ensure!(
        low == 1 && high == HIST_HIGH,
        "unsupported .hist range [{low},{high}]"
    );
    let max_inst = u64::from_le_bytes(bytes[20..28].try_into().unwrap());
    let mut hist = vec![0u64; HIST_HIGH as usize];
    for (i, chunk) in bytes[28..].chunks_exact(8).enumerate() {
        hist[i] = u64::from_le_bytes(chunk.try_into().unwrap());
    }
    Ok(Histogram { k, hist, max_inst })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_table_accumulates_frequency_bins() {
        let table = KmerTable {
            k: 4,
            keys: vec![0, 1, 2, 3, 4, 5, 6, 7],
            counts: vec![1, 2, 2, 3, 5, HIST_HIGH, HIST_HIGH + 10, 0],
        };
        let h = from_table(&table);
        assert_eq!(h.hist[0], 1); // count 1
        assert_eq!(h.hist[1], 2); // count 2
        assert_eq!(h.hist[2], 1); // count 3
        assert_eq!(h.hist[4], 1); // count 5
        assert_eq!(h.hist[HIST_HIGH as usize - 1], 2); // both folded
        assert_eq!(h.max_inst, (HIST_HIGH + HIST_HIGH + 10) as u64);
        // count 0 is skipped (does not occur in real tables)
    }

    #[test]
    fn write_matches_fastk_layout() {
        let table = KmerTable {
            k: 17,
            keys: vec![0, 1, 2],
            counts: vec![1, 2, 32767],
        };
        let h = from_table(&table);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.hist");
        write(&path, &h).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(bytes.len(), HIST_FILE_LEN);
        assert_eq!(&bytes[0..4], &17u32.to_le_bytes());
        assert_eq!(&bytes[4..8], &1u32.to_le_bytes()); // low
        assert_eq!(&bytes[8..12], &HIST_HIGH.to_le_bytes()); // high
        assert_eq!(&bytes[12..20], &1u64.to_le_bytes()); // ilowcnt = bin 1
        assert_eq!(&bytes[20..28], &32767u64.to_le_bytes()); // max_inst
                                                             // payload bins: index 0 (count 1), 1 (count 2), 32766 (folded)
        let bin = |i: usize| -> u64 {
            let off = 28 + i * 8;
            u64::from_le_bytes(bytes[off..off + 8].try_into().unwrap())
        };
        assert_eq!(bin(0), 1);
        assert_eq!(bin(1), 1);
        assert_eq!(bin(HIST_HIGH as usize - 1), 1);
    }

    #[test]
    fn load_roundtrip() {
        let table = KmerTable {
            k: 9,
            keys: vec![0, 1],
            counts: vec![4, 500],
        };
        let h = from_table(&table);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.hist");
        write(&path, &h).unwrap();
        let loaded = load(&path).unwrap();
        assert_eq!(loaded.k, 9);
        assert_eq!(loaded.hist, h.hist);
        assert_eq!(loaded.max_inst, h.max_inst);
    }

    #[test]
    fn load_rejects_foreign_layout() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.hist");
        std::fs::write(&path, vec![0u8; 100]).unwrap();
        assert!(load(&path).is_err(), "wrong size must be rejected");

        let mut bytes = vec![0u8; HIST_FILE_LEN];
        bytes[8..12].copy_from_slice(&1000u32.to_le_bytes()); // high != 32767
        std::fs::write(&path, &bytes).unwrap();
        assert!(load(&path).is_err(), "foreign range must be rejected");
    }

    #[test]
    fn estimate_peak_and_genome_size() {
        // Synthetic 30x coverage: 300 reads of 100 bp sampled from a
        // deterministic 1 kb genome. The peak should sit near 30x and the
        // estimated genome size near 1000 bp.
        let genome = random_genome(1000, 42);
        let reads = sample_reads(&genome, 100, 300, 7);
        let table = crate::libs::kmer::count::build_table(&reads, 17).unwrap();
        let est = estimate(&table);
        assert!(
            (8..=60).contains(&est.peak_cov),
            "peak coverage {} far from 30x",
            est.peak_cov
        );
        assert_eq!(est.total_distinct, table.counts.len() as u64);
        let ratio = est.genome_size / 1000.0;
        assert!(
            (0.5..=2.0).contains(&ratio),
            "genome size {} far from 1000 bp",
            est.genome_size
        );
    }

    /// Deterministic pseudo-random genome.
    fn random_genome(len: usize, seed: u64) -> Vec<u8> {
        let bases = *b"ACGT";
        let mut x = seed;
        (0..len)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                bases[(x >> 33) as usize & 3]
            })
            .collect()
    }

    /// Sample fixed-length reads at random positions (with replacement).
    fn sample_reads(genome: &[u8], read_len: usize, n: usize, seed: u64) -> Vec<Vec<u8>> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let start = (x >> 33) as usize % (genome.len() - read_len + 1);
                genome[start..start + read_len].to_vec()
            })
            .collect()
    }
}
