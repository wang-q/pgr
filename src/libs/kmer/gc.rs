//! GC-content vs k-mer frequency matrix (MerquryFK KatGC equivalent).

use super::KmerTable;
use std::io::Write;

/// Per-byte GC counts: 4 packed 2-bit bases per byte (C/G are GC).
fn gc_tables(k: usize) -> ([u8; 256], [u8; 256]) {
    let is_gc = |code: u8| usize::from(code == 1 || code == 2);
    let mut gc = [0u8; 256];
    for (x, g) in gc.iter_mut().enumerate() {
        let x = x as u8;
        *g = (is_gc(x >> 6) + is_gc((x >> 4) & 3) + is_gc((x >> 2) & 3) + is_gc(x & 3)) as u8;
    }
    // Trailing-byte table: only the leading k%4 bases of the last byte
    // belong to the k-mer (bases are high-aligned within the byte).
    let mut gcr = [0u8; 256];
    for (x, g) in gcr.iter_mut().enumerate() {
        let x = x as u8;
        *g = match k % 4 {
            0 => gc[x as usize],
            1 => is_gc(x >> 6) as u8,
            2 => (is_gc(x >> 6) + is_gc((x >> 4) & 3)) as u8,
            _ => (is_gc(x >> 6) + is_gc((x >> 4) & 3) + is_gc((x >> 2) & 3)) as u8,
        };
    }
    (gc, gcr)
}

/// Number of G/C bases in a packed 2-bit k-mer (A/T non-GC, C/G GC).
fn gc_count(key: u128, k: usize, gc: &[u8; 256], gcr: &[u8; 256]) -> usize {
    let kbyte = k.div_ceil(4);
    let mut packed = [0u8; 16];
    crate::libs::pgi::pack_kmer(key, k, &mut packed);
    let mut cnt = 0usize;
    for &b in &packed[..kbyte - 1] {
        cnt += gc[b as usize] as usize;
    }
    cnt + gcr[packed[kbyte - 1] as usize] as usize
}

/// GC (0..=k) × count (0..=hmax) matrix: `plot[gc][min(count, hmax)] += 1`.
///
/// `hmax` caps the count axis exactly like KatGC's `HMAX` (default 1000).
pub fn gc_matrix(table: &KmerTable, hmax: usize) -> Vec<Vec<u64>> {
    let (gc, gcr) = gc_tables(table.k);
    let mut plot = vec![vec![0u64; hmax + 1]; table.k + 1];
    for (key, &count) in table.keys.iter().zip(&table.counts) {
        let g = gc_count(*key, table.k, &gc, &gcr);
        plot[g][(count as usize).min(hmax)] += 1;
    }
    plot
}

/// Count-axis peak away from 0 (KatGC `xmax` / `ZMAX`).
#[derive(Debug, Clone, Copy)]
pub struct GcPeak {
    /// Count value at the global peak.
    pub xmax: usize,
    /// Highest bin value found at the peak.
    pub zmax: u64,
}

/// Find the count-axis peak away from zero, mirroring KatGC's search.
pub fn find_peak(plot: &[Vec<u64>], hmax: usize) -> anyhow::Result<GcPeak> {
    anyhow::ensure!(hmax >= 2, "hmax must be at least 2, got {hmax}");
    let mut zmax = 0u64;
    let mut xmax = 0usize;
    for row in plot {
        let mut k = 2usize;
        while k < hmax && row[k] < row[k - 1] {
            k += 1;
        }
        let mut ym = row[k];
        let mut xm = k;
        for (kk, &v) in row.iter().enumerate().skip(k) {
            if v >= ym {
                ym = v;
                xm = kk;
            }
        }
        if ym > zmax {
            zmax = ym;
            xmax = xm;
        }
    }
    if xmax == 0 || xmax >= hmax {
        anyhow::bail!("no maximal peak away from 0 in histogram interval [1,{hmax}]");
    }
    Ok(GcPeak { xmax, zmax })
}

/// Output x-limit: `xrel` times the peak (KatGC `XMAX = XREL * xmax`),
/// clamped to `hmax`.
pub fn x_limit(peak: GcPeak, xrel: f64, hmax: usize) -> usize {
    ((xrel * peak.xmax as f64) as usize).min(hmax)
}

/// Write the 4-neighbor-averaged matrix in KatGC `.kgc` format.
///
/// Rows are GC counts `0..k-1`, columns are count bins `0..xmax-1`; each
/// value averages the 2×2 block and is clamped to the peak value `zmax`.
pub fn write_kgc(
    w: &mut impl Write,
    plot: &[Vec<u64>],
    xmax: usize,
    zmax: u64,
) -> anyhow::Result<()> {
    writeln!(w, "GCP\tKF\tCount")?;
    for i in 0..plot.len() - 1 {
        let row0 = &plot[i];
        let row1 = &plot[i + 1];
        for a in 0..xmax {
            let val = (row0[a] + row0[a + 1] + row1[a] + row1[a + 1]) / 4;
            writeln!(w, "{}.5\t{}.5\t{}", i, a, val.min(zmax))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gc_naive(key: u128, k: usize) -> usize {
        let mut n = 0;
        for i in 0..k {
            let code = ((key >> (2 * (k - 1 - i))) & 3) as u8;
            if code == 1 || code == 2 {
                n += 1;
            }
        }
        n
    }

    #[test]
    fn gc_count_matches_naive() {
        for k in [1usize, 4, 7, 8, 17, 21, 51, 64] {
            let (gc, gcr) = gc_tables(k);
            for seed in 0..200u128 {
                // Pseudo-random key confined to 2k bits.
                let mask = if 2 * k >= 128 {
                    u128::MAX
                } else {
                    (1u128 << (2 * k)) - 1
                };
                let key = (seed * 0x9e3779b97f4a7c15 + 0xbf58476d1ce4e5b9) & mask;
                assert_eq!(
                    gc_count(key, k, &gc, &gcr),
                    gc_naive(key, k),
                    "k={k} seed={seed}"
                );
            }
        }
    }

    #[test]
    fn gc_matrix_accumulates_gc_and_count() {
        // k=4 keys: 0b00_00_00_00 (AAAA, gc=0), 0b00_01_10_00 (ACGT, gc=2),
        // 0b01_01_10_10 (CCGG, gc=4).
        let table = KmerTable {
            k: 4,
            keys: vec![0, 0b00011000, 0b01011010],
            counts: vec![1, 3, 5000],
        };
        let hmax = 100;
        let plot = gc_matrix(&table, hmax);
        assert_eq!(plot.len(), 5); // gc 0..=4
        assert_eq!(plot[0][1], 1); // AAAA, count 1
        assert_eq!(plot[2][3], 1); // ACGT, count 3
        assert_eq!(plot[4][hmax], 1); // CCGG, count capped at hmax
    }

    #[test]
    fn peak_and_x_limit_match_katgc() {
        // A single GC row with a peak at count 20; the other rows are empty
        // (empty rows contribute no peak).
        let hmax = 1000;
        let mut plot = vec![vec![0u64; hmax + 1]; 5];
        for (c, v) in plot[2].iter_mut().enumerate() {
            *v = if c == 20 { 1000 } else { 1 };
        }
        let peak = find_peak(&plot, hmax).unwrap();
        assert_eq!(peak.xmax, 20);
        assert_eq!(peak.zmax, 1000);
        assert_eq!(x_limit(peak, 2.1, hmax), 42);
        assert_eq!(x_limit(peak, 2.1, 30), 30); // clamped to hmax
    }

    #[test]
    fn write_kgc_averages_and_clamps() {
        let mut plot = vec![vec![0u64; 6]; 3];
        plot[0] = vec![0, 10, 10, 0, 0, 0];
        plot[1] = vec![0, 0, 0, 0, 0, 0];
        plot[2] = vec![0, 0, 0, 0, 0, 0];
        let mut out = Vec::new();
        write_kgc(&mut out, &plot, 3, 5).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "GCP\tKF\tCount");
        assert_eq!(lines.len(), 1 + 2 * 3); // rows 0..=1 × cols 0..=2
                                            // a=0: (0+10+0+0)/4 = 2.
        assert_eq!(lines[1], "0.5\t0.5\t2");
        // a=1: (10+10+0+0)/4 = 5, clamped to zmax 5.
        assert_eq!(lines[2], "0.5\t1.5\t5");
        // a=2: (10+0+0+0)/4 = 2.
        assert_eq!(lines[3], "0.5\t2.5\t2");
    }
}
