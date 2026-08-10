//! BBTools `kmercountexact`-compatible histogram text and peak calling.
//!
//! Ports `jgi.KmerCountExact` output: the `#Depth Count logScale` k-mer
//! histogram text and `CallPeaks` peak detection/summary, byte-identical to
//! `kmercountexact.sh` with default parameters (see
//! tests/bbtools/Lambda/README.md).

use crate::libs::kmer::KmerTable;
use anyhow::Result;
use std::io::Write;

/// Default histogram top bin (kmercountexact `histmax`).
pub const HIST_MAX: usize = 100000;
/// `logwidth` of the k-mer count log scaling.
const LOG_WIDTH: f64 = 0.1;
/// Peak detection limits (kmercountexact defaults).
const MIN_HEIGHT: u64 = 2;
const MIN_VOLUME: u64 = 5;
const MIN_WIDTH: usize = 3;
const MIN_PEAK: usize = 2;
const MAX_PEAK_COUNT: usize = 12;
const MAX_WIDTH_MULT: f64 = 2.5;

/// A called peak over a k-mer count histogram.
#[derive(Debug, Clone)]
pub struct Peak {
    /// Inclusive start depth.
    pub start: usize,
    /// Peak center depth.
    pub center: usize,
    /// Exclusive stop depth.
    pub stop: usize,
    /// Depth of the maximum count.
    pub max_pos: usize,
    /// Maximum count in the peak.
    pub max_height: u64,
    /// Count at `start`.
    pub start_height: u64,
    /// Count at `stop`.
    pub stop_height: u64,
    /// Minimum count left of center.
    pub left_min: u64,
    /// Minimum count right of center.
    pub right_min: u64,
    /// Sum of counts in `start..stop`.
    pub volume: u64,
    /// Sum of `count * depth` in `start..stop`.
    pub volume2: u64,
}

/// Builds the depth histogram of a k-mer table (counts capped at `max`).
pub fn histogram(table: &KmerTable, max: usize) -> Vec<u64> {
    let mut h = vec![0u64; max + 1];
    for &c in &table.counts {
        h[(c as usize).min(max)] += 1;
    }
    h
}

/// BBTools `CallPeaks.logScale`: weighted window smoothing over `passes`.
pub fn log_scale(array: &[u64], width: f64, scale: f64, passes: usize) -> Vec<u64> {
    let mut log = array.to_vec();
    for _ in 0..passes {
        let half_width = width / 2.0;
        let limit = array.len() as f64 - 0.00001;
        let mut out = vec![0u64; array.len()];
        for (pos, out_v) in out.iter_mut().enumerate().skip(1) {
            let center = pos as f64 + 0.5;
            let min = (center - half_width * pos as f64).max(0.0);
            let max = (center + half_width * pos as f64).min(limit);
            let mini = min as usize;
            let maxi = max as usize;
            if mini == maxi {
                *out_v = ((max - min) * log[mini] as f64 * scale).round() as u64;
            } else {
                let mut sum = log[mini] as f64 * (mini as f64 + 1.0 - min);
                sum += log[maxi] as f64 * (max - maxi as f64);
                for &v in &log[mini + 1..maxi] {
                    sum += v as f64;
                }
                *out_v = (sum * scale).round() as u64;
            }
        }
        log = out;
    }
    log
}

/// Writes the `#Depth Count logScale` text histogram (kmercountexact khist).
pub fn write_khist_text<W: Write>(out: &mut W, hist: &[u64], max: usize) -> Result<()> {
    let scaled = log_scale(hist, LOG_WIDTH, 1.0, 1);
    writeln!(out, "#Depth\tCount\tlogScale")?;
    for i in 1..=max.min(hist.len() - 1) {
        if hist[i] > 0 {
            writeln!(out, "{}\t{}\t{}", i, hist[i], scaled[i])?;
        }
    }
    Ok(())
}

/// Calls peaks over the histogram (BBTools CallPeaks.callPeaks).
pub fn call_peaks(original: &[u64]) -> Vec<Peak> {
    let scaled = log_scale(original, LOG_WIDTH, 1.0, 1);
    let mut peaks = detect(&scaled, original.len());
    cap_width(&mut peaks, MAX_WIDTH_MULT, &scaled);
    if peaks.len() > MAX_PEAK_COUNT {
        peaks = condense(peaks, MAX_PEAK_COUNT);
    }
    cap_width(&mut peaks, MAX_WIDTH_MULT, &scaled);
    if peaks.len() > 1 {
        let biggest = peaks[biggest_peak(&peaks)].volume;
        while peaks.len() > 1 && peaks[0].volume * 10000 < biggest {
            peaks.remove(0);
        }
    }
    let mut out: Vec<Peak> = Vec::new();
    recalculate(&mut peaks, original);
    for p in peaks {
        if p.volume >= MIN_VOLUME {
            out.push(p);
        }
    }
    out
}

/// Peak detection over a (possibly log-scaled) histogram.
fn detect(array: &[u64], length: usize) -> Vec<Peak> {
    let mut peaks = Vec::new();
    let mut dip0 = -1i64;
    for i in 1..length {
        if array[i - 1] < array[i] {
            dip0 = i as i64 - 1;
            break;
        }
    }
    if dip0 < 0 {
        return peaks;
    }
    let mut mode = 0usize; // 0 = UP, 1 = DOWN
    let mut start = dip0 as usize;
    let mut center = -1i64;
    let mut prev = array[dip0 as usize];
    let mut sum = prev;
    let mut sum2 = prev * dip0 as u64;
    let mut i = dip0 + 1;
    while i < length as i64 {
        let x = array[i as usize];
        if mode == 0 {
            if x < prev {
                mode = 1;
                center = i - 1;
            }
        } else if x > prev {
            mode = 0;
            let stop = i - 1;
            let max = array[center as usize];
            if (center as usize) >= MIN_PEAK
                && max >= MIN_HEIGHT
                && (stop as usize - start) >= MIN_WIDTH
                && sum >= MIN_VOLUME
            {
                // Middle of mesa.
                let mut c = center;
                let mut j = center - 1;
                while j >= 0 {
                    if array[j as usize] != max {
                        c = (center + j + 2) / 2;
                        break;
                    }
                    j -= 1;
                }
                // Middle of valley.
                let valley = array[stop as usize];
                let mut s = stop;
                let mut j = stop;
                while j >= 0 {
                    if array[j as usize] != valley {
                        if valley == 0 {
                            s = j + 1;
                        } else {
                            s = (stop + j + 2) / 2;
                        }
                        break;
                    }
                    j -= 1;
                }
                let start_h = array[start];
                let stop_h = array[s as usize];
                peaks.push(Peak {
                    start,
                    center: c as usize,
                    stop: s as usize,
                    max_pos: c as usize,
                    max_height: max,
                    start_height: start_h,
                    stop_height: stop_h,
                    left_min: start_h,
                    right_min: stop_h,
                    volume: sum,
                    volume2: sum2,
                });
            }
            start = stop as usize;
            sum = 0;
            sum2 = 0;
            center = -1;
            while i < length as i64 && array[i as usize] == 0 {
                i += 1;
            }
        }
        sum += x;
        sum2 += x * i as u64;
        prev = x;
        i += 1;
    }
    if mode == 1 {
        let stop = length as i64;
        let max = array[center as usize];
        let mut c = center;
        let mut j = center - 1;
        while j >= 0 {
            if array[j as usize] != max {
                c = (center + j + 2) / 2;
                break;
            }
            j -= 1;
        }
        let valley = array[(stop - 1) as usize];
        let mut s = stop - 1;
        let mut j = stop - 1;
        while j >= 0 {
            if array[j as usize] != valley {
                if valley == 0 {
                    s = j + 1;
                } else {
                    s = (stop - 1 + j + 2) / 2;
                }
                break;
            }
            j -= 1;
        }
        let stop_final = s.min(length as i64 - 1);
        let start_h = array[start];
        let stop_h = array[stop_final as usize];
        if (c as usize) >= MIN_PEAK
            && max >= MIN_HEIGHT
            && (stop_final as usize - start) >= MIN_WIDTH
            && sum >= MIN_VOLUME
        {
            peaks.push(Peak {
                start,
                center: c as usize,
                stop: stop_final as usize,
                max_pos: c as usize,
                max_height: max,
                start_height: start_h,
                stop_height: stop_h,
                left_min: start_h,
                right_min: stop_h,
                volume: sum,
                volume2: sum2,
            });
        }
    }
    peaks
}

/// Narrows peaks to `center/mult .. center*mult` (CallPeaks.capWidth).
fn cap_width(peaks: &mut [Peak], max_width_mult: f64, counts: &[u64]) {
    let mult = 1.0 / max_width_mult;
    for p in peaks.iter_mut() {
        let start = (p.start as f64).max(p.center as f64 * mult);
        let stop = (p.stop as f64).min(p.center as f64 * max_width_mult);
        p.start = start.round() as usize;
        p.stop = stop.round() as usize;
    }
    recalculate(peaks, counts);
}

/// Keeps the top peaks by height/volume thresholds (CallPeaks.condense).
fn condense(mut peaks: Vec<Peak>, max_count: usize) -> Vec<Peak> {
    let max_count = max_count.min(peaks.len()).max(1);
    let mut heights: Vec<u64> = peaks.iter().map(|p| p.max_height).collect();
    heights.sort_unstable();
    let hlimit = heights[heights.len() - max_count];
    let mc2 = (max_count + 1).div_ceil(2);
    let mut volumes: Vec<u64> = peaks.iter().map(|p| p.volume).collect();
    volumes.sort_unstable();
    let vlimit = volumes[volumes.len() - mc2];
    peaks.retain(|p| p.volume >= vlimit || p.max_height >= hlimit);
    peaks
}

/// Recomputes peak fields from a histogram (`Peak.recalculate`).
fn recalculate(peaks: &mut [Peak], array: &[u64]) {
    for p in peaks.iter_mut() {
        p.max_height = array[p.center];
        p.start_height = array[p.start];
        p.stop_height = array[p.stop];
        p.left_min = p.start_height;
        p.right_min = p.stop_height;
        p.max_pos = p.center;
        p.volume = 0;
        p.volume2 = 0;
        for (i, &x) in array.iter().enumerate().take(p.stop).skip(p.start) {
            if x > p.max_height {
                p.max_pos = i;
                p.max_height = x;
            }
            if i < p.center {
                p.left_min = p.left_min.min(x);
            } else if i > p.center {
                p.right_min = p.right_min.min(x);
            }
            p.volume += x;
            p.volume2 += x * i as u64;
        }
    }
}

fn biggest_peak(peaks: &[Peak]) -> usize {
    if peaks.len() < 2 {
        return peaks.len() - 1;
    }
    let mut loc = 0;
    for (i, p) in peaks.iter().enumerate() {
        if p.volume > peaks[loc].volume {
            loc = i;
        }
    }
    loc
}

fn second_biggest_peak(peaks: &[Peak]) -> usize {
    if peaks.len() < 2 {
        return peaks.len() - 1;
    }
    let mut biggest = &peaks[0];
    let mut second = &peaks[1];
    let mut bloc = 0usize;
    let mut sloc = 1usize;
    if second.volume > biggest.volume {
        std::mem::swap(&mut biggest, &mut second);
        bloc = 1;
        sloc = 0;
    }
    for (i, p) in peaks.iter().enumerate().skip(2) {
        if p.volume > second.volume {
            sloc = i;
            second = p;
            if second.volume > biggest.volume {
                std::mem::swap(&mut second, &mut biggest);
                sloc = bloc;
                bloc = i;
            }
        }
    }
    sloc
}

fn homozygous_peak(peaks: &[Peak], ploidy: usize, haploid_center: f64) -> usize {
    if peaks.len() < 2 {
        return peaks.len() - 1;
    }
    let target = haploid_center * ploidy as f64;
    let mut best = f64::MAX;
    let mut loc = 0;
    for (i, p) in peaks.iter().enumerate() {
        let dif = (target - p.center as f64).abs();
        if dif < best {
            best = dif;
            loc = i;
        }
    }
    loc
}

fn haploid_peak_center(peaks: &[Peak], ploidy: usize) -> f64 {
    let biggest = &peaks[biggest_peak(peaks)];
    let second = &peaks[second_biggest_peak(peaks)];
    if second.volume * 4 >= biggest.volume {
        biggest.center.min(second.center) as f64
    } else {
        biggest.center as f64 / ploidy as f64
    }
}

fn single_copy_kmer_fraction(het_rate: f32, k: usize, ploidy: usize) -> f32 {
    if ploidy < 2 {
        return 1.0;
    }
    let kmers_per_snp = k as f32;
    let single_copy = het_rate * kmers_per_snp;
    let asymptote = single_copy / (1.0 + single_copy);
    asymptote * 2.0
}

fn calc_ploidy(peaks: &[Peak], min_volume_fraction: f32) -> usize {
    if peaks.len() < 2 {
        return 1;
    }
    let biggest = &peaks[biggest_peak(peaks)];
    let second = &peaks[second_biggest_peak(peaks)];
    // ploidyLogic == 2
    if std::ptr::eq(second, biggest) {
        return 1;
    }
    if second.center < biggest.center {
        if (second.volume as f32) < biggest.volume as f32 * min_volume_fraction {
            return 1;
        }
    } else if second.volume * 4 < biggest.volume {
        return 1;
    }
    let max = biggest.center.max(second.center);
    let min = biggest.center.min(second.center);
    ((max as f64 / min as f64).round()) as usize
}

fn error_kmers(peaks: &[Peak], hist: &[u64], min_volume_fraction: f32) -> u64 {
    if peaks.is_empty() {
        return 0;
    }
    let first = first_genomic_peak(peaks, min_volume_fraction);
    if first.is_none() {
        return 0;
    }
    let start = peaks[first.unwrap()].start;
    hist[..start].iter().sum()
}

fn first_genomic_peak(peaks: &[Peak], min_fraction: f32) -> Option<usize> {
    let biggest = peaks[biggest_peak(peaks)].volume;
    let min_volume = (biggest as f32 * min_fraction) as u64;
    peaks.iter().position(|p| p.volume >= min_volume)
}

fn genome_size_in_peaks(peaks: &[Peak], haploid_center: f64) -> u64 {
    if peaks.is_empty() {
        return 0;
    }
    let mult = 1.0 / haploid_center.max(1.0);
    peaks
        .iter()
        .map(|p| p.volume * ((p.center as f64 * mult).round() as u64))
        .sum()
}

fn genome_size2(peaks: &[Peak], haploid_center: f64, hist: &[u64]) -> u64 {
    if peaks.is_empty() {
        return 0;
    }
    let mult = 1.0 / haploid_center.max(1.0);
    let start = peaks[0].start;
    let mut sum = 0u64;
    for (i, &v) in hist.iter().enumerate().skip(start) {
        sum += v * ((i as f64 * mult).round().max(1.0) as u64);
    }
    sum
}

fn repeat_size(peaks: &[Peak], ploidy: usize, haploid_center: f64) -> u64 {
    if peaks.len() < 2 {
        return 0;
    }
    let hom = homozygous_peak(peaks, ploidy, haploid_center);
    let mult = 1.0 / haploid_center.max(1.0);
    let mut sum = 0u64;
    for p in &peaks[hom + 1..] {
        sum += p.volume * (((p.center as f64 * mult).round() as u64).saturating_sub(1));
    }
    sum
}

fn repeat_size2(peaks: &[Peak], ploidy: usize, haploid_center: f64, hist: &[u64]) -> u64 {
    let mult = 1.0 / haploid_center.max(1.0);
    let mut valley =
        (haploid_center * ploidy as f64 * (1.2 + 1.0 / ploidy.max(2) as f64)).ceil() as usize;
    let hom = homozygous_peak(peaks, ploidy, haploid_center);
    if ploidy > 1 && hom < peaks.len() {
        valley = peaks[hom].stop + 1;
    }
    let mut sum = 0u64;
    for (i, &v) in hist.iter().enumerate().skip(valley) {
        let copies = ((i as f64 * mult).round() as u64).saturating_sub(1);
        sum += v * copies;
    }
    sum
}

fn calc_het_locations(peaks: &[Peak], ploidy: usize, haploid_center: f64, k: usize) -> u64 {
    if peaks.len() < 2 {
        return 0;
    }
    let hom = homozygous_peak(peaks, ploidy, haploid_center);
    let homo_center = peaks[hom].center as f64;
    let lim = ploidy / 2;
    let mut sum = 0u64;
    for p in &peaks[..hom] {
        let copy_count = ((p.center as f64 * ploidy as f64) / homo_center).round() as usize;
        if copy_count > lim {
            break;
        }
        sum += p.volume;
    }
    sum / k as u64
}

/// Writes the kmercountexact `peaks.txt` summary and peak table.
pub fn write_peaks_text<W: Write>(
    out: &mut W,
    peaks: &[Peak],
    k: usize,
    unique_kmers: u64,
    hist: &[u64],
) -> Result<()> {
    let min_het_rate = 0.0003f32;
    let min_volume_fraction = single_copy_kmer_fraction(min_het_rate, k, 2).min(1.0);
    let ploidy_est = calc_ploidy(peaks, min_volume_fraction);
    let ploidy = ploidy_est;
    let hap = haploid_peak_center(peaks, ploidy);
    let err = error_kmers(peaks, hist, min_volume_fraction);
    let gs_peaks = genome_size_in_peaks(peaks, hap);
    let gs_total = genome_size2(peaks, hap, hist);
    let rep = repeat_size(peaks, ploidy, hap);
    let rep2 = repeat_size2(peaks, ploidy, hap, hist);
    let hap_size = gs_total / ploidy as u64;
    let het = calc_het_locations(peaks, ploidy, hap, k);
    let het_rate = (het as f64 / hap_size as f64) / 2.0;
    let repeat_rate = rep as f64 / gs_peaks as f64;
    let repeat_rate2 = rep2 as f64 / gs_total as f64;

    let mut main_peak = &peaks[0];
    let mut ploidy_peak = &peaks[0];
    let target = hap * ploidy as f64;
    for p in peaks {
        if p.volume > main_peak.volume {
            main_peak = p;
        }
        if (p.center as f64 - target).abs() < (ploidy_peak.center as f64 - target).abs() {
            ploidy_peak = p;
        }
    }
    let haploid_cov =
        if target.max(ploidy_peak.center as f64) / target.min(ploidy_peak.center as f64) < 1.3 {
            ploidy_peak.center as i64
        } else {
            target.round() as i64
        };

    writeln!(out, "#k\t{}", k)?;
    writeln!(out, "#unique_kmers\t{}", unique_kmers)?;
    writeln!(out, "#error_kmers\t{}", err)?;
    writeln!(out, "#genomic_kmers\t{}", unique_kmers - err)?;
    writeln!(out, "#main_peak\t{}", main_peak.center)?;
    writeln!(out, "#genome_size_in_peaks\t{}", gs_peaks)?;
    writeln!(out, "#genome_size\t{}", gs_total)?;
    writeln!(out, "#haploid_genome_size\t{}", hap_size)?;
    writeln!(out, "#fold_coverage\t{}", hap.round())?;
    writeln!(out, "#haploid_fold_coverage\t{}", haploid_cov)?;
    writeln!(out, "#ploidy\t{}", ploidy)?;
    writeln!(out, "#percent_repeat_in_peaks\t{:.3}", 100.0 * repeat_rate)?;
    writeln!(out, "#percent_repeat\t{:.3}", 100.0 * repeat_rate2)?;
    let _ = het_rate;
    writeln!(out, "#start\tcenter\tstop\tmax\tvolume")?;
    for p in peaks {
        if p.volume >= MIN_VOLUME {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                p.start, p.center, p.stop, p.max_height, p.volume
            )?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_scale_matches_bbtools_rounding() {
        let mut h = vec![0u64; 100];
        h[1] = 38961;
        h[56] = 2077;
        let scaled = log_scale(&h, 0.1, 1.0, 1);
        assert_eq!(scaled[1], 3896);
        assert_eq!(scaled[56], 2077);
    }
}
