//! Paired-read overlap merging and overlap-based error correction
//! (BBMerge-compatible).

use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use crate::libs::fq::bbnet::CellNet;
use crate::libs::fq::overlap;
use crate::libs::fq::qual::{base_to_number, from_phred, to_phred};
use crate::libs::fq::tadpole::{extend_read_right, TadpoleOptions, TadpoleTable};
use crate::libs::nt::rev_comp;
use anyhow::Result;
use std::io::Write;

/// Histogram capacity (BBMerge `histlen`).
pub const HIST_LEN: usize = 2000;

/// Overlap-detection presets, mirroring `bbmerge.sh` `strict`/`vstrict`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Preset {
    /// Default (no `strict`/`vstrict` flag).
    Normal,
    /// `strict` parameter set.
    Strict,
    /// `vstrict` parameter set.
    VStrict,
}

/// Options for `fq merge`, with BBMerge defaults.
#[derive(Debug, Clone)]
pub struct MergeOptions {
    /// `minoverlap` (MIN_OVERLAPPING_BASES).
    pub min_overlap: usize,
    /// `minoverlap0` (MIN_OVERLAPPING_BASES_0).
    pub min_overlap0: usize,
    /// `mininsert`.
    pub min_insert: usize,
    /// `mininsert0`; `None` selects the BBMerge auto value.
    pub min_insert0: Option<usize>,
    /// `maxratio`.
    pub max_ratio: f32,
    /// `ratiomargin`.
    pub ratio_margin: f32,
    /// `ratiooffset`.
    pub ratio_offset: f32,
    /// `minsecondratio`.
    pub min_second_ratio: f32,
    /// `ratiominoverlapreduction`.
    pub ratio_reduction: usize,
    /// `minentropy`.
    pub min_entropy: usize,
    /// `efilter` ratio; `None` disables the filter.
    pub efilter: Option<f32>,
    /// `efilteroffset`.
    pub efilter_offset: f32,
    /// `pfilter` ratio; 0 disables the filter.
    pub pfilter: f32,
    /// `maxbad` bound for the ratio mode (MAX_MISMATCHES_R).
    pub max_bad: usize,
    /// `ecco`: correct by overlap without joining.
    pub ecco: bool,
    /// `mix`: also write unmerged reads to the main output.
    pub mix: bool,
    /// BBMerge `MAKE_VECTOR` state: bbmerge.sh always runs with it true,
    /// which forces the ratio pre-screen `maxratio` to 0.7 and disables the
    /// ambiguity/pfilter rejections. bbmerge-auto with a tadpole resets it.
    pub make_vector: bool,
    /// Optional BBMerge overlap-filter net (bbmerge.bbnet); required when
    /// `make_vector` is true.
    pub net: Option<CellNet>,
    /// `extend2`: extend unmerged reads by up to this many bases per attempt
    /// (via the tadpole k-mer graph) and retry the overlap, mirroring
    /// `bbmerge-auto.sh ... extend2=N`. Forces `make_vector=false`.
    pub extend2: usize,
    /// `rem` / `requireExtensionMatch`: require the extended overlap to match
    /// the unextended one before accepting an extended merge.
    pub rem: bool,
}

impl MergeOptions {
    /// Options for a preset, then overridden by explicit CLI values.
    pub fn from_preset(preset: Preset) -> Self {
        match preset {
            Preset::Normal => Self {
                min_overlap: 11,
                min_overlap0: 8,
                min_insert: 15,
                min_insert0: None,
                max_ratio: 0.09,
                ratio_margin: 5.5,
                ratio_offset: 0.55,
                min_second_ratio: 0.1,
                ratio_reduction: 3,
                min_entropy: 39,
                efilter: Some(6.0),
                efilter_offset: 0.05,
                pfilter: 0.00004,
                max_bad: 20,
                ecco: false,
                mix: false,
                make_vector: true,
                net: None,
                extend2: 0,
                rem: false,
            },
            Preset::Strict => Self {
                min_overlap0: 7,
                min_entropy: 42,
                max_ratio: 0.075,
                ratio_margin: 7.5,
                ratio_reduction: 4,
                efilter: Some(4.0),
                pfilter: 0.0008,
                min_second_ratio: 0.12,
                ..Self::from_preset(Preset::Normal)
            },
            Preset::VStrict => Self {
                min_overlap: 12,
                min_overlap0: 4,
                min_entropy: 52,
                max_ratio: 0.05,
                ratio_margin: 12.0,
                ratio_offset: 0.5,
                ratio_reduction: 4,
                efilter: Some(2.0),
                pfilter: 0.008,
                min_second_ratio: 0.16,
                ..Self::from_preset(Preset::Normal)
            },
        }
    }

    /// Resolves derived parameters (`min_insert`/`min_insert0`) like the
    /// BBMerge constructor does.
    pub fn resolved(&self) -> (usize, usize) {
        let min_insert = self.min_insert.max(self.min_overlap.max(self.min_overlap0));
        let min_insert0 = match self.min_insert0 {
            Some(x) => x,
            None => {
                let auto = (min_insert as f32 * 0.75).ceil() as usize;
                auto.max(5).max(self.min_overlap0).min(35)
            }
        };
        (min_insert, min_insert0.min(min_insert))
    }
}

/// Per-pair result of overlap processing.
#[derive(Debug, Clone)]
pub enum Outcome {
    /// Joined read (bases, qualities) and its insert size.
    Merged(Vec<u8>, Vec<u8>, usize),
    /// ecco: both reads were corrected in place; insert size.
    Corrected(usize),
    /// No overlapping solution.
    NoSolution,
    /// Overlap was ambiguous.
    Ambiguous,
    /// Overlap was shorter than `mininsert`.
    TooShort(usize),
}

impl Outcome {
    fn insert(&self) -> Option<usize> {
        match self {
            Outcome::Merged(_, _, x) | Outcome::Corrected(x) | Outcome::TooShort(x) => Some(*x),
            Outcome::NoSolution | Outcome::Ambiguous => None,
        }
    }
}

/// Accumulated merge statistics and insert-size histogram.
#[derive(Debug, Clone)]
pub struct MergeStats {
    /// Pairs processed.
    pub pairs: u64,
    /// Pairs with a positive insert (mated).
    pub joined: u64,
    /// Ambiguous overlaps.
    pub ambiguous: u64,
    /// No-solution pairs.
    pub no_solution: u64,
    /// Overlaps below `mininsert`.
    pub too_short: u64,
    /// Bases corrected in ecco mode.
    pub errors_corrected: u64,
    /// Smallest observed insert.
    pub insert_min: usize,
    /// Largest observed insert.
    pub insert_max: usize,
    /// Insert-size histogram (index = insert size).
    pub hist: [u64; HIST_LEN],
}

impl Default for MergeStats {
    fn default() -> Self {
        Self {
            pairs: 0,
            joined: 0,
            ambiguous: 0,
            no_solution: 0,
            too_short: 0,
            errors_corrected: 0,
            insert_min: 0,
            insert_max: 0,
            hist: [0; HIST_LEN],
        }
    }
}

/// Processes one read pair, mirroring BBMerge's per-pair pipeline.
///
/// `r2` is left in its original orientation. In ecco mode with a positive
/// insert, both reads are corrected in place.
#[allow(clippy::too_many_arguments)]
pub fn process_pair(
    r1: &mut SeqRecord,
    r2: &mut SeqRecord,
    opts: &MergeOptions,
    stats: &mut MergeStats,
) -> Outcome {
    let (min_insert, min_insert0) = opts.resolved();
    let seq1 = r1.sequence().to_vec();
    let seq2 = r2.sequence().to_vec();
    let raw1 = r1.quality_scores().to_vec();
    let raw2 = r2.quality_scores().to_vec();
    let has_qual = !raw1.is_empty() && !raw2.is_empty();

    if seq1.len() < 2 || seq2.len() < 2 {
        return Outcome::Ambiguous;
    }

    // r2 reverse-complemented for the overlap layout; its quality must be
    // reversed along with the bases to stay aligned (Read.reverseComplementFast
    // reverses both in BBTools).
    let rc2: Vec<u8> = rev_comp(&seq2).collect();
    let qual2: Vec<u8> = to_phred(&seq2, &raw2).into_iter().rev().collect();

    let min_overlap_entropy =
        min_overlap_from_entropy(&seq1, &rc2, opts.min_entropy).max(opts.min_overlap);
    let min0 = opts.min_overlap0 as isize - opts.ratio_reduction as isize;
    let min = min_overlap_entropy as isize - opts.ratio_reduction as isize;
    // BBMerge sets MAKE_VECTOR=true in main(), forcing maxRatio=0.7 in the
    // ratio pre-screen (BBMergeOverlapper.mateByOverlapRatioJava).
    let max_ratio = if opts.make_vector {
        0.7
    } else {
        opts.max_ratio
    };
    let res = overlap::mate_by_overlap_ratio(
        &seq1,
        &rc2,
        min0.max(0) as usize,
        min.max(0) as usize,
        min_insert0,
        min_insert,
        max_ratio,
        opts.min_second_ratio,
        opts.ratio_margin,
        opts.ratio_offset,
        0.95,
        0.95,
    );
    let mut best_insert = res.insert;
    let best_bad = res.bad;

    if opts.make_vector {
        // With MAKE_VECTOR=true BBMerge skips the ambiguity return, restores
        // the ratio insert even when the selection fell through, and guards
        // the pfilter rejection behind !MAKE_VECTOR, so every positive ratio
        // insert is accepted -- unless the overlap-filter net rejects it.
        if best_insert > 0 {
            let qual1 = to_phred(&seq1, &raw1);
            let mut v = Vec::with_capacity(23);
            v.push(min_overlap_entropy as f32 * 0.1);
            let max_bases = (seq1.len().max(seq2.len())).min(seq1.len() + seq2.len() - min_insert);
            v.push(expected_tip_errors(&seq1, &qual1, max_bases));
            v.push(expected_tip_errors(&rc2, &qual2, max_bases));
            v.push((seq1.len() as f32 - 100.0) * 0.01);
            v.push((seq2.len() as f32 - 100.0) * 0.01);
            v.push(res.insert as f32 * 0.004);
            v.push(res.best_overlap as f32 / (res.best_overlap as f32 + 50.0));
            v.push((res.best_bad + 1.0) / (res.best_bad + 5.0));
            v.push((res.best_good + 1.0) / (res.best_good + 5.0));
            v.push(res.best_ratio);
            v.push((res.bad as f32 + 1.0) / (res.bad as f32 + 5.0));
            v.push(res.second_insert as f32 * 0.004);
            v.push(res.second_overlap as f32 / (res.second_overlap as f32 + 50.0));
            v.push((res.second_bad + 1.0) / (res.second_bad + 5.0));
            v.push((res.second_good + 1.0) / (res.second_good + 5.0));
            v.push(res.second_ratio);
            v.push(res.second_bad_int as f32 / (res.second_bad_int as f32 + 5.0));
            v.push((res.second_ratio + 1.0) / (res.best_ratio + 1.0));
            v.push(res.second_bad / (res.best_bad + 8.0));
            v.push(res.second_good / (res.best_good + 8.0));
            v.push(
                (res.best_overlap as f32 + 1.0)
                    / (res.second_overlap as f32 + res.best_overlap as f32 + 1.0),
            );
            v.push(expected_mismatches(
                &seq1,
                &rc2,
                &qual1,
                &qual2,
                best_insert,
            ));
            v.push(probability(&seq1, &rc2, &qual1, &qual2, best_insert).sqrt() + 0.0000015);
            debug_assert_eq!(v.len(), 23);
            let net = opts.net.as_ref().expect("net required in make-vector mode");
            let score = net.feed_forward(&v);
            if score < net.cutoff {
                stats.no_solution += 1;
                return Outcome::NoSolution;
            }
        }
    } else {
        let mut ambig = res.ambig;
        // BBMerge trips ambig on a bad mismatch count before the filters.
        if best_bad > opts.max_bad as i32 {
            ambig = true;
        }
        if !ambig && best_insert > 0 && has_qual {
            // BBMerge converts ASCII quality to phred (applyQualOffset, -33)
            // at parse time; the probability filter works on phred values.
            let qual1 = to_phred(&seq1, &raw1);
            if opts.pfilter > 0.0 {
                let prob = probability(&seq1, &rc2, &qual1, &qual2, best_insert);
                if prob < opts.pfilter {
                    best_insert = -1;
                }
            }
        }
        // BBMerge routes a failed ratio scan through the "else" branch with
        // bestBad=99999, which trips the MAX_MISMATCHES_R check -> RET_AMBIG.
        if ambig || res.insert <= 0 {
            stats.ambiguous += 1;
            return Outcome::Ambiguous;
        }
    }
    if best_insert <= 0 {
        stats.no_solution += 1;
        return Outcome::NoSolution;
    }
    let insert = best_insert as usize;
    if insert < min_insert {
        stats.too_short += 1;
        return Outcome::TooShort(insert);
    }

    stats.joined += 1;
    stats.insert_min = stats.insert_min.min(insert);
    stats.insert_max = stats.insert_max.max(insert);
    stats.hist[insert.min(HIST_LEN - 1)] += 1;

    if opts.ecco {
        let qual1 = to_phred(&seq1, &raw1);
        let (c1, c2) = corrected_pair(&seq1, &rc2, &qual1, &qual2, insert);
        let mut errors = 0usize;
        for (i, &b) in c1.0.iter().enumerate().take(seq1.len()) {
            if seq1[i] != b && is_fully_defined(b) {
                errors += 1;
            }
        }
        let c2len = c2.0.len();
        for i in 0..c2len {
            let j = rc2.len() - c2len + i;
            if rc2[j] != c2.0[i] && is_fully_defined(c2.0[i]) {
                errors += 1;
            }
        }
        stats.errors_corrected += errors as u64;
        r1.set_sequence(c1.0);
        r1.set_quality(from_phred(&c1.1));
        // r2 back to original orientation.
        let rc: Vec<u8> = rev_comp(&c2.0).collect();
        r2.set_sequence(rc);
        let q: Vec<u8> = from_phred(&c2.1).into_iter().rev().collect();
        r2.set_quality(q);
        return Outcome::Corrected(insert);
    }

    let qual1 = to_phred(&seq1, &raw1);
    let (jbases, jqual) = join_reads(&seq1, &rc2, &qual1, &qual2, insert);
    Outcome::Merged(jbases, from_phred(&jqual), insert)
}

/// `calcMinOverlapByEntropyTail`/`Head` merged into one call: the tail of
/// `r1` (3' end) and the head of the reverse-complemented `r2`.
fn min_overlap_from_entropy(r1: &[u8], r2_rc: &[u8], minscore: usize) -> usize {
    calc_min_overlap_by_entropy_tail(r1, minscore)
        .max(calc_min_overlap_by_entropy_head(r2_rc, minscore))
}

/// `calcMinOverlapByEntropyTail`: walks from the 3' end, scoring unique and
/// twice-seen 3-mers until `ones*4+twos >= minscore`.
fn calc_min_overlap_by_entropy_tail(bases: &[u8], minscore: usize) -> usize {
    let mut counts = [0u16; 64];
    let mut kmer = 0usize;
    let mut len = 0usize;
    let mut ones = 0usize;
    let mut twos = 0usize;
    for (i, &b) in bases.iter().rev().enumerate() {
        let Some(n) = base_to_number(b).map(|x| x as usize) else {
            len = 0;
            kmer = 0;
            continue;
        };
        len += 1;
        kmer = ((kmer << 2) | n) & 0x3f;
        if len >= 3 {
            let c = &mut counts[kmer];
            *c += 1;
            if *c == 1 {
                ones += 1;
            } else if *c == 2 {
                twos += 1;
            }
            if ones * 4 + twos >= minscore {
                return i;
            }
        }
    }
    bases.len() + 1
}

/// `calcMinOverlapByEntropyHead`: walks from the 5' end.
fn calc_min_overlap_by_entropy_head(bases: &[u8], minscore: usize) -> usize {
    let mut counts = [0u16; 64];
    let mut kmer = 0usize;
    let mut len = 0usize;
    let mut ones = 0usize;
    let mut twos = 0usize;
    for (i, &b) in bases.iter().enumerate() {
        let Some(n) = base_to_number(b).map(|x| x as usize) else {
            len = 0;
            kmer = 0;
            continue;
        };
        len += 1;
        kmer = ((kmer << 2) | n) & 0x3f;
        if len >= 3 {
            let c = &mut counts[kmer];
            *c += 1;
            if *c == 1 {
                ones += 1;
            } else if *c == 2 {
                twos += 1;
            }
            if ones * 4 + twos >= minscore {
                return i;
            }
        }
    }
    bases.len() + 1
}

/// `probCorrect4`: quality ASCII value to base-correctness probability.
const PROB_CORRECT4: [f32; 60] = [
    0.0000, 0.2501, 0.3690, 0.4988, 0.6019, 0.6838, 0.7488, 0.8005, 0.8415, 0.8741, 0.9000, 0.9206,
    0.9369, 0.9499, 0.9602, 0.9684, 0.9749, 0.9800, 0.9842, 0.9874, 0.9900, 0.9921, 0.9937, 0.9950,
    0.9960, 0.9968, 0.9975, 0.9980, 0.9984, 0.9987, 0.9990, 0.9992, 0.9994, 0.9995, 0.9996, 0.9997,
    0.9997, 0.9998, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,
    0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,
];

/// `BBMergeOverlapper.probability`.
fn probability(a: &[u8], b: &[u8], aqual: &[u8], bqual: &[u8], insert: i32) -> f32 {
    let alen = a.len() as i32;
    let blen = b.len() as i32;
    if aqual.is_empty() || bqual.is_empty() {
        return 1.0;
    }
    let istart = if insert <= blen { 0 } else { insert - blen };
    let jstart = if insert >= blen { 0 } else { blen - insert };
    let mut prob_actual = 1f32;
    let mut prob_common = 1f32;
    let mut i = istart;
    let mut j = jstart;
    while i < insert && i < alen && j < blen {
        let ca = a[i as usize];
        let cb = b[j as usize];
        if ca != b'N' && cb != b'N' {
            // Java indexes with raw phred; clamp so pgr never panics on
            // qualities above the table (Java would crash on such input).
            let qi = (aqual[i as usize] as usize).min(PROB_CORRECT4.len() - 1);
            let qj = (bqual[j as usize] as usize).min(PROB_CORRECT4.len() - 1);
            let prob_c = PROB_CORRECT4[qi] * PROB_CORRECT4[qj];
            let prob_m = prob_c + (1.0 - prob_c) * 0.25;
            let prob_e = 1.0 - prob_m;
            prob_common *= prob_m.max(prob_e);
            prob_actual *= if ca == cb { prob_m } else { prob_e };
        }
        i += 1;
        j += 1;
    }
    (prob_actual / prob_common).sqrt()
}

/// `BBMergeOverlapper.expectedMismatches` (faithful, including its handling
/// of unequal read lengths).
fn expected_mismatches(a: &[u8], b: &[u8], aqual: &[u8], bqual: &[u8], insert: i32) -> f32 {
    let alen = a.len() as i32;
    let blen = b.len() as i32;
    if aqual.is_empty() || bqual.is_empty() {
        return (insert as f32 + 0.0) / 16.0;
    }
    let istart = if insert <= blen { 0 } else { insert - blen };
    let jstart = if insert <= alen { alen - insert } else { 0 };
    let mut expected = 0f32;
    let mut i = istart;
    let mut j = jstart;
    while i < insert && i < alen && j < blen {
        let ca = a[i as usize];
        let cb = b[j as usize];
        if ca != b'N' && cb != b'N' {
            let qi = (aqual[i as usize] as usize).min(PROB_CORRECT4.len() - 1);
            let qj = (bqual[j as usize] as usize).min(PROB_CORRECT4.len() - 1);
            let prob_c = PROB_CORRECT4[qi] * PROB_CORRECT4[qj];
            expected += 1.0 - prob_c;
        }
        i += 1;
        j += 1;
    }
    expected
}

/// `Read.expectedTipErrors`: expected errors in the 3' tail.
fn expected_tip_errors(bases: &[u8], quals: &[u8], max_bases: usize) -> f32 {
    if quals.is_empty() {
        return 0.0;
    }
    let limit0 = max_bases.min(quals.len());
    let limit = quals.len() - limit0;
    let mut sum = 0f32;
    let mut i = quals.len() as isize - 1;
    while i >= limit as isize {
        let b = bases[i as usize];
        let q = quals[i as usize];
        if is_fully_defined(b) {
            sum += prob_error(q);
        }
        i -= 1;
    }
    sum
}

/// `QualityTools.PROB_ERROR`: phred quality to error probability.
fn prob_error(q: u8) -> f32 {
    match q {
        0 => 0.75,
        1 => 0.7,
        _ => 10f32.powf(-0.1 * q as f32),
    }
}

/// `Read.joinRead`: builds the merged read by walking the overlap from the
/// 3' end. `b` is reverse-complemented.
fn join_reads(a: &[u8], b: &[u8], aqual: &[u8], bqual: &[u8], insert: usize) -> (Vec<u8>, Vec<u8>) {
    let length_sum = a.len() + b.len();
    let overlap = insert.min(length_sum.saturating_sub(insert));
    let has_qual = !aqual.is_empty() && !bqual.is_empty();
    let mut bases = vec![0u8; insert];
    let mut quals = vec![0u8; insert];

    if overlap == 0 {
        // Simple join with an N gap.
        let lim = insert.saturating_sub(b.len());
        bases[..a.len().min(insert)].copy_from_slice(&a[..a.len().min(insert)]);
        for x in bases[a.len().min(insert)..lim.min(insert)].iter_mut() {
            *x = b'N';
        }
        let bstart = lim;
        if bstart < insert {
            let n = (insert - bstart).min(b.len());
            bases[bstart..bstart + n].copy_from_slice(&b[..n]);
        }
        if has_qual {
            quals[..a.len().min(insert)].copy_from_slice(&aqual[..a.len().min(insert)]);
            if bstart < insert {
                let n = (insert - bstart).min(b.len());
                quals[bstart..bstart + n].copy_from_slice(&bqual[..n]);
            }
        }
    } else {
        let n = a.len().min(insert);
        bases[..n].copy_from_slice(&a[..n]);
        if has_qual {
            quals[..n].copy_from_slice(&aqual[..n]);
        }
        let mut i = insert as isize - 1;
        let mut j = b.len() as isize - 1;
        while i >= 0 && j >= 0 {
            let (ii, jj) = (i as usize, j as usize);
            let ca = bases[ii];
            let cb = b[jj];
            if has_qual {
                let qa = quals[ii];
                let qb = bqual[jj];
                if ca == 0 || ca == b'N' {
                    bases[ii] = cb;
                    quals[ii] = qb;
                } else if cb == 0 || cb == b'N' {
                    // keep a's base
                } else if ca == cb {
                    let q = (qa.max(qb) + qa.min(qb) / 4).min(50);
                    quals[ii] = q;
                } else {
                    bases[ii] = if qa > qb {
                        ca
                    } else if qa < qb {
                        cb
                    } else {
                        b'N'
                    };
                    quals[ii] = qa.max(qb) - qa.min(qb);
                }
            } else if ca == 0 || ca == b'N' {
                bases[ii] = cb;
            } else if cb != 0 && cb != b'N' && ca != cb {
                bases[ii] = ca.max(cb);
            }
            i -= 1;
            j -= 1;
        }
    }
    if has_qual {
        (bases, quals)
    } else {
        (bases, Vec::new())
    }
}

/// Returns the corrected pair (r1 consensus, r2-rc consensus).
type CorrectedPair = ((Vec<u8>, Vec<u8>), (Vec<u8>, Vec<u8>));

fn corrected_pair(a: &[u8], b: &[u8], aqual: &[u8], bqual: &[u8], insert: usize) -> CorrectedPair {
    let (jbases, jqual) = join_reads(a, b, aqual, bqual, insert);
    let lenj = jbases.len();
    let lim1 = lenj.min(a.len());
    let lim2 = lenj - lenj.min(b.len());
    let mut c1 = (jbases[..lim1].to_vec(), Vec::new());
    if !jqual.is_empty() {
        c1.1 = jqual[..lim1].to_vec();
    }
    let mut c2 = (jbases[lim2..].to_vec(), Vec::new());
    if !jqual.is_empty() {
        c2.1 = jqual[lim2..].to_vec();
    }
    (c1, c2)
}

fn is_fully_defined(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Runs overlap merge / ecco correction over a FASTQ input (1 interleaved
/// file or 2 paired files).
pub fn merge<W: Write>(
    infiles: &[String],
    out: &mut W,
    mut outu: Option<&mut W>,
    opts: &MergeOptions,
) -> Result<MergeStats> {
    if opts.make_vector && opts.net.is_none() {
        anyhow::bail!(
            "make-vector mode requires a BBMerge overlap net (bbmerge.bbnet); \
             pass --net or use --no-make-vector"
        );
    }
    let mut stats = MergeStats {
        insert_min: usize::MAX,
        ..MergeStats::default()
    };
    let mut reader1 = SeqReader::new(&infiles[0])?;
    let mut reader2 = if infiles.len() > 1 {
        Some(SeqReader::new(&infiles[1])?)
    } else {
        None
    };
    let mut rec1 = SeqRecord::new();
    let mut rec2 = SeqRecord::new();

    // `extend2` (bbmerge-auto / tadpole mode) builds a k-mer table from the
    // input reads and extends unmerged pairs, mirroring BBMerge's
    // `extendAndMerge` retry. BBMerge forces MAKE_VECTOR=false in this mode.
    let table = if opts.extend2 > 0 {
        if opts.make_vector {
            anyhow::bail!("--extend2 requires --no-make-vector (bbmerge-auto forces it)");
        }
        let reads: Vec<(Vec<u8>, Vec<u8>)> = {
            let mut r1 = SeqRecord::new();
            let mut r2 = SeqRecord::new();
            let mut reader1 = SeqReader::new(&infiles[0])?;
            let mut reader2 = if infiles.len() > 1 {
                Some(SeqReader::new(&infiles[1])?)
            } else {
                None
            };
            let mut out = Vec::new();
            loop {
                if !reader1.read_record(&mut r1)? {
                    break;
                }
                canonicalize_quality(&mut r1);
                out.push((
                    r1.sequence().to_vec(),
                    to_phred(r1.sequence(), r1.quality_scores()),
                ));
                let has2 = match reader2.as_mut() {
                    Some(r) => r.read_record(&mut r2)?,
                    None => reader1.read_record(&mut r2)?,
                };
                if !has2 {
                    anyhow::bail!("unpaired trailing read: {}", r1.name());
                }
                canonicalize_quality(&mut r2);
                out.push((
                    r2.sequence().to_vec(),
                    to_phred(r2.sequence(), r2.quality_scores()),
                ));
            }
            out
        };
        let t = TadpoleOptions {
            k: 81,
            ..TadpoleOptions::default()
        };
        Some(TadpoleTable::build(&reads, t.k, t.min_prob))
    } else {
        None
    };

    loop {
        if !reader1.read_record(&mut rec1)? {
            break;
        }
        let has2 = match reader2.as_mut() {
            Some(r) => r.read_record(&mut rec2)?,
            None => reader1.read_record(&mut rec2)?,
        };
        if !has2 {
            anyhow::bail!("unpaired trailing read: {}", rec1.name());
        }
        stats.pairs += 1;

        // BBTools converts quality to phred at parse time (applyQualOffset:
        // non-ACGT -> 0, ACGT -> max(2, raw-33)) and writes back +33, so the
        // output quality of every record (merged or not) is canonicalized.
        canonicalize_quality(&mut rec1);
        canonicalize_quality(&mut rec2);
        // BBMerge snapshots the pre-extension reads and restores them when a
        // pair is written unmerged (`originals` in `findOverlapInThread`).
        let orig1 = (rec1.sequence().to_vec(), rec1.quality_scores().to_vec());
        let orig2 = (rec2.sequence().to_vec(), rec2.quality_scores().to_vec());

        // BBMerge counts each pair exactly once, using the final
        // `processReadPair` result. The extend-retry below re-runs the
        // overlap detection, so snapshot the stats and only keep the retry's
        // contribution when extension is attempted.
        let stats_before = stats.clone();
        let mut outcome = process_pair(&mut rec1, &mut rec2, opts, &mut stats);

        // bbmerge-auto `extend2=N rem`: BBMerge extends every pair when
        // `rem` is set (`requireExtensionMatch || AMBIG || NO_SOLUTION`),
        // then re-checks the overlap, mirroring `processReadPair` ->
        // `extendAndMerge`.
        if opts.extend2 > 0
            && (opts.rem || matches!(outcome, Outcome::Ambiguous | Outcome::NoSolution))
        {
            stats.clone_from(&stats_before);
            let original_insert = outcome.insert();
            let t = TadpoleOptions {
                k: 81,
                ..TadpoleOptions::default()
            };
            let table = table.as_ref().unwrap();
            // BBMerge computes `lengthSum` from the *unextended* reads before
            // the extension block, and the rem acceptance rule uses it.
            let pre_len_sum = rec1.sequence().len() + rec2.sequence().len();
            // `extendAndMerge`: BBMerge runs `extendIterations` passes
            // (default 1) over the reads; each pass extends by `extend2` and
            // re-checks the overlap.
            let mut e1 = 0usize;
            let mut e2 = 0usize;
            let mut iterations = 0usize;
            loop {
                if iterations >= 1 {
                    break;
                }
                iterations += 1;
                let mut b1 = rec1.sequence().to_vec();
                let mut q1 = to_phred(&b1, rec1.quality_scores());
                let mut b2 = rec2.sequence().to_vec();
                let mut q2 = to_phred(&b2, rec2.quality_scores());
                e1 = extend_read_right(&mut b1, &mut q1, table, &t, opts.extend2);
                e2 = extend_read_right(&mut b2, &mut q2, table, &t, opts.extend2);
                rec1.set_sequence(b1);
                rec1.set_quality(from_phred(&q1));
                rec2.set_sequence(b2);
                rec2.set_quality(from_phred(&q2));
                // BBMerge always re-checks the overlap after the extension
                // pass (even when neither read extended), and that final
                // result is what gets counted in the histogram.
                outcome = process_pair(&mut rec1, &mut rec2, opts, &mut stats);
                if !matches!(outcome, Outcome::Ambiguous | Outcome::NoSolution) {
                    break;
                }
            }

            // `rem` acceptance rule: an extended overlap must match the
            // unextended one, unless it is a genuinely new overlap too long
            // to have been detectable before extension.
            if opts.rem && (e1 > 0 || e2 > 0) {
                let ext_insert = outcome.insert();
                if ext_insert != original_insert {
                    if original_insert.is_none() {
                        if let Some(insert) = ext_insert {
                            let approx_max = pre_len_sum.saturating_sub(26);
                            if insert <= approx_max || e1 + e2 < 12 {
                                outcome = Outcome::Ambiguous;
                            }
                        }
                    } else {
                        outcome = Outcome::Ambiguous;
                    }
                }
            }
        }

        match outcome {
            Outcome::Merged(bases, quals, _) => {
                write_record(out, rec1.name(), rec1.comment(), &bases, &quals)?;
            }
            Outcome::Corrected(_) => {
                write_record(
                    out,
                    rec1.name(),
                    rec1.comment(),
                    rec1.sequence(),
                    rec1.quality_scores(),
                )?;
                write_record(
                    out,
                    rec2.name(),
                    rec2.comment(),
                    rec2.sequence(),
                    rec2.quality_scores(),
                )?;
            }
            _ => {
                if opts.ecco && opts.mix {
                    write_record(
                        out,
                        rec1.name(),
                        rec1.comment(),
                        rec1.sequence(),
                        rec1.quality_scores(),
                    )?;
                    write_record(
                        out,
                        rec2.name(),
                        rec2.comment(),
                        rec2.sequence(),
                        rec2.quality_scores(),
                    )?;
                } else if let Some(w) = outu.as_deref_mut() {
                    // Write the pre-extension reads, like BBMerge.
                    rec1.set_sequence(orig1.0.clone());
                    rec1.set_quality(orig1.1.clone());
                    rec2.set_sequence(orig2.0.clone());
                    rec2.set_quality(orig2.1.clone());
                    write_record(
                        w,
                        rec1.name(),
                        rec1.comment(),
                        rec1.sequence(),
                        rec1.quality_scores(),
                    )?;
                    write_record(
                        w,
                        rec2.name(),
                        rec2.comment(),
                        rec2.sequence(),
                        rec2.quality_scores(),
                    )?;
                } else if opts.mix {
                    write_record(
                        out,
                        rec1.name(),
                        rec1.comment(),
                        rec1.sequence(),
                        rec1.quality_scores(),
                    )?;
                }
            }
        }
    }
    Ok(stats)
}

/// Applies the BBTools phred round-trip to a record's quality scores.
fn canonicalize_quality(rec: &mut SeqRecord) {
    if rec.quality_scores().is_empty() {
        return;
    }
    let seq = rec.sequence().to_vec();
    let raw = rec.quality_scores().to_vec();
    let phred = to_phred(&seq, &raw);
    rec.set_quality(from_phred(&phred));
}

/// Writes one record, keeping the header/comment and format of the input.
fn write_record<W: Write>(
    w: &mut W,
    name: &bstr::BStr,
    comment: &bstr::BStr,
    seq: &[u8],
    qual: &[u8],
) -> Result<()> {
    let header = if comment.is_empty() {
        name.to_string()
    } else {
        format!("{} {}", name, comment)
    };
    if qual.is_empty() {
        crate::libs::fmt::fq::write_fa(w, &header, seq)?;
    } else {
        crate::libs::fmt::fq::write_fq(w, &header, seq, qual)?;
    }
    Ok(())
}

/// Writes the BBMerge `ihist` format.
pub fn write_ihist<W: Write>(w: &mut W, stats: &MergeStats) -> Result<()> {
    let sum: u64 = stats.hist.iter().sum();
    let sum = sum.max(1);
    let sum2: u64 = stats
        .hist
        .iter()
        .enumerate()
        .map(|(i, &c)| i as u64 * c)
        .sum();
    let mean = sum2 as f64 / sum as f64;
    let median = percentile_histogram(&stats.hist, 0.5);
    let mode = calc_mode_histogram(&stats.hist);
    let stdev = {
        let sumdev2: f64 = stats
            .hist
            .iter()
            .enumerate()
            .map(|(i, &c)| c as f64 * (mean - i as f64).powi(2))
            .sum();
        (sumdev2 / sum as f64).sqrt()
    };
    let percent = if stats.pairs > 0 {
        stats.joined as f64 * 100.0 / stats.pairs as f64
    } else {
        0.0
    };
    writeln!(w, "#Mean\t{}", fmt_3f(mean))?;
    writeln!(w, "#Median\t{}", median)?;
    writeln!(w, "#Mode\t{}", mode)?;
    writeln!(w, "#STDev\t{}", fmt_3f(stdev))?;
    writeln!(w, "#PercentOfPairs\t{}", fmt_3f(percent))?;
    writeln!(w, "#InsertSize\tCount")?;
    for (i, &c) in stats.hist.iter().enumerate() {
        if c > 0 && i <= stats.insert_max {
            writeln!(w, "{}\t{}", i, c)?;
        }
    }
    Ok(())
}

fn percentile_histogram(hist: &[u64; HIST_LEN], fraction: f64) -> usize {
    let sum: u64 = hist.iter().sum();
    let target = (sum as f64 * fraction) as u64;
    let mut acc = 0u64;
    for (i, &c) in hist.iter().enumerate() {
        acc += c;
        if acc >= target {
            return i;
        }
    }
    HIST_LEN - 1
}

fn calc_mode_histogram(hist: &[u64; HIST_LEN]) -> usize {
    let median = percentile_histogram(hist, 0.5);
    let mut mode = 0usize;
    let mut mode_count = hist[0];
    for (i, &c) in hist.iter().enumerate().skip(1) {
        if c > mode_count || (c == mode_count && i.abs_diff(median) < mode.abs_diff(median)) {
            mode = i;
            mode_count = c;
        }
    }
    mode
}

/// Java `String.format("%.3f")`-equivalent (half-up) formatting.
fn fmt_3f(x: f64) -> String {
    let scaled = (x * 1000.0).round() as i64;
    format!("{}.{:03}", scaled / 1000, (scaled % 1000).abs())
}

#[cfg(test)]
mod tests_hist {
    use super::*;

    #[test]
    fn entropy_returns_high_for_low_complexity() {
        // Poly-A tail: few unique 3-mers -> never reaches a high score.
        let poly = b"AAAAAAAAAAAAAAAAAAAA";
        let r = calc_min_overlap_by_entropy_tail(poly, 39);
        assert_eq!(r, poly.len() + 1);
    }

    #[test]
    fn join_consensus_matches_expected() {
        let a = b"ACGTACGTACGT";
        let b = b"ACGTACGTACGT"; // rc of itself
        let qa = b"IIIIIIIIIIII";
        let qb = b"IIIIIIIIIIII";
        let (bases, quals) = join_reads(a, b, qa, qb, 12);
        assert_eq!(&bases, b"ACGTACGTACGT");
        assert_eq!(quals.len(), 12);
    }

    #[test]
    fn merge_stats_default() {
        let s = MergeStats {
            insert_min: usize::MAX,
            ..MergeStats::default()
        };
        assert_eq!(s.pairs, 0);
        assert_eq!(s.insert_min, usize::MAX);
    }

    #[test]
    fn options_resolve_like_bbmerge() {
        let o = MergeOptions::from_preset(Preset::Strict);
        let (min_insert, min_insert0) = o.resolved();
        assert_eq!(min_insert, 15);
        assert_eq!(min_insert0, 12);
        let v = MergeOptions::from_preset(Preset::VStrict);
        let (mi, mi0) = v.resolved();
        assert_eq!(mi, 15);
        assert_eq!(mi0, 12);
    }
}
