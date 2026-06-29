use intspan::{IntSpan, Range};

use super::intspan_util::indel_intspan;

/// Trims pure dash regions
///
/// ```
/// let mut seqs = vec![
///     "AAAATTTTTG".to_string(),
///     "AAAATTTTTG".to_string(),
///     "AAAATTTTTG".to_string(),
/// ];
/// pgr::libs::alignment::trim_pure_dash(&mut seqs);
/// assert_eq!(seqs[0].len(), 10);
///
/// let mut seqs = vec![
///     "-AA--TTTGG".to_string(),
///     "-AA--TTTGG".to_string(),
///     "-AA--TTTGG".to_string(),
/// ];
/// pgr::libs::alignment::trim_pure_dash(&mut seqs);
/// assert_eq!(seqs[0].len(), 7);
///
/// let mut seqs = vec![
///     "-AA--TTTGG".to_string(),
///     "-AAA-TTTGG".to_string(),
///     "AAA--TTTTG".to_string(),
/// ];
/// pgr::libs::alignment::trim_pure_dash(&mut seqs);
/// assert_eq!(seqs[0].len(), 9);
///
/// ```
pub fn trim_pure_dash(seqs: &mut [String]) {
    let mut trim_region = IntSpan::new();
    let seq_count = seqs.len();

    for (i, seq) in seqs.iter().enumerate() {
        let ints = indel_intspan(seq.as_bytes().to_vec().as_ref());
        if i == 0 {
            trim_region.merge(&ints);
        } else {
            trim_region = trim_region.intersect(&ints);
        }
    }

    // eprintln!("trim_region = {:#?}", trim_region.to_string());

    // trim all segments in trim_region
    for (lower, upper) in trim_region.spans().iter().rev() {
        for seq in seqs.iter_mut().take(seq_count) {
            seq.replace_range((*lower as usize - 1)..*upper as usize, "");
        }
    }
}

pub(super) fn align_indel_ints(seqs: &mut [String], count: usize) -> (IntSpan, IntSpan) {
    let mut union_ints = IntSpan::new();
    let mut intersect_ints = IntSpan::new();

    for (i, seq) in seqs.iter().enumerate().take(count) {
        let ints = indel_intspan(seq.as_bytes().to_vec().as_ref());

        if i == 0 {
            union_ints.merge(&ints);
            intersect_ints.merge(&ints);
        } else {
            union_ints = union_ints.union(&ints);
            intersect_ints = intersect_ints.intersect(&ints);
        }
    }

    (union_ints, intersect_ints)
}

/// Trims outgroup-only regions
///
/// Iff. intersect is superset of union
///     T G----C
///     Q G----C
///     O GAAAAC
///
/// ```
/// let mut seqs = vec![
///     "AAAATTTTTG".to_string(),
///     "AAAATTTTTG".to_string(),
///     "AAAATTTTTG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// assert_eq!(seqs[0].len(), 10);
///
/// let mut seqs = vec![
///     "-AA--TTTGG".to_string(),
///     "-AA--TTTGG".to_string(),
///     "-AA--TTTGG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// assert_eq!(seqs[0].len(), 7);
///
/// let mut seqs = vec![
///     "-AA--TTTGG".to_string(),
///     "-AAA-TTTGG".to_string(),
///     "AAA--TTTTG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// assert_eq!(seqs[0].len(), 9);
///
/// let mut seqs = vec![
///     "AAA--TT-GG".to_string(),
///     "AAAATTT-GG".to_string(),
///     "AAA--TTTTG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// assert_eq!(seqs[0].len(), 9);
///
/// let mut seqs = vec![
///     "-AA--TT-GG".to_string(),
///     "-AAA-TT-GG".to_string(),
///     "AAA--TTTTG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// assert_eq!(seqs[0].len(), 8);
///
/// ```
pub fn trim_outgroup(seqs: &mut [String]) {
    let seq_count = seqs.len();

    assert!(seq_count >= 3, "Need three or more sequences");

    // Last seq is the outgroup
    let (union_ints, intersect_ints) = align_indel_ints(seqs, seq_count - 1);

    // find trim_region
    let mut trim_region = IntSpan::new();
    for (lower, upper) in union_ints.spans() {
        let ints = IntSpan::from_pair(lower, upper);
        if intersect_ints.superset(&ints) {
            trim_region.merge(&ints);
        }
    }

    // trim all segments in trim_region
    for (lower, upper) in trim_region.spans().iter().rev() {
        for seq in seqs.iter_mut().take(seq_count) {
            seq.replace_range((*lower as usize - 1)..*upper as usize, "");
        }
    }
}

/// Records complex ingroup indels (ingroup-outgroup complex indels are not identified here)
///
/// After trim_outgroup(), All ingroup intersect ints are parts of complex indels
/// intersect 4-5, union 2-5
///     T GGA--C
///     Q G----C
///     O GGAGAC
/// result, complex_region 2-3
///     T GGAC
///     Q G--C
///     O GGAC
///
/// ```
/// let mut seqs = vec![
///     "AAAATTTTTG".to_string(),
///     "AAAATTTTTG".to_string(),
///     "AAAATTTTTG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// let complex = pgr::libs::alignment::trim_complex_indel(&mut seqs);
/// assert_eq!(seqs[0].len(), 10);
/// assert_eq!(complex.to_string(), "-");
///
/// let mut seqs = vec![
///     "-AA--TTTGG".to_string(),
///     "-AA--TTTGG".to_string(),
///     "-AA--TTTGG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// let complex = pgr::libs::alignment::trim_complex_indel(&mut seqs);
/// assert_eq!(seqs[0].len(), 7);
/// assert_eq!(complex.to_string(), "-");
///
/// let mut seqs = vec![
///     "-AA--TTTGG".to_string(),
///     "-AAA-TTTGG".to_string(),
///     "AAA--TTTGG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// let complex = pgr::libs::alignment::trim_complex_indel(&mut seqs);
/// assert_eq!(seqs[0].len(), 8);
/// assert_eq!(complex.to_string(), "3");
///
/// let mut seqs = vec![
///     "AAA--TT-GG".to_string(),
///     "AAAATTT-GG".to_string(),
///     "AAA--TTTTG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// let complex = pgr::libs::alignment::trim_complex_indel(&mut seqs);
/// assert_eq!(seqs[0].len(), 9);
/// assert_eq!(complex.to_string(), "-");
///
/// let mut seqs = vec![
///     "-AA--TT-GG".to_string(),
///     "-AAA-TT-GG".to_string(),
///     "AAA--TTTTG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// let complex = pgr::libs::alignment::trim_complex_indel(&mut seqs);
/// assert_eq!(seqs[0].len(), 7);
/// assert_eq!(complex.to_string(), "3");
///
/// let mut seqs = vec![
///     "-AA--TTTGG".to_string(),
///     "-AAA-TT-GG".to_string(),
///     "AAA--TTTTG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// let complex = pgr::libs::alignment::trim_complex_indel(&mut seqs);
/// assert_eq!(seqs[0].len(), 8);
/// assert_eq!(complex.to_string(), "3");
///
/// let mut seqs = vec![
///     "-AA--TTTGG".to_string(),
///     "-AAA-TT--G".to_string(),
///     "AAA--TT--G".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// let complex = pgr::libs::alignment::trim_complex_indel(&mut seqs);
/// assert_eq!(seqs[0].len(), 8);
/// assert_eq!(complex.to_string(), "3");
///
/// let mut seqs = vec![
///     "-AA--TTTGG".to_string(),
///     "-AAA-TT--G".to_string(),
///     "AAA--TT-GG".to_string(),
/// ];
/// pgr::libs::alignment::trim_outgroup(&mut seqs);
/// let complex = pgr::libs::alignment::trim_complex_indel(&mut seqs);
/// assert_eq!(seqs[0].len(), 8);
/// assert_eq!(complex.to_string(), "3");
///
/// ```
pub fn trim_complex_indel(seqs: &mut [String]) -> String {
    let seq_count = seqs.len();

    assert!(seq_count >= 3, "Need three or more sequences");

    // Last seq is the outgroup
    let (union_ints, intersect_ints) = align_indel_ints(seqs, seq_count - 1);

    // find ingroup complex_region
    let mut complex_region = IntSpan::new();
    for (lower, upper) in intersect_ints.spans().iter().rev() {
        let sub_intersect_ints = IntSpan::from_pair(*lower, *upper);

        // trim sequences, including outgroup
        for seq in seqs.iter_mut().take(seq_count) {
            seq.replace_range((*lower as usize - 1)..*upper as usize, "");
        }

        // add to complex_region
        for sub_union_ints in union_ints.intses() {
            if sub_union_ints.superset(&sub_intersect_ints) {
                complex_region.merge(&sub_union_ints);
            }
        }

        // modify all related set
        // intersect_ints is not affected
        // union_ints is affected
        complex_region = complex_region.banish(*lower, *upper);
    }

    complex_region.to_string()
}

/// Trims head and tail indels.
/// Returns a Vecter of Tuple(head, tail) corresponding to the bases deleted from each sequence
///
/// If chop length set to 1, the first indel will be trimmed.
/// Length set to 5 and the second indel will also be trimmed.
/// GAAA--C...
/// --AAAGC...
/// GAAAAGC...
///
/// ```
/// let seqs = vec![
///     "-AA--TTTGGCGCGCGCGCGCGCGCGC".to_string(),
///     "-AAAATT--GCGCGCGCGCGCGCGC-C".to_string(),
///     "AAA--TT-GGCGCGCGCGCGCGCGCGC".to_string(),
/// ];
/// let ranges = vec![
///     intspan::Range::from_str("I(+):101-124"),
///     intspan::Range::from_str("1:1-23"),
///     intspan::Range::from_str("a(-):101-124"),
/// ];
///
/// let mut seqc = seqs.clone();
/// let mut rangec = ranges.clone();
/// pgr::libs::alignment::trim_head_tail(&mut seqc, &mut rangec, 0);
/// assert_eq!(seqc[0].len(), 27);
/// assert_eq!(rangec[0].start, 101);
/// assert_eq!(rangec[1].start, 1);
///
/// let mut seqc = seqs.clone();
/// let mut rangec = ranges.clone();
/// pgr::libs::alignment::trim_head_tail(&mut seqc, &mut rangec, 1); // head 1
/// assert_eq!(seqc[0].len(), 26);
/// assert_eq!(rangec[0].start, 101);
/// assert_eq!(rangec[1].start, 1);
/// assert_eq!(rangec[2].start, 101);
/// assert_eq!(rangec[2].end, 123);
///
/// let mut seqc = seqs.clone();
/// let mut rangec = ranges.clone();
/// pgr::libs::alignment::trim_head_tail(&mut seqc, &mut rangec, 2); // head 1, tail 2
/// assert_eq!(seqc[0].len(), 24);
/// assert_eq!(rangec[0].start, 101);
/// assert_eq!(rangec[0].end, 122);
/// assert_eq!(rangec[1].start, 1);
/// assert_eq!(rangec[1].end, 22);
/// assert_eq!(rangec[2].start, 103);
/// assert_eq!(rangec[2].end, 123);
///
/// let mut seqc = seqs.clone();
/// let mut rangec = ranges.clone();
/// pgr::libs::alignment::trim_head_tail(&mut seqc, &mut rangec, 4); // head 5, tail 2
/// assert_eq!(seqc[0].len(), 20);
/// assert_eq!(rangec[0].start, 103);
/// assert_eq!(rangec[0].end, 122);
/// assert_eq!(rangec[1].start, 5);
/// assert_eq!(rangec[1].end, 22);
/// assert_eq!(rangec[2].start, 103);
/// assert_eq!(rangec[2].end, 121);
/// ```
///
/// ```
/// let seqs = vec![
///     "-AA--TTTGGCATGCATG123456789".to_string(),
///     "-AAAATT--GCATGCATG1234567-9".to_string(),
///     "AAA--TT-GGCATGCATG123456789".to_string(),
///     "AAA--TT-GGCATGCATG1234567--".to_string(),
/// ];
/// let ranges = vec![
///     intspan::Range::from_str("I(+):101-124"),
///     intspan::Range::from_str("1:1-23"),
///     intspan::Range::from_str("a(-):101-124"),
///     intspan::Range::from_str("b(-):1-22"),
/// ];
///
/// let mut seqc = seqs.clone();
/// let mut rangec = ranges.clone();
/// pgr::libs::alignment::trim_head_tail(&mut seqc, &mut rangec, 4); // head 5, tail 2
/// assert_eq!(seqc[0].len(), 20);
/// assert_eq!(seqc[0], "TTTGGCATGCATG1234567".to_string());
/// assert_eq!(rangec[0].start, 103);
/// assert_eq!(rangec[0].end, 122);
/// assert_eq!(seqc[1], "TT--GCATGCATG1234567".to_string());
/// assert_eq!(rangec[1].start, 5);
/// assert_eq!(rangec[1].end, 22);
/// assert_eq!(seqc[2], "TT-GGCATGCATG1234567".to_string());
/// assert_eq!(rangec[2].start, 103);
/// assert_eq!(rangec[2].end, 121);
/// assert_eq!(seqc[3], "TT-GGCATGCATG1234567".to_string());
/// assert_eq!(rangec[3].start, 1, "negative strand");
/// assert_eq!(rangec[3].end, 19);
///
/// ```
pub fn trim_head_tail(seqs: &mut [String], ranges: &mut [Range], chop: usize) {
    let seq_count = seqs.len();

    if chop == 0 {
        return;
    }

    // chop region covers all
    let align_len = seqs.first().unwrap().len();
    if chop * 2 >= align_len {
        return;
    }

    // include all seqs
    let (indel_ints, _) = align_indel_ints(seqs, seq_count);

    // There're no indels at all
    if indel_ints.is_empty() {
        return;
    }

    // head indels to be trimmed
    {
        let head_ints = IntSpan::from_pair(1, chop as i32);
        let head_indel_ints = indel_ints.find_islands_ints(&head_ints);

        if !head_indel_ints.is_empty() {
            for _ in 1..=(head_indel_ints.max() as usize) {
                for i in 0..seq_count {
                    let base = seqs[i].remove(0);
                    if base != '-' {
                        if ranges[i].strand == "+" || ranges[i].strand.is_empty() {
                            ranges[i].start += 1;
                        } else {
                            ranges[i].end -= 1;
                        }
                    }
                }
            }
        }
    }

    // tail indels to be trimmed
    {
        let tail_ints = IntSpan::from_pair((align_len - chop + 1) as i32, align_len as i32);
        let tail_indel_ints = indel_ints.find_islands_ints(&tail_ints);

        if !tail_indel_ints.is_empty() {
            for _ in (tail_indel_ints.min() as usize)..=align_len {
                // record current length
                let cur_len = seqs.first().unwrap().len();
                for i in 0..seq_count {
                    let base = seqs[i].remove(cur_len - 1);
                    if base != '-' {
                        if ranges[i].strand == "+" || ranges[i].strand.is_empty() {
                            ranges[i].end -= 1;
                        } else {
                            ranges[i].start += 1;
                        }
                    }
                }
            }
        }
    }
}
