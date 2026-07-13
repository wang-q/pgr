use crate::libs::nt::NT_VAL;
use anyhow::bail;
use itertools::Itertools;

/// Divergence (D) between two sequences
///
/// ```ignore
/// //           * **  **
/// let seq1 = b"GTCTGCATGCN";
/// let seq2 = b"TTTAGCTAgc-";
/// // difference 5
/// // comparable 10
/// assert_eq!(pgr::libs::alignment::pair_d(seq1, seq2).unwrap(), 0.5);
/// ```
pub fn pair_d(seq1: &[u8], seq2: &[u8]) -> anyhow::Result<f32> {
    if seq1.len() != seq2.len() {
        bail!(
            "Two sequences of different length ({}!={})",
            seq1.len(),
            seq2.len()
        );
    }

    let mut comparable = 0;
    let mut difference = 0;

    for (base1, base2) in seq1.iter().zip(seq2) {
        if NT_VAL[*base1 as usize] <= 3 && NT_VAL[*base2 as usize] <= 3 {
            comparable += 1;
            if !base1.eq_ignore_ascii_case(base2) {
                difference += 1;
            }
        }
    }

    if comparable == 0 {
        bail!("Comparable bases shouldn't be zero");
    }

    Ok(difference as f32 / comparable as f32)
}

/// Basic stats on alignments
///
/// ```ignore
/// let seqs = vec![
///     //        *
///     b"AAAATTTTGG".as_ref(),
///     b"aaaatttttg".as_ref(),
/// ];
/// assert_eq!(pgr::libs::alignment::alignment_stat(&seqs).unwrap(), (10, 10, 1, 0, 0, 0.1,));
///
/// let seqs = vec![
///     //*          * *
///     b"TTAGCCGCTGAGAAGCC".as_ref(),
///     b"GTAGCCGCTGA-AGGCC".as_ref(),
/// ];
/// assert_eq!(pgr::libs::alignment::alignment_stat(&seqs).unwrap(), (17, 16, 2, 1, 0, 0.125,));
///
/// let seqs = vec![
///     //    * **    *   ** *   *
///     b"GATTATCATCACCCCAGCCACATW".as_ref(),
///     b"GATTTT--TCACTCCATTCGCATA".as_ref(),
/// ];
/// assert_eq!(pgr::libs::alignment::alignment_stat(&seqs).unwrap(), (24, 21, 5, 2, 1, 0.238,));
///
/// ```
pub fn alignment_stat(seqs: &[&[u8]]) -> anyhow::Result<(i32, i32, i32, i32, i32, f32)> {
    let seq_count = seqs.len();
    if seq_count == 0 {
        bail!("Need sequences");
    }

    let length = seqs[0].len();

    let mut comparable = 0;
    let mut difference = 0;
    let mut gap = 0;
    let mut ambiguous = 0;

    // For each position, search for polymorphic sites
    #[allow(clippy::needless_range_loop)]
    for pos in 0..length {
        let mut column = vec![];
        for seq in seqs.iter().take(seq_count) {
            column.push(seq[pos].to_ascii_uppercase());
        }
        column = column.into_iter().unique().collect();

        if column.iter().all(|e| NT_VAL[*e as usize] <= 3) {
            comparable += 1;
            if column.iter().any(|e| *e != column[0]) {
                difference += 1;
            }
        } else if column.contains(&b'-') {
            gap += 1;
        } else {
            ambiguous += 1;
        }
    }

    if comparable == 0 {
        bail!("Comparable bases shouldn't be zero");
    }

    if seq_count < 2 {
        return Ok((length as i32, comparable, difference, gap, ambiguous, 0.0));
    }

    let mut dists = vec![];
    for i in 0..seq_count {
        for j in i + 1..seq_count {
            let dist = pair_d(seqs[i], seqs[j])?;
            dists.push(dist);
        }
    }

    let mean_d = f32::trunc(dists.iter().sum::<f32>() / dists.len() as f32 * 1000.0) / 1000.0;

    Ok((
        length as i32,
        comparable,
        difference,
        gap,
        ambiguous,
        mean_d,
    ))
}
