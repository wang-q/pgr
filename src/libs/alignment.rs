use crate::*;
use anyhow::anyhow;
use bio::io::fasta;
use itertools::Itertools;
use std::collections::HashSet;
use std::io::Write;
use std::process::Command;
use std::str;

lazy_static! {
    static ref BASES: HashSet<u8> = vec![b'a', b'g', b'c', b't', b'A', b'G', b'C', b'T',]
        .into_iter()
        .collect();
}

/// Divergence (D) between two sequences
///
/// ```
/// //           * **  **
/// let seq1 = b"GTCTGCATGCN";
/// let seq2 = b"TTTAGCTAgc-";
/// // difference 5
/// // comparable 10
/// assert_eq!(intspan::pair_d(seq1, seq2), 0.5);
/// ```
pub fn pair_d(seq1: &[u8], seq2: &[u8]) -> f32 {
    assert_eq!(
        seq1.len(),
        seq2.len(),
        "Two sequences of different length ({}!={})",
        seq1.len(),
        seq2.len()
    );

    let mut comparable = 0;
    let mut difference = 0;

    for (base1, base2) in seq1.iter().zip(seq2) {
        if BASES.contains(base1) && BASES.contains(base2) {
            comparable += 1;
            if base1.to_ascii_uppercase() != base2.to_ascii_uppercase() {
                difference += 1;
            }
        }
    }

    assert_ne!(comparable, 0, "Comparable bases shouldn't be zero");

    // eprintln!("{} {}", difference, comparable);

    difference as f32 / comparable as f32
}

/// Basic stats on alignments
///
/// ```
/// let seqs = vec![
///     //        *
///     b"AAAATTTTGG".as_ref(),
///     b"aaaatttttg".as_ref(),
/// ];
/// assert_eq!(intspan::alignment_stat(&seqs), (10, 10, 1, 0, 0, 0.1,));
///
/// let seqs = vec![
///     //*          * *
///     b"TTAGCCGCTGAGAAGCC".as_ref(),
///     b"GTAGCCGCTGA-AGGCC".as_ref(),
/// ];
/// assert_eq!(intspan::alignment_stat(&seqs), (17, 16, 2, 1, 0, 0.125,));
///
/// let seqs = vec![
///     //    * **    *   ** *   *
///     b"GATTATCATCACCCCAGCCACATW".as_ref(),
///     b"GATTTT--TCACTCCATTCGCATA".as_ref(),
/// ];
/// assert_eq!(intspan::alignment_stat(&seqs), (24, 21, 5, 2, 1, 0.238,));
///
/// ```
pub fn alignment_stat(seqs: &[&[u8]]) -> (i32, i32, i32, i32, i32, f32) {
    let seq_count = seqs.len();
    assert_ne!(seq_count, 0, "Need sequences");

    let length = seqs[0].len();

    let mut comparable = 0;
    let mut difference = 0;
    let mut gap = 0;
    let mut ambiguous = 0;

    // For each position, search for polymorphic sites
    for pos in 0..length {
        let mut column = vec![];
        for i in 0..seq_count {
            column.push(seqs[i][pos].to_ascii_uppercase());
        }
        column = column.into_iter().unique().collect();

        if column.clone().into_iter().all(|e| BASES.contains(&e)) {
            comparable += 1;
            if column.clone().into_iter().any(|e| e != column[0]) {
                difference += 1;
            }
        } else if column.clone().into_iter().any(|e| e == b'-') {
            gap += 1;
        } else {
            ambiguous += 1;
        }
    }

    assert_ne!(comparable, 0, "Comparable bases shouldn't be zero");

    let mut dists = vec![];
    for i in 0..seq_count {
        for j in i + 1..seq_count {
            let dist = pair_d(seqs[i], seqs[j]);
            dists.push(dist);
        }
    }

    let mean_d = f32::trunc(dists.iter().sum::<f32>() / dists.len() as f32 * 1000.0) / 1000.0;

    (
        length as i32,
        comparable,
        difference,
        gap,
        ambiguous,
        mean_d,
    )
}

pub fn indel_intspan(seqs: &[u8]) -> IntSpan {
    let mut positions = vec![];

    for (i, base) in seqs.iter().enumerate() {
        if *base == b'-' {
            positions.push(i as i32);
        }
    }

    let mut ints = IntSpan::new();
    ints.add_vec(&positions);

    ints
}

/// ```
/// # // scoop install clustalw
/// match which::which("clustalw") {
///     Ok(_) => {
///         let seqs = vec![
///             //           *
///             b"TTAGCCGCTGAGAAGC".as_ref(),
///             b"TTAGCCGCTGAGAAGC".as_ref(),
///             b"TTAGCCGCTGAAAGC".as_ref(),
///         ];
///         let alns = intspan::align_seqs(&seqs, "clustalw").unwrap();
///         assert_eq!(alns[2], "TTAGCCGCTGA-AAGC".to_string());
///
///     }
///     Err(_) => {}
/// }
/// ```
// cargo test --doc utils::get_consensus_poa
pub fn align_seqs(seqs: &[&[u8]], aligner: &str) -> anyhow::Result<Vec<String>> {
    let mut bin = String::new();

    match aligner {
        "clustalw" => {
            for e in &["clustalw", "clustal-w", "clustalw2"] {
                if let Ok(pth) = which::which(e) {
                    bin = pth.to_string_lossy().to_string();
                    break;
                }
            }
        }
        "muscle" => {
            for e in &["muscle"] {
                if let Ok(pth) = which::which(e) {
                    bin = pth.to_string_lossy().to_string();
                    break;
                }
            }
        }
        "mafft" => {
            for e in &["mafft"] {
                if let Ok(pth) = which::which(e) {
                    bin = pth.to_string_lossy().to_string();
                    break;
                }
            }
        }
        _ => {
            return Err(anyhow!("Unrecognized aligner: {}", aligner));
        }
    };

    eprintln!("bin = {:#?}", bin);

    if bin.is_empty() {
        return Err(anyhow!("Can't find the external command: {}", aligner));
    }

    let mut seq_in = tempfile::Builder::new()
        .prefix("seq-in-")
        .suffix(".fasta")
        .rand_bytes(8)
        .tempfile()?;
    for (i, seq) in seqs.iter().enumerate() {
        write!(seq_in, ">seq-{}\n{:?}\n", i, str::from_utf8(seq).unwrap())?;
    }
    let seq_in_path = seq_in.into_temp_path();

    let seq_out = tempfile::Builder::new()
        .prefix("seq-out-")
        .suffix(".fasta")
        .rand_bytes(8)
        .tempfile()?;
    let seq_out_path = seq_out.into_temp_path();

    eprintln!("seq_in_path = {:#?}", seq_in_path);

    let output = match aligner {
        "clustalw" => Command::new(bin)
            .arg("-align")
            .arg("-type=dna")
            .arg("-output=fasta")
            .arg("-outorder=input")
            .arg("-quiet")
            .arg(format!(
                "-infile={}",
                seq_in_path.to_string_lossy().to_string()
            ))
            .arg(format!(
                "-outfile={}",
                seq_out_path.to_string_lossy().to_string()
            ))
            .output()?,
        "muscle" => Command::new(bin)
            .arg("-quiet")
            .arg("-in")
            .arg(seq_in_path.to_string_lossy().to_string())
            .arg("-out")
            .arg(seq_out_path.to_string_lossy().to_string())
            .output()?,
        "mafft" => Command::new(bin)
            .arg("-quiet")
            .arg("-auto")
            .arg(seq_in_path.to_string_lossy().to_string())
            .arg(">")
            .arg(seq_out_path.to_string_lossy().to_string())
            .output()?,
        _ => unreachable!(),
    };

    eprintln!("output = {:#?}", output);

    if !output.status.success() {
        return Err(anyhow!("Command executed with failing error code"));
    }

    let mut out_seq = vec![];
    let reader = reader(seq_out_path.to_string_lossy().as_ref());
    let fa_in = fasta::Reader::new(reader);
    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result.unwrap();
        out_seq.push(String::from_utf8(record.seq().to_vec()).unwrap());
    }

    // closing the `TempPath` explicitly
    seq_in_path.close()?;
    seq_out_path.close()?;

    Ok(out_seq)
}
