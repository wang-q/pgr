use crate::libs::poa::{AlignmentParams, AlignmentType, Poa};
use anyhow::anyhow;
use bio::io::fasta;
use intspan::IntSpan;
use std::fs;
use std::io::{BufRead, Write};
use std::path::Path;
use std::process::Command;
use std::str;

use super::coords::indel_intspan;
use crate::reader;

/// ```
/// match which::which("spoa") {
///     Ok(_) => {
///         let seqs = vec![
///         //              *
///             b"TTAGCCGCTGAGAAGC".as_ref(),
///             b"TTAGCCGCTGAGAAGC".as_ref(),
///             b"TTAGCCGCTGA-AAGC".as_ref(),
///         ];
///         let cons = pgr::libs::alignment::get_consensus_poa_external(&seqs, 5, -4, -8, -6, 0).unwrap();
///         assert_eq!(cons, "TTAGCCGCTGAGAAGC".to_string());
///
///         let seqs = vec![
///         //      *   **
///             b"AAATTTTGG".as_ref(),
///             b"AAAATTTTT".as_ref(),
///         ];
///         let cons = pgr::libs::alignment::get_consensus_poa_external(&seqs, 5, -4, -8, -6, 0).unwrap();
///         assert_eq!(cons, "AAAATTTTGG".to_string());
///
///         let seqs = vec![
///         //           **
///             b"AAAATTTTGG".as_ref(),
///             b"AAAATTTTTG".as_ref(),
///         ];
///         let cons = pgr::libs::alignment::get_consensus_poa_external(&seqs, 5, -4, -8, -6, 0).unwrap();
///         assert_eq!(cons, "AAAATTTTTG".to_string());
///
///         let seqs = vec![
///         //
///             b"AAAATTTTGG".as_ref(),
///         ];
///         let cons = pgr::libs::alignment::get_consensus_poa_external(&seqs, 5, -4, -8, -6, 0).unwrap();
///         assert_eq!(cons, "AAAATTTTGG".to_string());
///
///     }
///     Err(_) => {}
/// }
/// ```
// cargo test --doc alignment::get_consensus_poa_external
pub fn get_consensus_poa_external(
    seqs: &[&[u8]],
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    algo_code: i32,
) -> anyhow::Result<String> {
    let mut bin = String::new();
    for e in &["spoa"] {
        if let Ok(pth) = which::which(e) {
            bin = pth.to_string_lossy().to_string();
            break;
        }
    }

    if bin.is_empty() {
        return Err(anyhow!("Can't find the external command"));
    }

    let mut seq_in = tempfile::Builder::new()
        .prefix("seq-in-")
        .suffix(".fasta")
        .rand_bytes(8)
        .tempfile()?;

    for (i, seq) in seqs.iter().enumerate() {
        write!(seq_in, ">seq-{}\n{}\n", i, str::from_utf8(seq).unwrap())?;
    }
    let seq_in_path = seq_in.into_temp_path();

    let mut seq = String::new();
    let output = Command::new(bin)
        .arg("--result")
        .arg("0")
        .arg("-m")
        .arg(match_score.to_string())
        .arg("-n")
        .arg(mismatch_score.to_string())
        .arg("-g")
        .arg(gap_open.to_string())
        .arg("-e")
        .arg(gap_extend.to_string())
        .arg("-l")
        .arg(algo_code.to_string())
        .arg(seq_in_path.to_string_lossy().to_string())
        .output()?;

    if !output.status.success() {
        return Err(anyhow!("Command executed with failing error code"));
    }

    // closing the `TempPath` explicitly
    seq_in_path.close()?;

    for line in output.stdout.lines().map_while(Result::ok) {
        // header
        if line.starts_with('>') {
            continue;
        }

        seq += line.as_str();
    }

    Ok(seq)
}

/// ```
/// let seqs = vec![
/// //              *
///     b"TTAGCCGCTGAGAAGC".as_ref(),
///     b"TTAGCCGCTGAGAAGC".as_ref(),
///     b"TTAGCCGCTGA-AAGC".as_ref(),
/// ];
/// let cons = pgr::libs::alignment::get_consensus_poa_builtin(&seqs, 5, -4, -8, -6, 0).unwrap();
/// assert_eq!(cons, "TTAGCCGCTGAGAAGC".to_string());
///
/// let seqs = vec![
/// //      *   **
///     b"AAATTTTGG".as_ref(),
///     b"AAAATTTTT".as_ref(),
/// ];
/// let cons = pgr::libs::alignment::get_consensus_poa_builtin(&seqs, 5, -4, -8, -6, 0).unwrap();
/// assert_eq!(cons, "AAAATTTTGG".to_string());
///
/// let seqs = vec![
/// //           **
///     b"AAAATTTTGG".as_ref(),
///     b"AAAATTTTTG".as_ref(),
/// ];
/// let cons = pgr::libs::alignment::get_consensus_poa_builtin(&seqs, 5, -4, -8, -6, 0).unwrap();
/// assert_eq!(cons, "AAAATTTTTG".to_string());
///
/// let seqs = vec![
/// //
///     b"AAAATTTTGG".as_ref(),
/// ];
/// let cons = pgr::libs::alignment::get_consensus_poa_builtin(&seqs, 5, -4, -8, -6, 0).unwrap();
/// assert_eq!(cons, "AAAATTTTGG".to_string());
/// ```
// cargo test --doc alignment::get_consensus_poa_builtin
pub fn get_consensus_poa_builtin(
    seqs: &[&[u8]],
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    algo_code: i32,
) -> anyhow::Result<String> {
    let params = AlignmentParams {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    };
    let align_type = match algo_code {
        0 => AlignmentType::Local,
        1 => AlignmentType::Global,
        2 => AlignmentType::SemiGlobal,
        _ => AlignmentType::Global,
    };

    let mut poa = Poa::new(params, align_type);

    for seq in seqs {
        poa.add_sequence(seq);
    }

    let consensus = poa.consensus();
    let consensus_str = String::from_utf8(consensus)?;
    Ok(consensus_str)
}

/// Returns Strings to avoid lifetime issues
///
/// ```
/// match which::which("clustalw") {
///     Ok(_) => {
///         let seqs = vec![
///            //           *
///             "TTAGCCGCTGAGAAGC".to_string(),
///             "TTAGCCGCTGAGAAGC".to_string(),
///             "TTAGCCGCTGAAAGC".to_string(),
///         ];
///         let alns = pgr::libs::alignment::align_seqs(&seqs, "clustalw").unwrap();
///         assert_eq!(alns[2], "TTAGCCGCTGA-AAGC".to_string());
///
///     }
///     Err(_) => {}
/// }
/// ```
// scoop install clustalw
pub fn align_seqs(seqs: &[String], aligner: &str) -> anyhow::Result<Vec<String>> {
    // find external aligner
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
        "spoa" => {
            for e in &["spoa"] {
                if let Ok(pth) = which::which(e) {
                    bin = pth.to_string_lossy().to_string();
                    break;
                }
            }
        }
        "builtin" => {
            let params = AlignmentParams::default();
            let mut poa = Poa::new(params, AlignmentType::Global);
            for seq in seqs {
                poa.add_sequence(seq.as_bytes());
            }
            return Ok(poa.msa());
        }
        _ => {
            return Err(anyhow!("Unrecognized aligner: {}", aligner));
        }
    };

    if bin.is_empty() {
        return Err(anyhow!("Can't find the external command: {}", aligner));
    }

    // Create temp in/out files
    let mut seq_in = tempfile::Builder::new()
        .prefix("seq-in-")
        .suffix(".fa")
        .rand_bytes(8)
        .tempfile()?;

    // muscle may alter the sequence positions in alignments
    // clustalw wouldn't do this
    for (i, seq) in seqs.iter().enumerate() {
        write!(seq_in, ">seq-{}\n{}\n", i, seq)?;
    }
    let seq_in_path = seq_in.into_temp_path();

    let seq_out = tempfile::Builder::new()
        .prefix("seq-out-")
        .suffix("")
        .rand_bytes(8)
        .tempfile()?;
    let seq_out_path = seq_out.into_temp_path();

    // Run
    let output = match aligner {
        "clustalw" => Command::new(bin)
            .arg("-align")
            .arg("-type=dna")
            .arg("-output=fasta")
            .arg("-outorder=input")
            .arg("-quiet")
            .arg(format!("-infile={}", seq_in_path.to_string_lossy()))
            .arg(format!("-outfile={}", seq_out_path.to_string_lossy()))
            .output()?,
        "muscle" => Command::new(bin)
            .arg("-quiet")
            .arg("-in")
            .arg(seq_in_path.to_string_lossy().to_string())
            .arg("-out")
            .arg(seq_out_path.to_string_lossy().to_string())
            .output()?,
        "mafft" => Command::new(bin)
            .arg("--quiet")
            .arg("--auto")
            .arg(seq_in_path.to_string_lossy().to_string())
            .output()?,
        "spoa" => Command::new(bin)
            .arg("-r")
            .arg("1")
            .arg(seq_in_path.to_string_lossy().to_string())
            .output()?,
        _ => anyhow::bail!("unsupported aligner: {aligner}"),
    };

    if !output.status.success() {
        return Err(anyhow!(
            "Command executed with failing error code: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    // can't use redirect in Command
    if aligner == "mafft" || aligner == "spoa" {
        if output.stdout.is_empty() {
            return Err(anyhow!(
                "Command executed but returned empty stdout: {}",
                String::from_utf8_lossy(&output.stderr)
            ));
        }
        let mut f = crate::libs::io::writer(&seq_out_path.to_string_lossy())?;
        f.write_all(&output.stdout)?;
    }
    // delete .dnd files created by clustalw
    if aligner == "clustalw" {
        let dnd = Path::new(&seq_in_path).with_extension("dnd");
        if dnd.exists() {
            let _ = fs::remove_file(dnd);
        }
    }

    // Load outputs
    let mut out_seqs = Vec::with_capacity(seqs.len());
    // init elements in the vector
    for _ in 0..seqs.len() {
        out_seqs.push("".to_string());
    }
    let reader = reader(seq_out_path.to_string_lossy().as_ref())?;
    let fa_in = fasta::Reader::new(reader);
    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result.unwrap();

        let idx = record.id().to_string();
        let idx = idx.replace("seq-", "");
        let idx = idx.parse::<usize>().unwrap();

        out_seqs[idx] = String::from_utf8(record.seq().to_vec().to_ascii_uppercase()).unwrap();
    }

    // closing the `TempPath` explicitly
    seq_in_path.close()?;
    seq_out_path.close()?;

    Ok(out_seqs)
}

pub fn align_seqs_quick(
    seqs: &[String],
    aligner: &str,
    pad: i32,
    fill: i32,
) -> anyhow::Result<Vec<String>> {
    let count = seqs.len();
    let align_len = seqs.first().unwrap().len() as i32;

    // realign regions
    let mut realign_ints = IntSpan::new();
    // Add head and tail
    realign_ints.add_pair(1, pad);
    realign_ints.add_pair(align_len - pad, align_len);

    for seq in seqs.iter().take(count) {
        let mut ints = indel_intspan(seq.as_bytes().to_vec().as_ref());
        ints = ints.pad(pad);
        realign_ints.merge(&ints);
    }
    // join adjacent realign regions
    realign_ints = realign_ints.fill(fill);
    realign_ints = realign_ints.intersect(&IntSpan::from_pair(1, align_len));

    // all segments
    let mut aligned: Vec<String> = seqs.to_owned();
    for (lower, upper) in realign_ints.spans().iter().rev() {
        let mut subseqs = vec![];
        let start = *lower as usize - 1;
        let end = *upper as usize;

        // extract subseqs
        for a in aligned.iter().take(count) {
            let subseq = &a[start..end];
            subseqs.push(subseq.to_string());
        }
        let subseqs = align_seqs(&subseqs, aligner)?;

        // put aligned subseqs back
        for (a, s) in aligned.iter_mut().take(count).zip(subseqs.iter()) {
            a.replace_range(start..end, s.as_str());
        }
    }

    Ok(aligned)
}
