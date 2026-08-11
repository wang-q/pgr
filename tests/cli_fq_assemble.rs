#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::Read;

/// Parses a FASTA into (header, sequence) pairs.
fn parse_fa(data: &[u8]) -> Vec<(String, String)> {
    let mut recs = Vec::new();
    let mut cur: Option<(String, String)> = None;
    for line in std::str::from_utf8(data).unwrap().lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(c) = cur.take() {
                recs.push(c);
            }
            cur = Some((rest.to_string(), String::new()));
        } else if let Some(c) = cur.as_mut() {
            c.1.push_str(line);
        }
    }
    if let Some(c) = cur {
        recs.push(c);
    }
    recs
}

fn read_gz(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    let mut dec = flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap());
    dec.read_to_end(&mut out).unwrap();
    out
}

/// The contig set is deterministic and statistically matches the BBTools
/// `tadpole.sh threads=1` golden. The pre-pop contig set is byte-identical;
/// bubble-popping resolutions can differ because BBTools' expand order
/// depends on its memory-dependent hash-table layout (see
/// notes/design/fq-assemble.md), so the popped output is compared by
/// sequence set rather than byte-for-byte.
#[test]
fn command_fq_assemble_matches_tadpole_contig_set() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("contigs.fa");
    PgrCmd::new()
        .args(&[
            "fq",
            "assemble",
            "tests/bbtools/Lambda/R1.2k.fq.gz",
            "tests/bbtools/Lambda/R2.2k.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
        ])
        .assert()
        .success();
    let pgr = parse_fa(&std::fs::read(&out).unwrap());
    let golden = parse_fa(&read_gz(
        "tests/bbtools/Lambda/golden/tadpole_contigs31.fasta.gz",
    ));

    // Same total assembled bases as the reference.
    let pgr_bases: usize = pgr.iter().map(|(_, s)| s.len()).sum();
    let golden_bases: usize = golden.iter().map(|(_, s)| s.len()).sum();
    // Bubble-resolution differences can shift a few contigs between the
    // kept and merged sets, so allow a small total-base delta.
    assert!(
        (pgr_bases as i64 - golden_bases as i64).abs() <= 100,
        "bases {pgr_bases} vs {golden_bases}"
    );

    // Contig count within 1 of the reference (pgr is deterministic; the
    // remaining bubble-resolution differences are documented).
    assert!(
        (pgr.len() as i64 - golden.len() as i64).abs() <= 1,
        "{} vs {}",
        pgr.len(),
        golden.len()
    );

    // At least 90% of the reference contigs are present verbatim.
    let golden_seqs: std::collections::HashSet<&str> =
        golden.iter().map(|(_, s)| s.as_str()).collect();
    let shared = pgr
        .iter()
        .filter(|(_, s)| golden_seqs.contains(s.as_str()))
        .count();
    assert!(
        shared * 10 >= golden.len() * 9,
        "shared {shared}/{}",
        golden.len()
    );
}

/// Repeated runs produce byte-identical output (deterministic scan order).
#[test]
fn command_fq_assemble_is_deterministic() {
    let out_dir = tempfile::tempdir().unwrap();
    let out1 = out_dir.path().join("a.fa");
    let out2 = out_dir.path().join("b.fa");
    for out in [&out1, &out2] {
        PgrCmd::new()
            .args(&[
                "fq",
                "assemble",
                "tests/bbtools/Lambda/R1.2k.fq.gz",
                "tests/bbtools/Lambda/R2.2k.fq.gz",
                "-o",
                out.to_str().unwrap(),
                "--kmer",
                "31",
            ])
            .assert()
            .success();
    }
    assert_eq!(std::fs::read(&out1).unwrap(), std::fs::read(&out2).unwrap());
}

/// A zero k-mer length must fail cleanly instead of panicking.
#[test]
fn command_fq_assemble_rejects_zero_kmer() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fa");
    PgrCmd::new()
        .args(&[
            "fq",
            "assemble",
            "tests/bbtools/Lambda/R1.2k.fq.gz",
            "tests/bbtools/Lambda/R2.2k.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "0",
        ])
        .assert()
        .failure();
}

/// Assembles a small synthetic repeat into contigs.
#[test]
fn command_fq_assemble_small_synthetic() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    // Two identical 60 bp reads: the k=31 graph assembles them into a
    // single 60 bp contig (below the default 124 bp output threshold, so
    // use --min-contig-len 1).
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    std::fs::write(&infile, format!(">r1\n{seq}\n>r2\n{seq}\n")).unwrap();
    PgrCmd::new()
        .args(&[
            "fq",
            "assemble",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    assert!(!recs.is_empty());
    assert!(recs.iter().any(|(_, s)| s.contains("ACGTACGT")));
}
