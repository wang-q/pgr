#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

const REF: &str = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATC\
GCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGAGTCCTAGGCTAACGGTTCGATCGATCGATC\
GATCGTAGGGAAACAGAACAGTCCTAGGCTAACGGTTCGATCGATCGATCGATCGTAGG";

fn rev_comp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => c,
        })
        .collect()
}

fn write_fastq(path: &std::path::Path, reads: &[(&str, &str)]) {
    let mut s = String::new();
    for (name, seq) in reads {
        s.push_str(&format!("@{name}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
    }
    std::fs::write(path, s).unwrap();
}

fn sam_lines(path: &std::path::Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| l.split('\t').map(|f| f.to_string()).collect())
        .collect()
}

/// Exact forward/reverse reads map with correct flags and positions.
#[test]
fn command_asm_map_exact_forward_reverse() {
    let out_dir = tempfile::tempdir().unwrap();
    let ref_file = out_dir.path().join("ut.fa");
    let reads_file = out_dir.path().join("reads.fq");
    let outm = out_dir.path().join("mapped.sam");
    let outu = out_dir.path().join("unmapped.sam");
    std::fs::write(&ref_file, format!(">ut\n{REF}\n")).unwrap();
    let fwd = &REF[10..60];
    let rev = rev_comp(&REF[30..80]);
    write_fastq(
        &reads_file,
        &[
            ("fwd", fwd),
            ("rev", &rev),
            ("mm", &mismatch(&REF[50..100])),
            ("gap", &with_gap(&REF[70..120])),
        ],
    );
    PgrCmd::new()
        .args(&[
            "asm",
            "map",
            ref_file.to_str().unwrap(),
            reads_file.to_str().unwrap(),
            "--outm",
            outm.to_str().unwrap(),
            "--outu",
            outu.to_str().unwrap(),
            "-k",
            "31",
        ])
        .assert()
        .success();
    let mapped = sam_lines(&outm);
    let unmapped = sam_lines(&outu);
    assert_eq!(mapped.len(), 2, "mapped: {mapped:?}");
    // fwd: flag 0, 1-based pos 11 (ref index 10).
    assert_eq!(mapped[0][0], "fwd");
    assert_eq!(mapped[0][1], "0");
    assert_eq!(mapped[0][2], "ut");
    assert_eq!(mapped[0][3], "11");
    assert_eq!(mapped[0][5], "50M");
    // rev: flag 16 (reverse strand), 1-based pos 31 (ref index 30).
    assert_eq!(mapped[1][0], "rev");
    assert_eq!(mapped[1][1], "16");
    assert_eq!(mapped[1][3], "31");
    // Mismatch and gap reads are unmapped (perfectmode semantics).
    assert_eq!(unmapped.len(), 2);
    assert!(unmapped.iter().any(|f| f[0] == "mm"));
    assert!(unmapped.iter().any(|f| f[0] == "gap"));
}

fn mismatch(seq: &str) -> String {
    let mut b = seq.as_bytes().to_vec();
    b[20] = if b[20] == b'A' { b'C' } else { b'A' };
    String::from_utf8(b).unwrap()
}

fn with_gap(seq: &str) -> String {
    let mut b = seq.as_bytes().to_vec();
    b.remove(20);
    String::from_utf8(b).unwrap()
}

/// A read matching a repeated region is reported at every position
/// (`ambiguous=all`).
#[test]
fn command_asm_map_ambiguous_all() {
    let out_dir = tempfile::tempdir().unwrap();
    let ref_file = out_dir.path().join("ut.fa");
    let reads_file = out_dir.path().join("reads.fq");
    let outm = out_dir.path().join("mapped.sam");
    std::fs::write(&ref_file, format!(">ut\n{REF}\n")).unwrap();
    let rep = &REF[100..140];
    assert_eq!(&REF[150..190], rep, "repeated region setup");
    write_fastq(&reads_file, &[("dup", rep)]);
    PgrCmd::new()
        .args(&[
            "asm",
            "map",
            ref_file.to_str().unwrap(),
            reads_file.to_str().unwrap(),
            "--outm",
            outm.to_str().unwrap(),
            "-k",
            "31",
        ])
        .assert()
        .success();
    let mapped = sam_lines(&outm);
    assert_eq!(mapped.len(), 2, "mapped: {mapped:?}");
    assert_eq!(mapped[0][3], "101");
    assert_eq!(mapped[1][3], "151");
}

/// A read shorter than k is unmapped.
#[test]
fn command_asm_map_short_read_unmapped() {
    let out_dir = tempfile::tempdir().unwrap();
    let ref_file = out_dir.path().join("ut.fa");
    let reads_file = out_dir.path().join("reads.fq");
    let outm = out_dir.path().join("mapped.sam");
    let outu = out_dir.path().join("unmapped.sam");
    std::fs::write(&ref_file, format!(">ut\n{REF}\n")).unwrap();
    write_fastq(&reads_file, &[("short", &REF[10..30])]);
    PgrCmd::new()
        .args(&[
            "asm",
            "map",
            ref_file.to_str().unwrap(),
            reads_file.to_str().unwrap(),
            "--outm",
            outm.to_str().unwrap(),
            "--outu",
            outu.to_str().unwrap(),
            "-k",
            "31",
        ])
        .assert()
        .success();
    assert!(sam_lines(&outm).is_empty());
    assert_eq!(sam_lines(&outu).len(), 1);
}

/// Per-base coverage is derived from the mapped SAM via `pgr sam to-rg` and
/// `pgr rg coverage` (the compositional anchors pipeline).
#[test]
fn command_asm_map_coverage_pipeline() {
    let out_dir = tempfile::tempdir().unwrap();
    let ref_file = out_dir.path().join("ut.fa");
    let reads_file = out_dir.path().join("reads.fq");
    let outm = out_dir.path().join("mapped.sam");
    std::fs::write(&ref_file, format!(">ut\n{REF}\n")).unwrap();
    write_fastq(&reads_file, &[("a", &REF[10..60]), ("b", &REF[20..70])]);
    PgrCmd::new()
        .args(&[
            "asm",
            "map",
            ref_file.to_str().unwrap(),
            reads_file.to_str().unwrap(),
            "--outm",
            outm.to_str().unwrap(),
            "-k",
            "31",
        ])
        .assert()
        .success();
    let (rg, _) = PgrCmd::new()
        .args(&["sam", "to-rg", outm.to_str().unwrap()])
        .run();
    // One range per mapped read, 1-based inclusive (ref 0-based 10..60 and
    // 20..70).
    assert_eq!(rg.lines().collect::<Vec<_>>(), ["ut:11-60", "ut:21-70"]);
    let rg_file = out_dir.path().join("mapped.rg");
    std::fs::write(&rg_file, &rg).unwrap();
    let cov_json = out_dir.path().join("cov.json");
    PgrCmd::new()
        .args(&[
            "rg",
            "coverage",
            rg_file.to_str().unwrap(),
            "-m",
            "2",
            "-o",
            cov_json.to_str().unwrap(),
        ])
        .assert()
        .success();
    // The 40-bp overlap (1-based 21..60) is the only region covered twice.
    let cov = std::fs::read_to_string(&cov_json).unwrap();
    assert!(cov.contains("21-60"), "cov: {cov}");
}

/// Output is deterministic across runs.
#[test]
fn command_asm_map_deterministic() {
    let out_dir = tempfile::tempdir().unwrap();
    let ref_file = out_dir.path().join("ut.fa");
    let reads_file = out_dir.path().join("reads.fq");
    std::fs::write(&ref_file, format!(">ut\n{REF}\n")).unwrap();
    write_fastq(
        &reads_file,
        &[("a", &REF[10..60]), ("b", &rev_comp(&REF[30..80]))],
    );
    let mut outs = Vec::new();
    for i in 0..2 {
        let outm = out_dir.path().join(format!("m{i}.sam"));
        let outu = out_dir.path().join(format!("u{i}.sam"));
        PgrCmd::new()
            .args(&[
                "asm",
                "map",
                ref_file.to_str().unwrap(),
                reads_file.to_str().unwrap(),
                "--outm",
                outm.to_str().unwrap(),
                "--outu",
                outu.to_str().unwrap(),
                "-k",
                "31",
            ])
            .assert()
            .success();
        outs.push((std::fs::read(&outm).unwrap(), std::fs::read(&outu).unwrap()));
    }
    assert_eq!(outs[0], outs[1]);
}

/// Lambda dataset sanity: deterministic, some reads map, coverage emitted.
#[test]
fn command_asm_map_lambda_sanity() {
    let out_dir = tempfile::tempdir().unwrap();
    let ref_file = out_dir.path().join("ut.fa");
    let reads_file = out_dir.path().join("reads.fq.gz");
    let outm = out_dir.path().join("mapped.sam");
    // Reference: assembled Lambda contigs (golden from the assemble tests).
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "asm",
            "contig",
            "tests/bbtools/Lambda/R1.2k.fq.gz",
            "tests/bbtools/Lambda/R2.2k.fq.gz",
            "-o",
            ref_file.to_str().unwrap(),
            "--kmer",
            "31",
            "--no-bubbles",
        ])
        .run();
    assert!(stdout.is_empty());
    // Reuse the same reads (a subset) as the reads to map.
    std::fs::copy("tests/bbtools/Lambda/R1.2k.fq.gz", &reads_file).unwrap();
    PgrCmd::new()
        .args(&[
            "asm",
            "map",
            ref_file.to_str().unwrap(),
            reads_file.to_str().unwrap(),
            "--outm",
            outm.to_str().unwrap(),
            "-k",
            "31",
        ])
        .assert()
        .success();
    let mapped = sam_lines(&outm);
    assert!(!mapped.is_empty(), "no reads mapped");
    let (rg, _) = PgrCmd::new()
        .args(&["sam", "to-rg", outm.to_str().unwrap()])
        .run();
    assert!(!rg.trim().is_empty());
}
