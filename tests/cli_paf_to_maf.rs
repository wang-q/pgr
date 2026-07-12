#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::path::PathBuf;

/// Return the absolute path to a fixture in `tests/paf/input`.
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/paf/input")
        .join(name)
}

// ── paf to-maf (pairwise MAF from CIGAR) ─────────────────────────

#[test]
fn command_paf_to_maf_strict_name_validation() {
    // PAF references A and B; TSV only has A.
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "-f",
            fixture("A_only.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run_fail();
    assert!(
        stderr.contains("FASTA TSV is missing") && stderr.contains("B"),
        "missing strict validation error for B"
    );
}

#[test]
fn command_paf_to_maf_perfect_match() {
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "-f",
            fixture("AB.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(stdout.contains("##maf version=1"), "missing MAF header");
    assert!(stdout.contains("a"), "missing alignment header");
    // target line first, query line second
    assert!(
        stdout.contains("s\tB\t0\t10\t+\t10\tACGTACGTAC"),
        "missing/incorrect target line"
    );
    assert!(
        stdout.contains("s\tA\t0\t10\t+\t10\tACGTACGTAC"),
        "missing/incorrect query line"
    );
}

#[test]
fn command_paf_to_maf_with_insertion() {
    // CIGAR: 4= 3I 3= → target 0-7, query 0-10
    // target: ACGT---ACG  (4 match + 3 gaps + 3 match)
    // query:  ACGTACGTAC  (4 match + 3 bases + 3 match)
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t7\t7\t10\t255\tcg:Z:4=3I3=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-7",
            "-f",
            fixture("AB.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(
        stdout.contains("ACGT---ACG"),
        "target alignment should contain gaps for insertion"
    );
    assert!(
        stdout.contains("ACGTACGTAC"),
        "query alignment should contain full sequence"
    );
    assert!(
        stdout.contains("s\tB\t0\t7\t+\t10\tACGT---ACG"),
        "target size should be 7"
    );
    assert!(
        stdout.contains("s\tA\t0\t10\t+\t10\tACGTACGTAC"),
        "query size should be 10"
    );
}

#[test]
fn command_paf_to_maf_with_deletion() {
    // CIGAR: 4= 3D 3= → target 0-10, query 0-7
    // target: ACGTACGTAC, query: ACGT---ACG
    let paf = "A\t7\t0\t7\t+\tB\t10\t0\t10\t7\t10\t255\tcg:Z:4=3D3=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "-f",
            fixture("AB_A7.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(
        stdout.contains("ACGTACGTAC"),
        "target alignment should contain full sequence"
    );
    assert!(
        stdout.contains("ACGT---ACG"),
        "query alignment should contain gaps for deletion"
    );
}

#[test]
fn command_paf_to_maf_trimmed_subregion() {
    // Full alignment 10= over B:0-10. Query B:2-8 should trim CIGAR to 6=.
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:2-8",
            "-f",
            fixture("AB.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // B[2..8) = GTACGT, A[2..8) = GTACGT
    assert!(
        stdout.contains("s\tB\t2\t6\t+\t10\tGTACGT"),
        "target should be trimmed to B:2-8"
    );
    assert!(
        stdout.contains("s\tA\t2\t6\t+\t10\tGTACGT"),
        "query should be trimmed to A:2-8"
    );
}

#[test]
fn command_paf_to_maf_reverse_strand_perfect_match() {
    // '-' strand perfect match: target B forward == RC(query A forward).
    // A forward = GTACGTACGT, RC = ACGTACGTAC = B forward.
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "-f",
            fixture("AB_rc.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(stdout.contains("##maf version=1"), "missing MAF header");
    assert!(
        stdout.contains("s\tB\t0\t10\t+\t10\tACGTACGTAC"),
        "missing/incorrect target line for '-' strand record"
    );
    // q_start_maf = srcSize - qe = 10 - 10 = 0; q_size = 10.
    assert!(
        stdout.contains("s\tA\t0\t10\t-\t10\tACGTACGTAC"),
        "missing/incorrect query line for '-' strand record (RC not applied)"
    );
}

#[test]
fn command_paf_to_maf_reverse_strand_with_insertion() {
    // '-' strand alignment with insertion: CIGAR 4=3I3= (7 target, 10 query cols).
    // A forward = GTACGTACGT, RC(A) = ACGTACGTAC.
    // target B = ACGT (RC(A)[0:4]) + TAC (RC(A)[7:10]) = ACGTTAC (7 bp).
    let paf = "A\t10\t0\t10\t-\tB\t7\t0\t7\t7\t7\t255\tcg:Z:4=3I3=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-7",
            "-f",
            fixture("AB_rc_B7_1.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(
        stdout.contains("ACGT---TAC"),
        "target alignment should contain gaps for insertion on '-' strand"
    );
    assert!(
        stdout.contains("ACGTACGTAC"),
        "query alignment should be RC of A forward on '-' strand"
    );
    assert!(
        stdout.contains("s\tA\t0\t10\t-\t10\tACGTACGTAC"),
        "missing/incorrect query s-line for '-' strand with insertion"
    );
}

#[test]
fn command_paf_to_maf_reverse_strand_subinterval_first_half() {
    // '-' strand perfect match, sub-interval query on the first half.
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-5",
            "-f",
            fixture("AB_rc.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(
        stdout.contains("s\tB\t0\t5\t+\t10\tACGTA"),
        "missing/incorrect target line for '-' strand sub-interval (first half)"
    );
    assert!(
        stdout.contains("s\tA\t0\t5\t-\t10\tACGTA"),
        "missing/incorrect query line for '-' strand sub-interval (first half)"
    );
    assert!(
        !stdout.contains("CGTAC"),
        "regression: query sequence looks like RC of forward A[0:5]"
    );
}

#[test]
fn command_paf_to_maf_reverse_strand_subinterval_second_half() {
    // Same setup, but query B:5-10.
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:5-10",
            "-f",
            fixture("AB_rc.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(
        stdout.contains("s\tB\t5\t5\t+\t10\tCGTAC"),
        "missing/incorrect target line for '-' strand sub-interval (second half)"
    );
    assert!(
        stdout.contains("s\tA\t5\t5\t-\t10\tCGTAC"),
        "missing/incorrect query line for '-' strand sub-interval (second half)"
    );
}

#[test]
fn command_paf_to_maf_reverse_strand_subinterval_with_insertion() {
    // '-' strand with insertion, sub-interval query on the trailing target segment.
    let paf = "A\t10\t0\t10\t-\tB\t7\t0\t7\t7\t7\t255\tcg:Z:5=3I2=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:5-7",
            "-f",
            fixture("AB_rc_B7_2.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(
        stdout.contains("s\tB\t5\t2\t+\t7\t---AC"),
        "missing/incorrect target line for '-' strand sub-interval with insertion"
    );
    assert!(
        stdout.contains("s\tA\t5\t5\t-\t10\tCGTAC"),
        "missing/incorrect query line for '-' strand sub-interval with insertion"
    );
}

// ── paf to-maf --msa (multi-way MSA via POA) ─────────────────────

#[test]
fn command_paf_to_maf_msa_three_genomes_transitive() {
    // Three genomes A/B/C, all 10 bp, identical sequence ACGTACGTAC.
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "--transitive",
            "--msa",
            "-f",
            fixture("ABC.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(
        stderr.contains("Total results:") && !stderr.contains("Total results: 0"),
        "expected non-zero results"
    );
    assert!(stdout.contains("##maf version=1"), "missing MAF header");
    let a_count = stdout.matches("\na\n").count() + if stdout.starts_with("a\n") { 1 } else { 0 };
    assert_eq!(a_count, 1, "expected exactly one MAF block, got {a_count}");
    let s_count = stdout.lines().filter(|l| l.starts_with("s\t")).count();
    assert_eq!(s_count, 3, "expected 3 s-lines, got {s_count}");
    for line in stdout.lines().filter(|l| l.starts_with("s\t")) {
        assert!(
            line.ends_with("ACGTACGTAC"),
            "expected gap-free ACGTACGTAC in s-line: {line}"
        );
    }
    let first_s = stdout.lines().find(|l| l.starts_with("s\t")).unwrap();
    assert!(
        first_s.starts_with("s\tB\t"),
        "expected B as the first s-line, got {first_s}"
    );
}

#[test]
fn command_paf_to_maf_msa_with_snp() {
    // B = ACGTACGTAC, A = ACGTACGTAC, C = ACGTTCGTAC (SNP at pos 4).
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "--transitive",
            "--msa",
            "-f",
            fixture("ABC_snp.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    let s_count = stdout.lines().filter(|l| l.starts_with("s\t")).count();
    assert_eq!(s_count, 3, "expected 3 s-lines, got {s_count}");
    for line in stdout.lines().filter(|l| l.starts_with("s\t")) {
        let aln = line.split('\t').next_back().unwrap();
        assert_eq!(aln.len(), 10, "expected 10-char alignment, got '{aln}'");
    }
    let b_line = stdout.lines().find(|l| l.starts_with("s\tB\t")).unwrap();
    let c_line = stdout.lines().find(|l| l.starts_with("s\tC\t")).unwrap();
    let b_aln = b_line.split('\t').next_back().unwrap();
    let c_aln = c_line.split('\t').next_back().unwrap();
    let diffs: Vec<usize> = b_aln
        .chars()
        .zip(c_aln.chars())
        .enumerate()
        .filter_map(|(i, (a, b))| if a != b { Some(i) } else { None })
        .collect();
    assert_eq!(
        diffs,
        vec![4],
        "expected single SNP at pos 4, got {diffs:?}"
    );
}

#[test]
fn command_paf_to_maf_msa_reverse_strand_query() {
    // B = ACGTACGTAC, A = GTACGTACGT on '-' strand, RC(A) = B.
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "--transitive",
            "--msa",
            "-f",
            fixture("AB_rc.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    let s_count = stdout.lines().filter(|l| l.starts_with("s\t")).count();
    assert_eq!(s_count, 2, "expected 2 s-lines (B + A), got {s_count}");
    let a_line = stdout.lines().find(|l| l.starts_with("s\tA\t")).unwrap();
    assert!(
        a_line.contains("\t-\t"),
        "A should be on '-' strand: {a_line}"
    );
    let a_aln = a_line.split('\t').next_back().unwrap();
    assert_eq!(a_aln, "ACGTACGTAC", "expected RC(A) gap-free: {a_aln}");
}
