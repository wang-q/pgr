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

// ── paf to-fas (pairwise FAS from CIGAR) ─────────────────────────

#[test]
fn command_paf_to_fas_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "to-fas", "--help"]).run();
    assert!(stdout.contains("block FASTA"));
    assert!(stdout.contains("--fasta-tsv"));
}

#[test]
fn command_paf_to_fas_strict_name_validation() {
    // PAF references A and B; TSV only has A.
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-fas",
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
fn command_paf_to_fas_perfect_match() {
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-fas",
            "stdin",
            "B:0-10",
            "-f",
            fixture("AB.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(
        stdout.contains(">B(+):1-10\nACGTACGTAC"),
        "missing/incorrect target record"
    );
    assert!(
        stdout.contains(">A(+):1-10\nACGTACGTAC"),
        "missing/incorrect query record"
    );
}

#[test]
fn command_paf_to_fas_with_insertion() {
    // CIGAR: 4= 3I 3= → target 0-7, query 0-10
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t7\t7\t10\t255\tcg:Z:4=3I3=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-fas",
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
}

#[test]
fn command_paf_to_fas_with_deletion() {
    // CIGAR: 4= 3D 3= → target 0-10, query 0-7
    let paf = "A\t7\t0\t7\t+\tB\t10\t0\t10\t7\t10\t255\tcg:Z:4=3D3=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-fas",
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
fn command_paf_to_fas_reverse_strand_perfect_match() {
    // '-' strand perfect match: target B forward == RC(query A forward).
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-fas",
            "stdin",
            "B:0-10",
            "-f",
            fixture("AB_rc.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(
        stdout.contains(">B(+):1-10\nACGTACGTAC"),
        "missing/incorrect target record for '-' strand"
    );
    assert!(
        stdout.contains(">A(-):1-10\nACGTACGTAC"),
        "missing/incorrect query record for '-' strand (RC not applied)"
    );
}

// ── paf to-fas --msa (multi-way MSA via POA) ─────────────────────

#[test]
fn command_paf_to_fas_msa_three_genomes_transitive() {
    // Three genomes A/B/C, all 10 bp, identical sequence ACGTACGTAC.
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-fas",
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
    let records: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(records.len(), 3, "expected 3 records, got {records:?}");
    assert!(
        records[0].contains(">B("),
        "expected B as the first record, got {records:?}"
    );
    let seqs: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('>'))
        .collect();
    for s in &seqs {
        assert_eq!(*s, "ACGTACGTAC", "expected gap-free ACGTACGTAC, got '{s}'");
    }
}

#[test]
fn command_paf_to_fas_msa_with_snp() {
    // B = ACGTACGTAC, A = ACGTACGTAC, C = ACGTTCGTAC (SNP at pos 4).
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-fas",
            "stdin",
            "B:0-10",
            "--transitive",
            "--msa",
            "-f",
            fixture("ABC_snp.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    let records: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(records.len(), 3, "expected 3 records");
    let seqs: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('>'))
        .collect();
    for s in &seqs {
        assert_eq!(s.len(), 10, "expected 10-char alignment, got '{s}'");
    }
    let b_idx = records.iter().position(|r| r.contains(">B(")).unwrap();
    let c_idx = records.iter().position(|r| r.contains(">C(")).unwrap();
    let b_aln = seqs[b_idx];
    let c_aln = seqs[c_idx];
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
fn command_paf_to_fas_msa_reverse_strand_query() {
    // B = ACGTACGTAC, A = GTACGTACGT on '-' strand, RC(A) = B.
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-fas",
            "stdin",
            "B:0-10",
            "--transitive",
            "--msa",
            "-f",
            fixture("AB_rc.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    let records: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(records.len(), 2, "expected 2 records (B + A)");
    assert!(
        records.iter().any(|r| r.contains(">A(-):")),
        "A should be on '-' strand: {records:?}"
    );
    let a_idx = records.iter().position(|r| r.contains(">A(")).unwrap();
    let seqs: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('>'))
        .collect();
    assert_eq!(seqs[a_idx], "ACGTACGTAC", "expected RC(A) gap-free");
}

// ── pipeline: paf to-fas | fas to-vcf ────────────────────────────

#[test]
fn command_paf_to_fas_pipeline_to_vcf() {
    // B = ACGTACGTAC (target, REF)
    // A = ACGTACGTAC (identical to B)
    // C = ACGTTCGTAC (SNP at position 4)
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";

    // Step 1: paf to-fas --msa → block FASTA
    let (fas_stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-fas",
            "stdin",
            "B:0-10",
            "--transitive",
            "--msa",
            "-f",
            fixture("ABC_snp.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(!fas_stdout.is_empty(), "to-fas output should not be empty");

    // Step 2: fas to-vcf ← block FASTA from step 1
    let (vcf_stdout, _stderr) = PgrCmd::new()
        .args(&["fas", "to-vcf", "stdin"])
        .stdin(&fas_stdout)
        .run();
    assert!(
        vcf_stdout.contains("##fileformat=VCFv4.2"),
        "missing VCF header: {vcf_stdout}"
    );
    assert!(
        vcf_stdout.contains("#CHROM"),
        "missing #CHROM header: {vcf_stdout}"
    );
    let body: Vec<&str> = vcf_stdout
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(body.len(), 1, "expected 1 variant row, got {body:?}");
    let fields: Vec<&str> = body[0].split('\t').collect();
    assert_eq!(fields[0], "B", "CHROM");
    assert_eq!(fields[1], "5", "POS (1-based)");
    assert_eq!(fields[3], "A", "REF");
    assert_eq!(fields[4], "T", "ALT");
}
