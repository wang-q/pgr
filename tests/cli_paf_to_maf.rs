#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── helper ───────────────────────────────────────────────────────

/// Write `content` to `<dir>/<name>.fa`, then BGZF-compress it via `pgr fa gz`
/// (which also creates the .gzi index required for random access). Returns
/// the path to the produced `.fa.gz` file.
fn write_bgzf_fa(dir: &std::path::Path, name: &str, content: &str) -> String {
    use std::fs;
    let fa_path = dir.join(format!("{name}.fa"));
    fs::write(&fa_path, content).unwrap();
    let fa_str = fa_path.to_string_lossy().into_owned();
    let (out, _) = PgrCmd::new().args(&["fa", "gz", &fa_str]).run();
    let _ = out;
    let gz_path = format!("{fa_str}.gz");
    assert!(
        std::path::Path::new(&gz_path).exists(),
        "pgr fa gz failed to produce {gz_path}"
    );
    gz_path
}

// ── paf to-maf (pairwise MAF from CIGAR) ─────────────────────────

#[test]
fn command_paf_to_maf_strict_name_validation() {
    use std::fs;
    // PAF references A and B; TSV only has A.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\n")).unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::fs;
    // CIGAR: 4= 3I 3= → target 0-7, query 0-10
    // target: ACGT---ACG  (4 match + 3 gaps + 3 match)
    // query:  ACGTACGTAC  (4 match + 3 bases + 3 match)
    // But query[7..10] should be "TAC" (query = ACGTACGTAC, idx 7,8,9 = T,A,C)
    // and query[4..7] = "ACG"
    // So query alignment = ACGT + ACG + TAC = ACGTACGTAC
    // target alignment = ACGT + --- + ACG = ACGT---ACG
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t7\t7\t10\t255\tcg:Z:4=3I3=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-7",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // target has gaps where query inserted
    assert!(
        stdout.contains("ACGT---ACG"),
        "target alignment should contain gaps for insertion"
    );
    assert!(
        stdout.contains("ACGTACGTAC"),
        "query alignment should contain full sequence"
    );
    // sizes: target 7 non-gap, query 10 non-gap
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
    use std::fs;
    // CIGAR: 4= 3D 3= → target 0-10, query 0-7
    // target: ACGTACGTAC (4 match + 3 bases + 3 match)
    // query:  ACGT---ACG (4 match + 3 gaps + 3 match)
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t7\t0\t7\t+\tB\t10\t0\t10\t7\t10\t255\tcg:Z:4=3D3=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACG\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::fs;
    // Full alignment 10= over B:0-10. Query B:2-8 should trim CIGAR to 6=.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:2-8",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::fs;
    // '-' strand perfect match: target B forward == RC(query A forward).
    // A forward = GTACGTACGT, RC = ACGTACGTAC = B forward.
    // CIGAR 10= describes 10 alignment columns of target vs RC(query).
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(stdout.contains("##maf version=1"), "missing MAF header");
    // Target line: forward strand, original sequence.
    assert!(
        stdout.contains("s\tB\t0\t10\t+\t10\tACGTACGTAC"),
        "missing/incorrect target line for '-' strand record"
    );
    // Query line: '-' strand, displayed sequence is RC of A forward.
    // q_start_maf = srcSize - qe = 10 - 10 = 0; q_size = 10.
    assert!(
        stdout.contains("s\tA\t0\t10\t-\t10\tACGTACGTAC"),
        "missing/incorrect query line for '-' strand record (RC not applied)"
    );
}

#[test]
fn command_paf_to_maf_reverse_strand_with_insertion() {
    use std::fs;
    // '-' strand alignment with insertion: CIGAR 4=3I3= (7 target, 10 query cols).
    // A forward = GTACGTACGT, RC(A) = ACGTACGTAC.
    // target B = ACGT (RC(A)[0:4]) + TAC (RC(A)[7:10]) = ACGTTAC (7 bp).
    // Expected alignment columns:
    //   target: ACGT---TAC  (4 match + 3 gaps + 3 match)
    //   query:  ACGTACGTAC  (RC of A forward, walked left-to-right)
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t-\tB\t7\t0\t7\t7\t7\t255\tcg:Z:4=3I3=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-7",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // Target has gaps where query inserted.
    assert!(
        stdout.contains("ACGT---TAC"),
        "target alignment should contain gaps for insertion on '-' strand"
    );
    // Query alignment is RC of A forward, walked left-to-right.
    assert!(
        stdout.contains("ACGTACGTAC"),
        "query alignment should be RC of A forward on '-' strand"
    );
    // Query line should be '-' strand with q_start = srcSize - qe = 0.
    assert!(
        stdout.contains("s\tA\t0\t10\t-\t10\tACGTACGTAC"),
        "missing/incorrect query s-line for '-' strand with insertion"
    );
}

#[test]
fn command_paf_to_maf_reverse_strand_subinterval_first_half() {
    use std::fs;
    // '-' strand perfect match, sub-interval query on the first half.
    // A forward = GTACGTACGT (G0 T1 A2 C3 G4 T5 A6 C7 G8 T9), RC = ACGTACGTAC = B.
    // PAF CIGAR 10= aligns RC(A) vs B left-to-right.
    // Query B:0-5 → CIGAR first 5 query bases = RC(A)[0:5] = ACGTA, which
    // corresponds to forward A[5:10) = TACGT (RC = ACGTA).
    // Before the project() fix this returned forward A[0:5) = GTACG
    // (RC = CGTAC), which is wrong.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-5",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // Target sub-interval B[0:5] = ACGTA, +strand, start 0, size 5, srcSize 10.
    assert!(
        stdout.contains("s\tB\t0\t5\t+\t10\tACGTA"),
        "missing/incorrect target line for '-' strand sub-interval (first half)"
    );
    // Query: forward A[5:10] RC'd = ACGTA. q_start_maf = srcSize - qe = 10-10 = 0.
    assert!(
        stdout.contains("s\tA\t0\t5\t-\t10\tACGTA"),
        "missing/incorrect query line for '-' strand sub-interval (first half); \
         this verifies project() maps RC offset back to forward [5,10)"
    );
    // Sanity: the buggy forward A[0:5] mapping must NOT appear.
    assert!(
        !stdout.contains("CGTAC"),
        "regression: query sequence looks like RC of forward A[0:5] — project() \
         did not convert RC offset to forward coordinates on '-' strand"
    );
}

#[test]
fn command_paf_to_maf_reverse_strand_subinterval_second_half() {
    use std::fs;
    // Same setup as the first-half test, but query B:5-10.
    // CIGAR last 5 query bases = RC(A)[5:10] = CGTAC, corresponding to
    // forward A[0:5) = GTACG (RC = CGTAC).
    // q_start_maf = srcSize - qe = 10 - 5 = 5.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:5-10",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // Target sub-interval B[5:10] = CGTAC.
    assert!(
        stdout.contains("s\tB\t5\t5\t+\t10\tCGTAC"),
        "missing/incorrect target line for '-' strand sub-interval (second half)"
    );
    // Query: forward A[0:5] RC'd = CGTAC. q_start_maf = srcSize - qe = 10-5 = 5.
    assert!(
        stdout.contains("s\tA\t5\t5\t-\t10\tCGTAC"),
        "missing/incorrect query line for '-' strand sub-interval (second half); \
         this verifies project() maps RC offset back to forward [0,5)"
    );
}

#[test]
fn command_paf_to_maf_reverse_strand_subinterval_with_insertion() {
    use std::fs;
    // '-' strand with insertion, sub-interval query on the trailing target
    // segment (which spans op2 3I + op3 2=).
    // A forward = GTACGTACGT (G0 T1 A2 C3 G4 T5 A6 C7 G8 T9), RC = ACGTACGTAC.
    // CIGAR 5=3I2= aligns RC(A) vs B (7 bp): B[0:5]=ACGTA, B[5:7]=AC, with
    // RC(A)[5:8]=CGT inserted between them.
    // Query B:5-7 → project returns forward A[0,5) (union of op2 RC[5,8)→fwd
    // [2,5) and op3 RC[8,10)→fwd [0,2)).
    // q_seq = RC(A[0:5]) = RC(GTACG) = CGTAC. qs_eff = rec_qe - qe = 10-5 = 5.
    // Alignment columns: op2 3I at ct=5 emits q_seq[0..3]=CGT with target
    // gaps; op3 2= at ct=[5,7) emits q_seq[3..5]=AC paired with B[5:7]=AC.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t-\tB\t7\t0\t7\t7\t7\t255\tcg:Z:5=3I2=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTAAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:5-7",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // Target sub-interval B[5:7] = AC, with 3 gap columns before it from the
    // insertion op (which sits at target position 5, inside [5,7)).
    assert!(
        stdout.contains("s\tB\t5\t2\t+\t7\t---AC"),
        "missing/incorrect target line for '-' strand sub-interval with insertion"
    );
    // Query: RC(A[0:5]) = CGTAC. q_start_maf = srcSize - qe = 10-5 = 5.
    assert!(
        stdout.contains("s\tA\t5\t5\t-\t10\tCGTAC"),
        "missing/incorrect query line for '-' strand sub-interval with insertion"
    );
}

// ── paf to-maf --msa (multi-way MSA via POA) ─────────────────────

#[test]
fn command_paf_to_maf_msa_three_genomes_transitive() {
    use std::fs;
    // Three genomes A/B/C, all 10 bp, identical sequence ACGTACGTAC.
    // A-B and A-C alignments → query B:0-10 with --transitive gathers
    // {B(target), A, C} into one region. --msa merges them into a single
    // 3-sequence MAF block via POA. Since all sequences are identical, the
    // MSA columns should be gap-free and all three `s` lines equal.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "--transitive",
            "--msa",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(
        stderr.contains("Total results:") && !stderr.contains("Total results: 0"),
        "expected non-zero results"
    );
    assert!(stdout.contains("##maf version=1"), "missing MAF header");
    // Exactly one `a` block (multi-way).
    let a_count = stdout.matches("\na\n").count() + if stdout.starts_with("a\n") { 1 } else { 0 };
    assert_eq!(a_count, 1, "expected exactly one MAF block, got {a_count}");
    // Three `s` lines (target B + queries A, C).
    let s_count = stdout.lines().filter(|l| l.starts_with("s\t")).count();
    assert_eq!(s_count, 3, "expected 3 s-lines, got {s_count}");
    // All identical → each s-line should end with ACGTACGTAC (no gaps).
    for line in stdout.lines().filter(|l| l.starts_with("s\t")) {
        assert!(
            line.ends_with("ACGTACGTAC"),
            "expected gap-free ACGTACGTAC in s-line: {line}"
        );
    }
    // Target B should appear first.
    let first_s = stdout.lines().find(|l| l.starts_with("s\t")).unwrap();
    assert!(
        first_s.starts_with("s\tB\t"),
        "expected B as the first s-line, got {first_s}"
    );
}

#[test]
fn command_paf_to_maf_msa_with_snp() {
    use std::fs;
    // B = ACGTACGTAC (target)
    // A = ACGTACGTAC (identical to B)
    // C = ACGTTCGTAC (SNP at position 4: A→T)
    // A-B and A-C alignments, query B:0-10 --transitive --msa.
    // POA should produce a 3-sequence MSA with one SNP column.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nACGTTCGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "--transitive",
            "--msa",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    let s_count = stdout.lines().filter(|l| l.starts_with("s\t")).count();
    assert_eq!(s_count, 3, "expected 3 s-lines, got {s_count}");
    // All three s-lines should have length 10 (no gaps introduced for a SNP).
    for line in stdout.lines().filter(|l| l.starts_with("s\t")) {
        let aln = line.split('\t').next_back().unwrap();
        assert_eq!(aln.len(), 10, "expected 10-char alignment, got '{aln}'");
    }
    // C should differ from B at position 4 (0-indexed).
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
    use std::fs;
    // B = ACGTACGTAC (target, forward)
    // A = GTACGTACGT, aligned to B on '-' strand. RC(A) = ACGTACGTAC = B,
    // so after reverse-complementation A's aligned sequence equals B.
    // Query B:0-10 --transitive --msa: A is RC'd before POA, then both
    // sequences are identical → gap-free MSA.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-maf",
            "stdin",
            "B:0-10",
            "--transitive",
            "--msa",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    let s_count = stdout.lines().filter(|l| l.starts_with("s\t")).count();
    assert_eq!(s_count, 2, "expected 2 s-lines (B + A), got {s_count}");
    // A should be emitted with strand '-'.
    let a_line = stdout.lines().find(|l| l.starts_with("s\tA\t")).unwrap();
    assert!(
        a_line.contains("\t-\t"),
        "A should be on '-' strand: {a_line}"
    );
    // A's aligned sequence should be RC(GTACGTACGT) = ACGTACGTAC, gap-free.
    let a_aln = a_line.split('\t').next_back().unwrap();
    assert_eq!(a_aln, "ACGTACGTAC", "expected RC(A) gap-free: {a_aln}");
}
