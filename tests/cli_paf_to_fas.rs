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

// ── paf to-fas (pairwise FAS from CIGAR) ─────────────────────────

#[test]
fn command_paf_to_fas_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "to-fas", "--help"]).run();
    assert!(stdout.contains("block FASTA"));
    assert!(stdout.contains("--fasta-tsv"));
}

#[test]
fn command_paf_to_fas_strict_name_validation() {
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
            "paf", "to-fas", "stdin", "B:0-10", "-f", tsv.to_str().unwrap(),
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
            "paf", "to-fas", "stdin", "B:0-10", "-f", tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // Target header first, query header second.
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
    use std::fs;
    // CIGAR: 4= 3I 3= → target 0-7, query 0-10
    // target: ACGT---ACG, query: ACGTACGTAC
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t7\t7\t10\t255\tcg:Z:4=3I3=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf", "to-fas", "stdin", "B:0-7", "-f", tsv.to_str().unwrap(),
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
}

#[test]
fn command_paf_to_fas_with_deletion() {
    use std::fs;
    // CIGAR: 4= 3D 3= → target 0-10, query 0-7
    // target: ACGTACGTAC, query: ACGT---ACG
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t7\t0\t7\t+\tB\t10\t0\t10\t7\t10\t255\tcg:Z:4=3D3=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACG\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf", "to-fas", "stdin", "B:0-10", "-f", tsv.to_str().unwrap(),
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
    use std::fs;
    // '-' strand perfect match: target B forward == RC(query A forward).
    // A forward = GTACGTACGT, RC = ACGTACGTAC = B forward.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf", "to-fas", "stdin", "B:0-10", "-f", tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // Target forward strand.
    assert!(
        stdout.contains(">B(+):1-10\nACGTACGTAC"),
        "missing/incorrect target record for '-' strand"
    );
    // Query '-' strand, displayed sequence is RC of A forward.
    assert!(
        stdout.contains(">A(-):1-10\nACGTACGTAC"),
        "missing/incorrect query record for '-' strand (RC not applied)"
    );
}

// ── paf to-fas --msa (multi-way MSA via POA) ─────────────────────

#[test]
fn command_paf_to_fas_msa_three_genomes_transitive() {
    use std::fs;
    // Three genomes A/B/C, all 10 bp, identical sequence ACGTACGTAC.
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
            "paf", "to-fas", "stdin", "B:0-10", "-t", "--msa", "-f", tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(
        stderr.contains("Total results:") && !stderr.contains("Total results: 0"),
        "expected non-zero results"
    );
    // Three records (target B + queries A, C), all identical.
    let records: Vec<&str> = stdout
        .lines()
        .filter(|l| l.starts_with('>'))
        .collect();
    assert_eq!(records.len(), 3, "expected 3 records, got {records:?}");
    // Target B should appear first.
    assert!(
        records[0].contains(">B("),
        "expected B as the first record, got {records:?}"
    );
    // All sequences should be gap-free ACGTACGTAC.
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
    use std::fs;
    // B = ACGTACGTAC (target), A = ACGTACGTAC, C = ACGTTCGTAC (SNP at pos 4).
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
            "paf", "to-fas", "stdin", "B:0-10", "-t", "--msa", "-f", tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    let records: Vec<&str> = stdout
        .lines()
        .filter(|l| l.starts_with('>'))
        .collect();
    assert_eq!(records.len(), 3, "expected 3 records");
    // All alignment lines should have length 10 (no gaps for a SNP).
    let seqs: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('>'))
        .collect();
    for s in &seqs {
        assert_eq!(s.len(), 10, "expected 10-char alignment, got '{s}'");
    }
    // C should differ from B at position 4 (0-indexed).
    let b_idx = records
        .iter()
        .position(|r| r.contains(">B("))
        .unwrap();
    let c_idx = records
        .iter()
        .position(|r| r.contains(">C("))
        .unwrap();
    let b_aln = seqs[b_idx];
    let c_aln = seqs[c_idx];
    let diffs: Vec<usize> = b_aln
        .chars()
        .zip(c_aln.chars())
        .enumerate()
        .filter_map(|(i, (a, b))| if a != b { Some(i) } else { None })
        .collect();
    assert_eq!(diffs, vec![4], "expected single SNP at pos 4, got {diffs:?}");
}

#[test]
fn command_paf_to_fas_msa_reverse_strand_query() {
    use std::fs;
    // B = ACGTACGTAC (target, forward)
    // A = GTACGTACGT, aligned to B on '-' strand. RC(A) = ACGTACGTAC = B.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf", "to-fas", "stdin", "B:0-10", "-t", "--msa", "-f", tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    let records: Vec<&str> = stdout
        .lines()
        .filter(|l| l.starts_with('>'))
        .collect();
    assert_eq!(records.len(), 2, "expected 2 records (B + A)");
    // A should be emitted with strand '-'.
    assert!(
        records.iter().any(|r| r.contains(">A(-):")),
        "A should be on '-' strand: {records:?}"
    );
    // A's aligned sequence should be RC(GTACGTACGT) = ACGTACGTAC, gap-free.
    let a_idx = records
        .iter()
        .position(|r| r.contains(">A("))
        .unwrap();
    let seqs: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('>'))
        .collect();
    assert_eq!(seqs[a_idx], "ACGTACGTAC", "expected RC(A) gap-free");
}

// ── pipeline: paf to-fas | fas to-vcf ────────────────────────────

#[test]
fn command_paf_to_fas_pipeline_to_vcf() {
    use std::fs;
    // B = ACGTACGTAC (target, REF)
    // A = ACGTACGTAC (identical to B)
    // C = ACGTTCGTAC (SNP at position 4: A->T)
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

    // Step 1: paf to-fas --msa → block FASTA
    let (fas_stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf", "to-fas", "stdin", "B:0-10", "-t", "--msa", "-f", tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(
        !fas_stdout.is_empty(),
        "to-fas output should not be empty"
    );

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
    // Should have one SNP row (A->T at position 5, 1-based).
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
