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

// ── paf to-vcf ───────────────────────────────────────────────────

#[test]
fn command_paf_to_vcf_with_snp() {
    // B = ACGTACGTAC (target, REF)
    // A = ACGTACGTAC (identical to B)
    // C = ACGTTCGTAC (SNP at position 4 (0-indexed): A->T)
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("ABC_snp.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    assert!(
        lines.iter().any(|l| l.starts_with("##fileformat=VCFv4.2")),
        "missing VCF fileformat header: {stdout}"
    );
    let header = lines
        .iter()
        .find(|l| l.starts_with("#CHROM"))
        .expect("missing #CHROM header");
    assert!(
        header.contains("\tB\tA\tC"),
        "sample columns should be B A C (target first): {header}"
    );

    let body: Vec<&str> = lines
        .iter()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .copied()
        .collect();
    assert_eq!(
        body.len(),
        1,
        "expected 1 variant row, got {}: {body:?}",
        body.len()
    );
    let fields: Vec<&str> = body[0].split('\t').collect();
    assert_eq!(fields[0], "B", "CHROM");
    assert_eq!(fields[1], "5", "POS (1-based)");
    assert_eq!(fields[3], "A", "REF");
    assert_eq!(fields[4], "T", "ALT");
    assert_eq!(fields[8], "GT", "FORMAT");
    assert_eq!(fields.len(), 12, "8 fixed + 3 samples = 12 columns");
    assert_eq!(fields[9], "0", "B (target=REF) -> GT 0");
    assert_eq!(fields[10], "0", "A (identical to REF) -> GT 0");
    assert_eq!(fields[11], "1", "C (ALT T) -> GT 1");
}

#[test]
fn command_paf_to_vcf_no_variant() {
    // All three genomes identical -> no substitution -> body empty (header only).
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("ABC.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let body: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert!(
        body.is_empty(),
        "expected no variants for identical sequences, got: {body:?}"
    );
}

#[test]
fn command_paf_to_vcf_with_del() {
    // B = ACGTACGTAC (target, REF, 10bp)
    // A = ACGTACGTAC (identical to B)
    // C = ACGTCGTAC  (9bp, deletion of "A" at B's position 4)
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t9\t0\t9\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:4=1D5=\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("ABC_del.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let body: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(body.len(), 1, "expected 1 DEL variant, got: {body:?}");
    let fields: Vec<&str> = body[0].split('\t').collect();
    assert_eq!(fields[0], "B", "CHROM");
    assert_eq!(fields[1], "4", "POS (1-based, anchor col)");
    assert_eq!(fields[3], "TA", "REF (anchor + deleted base)");
    assert_eq!(fields[4], "T", "ALT (anchor only = deletion)");
    assert_eq!(fields[8], "GT", "FORMAT");
    assert_eq!(fields[9], "0", "B (target=REF) -> GT 0");
    assert_eq!(fields[10], "0", "A (identical to REF) -> GT 0");
    assert_eq!(fields[11], "1", "C (deletion) -> GT 1");
}

#[test]
fn command_paf_to_vcf_with_ins() {
    // B = ACGTACGTAC  (target, REF, 10bp)
    // A = ACGTACGTAC  (identical to B)
    // C = ACGTAGCGTAC (11bp, insertion of "G" after B's position 4)
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t11\t0\t11\t+\tA\t10\t0\t10\t10\t11\t255\tcg:Z:5=1I5=\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("ABC_ins.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let body: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(body.len(), 1, "expected 1 INS variant, got: {body:?}");
    let fields: Vec<&str> = body[0].split('\t').collect();
    assert_eq!(fields[0], "B", "CHROM");
    assert_eq!(fields[1], "5", "POS (1-based, anchor col)");
    assert_eq!(fields[3], "A", "REF (anchor only)");
    assert_eq!(fields[4], "AG", "ALT (anchor + inserted base)");
    assert_eq!(fields[8], "GT", "FORMAT");
    assert_eq!(fields[9], "0", "B (target=REF) -> GT 0");
    assert_eq!(fields[10], "0", "A (identical to REF) -> GT 0");
    assert_eq!(fields[11], "1", "C (insertion) -> GT 1");
}

#[test]
fn command_paf_to_vcf_left_align_ins() {
    // B = GACTTTTTTTTCAC  (target, REF, 14bp)
    // A = GACTTTTTTTTCAC  (identical to B)
    // C = GACTTTTTTTTTCAC (15bp, extra T inside the T run)
    let paf = "\
A\t14\t0\t14\t+\tB\t14\t0\t14\t14\t14\t255\tcg:Z:14=\n\
C\t15\t0\t15\t+\tA\t14\t0\t14\t14\t15\t255\tcg:Z:11=1I3=\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-14",
            "--transitive",
            "-f",
            fixture("ABC_14ins.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let body: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(body.len(), 1, "expected 1 INS variant, got: {body:?}");
    let fields: Vec<&str> = body[0].split('\t').collect();
    assert_eq!(fields[0], "B", "CHROM");
    assert_eq!(fields[1], "3", "POS left-aligned to T-run boundary");
    assert_eq!(fields[3], "C", "REF (anchor = base before T run)");
    assert_eq!(fields[4], "CT", "ALT (anchor + inserted T)");
    assert_eq!(fields[8], "GT", "FORMAT");
    assert_eq!(fields[9], "0", "B (target=REF) -> GT 0");
    assert_eq!(fields[10], "0", "A (identical to REF) -> GT 0");
    assert_eq!(fields[11], "1", "C (insertion) -> GT 1");
}

#[test]
fn command_paf_to_vcf_left_align_del() {
    // B = GACTTTTTTTTTCAC (target, REF, 15bp)
    // A = GACTTTTTTTTTCAC (identical to B)
    // C = GACTTTTTTTTCAC  (14bp, one fewer T in the T run)
    let paf = "\
A\t15\t0\t15\t+\tB\t15\t0\t15\t15\t15\t255\tcg:Z:15=\n\
C\t14\t0\t14\t+\tA\t15\t0\t15\t14\t15\t255\tcg:Z:11=1D3=\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-15",
            "--transitive",
            "-f",
            fixture("ABC_15del.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let body: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(body.len(), 1, "expected 1 DEL variant, got: {body:?}");
    let fields: Vec<&str> = body[0].split('\t').collect();
    assert_eq!(fields[0], "B", "CHROM");
    assert_eq!(fields[1], "3", "POS left-aligned to T-run boundary");
    assert_eq!(fields[3], "CT", "REF (anchor + deleted T)");
    assert_eq!(fields[4], "C", "ALT (anchor only = deletion)");
    assert_eq!(fields[8], "GT", "FORMAT");
    assert_eq!(fields[9], "0", "B (target=REF) -> GT 0");
    assert_eq!(fields[10], "0", "A (identical to REF) -> GT 0");
    assert_eq!(fields[11], "1", "C (deletion) -> GT 1");
}
