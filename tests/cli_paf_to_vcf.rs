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

// ── paf to-vcf ───────────────────────────────────────────────────

#[test]
fn command_paf_to_vcf_with_snp() {
    use std::fs;
    // B = ACGTACGTAC (target, REF)
    // A = ACGTACGTAC (identical to B)
    // C = ACGTTCGTAC (SNP at position 4 (0-indexed): A->T)
    // A-B and A-C alignments, query B:0-10 --transitive.
    // VCF should emit one row: CHROM=B, POS=5 (1-based), REF=A, ALT=T,
    // GT: B=0, A=0, C=1.
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
            "to-vcf",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    // Header lines.
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

    // Body rows: exactly one SNP at pos 5 (1-based), REF=A, ALT=T.
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
    // FORMAT = GT, then 3 samples in order B, A, C.
    assert_eq!(fields[8], "GT", "FORMAT");
    assert_eq!(fields.len(), 12, "8 fixed + 3 samples = 12 columns");
    assert_eq!(fields[9], "0", "B (target=REF) -> GT 0");
    assert_eq!(fields[10], "0", "A (identical to REF) -> GT 0");
    assert_eq!(fields[11], "1", "C (ALT T) -> GT 1");
}

#[test]
fn command_paf_to_vcf_no_variant() {
    use std::fs;
    // All three genomes identical -> no substitution -> body empty (header only).
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

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::fs;
    // B = ACGTACGTAC (target, REF, 10bp)
    // A = ACGTACGTAC (identical to B)
    // C = ACGTCGTAC  (9bp, deletion of "A" at B's position 4)
    // A-B full align; C-A aligns with 1bp deletion in C.
    // POA MSA: B/A = ACGTACGTAC, C = ACGT-CGTAC (gap at col 4).
    // DEL variant: anchor=T (col 3), REF="TA", ALT="T", POS=4 (1-based).
    // GT: B=0, A=0, C=1.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t9\t0\t9\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:4=1D5=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nACGTCGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::fs;
    // B = ACGTACGTAC  (target, REF, 10bp)
    // A = ACGTACGTAC  (identical to B)
    // C = ACGTAGCGTAC (11bp, insertion of "G" after B's position 4)
    // A-B full align; C-A aligns with 1bp insertion in C.
    // POA MSA: B/A = ACGTA-CGTAC, C = ACGTAGCGTAC (gap in target at col 5).
    // INS variant: anchor=A (col 4), REF="A", ALT="AG", POS=5 (1-based).
    // GT: B=0, A=0, C=1.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t11\t0\t11\t+\tA\t10\t0\t10\t10\t11\t255\tcg:Z:5=1I5=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nACGTAGCGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::fs;
    // B = GACTTTTTTTTCAC  (target, REF, 14bp; GAC + TTTTTTTT + CAC)
    // A = GACTTTTTTTTCAC  (identical to B)
    // C = GACTTTTTTTTTCAC (15bp, extra T inside the T run)
    // POA MSA left-aligns the gap to the T-run boundary:
    //   B/A = GAC-TTTTTTTTCAC, C = GACTTTTTTTTTCAC (target gap at col 3).
    // Anchor = C (col 2). left_align_indels cannot shift further left
    //   (preceding base A != inserted base T). So POS=3, REF=C, ALT=CT.
    // GT: B=0, A=0, C=1.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t14\t0\t14\t+\tB\t14\t0\t14\t14\t14\t255\tcg:Z:14=\n\
C\t15\t0\t15\t+\tA\t14\t0\t14\t14\t15\t255\tcg:Z:11=1I3=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nGACTTTTTTTTCAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nGACTTTTTTTTCAC\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nGACTTTTTTTTTCAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-14",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::fs;
    // B = GACTTTTTTTTTCAC (target, REF, 15bp; GAC + TTTTTTTTT + CAC)
    // A = GACTTTTTTTTTCAC (identical to B)
    // C = GACTTTTTTTTCAC  (14bp, one fewer T in the T run)
    // POA MSA left-aligns the gap to the T-run boundary:
    //   B/A = GACTTTTTTTTTCAC, C = GAC-TTTTTTTTCAC (gap in C at col 3).
    // Anchor = C (col 2). left_align_indels cannot shift further left
    //   (preceding base A != deleted base T). So POS=3, REF=CT, ALT=C.
    // GT: B=0, A=0, C=1.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t15\t0\t15\t+\tB\t15\t0\t15\t15\t15\t255\tcg:Z:15=\n\
C\t14\t0\t14\t+\tA\t15\t0\t15\t14\t15\t255\tcg:Z:11=1D3=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nGACTTTTTTTTTCAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nGACTTTTTTTTTCAC\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nGACTTTTTTTTCAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-15",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
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
