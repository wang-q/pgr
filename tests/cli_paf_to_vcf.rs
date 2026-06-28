#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── helper ───────────────────────────────────────────────────────

/// Write `content` to a plain `path`, then BGZF-compress it via `pgr fa gz`
/// (which also creates the .gzi index required for random access).
fn write_bgzf_fa(path_no_gz: &str, content: &str) -> String {
    use std::fs;
    fs::write(path_no_gz, content).unwrap();
    let (out, _) = PgrCmd::new().args(&["fa", "gz", path_no_gz]).run();
    let _ = out;
    let gz_path = format!("{path_no_gz}.gz");
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
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_vcf_snp_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_vcf_snp_B.fa", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_vcf_snp_C.fa", ">C\nACGTTCGTAC\n");
    let tsv = "/tmp/pgr_vcf_snp.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-vcf", "stdin", "B:0-10", "-t", "-f", tsv])
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

    let _ = fs::remove_file("/tmp/pgr_vcf_snp_A.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_snp_B.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_snp_C.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(&c_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(format!("{c_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
    let _ = fs::remove_file(format!("{c_fa}.loc"));
}

#[test]
fn command_paf_to_vcf_no_variant() {
    use std::fs;
    // All three genomes identical -> no substitution -> body empty (header only).
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_vcf_novar_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_vcf_novar_B.fa", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_vcf_novar_C.fa", ">C\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_vcf_novar.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-vcf", "stdin", "B:0-10", "-t", "-f", tsv])
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

    let _ = fs::remove_file("/tmp/pgr_vcf_novar_A.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_novar_B.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_novar_C.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(&c_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(format!("{c_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
    let _ = fs::remove_file(format!("{c_fa}.loc"));
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
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t9\t0\t9\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:4=1D5=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_vcf_del_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_vcf_del_B.fa", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_vcf_del_C.fa", ">C\nACGTCGTAC\n");
    let tsv = "/tmp/pgr_vcf_del.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-vcf", "stdin", "B:0-10", "-t", "-f", tsv])
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

    let _ = fs::remove_file("/tmp/pgr_vcf_del_A.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_del_B.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_del_C.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(&c_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(format!("{c_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
    let _ = fs::remove_file(format!("{c_fa}.loc"));
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
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t11\t0\t11\t+\tA\t10\t0\t10\t10\t11\t255\tcg:Z:5=1I5=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_vcf_ins_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_vcf_ins_B.fa", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_vcf_ins_C.fa", ">C\nACGTAGCGTAC\n");
    let tsv = "/tmp/pgr_vcf_ins.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-vcf", "stdin", "B:0-10", "-t", "-f", tsv])
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

    let _ = fs::remove_file("/tmp/pgr_vcf_ins_A.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_ins_B.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_ins_C.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(&c_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(format!("{c_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
    let _ = fs::remove_file(format!("{c_fa}.loc"));
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
    let paf = "\
A\t14\t0\t14\t+\tB\t14\t0\t14\t14\t14\t255\tcg:Z:14=\n\
C\t15\t0\t15\t+\tA\t14\t0\t14\t14\t15\t255\tcg:Z:11=1I3=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_vcf_la_ins_A.fa", ">A\nGACTTTTTTTTCAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_vcf_la_ins_B.fa", ">B\nGACTTTTTTTTCAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_vcf_la_ins_C.fa", ">C\nGACTTTTTTTTTCAC\n");
    let tsv = "/tmp/pgr_vcf_la_ins.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-vcf", "stdin", "B:0-14", "-t", "-f", tsv])
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

    let _ = fs::remove_file("/tmp/pgr_vcf_la_ins_A.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_la_ins_B.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_la_ins_C.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(&c_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(format!("{c_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
    let _ = fs::remove_file(format!("{c_fa}.loc"));
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
    let paf = "\
A\t15\t0\t15\t+\tB\t15\t0\t15\t15\t15\t255\tcg:Z:15=\n\
C\t14\t0\t14\t+\tA\t15\t0\t15\t14\t15\t255\tcg:Z:11=1D3=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_vcf_la_del_A.fa", ">A\nGACTTTTTTTTTCAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_vcf_la_del_B.fa", ">B\nGACTTTTTTTTTCAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_vcf_la_del_C.fa", ">C\nGACTTTTTTTTCAC\n");
    let tsv = "/tmp/pgr_vcf_la_del.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-vcf", "stdin", "B:0-15", "-t", "-f", tsv])
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

    let _ = fs::remove_file("/tmp/pgr_vcf_la_del_A.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_la_del_B.fa");
    let _ = fs::remove_file("/tmp/pgr_vcf_la_del_C.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(&c_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(format!("{c_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
    let _ = fs::remove_file(format!("{c_fa}.loc"));
}
