#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

/// Return a fixture path under `tests/genome`.
fn genome(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/genome")
        .join(name)
}

/// Return a fixture path under `tests/paf/input`.
fn paf_input(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/paf/input")
        .join(name)
}

/// Extract the invariant PAF columns — `qname qlen qstart qend strand tname
/// tlen tstart tend mapq` (indices 0-8 and 11).
///
/// These fields come directly from the `.1aln` `A`/`R` lines and the GDB
/// skeleton, so pgr reproduces them identically to FastGA `ALNtoPAF`. The
/// CIGAR-derived fields (`matches`, `block`, `dv:f`, `df:i`, `cg:Z`) are
/// *excluded* because FastGA recomputes them via its own GREEDIEST box DP +
/// `Gap_Improver`, which is not uniquely determined by the trace points.
fn invariant_cols(line: &str) -> Vec<&str> {
    let f: Vec<&str> = line.split('\t').collect();
    vec![f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[11]]
}

// ── read side: pgr 1aln stat / to-paf / to-psl ───────────────────

#[test]
fn command_1aln_stat_reports_header_and_records() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "1aln",
            "stat",
            genome("mg1655-sakai.1aln").to_str().unwrap(),
        ])
        .run();
    for expected in [
        "tspace\t100",
        "records\t700",
        "skeletons\t2",
        "scaffolds\t4",
        "contigs\t4",
        "refs\t3",
    ] {
        assert!(
            stdout.contains(expected),
            "missing {expected:?} in stat:\n{stdout}"
        );
    }
}

#[test]
fn command_1aln_to_paf_matches_fastga_coordinates() {
    let temp = TempDir::new().unwrap();
    let out = temp.path().join("out.paf");
    PgrCmd::new()
        .args(&[
            "1aln",
            "to-paf",
            genome("mg1655-sakai.1aln").to_str().unwrap(),
            "--ref-seq",
            genome("mg1655.fa.gz").to_str().unwrap(),
            "--query-seq",
            genome("sakai.fa.gz").to_str().unwrap(),
            "--cigar",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();

    let golden = fs::read_to_string(genome("mg1655-sakai.expected.paf")).unwrap();
    let pgr = fs::read_to_string(&out).unwrap();
    let golden_lines: Vec<&str> = golden.lines().filter(|l| !l.is_empty()).collect();
    let pgr_lines: Vec<String> = pgr
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| invariant_cols(l).join("\t"))
        .collect();
    assert_eq!(
        pgr_lines.len(),
        golden_lines.len(),
        "record count differs from FastGA golden"
    );
    for (i, (g, p)) in golden_lines.iter().zip(pgr_lines.iter()).enumerate() {
        assert_eq!(p, g, "record {i} invariant fields differ from FastGA");
    }
}

#[test]
fn command_1aln_to_psl_produces_one_record_per_alignment() {
    let temp = TempDir::new().unwrap();
    let out = temp.path().join("out.psl");
    PgrCmd::new()
        .args(&[
            "1aln",
            "to-psl",
            genome("mg1655-sakai.1aln").to_str().unwrap(),
            "--ref-seq",
            genome("mg1655.fa.gz").to_str().unwrap(),
            "--query-seq",
            genome("sakai.fa.gz").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    let content = fs::read_to_string(&out).unwrap();
    let n = content.lines().filter(|l| !l.is_empty()).count();
    assert_eq!(n, 700, "expected one PSL record per alignment");
    assert!(
        content.contains("NC_000913"),
        "missing reference name in PSL"
    );
}

// ── write side: pgr paf to-1aln / pgr maf to-1aln round trips ────

#[test]
fn command_paf_to_1aln_round_trips_coordinates() {
    let temp = TempDir::new().unwrap();
    let paf = temp.path().join("in.paf");
    let aln = temp.path().join("out.1aln");
    let back = temp.path().join("back.paf");
    let input = concat!(
        "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n",
        "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n",
        "A\t10\t0\t10\t+\tB\t10\t0\t7\t7\t10\t255\tcg:Z:4=3I3=\n",
    );
    fs::write(&paf, input).unwrap();

    PgrCmd::new()
        .args(&[
            "paf",
            "to-1aln",
            paf.to_str().unwrap(),
            "-o",
            aln.to_str().unwrap(),
        ])
        .run();
    assert!(aln.exists(), "paf to-1aln did not produce an output file");

    PgrCmd::new()
        .args(&[
            "1aln",
            "to-paf",
            aln.to_str().unwrap(),
            "--ref-seq",
            paf_input("A.fa.gz").to_str().unwrap(),
            "--query-seq",
            paf_input("B.fa.gz").to_str().unwrap(),
            "-o",
            back.to_str().unwrap(),
        ])
        .run();

    let lines: Vec<String> = fs::read_to_string(&back)
        .unwrap()
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| invariant_cols(l).join("\t"))
        .collect();
    let expected = vec![
        "A\t10\t0\t10\t+\tB\t10\t0\t10\t255",
        "A\t10\t0\t10\t-\tB\t10\t0\t10\t255",
        "A\t10\t0\t10\t+\tB\t10\t0\t7\t255",
    ];
    assert_eq!(lines, expected);
}

#[test]
fn command_maf_to_1aln_round_trips_coordinates() {
    let temp = TempDir::new().unwrap();
    let maf = temp.path().join("in.maf");
    let aln = temp.path().join("out.1aln");
    let back = temp.path().join("back.paf");
    let input = concat!(
        "##maf version=1\n",
        "\n",
        "a score=1\n",
        "s A 0 10 + 10 ACGTACGTAC\n",
        "s B 0 10 + 10 ACGTACGTAC\n",
        "\n",
        "a score=1\n",
        "s A 0 10 + 10 ACGTACGTAC\n",
        "s B 0 8 + 10 ACGTAC--AC\n",
        "\n",
        "a score=1\n",
        "s A 0 10 + 10 ACGTACGTAC\n",
        "s B 0 10 - 10 ACGTACGTAC\n",
    );
    fs::write(&maf, input).unwrap();

    PgrCmd::new()
        .args(&[
            "maf",
            "to-1aln",
            maf.to_str().unwrap(),
            "-o",
            aln.to_str().unwrap(),
        ])
        .run();
    assert!(aln.exists(), "maf to-1aln did not produce an output file");

    PgrCmd::new()
        .args(&[
            "1aln",
            "to-paf",
            aln.to_str().unwrap(),
            "--ref-seq",
            paf_input("A.fa.gz").to_str().unwrap(),
            "--query-seq",
            paf_input("B.fa.gz").to_str().unwrap(),
            "-o",
            back.to_str().unwrap(),
        ])
        .run();

    let lines: Vec<String> = fs::read_to_string(&back)
        .unwrap()
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| invariant_cols(l).join("\t"))
        .collect();
    let expected = vec![
        "A\t10\t0\t10\t+\tB\t10\t0\t10\t255",
        "A\t10\t0\t10\t+\tB\t10\t0\t8\t255",
        "A\t10\t0\t10\t-\tB\t10\t0\t10\t255",
    ];
    assert_eq!(lines, expected);
}

// ── write side error handling ─────────────────────────────────────

#[test]
fn command_paf_to_1aln_requires_cigar() {
    // A PAF record without cg:Z must fail with a friendly error, not panic.
    let temp = TempDir::new().unwrap();
    let paf = temp.path().join("in.paf");
    let aln = temp.path().join("out.1aln");
    fs::write(&paf, "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\n").unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-1aln",
            paf.to_str().unwrap(),
            "-o",
            aln.to_str().unwrap(),
        ])
        .run_fail();
    assert!(stderr.contains("missing a cg:Z CIGAR"), "stderr: {stderr}");
}

#[test]
fn command_paf_to_1aln_rejects_stdout() {
    // The ONEcode container is binary, so stdout is rejected.
    let temp = TempDir::new().unwrap();
    let paf = temp.path().join("in.paf");
    fs::write(
        &paf,
        "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n",
    )
    .unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "to-1aln", paf.to_str().unwrap(), "-o", "stdout"])
        .run_fail();
    assert!(
        stderr.contains("requires a real output file"),
        "stderr: {stderr}"
    );
}
