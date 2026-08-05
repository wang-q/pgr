#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── paf top-level help ──────────────────────────────────────────

#[test]
fn command_paf_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "--help"]).run();
    assert!(stdout.contains("Manipulates PAF"));
    assert!(stdout.contains("index"));
}

// ── paf index ───────────────────────────────────────────────────

#[test]
fn command_paf_index_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "index", "--help"]).run();
    assert!(stdout.contains("Builds interval-tree index"));
    assert!(stdout.contains("infiles"));
}

#[test]
fn command_paf_index_single_file() {
    let paf = "\
q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255\tcg:Z:50M\tgi:f:0.9
q2\t300\t10\t60\t-\tt1\t200\t10\t60\t45\t50\t255\tcg:Z:50M
";
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", "stdin"])
        .stdin(paf)
        .run();
    assert!(stderr.contains("sequences: 3"));
    assert!(stderr.contains("targets:   1"));
}

#[test]
fn command_paf_index_no_cigar() {
    let paf = "\
q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255
q2\t300\t10\t60\t+\tt2\t400\t10\t60\t45\t50\t255
";
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", "stdin"])
        .stdin(paf)
        .run();
    assert!(stderr.contains("sequences: 4"));
    assert!(stderr.contains("targets:   2"));
}

#[test]
fn command_paf_index_empty() {
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", "stdin"])
        .stdin("")
        .run();
    assert!(stderr.contains("sequences: 0"));
    assert!(stderr.contains("targets:   0"));
}

#[test]
fn command_paf_index_comments_and_blanks() {
    let paf = "\
# header comment

q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255\tcg:Z:50M

# another comment
q2\t300\t10\t60\t-\tt1\t200\t10\t60\t45\t50\t255\tcg:Z:50M
";
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", "stdin"])
        .stdin(paf)
        .run();
    assert!(stderr.contains("sequences: 3"));
    assert!(stderr.contains("targets:   1"));
}

#[test]
fn command_paf_index_invalid() {
    PgrCmd::new()
        .args(&["paf", "index", "stdin"])
        .stdin("invalid line\n")
        .run_fail();
}

#[test]
fn command_paf_index_multiple_files() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let p1 = temp.path().join("multi_a.paf");
    let p2 = temp.path().join("multi_b.paf");
    let idx = temp.path().join("multi.paf.idx");
    fs::write(
        &p1,
        "A\t100\t0\t50\t+\tX\t200\t0\t50\t45\t50\t255\tcg:Z:50M\n",
    )
    .unwrap();
    fs::write(
        &p2,
        "B\t100\t0\t50\t+\tX\t200\t50\t100\t45\t50\t255\tcg:Z:50M\n",
    )
    .unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "index",
            p1.to_str().unwrap(),
            p2.to_str().unwrap(),
            "-o",
            idx.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("Building PAF index from 2 file"));
    assert!(stderr.contains("saved to"));
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", idx.to_str().unwrap(), "X:0-100"])
        .run();
    assert!(stdout.contains("A\t0\t0\t50\t+\tX"), "A not found");
    assert!(stdout.contains("B\t0\t0\t50\t+\tX"), "B not found");
}

// ── persist roundtrip (index save → query load) ─────────────────

#[test]
fn command_paf_index_save_and_query() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("persist.paf");
    let idx_path = temp.path().join("persist.paf.idx");
    fs::write(
        &paf_path,
        "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M
",
    )
    .unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "index",
            paf_path.to_str().unwrap(),
            "-o",
            idx_path.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("saved to"));
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "query", idx_path.to_str().unwrap(), "B:0-100"])
        .run();
    assert!(stderr.contains("Loading index"));
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
    assert!(stdout.contains("C\t0\t0\t50\t+\tB"), "C not found");
}

#[test]
fn command_paf_query_bad_idx_magic() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let bad_path = temp.path().join("bad.paf.idx");
    fs::write(&bad_path, "garbage data\n").unwrap();
    PgrCmd::new()
        .args(&["paf", "query", bad_path.to_str().unwrap(), "B:0-100"])
        .run_fail();
}

#[test]
fn command_paf_query_direct_vs_idx_same_result() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("compare.paf");
    let idx_path = temp.path().join("compare.paf.idx");
    fs::write(
        &paf_path,
        "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M
",
    )
    .unwrap();
    let (direct_out, _) = PgrCmd::new()
        .args(&["paf", "query", paf_path.to_str().unwrap(), "B:0-100"])
        .run();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "index",
            paf_path.to_str().unwrap(),
            "-o",
            idx_path.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("saved to"));
    let (idx_out, stderr) = PgrCmd::new()
        .args(&["paf", "query", idx_path.to_str().unwrap(), "B:0-100"])
        .run();
    assert!(stderr.contains("Loading index"));
    assert_eq!(direct_out, idx_out, "PAF direct vs .idx results differ");
}

#[test]
fn command_paf_query_transitive_from_idx() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("bfs_idx.paf");
    let idx_path = temp.path().join("bfs_idx.paf.idx");
    fs::write(
        &paf_path,
        "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t100\t+\tA\t100\t0\t100\t90\t100\t255\tcg:Z:100M
",
    )
    .unwrap();
    let _ = PgrCmd::new()
        .args(&[
            "paf",
            "index",
            paf_path.to_str().unwrap(),
            "-o",
            idx_path.to_str().unwrap(),
        ])
        .run();
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            idx_path.to_str().unwrap(),
            "B:0-100",
            "--transitive",
        ])
        .run();
    assert!(stderr.contains("Loading index"));
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A (1-hop) not found");
    assert!(stdout.contains("C\t0\t0\t100\t+\tA"), "C (2-hop) not found");
}

// ── -o overwrite protection (data safety) ────────────────────────

#[test]
fn command_paf_index_output_same_as_input_rejected() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("overwrite.paf");
    let original = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    fs::write(&paf_path, original).unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "index",
            paf_path.to_str().unwrap(),
            "-o",
            paf_path.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("is also an input file"),
        "should reject -o == input: {stderr}"
    );
    // Input file must remain untouched.
    let after = fs::read_to_string(&paf_path).unwrap();
    assert_eq!(after, original, "input file must not be overwritten");
}

#[test]
fn command_paf_query_output_same_as_input_rejected() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("overwrite_q.paf");
    let original = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    fs::write(&paf_path, original).unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            paf_path.to_str().unwrap(),
            "B:0-100",
            "-o",
            paf_path.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("is also an input file"),
        "should reject -o == input: {stderr}"
    );
    let after = fs::read_to_string(&paf_path).unwrap();
    assert_eq!(after, original, "input file must not be overwritten");
}

#[test]
fn command_paf_query_min_identity_out_of_range_rejected() {
    // --min-identity is a ratio and must be in [0.0, 1.0]; out-of-range
    // values are rejected at parse time (before any input is read).
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:0-10", "--min-identity", "1.5"])
        .stdin("")
        .run_fail();
    assert!(
        stderr.contains("must be in 0.0..=1.0"),
        "should reject identity > 1.0: {stderr}"
    );
    PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:0-10", "--min-identity", "-0.1"])
        .stdin("")
        .run_fail();
    // In-range values are accepted.
    PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:0-10", "--min-identity", "0.5"])
        .stdin("")
        .run();
}

#[test]
fn command_paf_poa_score_params_validated() {
    // --match must be a positive reward; penalties must be non-positive.
    // POA args live on the to-* commands (here to-vcf); validated at parse time.
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "to-vcf", "stdin", "B:0-10", "--match", "0"])
        .stdin("")
        .run_fail();
    assert!(
        stderr.contains("--match must be > 0"),
        "should reject match=0: {stderr}"
    );
    PgrCmd::new()
        .args(&["paf", "to-vcf", "stdin", "B:0-10", "--match", "-5"])
        .stdin("")
        .run_fail();
    for (flag, bad) in [
        ("--mismatch", "1"),
        ("--gap-open", "2"),
        ("--gap-extend", "3"),
    ] {
        let (_, stderr) = PgrCmd::new()
            .args(&["paf", "to-vcf", "stdin", "B:0-10", flag, bad])
            .stdin("")
            .run_fail();
        assert!(
            stderr.contains("must be <= 0"),
            "should reject {flag}={bad}: {stderr}"
        );
    }
    // Valid values are accepted (the nonexistent fasta-tsv then fails later,
    // proving the score params passed parse-time validation).
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-vcf",
            "stdin",
            "B:0-10",
            "--match",
            "3",
            "--mismatch",
            "-2",
            "--gap-open",
            "-4",
            "--gap-extend",
            "-1",
            "-f",
            "tests/paf/input/nonexistent.tsv",
        ])
        .stdin("")
        .run_fail();
    assert!(
        !stderr.contains("--match must be > 0") && !stderr.contains("must be <= 0"),
        "valid score params must pass parse-time validation: {stderr}"
    );
}

#[test]
fn command_paf_negative_query_filters_rejected() {
    // Distance/length filters are non-negative. The space-separated `-5` is
    // rejected by clap as an unknown flag; the `=` form reaches the value
    // parser and would otherwise be silently treated as "off".
    for flag in [
        "--min-dist",
        "--min-output-len",
        "--merge-distance",
        "--min-chain-length",
    ] {
        let (_, stderr) = PgrCmd::new()
            .args(&["paf", "query", "stdin", "B:0-10", flag, "-5"])
            .stdin("")
            .run_fail();
        assert!(
            stderr.contains("must be >= 0") || stderr.contains("unexpected argument"),
            "should reject {flag}=-5: {stderr}"
        );
        // The `=` form reaches the value parser and must be rejected there.
        let (_, stderr) = PgrCmd::new()
            .args(&["paf", "query", "stdin", "B:0-10", &format!("{flag}=-5")])
            .stdin("")
            .run_fail();
        assert!(
            stderr.contains("must be >= 0"),
            "should reject {flag}=-5 via = syntax: {stderr}"
        );
    }
}

#[test]
fn command_paf_to_bed_output_same_as_input_rejected() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("overwrite_bed.paf");
    let original = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    fs::write(&paf_path, original).unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-bed",
            paf_path.to_str().unwrap(),
            "B:0-100",
            "-o",
            paf_path.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("is also an input file"),
        "should reject -o == input: {stderr}"
    );
    let after = fs::read_to_string(&paf_path).unwrap();
    assert_eq!(after, original, "input file must not be overwritten");
}

#[test]
fn command_paf_graph_output_same_as_input_rejected() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("overwrite_graph.paf");
    let original = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    fs::write(&paf_path, original).unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "graph",
            paf_path.to_str().unwrap(),
            "-o",
            paf_path.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("is also an input file"),
        "should reject -o == input: {stderr}"
    );
    let after = fs::read_to_string(&paf_path).unwrap();
    assert_eq!(after, original, "input file must not be overwritten");
}
