#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── paf help ─────────────────────────────────────────────────────

#[test]
fn command_paf_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "--help"]).run();
    assert!(stdout.contains("Manipulate PAF"));
    assert!(stdout.contains("index"));
}

#[test]
fn command_paf_index_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "index", "--help"]).run();
    assert!(stdout.contains("Build interval-tree index"));
    assert!(stdout.contains("infiles"));
}

// ── paf index ────────────────────────────────────────────────────

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

// ── paf query ────────────────────────────────────────────────────

#[test]
fn command_paf_query_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "query", "--help"]).run();
    assert!(stdout.contains("Query PAF index"));
    assert!(stdout.contains("--transitive"));
}

#[test]
fn command_paf_query_basic() {
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M
";
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:0-100"])
        .stdin(paf)
        .run();
    assert!(stdout.contains("A\t0\t100\tB\t0\t100"));
    assert!(stdout.contains("C\t0\t50\tB\t50\t100"));
}

#[test]
fn command_paf_query_transitive() {
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t100\t+\tA\t100\t0\t100\t90\t100\t255\tcg:Z:100M
";
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:0-100", "--transitive"])
        .stdin(paf)
        .run();
    assert!(stdout.contains("A\t0\t100\tB"), "A (1-hop) not found");
    assert!(stdout.contains("C\t0\t100\tA"), "C (2-hop) not found");
}

#[test]
fn command_paf_query_not_found() {
    let paf = "\
A\t100\t0\t50\t+\tB\t100\t0\t50\t45\t50\t255\tcg:Z:50M
";
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:100-200"])
        .stdin(paf)
        .run();
    assert!(stderr.contains("No results found"));
}

#[test]
fn command_paf_query_bad_region() {
    PgrCmd::new()
        .args(&["paf", "query", "stdin", "bad_region"])
        .stdin("A\t100\t0\t50\t+\tB\t100\t0\t50\t45\t50\t255\n")
        .run_fail();
}

#[test]
fn command_paf_query_missing_target() {
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "Z:0-100"])
        .stdin("A\t100\t0\t50\t+\tB\t100\t0\t50\t45\t50\t255\n")
        .run();
    assert!(stderr.contains("not found"));
}

// ── persist roundtrip ────────────────────────────────────────────

#[test]
fn command_paf_index_save_and_query() {
    use std::fs;
    let paf_path = "/tmp/pgr_cli_test_persist.paf";
    let idx_path = "/tmp/pgr_cli_test_persist.paf.idx";

    fs::write(
        paf_path,
        "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M
",
    )
    .unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", paf_path, "-o", idx_path])
        .run();
    assert!(stderr.contains("saved to"), "index save failed");

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "query", idx_path, "B:0-100"])
        .run();
    assert!(stderr.contains("Loading index"), "should load from idx");
    assert!(stdout.contains("A\t0\t100\tB\t0\t100"), "A not found");
    assert!(stdout.contains("C\t0\t50\tB\t50\t100"), "C not found");

    let _ = fs::remove_file(paf_path);
    let _ = fs::remove_file(idx_path);
}

#[test]
fn command_paf_query_bad_idx_magic() {
    use std::fs;
    let bad_path = "/tmp/pgr_cli_test_bad.paf.idx";
    fs::write(bad_path, "garbage data\n").unwrap();

    PgrCmd::new()
        .args(&["paf", "query", bad_path, "B:0-100"])
        .run_fail();

    let _ = fs::remove_file(bad_path);
}

#[test]
fn command_paf_query_direct_vs_idx_same_result() {
    use std::fs;
    let paf_path = "/tmp/pgr_cli_test_compare.paf";
    let idx_path = "/tmp/pgr_cli_test_compare.paf.idx";

    fs::write(
        paf_path,
        "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M
",
    )
    .unwrap();

    // Query directly from PAF
    let (direct_out, _) = PgrCmd::new()
        .args(&["paf", "query", paf_path, "B:0-100"])
        .run();

    // Query from saved index
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", paf_path, "-o", idx_path])
        .run();
    assert!(stderr.contains("saved to"));

    let (idx_out, stderr) = PgrCmd::new()
        .args(&["paf", "query", idx_path, "B:0-100"])
        .run();
    assert!(stderr.contains("Loading index"));

    // Results must be identical
    assert_eq!(direct_out, idx_out, "PAF direct vs .idx results differ");

    let _ = fs::remove_file(paf_path);
    let _ = fs::remove_file(idx_path);
}

#[test]
fn command_paf_query_transitive_from_idx() {
    use std::fs;
    let paf_path = "/tmp/pgr_cli_test_bfs_idx.paf";
    let idx_path = "/tmp/pgr_cli_test_bfs_idx.paf.idx";

    fs::write(
        paf_path,
        "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t100\t+\tA\t100\t0\t100\t90\t100\t255\tcg:Z:100M
",
    )
    .unwrap();

    let _ = PgrCmd::new()
        .args(&["paf", "index", paf_path, "-o", idx_path])
        .run();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "query", idx_path, "B:0-100", "--transitive"])
        .run();
    assert!(stderr.contains("Loading index"));
    assert!(stdout.contains("A\t0\t100\tB"), "A (1-hop) not found");
    assert!(stdout.contains("C\t0\t100\tA"), "C (2-hop) not found");

    let _ = fs::remove_file(paf_path);
    let _ = fs::remove_file(idx_path);
}

// ── real-data validation ────────────────────────────────────────

#[test]
fn command_maf_to_paf_real_multiz_spar() {
    let (stdout, _) = PgrCmd::new()
        .args(&["maf", "to-paf", "tests/multiz/S288cvsSpar.maf"])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields[0], "Spar.gi_29362594", "query name");
    assert_eq!(fields[4], "-", "reverse strand");
    assert_eq!(fields[5], "S288c.I", "target name");
    assert_eq!(fields[11], "255", "mapq");
}

#[test]
fn command_maf_to_paf_real_multiz_rm11() {
    let (stdout, _) = PgrCmd::new()
        .args(&["maf", "to-paf", "tests/multiz/S288cvsRM11_1a.maf"])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields[0], "RM11_1a.scaffold_17", "query name");
    assert_eq!(fields[4], "+", "forward strand");
    assert_eq!(fields[5], "S288c.I", "target name");
    assert_eq!(fields[9], "456", "all 456 bases match");
}

#[test]
fn command_paf_query_real_multiz_transitive() {
    use std::fs;
    use std::process::Command;
    let paf_path = "/tmp/pgr_real_test_merged.paf";
    let idx_path = "/tmp/pgr_real_test_merged.paf.idx";

    // Build merged PAF from two MAF files
    let spar_out = Command::new(std::env::current_dir().unwrap().join("target/debug/pgr"))
        .args(["maf", "to-paf", "tests/multiz/S288cvsSpar.maf"])
        .output()
        .unwrap();
    let rm11_out = Command::new(std::env::current_dir().unwrap().join("target/debug/pgr"))
        .args(["maf", "to-paf", "tests/multiz/S288cvsRM11_1a.maf"])
        .output()
        .unwrap();
    let mut merged = spar_out.stdout.clone();
    merged.extend_from_slice(&rm11_out.stdout);
    fs::write(paf_path, &merged).unwrap();

    // Index
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", paf_path, "-o", idx_path])
        .run();
    assert!(stderr.contains("saved to"));

    // Transitive query
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            idx_path,
            "S288c.I:26000-30000",
            "--transitive",
        ])
        .run();
    assert!(stderr.contains("Loading index"));
    assert!(stdout.contains("Spar.gi_29362594"), "Spar not found");
    assert!(stdout.contains("RM11_1a.scaffold_17"), "RM11 not found");

    let _ = fs::remove_file(paf_path);
    let _ = fs::remove_file(idx_path);
}

#[test]
fn command_paf_index_multiple_files() {
    use std::fs;
    let p1 = "/tmp/pgr_multi_a.paf";
    let p2 = "/tmp/pgr_multi_b.paf";
    let idx = "/tmp/pgr_multi.paf.idx";
    fs::write(
        p1,
        "A\t100\t0\t50\t+\tX\t200\t0\t50\t45\t50\t255\tcg:Z:50M\n",
    )
    .unwrap();
    fs::write(
        p2,
        "B\t100\t0\t50\t+\tY\t200\t0\t50\t45\t50\t255\tcg:Z:50M\n",
    )
    .unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", p1, p2, "-o", idx])
        .run();
    assert!(stderr.contains("saved to"));
    // V1: multi-file processes sequentially; last file wins when -o is used.
    // Query Y (from p2, the last file) — must be in the saved index.
    let (stdout, _) = PgrCmd::new().args(&["paf", "query", idx, "Y:0-50"]).run();
    assert!(
        stdout.contains("B\t0\t50\tY"),
        "B (from p2) should be in saved index"
    );
    let _ = fs::remove_file(p1);
    let _ = fs::remove_file(p2);
    let _ = fs::remove_file(idx);
}

#[test]
fn command_paf_query_max_depth_1() {
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\nC\t100\t0\t100\t+\tA\t100\t0\t100\t90\t100\t255\tcg:Z:100M\n";
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "--transitive",
            "--max-depth",
            "1",
        ])
        .stdin(paf)
        .run();
    assert!(stdout.contains("A\t0\t100\tB"), "A (1-hop) not found");
    assert!(
        !stdout.contains("C\t"),
        "C should NOT appear: max-depth=1 stops before 2nd hop"
    );
}
