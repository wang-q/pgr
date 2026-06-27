#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── paf top-level help ──────────────────────────────────────────

#[test]
fn command_paf_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "--help"]).run();
    assert!(stdout.contains("Manipulate PAF"));
    assert!(stdout.contains("index"));
}

// ── paf index ───────────────────────────────────────────────────

#[test]
fn command_paf_index_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "index", "--help"]).run();
    assert!(stdout.contains("Build interval-tree index"));
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
        "B\t100\t0\t50\t+\tX\t200\t50\t100\t45\t50\t255\tcg:Z:50M\n",
    )
    .unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", p1, p2, "-o", idx])
        .run();
    assert!(stderr.contains("Building PAF index from 2 file"));
    assert!(stderr.contains("saved to"));
    let (stdout, _) = PgrCmd::new().args(&["paf", "query", idx, "X:0-100"]).run();
    assert!(stdout.contains("A\t0\t0\t50\t+\tX"), "A not found");
    assert!(stdout.contains("B\t0\t0\t50\t+\tX"), "B not found");
    let _ = fs::remove_file(p1);
    let _ = fs::remove_file(p2);
    let _ = fs::remove_file(idx);
}

// ── persist roundtrip (index save → query load) ─────────────────

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
    assert!(stderr.contains("saved to"));
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "query", idx_path, "B:0-100"])
        .run();
    assert!(stderr.contains("Loading index"));
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
    assert!(stdout.contains("C\t0\t0\t50\t+\tB"), "C not found");
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
    let (direct_out, _) = PgrCmd::new()
        .args(&["paf", "query", paf_path, "B:0-100"])
        .run();
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", paf_path, "-o", idx_path])
        .run();
    assert!(stderr.contains("saved to"));
    let (idx_out, stderr) = PgrCmd::new()
        .args(&["paf", "query", idx_path, "B:0-100"])
        .run();
    assert!(stderr.contains("Loading index"));
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
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A (1-hop) not found");
    assert!(stdout.contains("C\t0\t0\t100\t+\tA"), "C (2-hop) not found");
    let _ = fs::remove_file(paf_path);
    let _ = fs::remove_file(idx_path);
}
