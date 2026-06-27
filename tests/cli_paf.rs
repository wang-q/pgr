#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

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

// ── paf query (PAF output) ───────────────────────────────────────

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
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
    assert!(stdout.contains("C\t0\t0\t50\t+\tB"), "C not found");
    assert!(stdout.contains("gi:f:"), "gi tag missing");
    assert!(stdout.contains("cg:Z:"), "cg tag missing");
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
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A (1-hop) not found");
    assert!(stdout.contains("C\t0\t0\t100\t+\tA"), "C (2-hop) not found");
}

#[test]
fn command_paf_query_not_found() {
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:100-200"])
        .stdin("A\t100\t0\t50\t+\tB\t100\t0\t50\t45\t50\t255\tcg:Z:50M\n")
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
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A (1-hop) not found");
    assert!(!stdout.contains("C\t"), "C should NOT appear: max-depth=1");
}

#[test]
fn command_paf_query_subset_filter() {
    use std::fs;
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\nC\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M\n";
    let list = "/tmp/pgr_subset.txt";
    fs::write(list, "A\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "--subset-sequence-list",
            list,
        ])
        .stdin(paf)
        .run();
    assert!(stdout.contains("A"), "A should be included");
    assert!(!stdout.contains("C"), "C should be excluded");
    let _ = fs::remove_file(list);
}

// ── paf query -o bed (BED3 output) ───────────────────────────────

#[test]
fn command_paf_query_bed_output() {
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M
";
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:0-100", "-o", "bed"])
        .stdin(paf)
        .run();
    // BED3: name start end (tab-separated), no strand/cigar/gi
    let lines: Vec<&str> = stdout.lines().filter(|l| !l.is_empty()).collect();
    assert!(
        lines.iter().all(|l| l.split('\t').count() == 3),
        "BED3 expected"
    );
    assert!(stdout.contains("A\t0\t100"), "A BED3 line missing");
    assert!(stdout.contains("C\t0\t50"), "C BED3 line missing");
    assert!(!stdout.contains("cg:Z:"), "BED should not contain cg tag");
    assert!(!stdout.contains("gi:f:"), "BED should not contain gi tag");
}

#[test]
fn command_paf_query_bed_output_reverse_strand() {
    // Reverse-strand alignment: query coords should still be emitted as (min, max)
    let paf = "A\t100\t0\t100\t-\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:0-100", "-o", "bed"])
        .stdin(paf)
        .run();
    assert!(
        stdout.contains("A\t0\t100"),
        "A BED3 line missing (reverse strand)"
    );
}

// ── paf query -b (batch BED regions) ─────────────────────────────

#[test]
fn command_paf_query_batch_bed() {
    use std::fs;
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tD\t100\t50\t100\t45\t50\t255\tcg:Z:50M
";
    let bed = "/tmp/pgr_batch_regions.bed";
    fs::write(bed, "B\t0\t100\nD\t50\t100\n# comment line\n\n").unwrap();
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "-b", bed, "-o", "bed"])
        .stdin(paf)
        .run();
    assert!(stdout.contains("A\t0\t100"), "A (from region B) missing");
    assert!(stdout.contains("C\t0\t50"), "C (from region D) missing");
    assert!(stderr.contains("Total results: 2"), "total count missing");
    let _ = fs::remove_file(bed);
}

#[test]
fn command_paf_query_batch_bed_skips_unknown_target() {
    use std::fs;
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let bed = "/tmp/pgr_batch_unknown.bed";
    fs::write(bed, "B\t0\t100\nZ\t0\t100\n").unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "-b", bed, "-o", "bed"])
        .stdin(paf)
        .run();
    assert!(
        stderr.contains("not found in index, skipping"),
        "missing skip warning for unknown target"
    );
    let _ = fs::remove_file(bed);
}

#[test]
fn command_paf_query_region_and_bed_mutually_exclusive() {
    use std::fs;
    let bed = "/tmp/pgr_mutex.bed";
    fs::write(bed, "B\t0\t100\n").unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:0-100", "-b", bed])
        .stdin("A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n")
        .run_fail();
    assert!(
        stderr.contains("mutually exclusive"),
        "missing mutual-exclusion error"
    );
    let _ = fs::remove_file(bed);
}

#[test]
fn command_paf_query_no_region_no_bed_fails() {
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "query", "stdin"])
        .stdin("A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n")
        .run_fail();
    assert!(
        stderr.contains("must be provided"),
        "missing required-region error"
    );
}

// ── paf query --min-degree / --min-chain-length ─────────────────

#[test]
fn command_paf_query_min_degree_passes() {
    // 2 distinct queries (A, C) align to B; --min-degree 2 keeps both
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M
";
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "--min-degree",
            "2",
            "-o",
            "bed",
        ])
        .stdin(paf)
        .run();
    assert!(
        stdout.contains("A\t0\t100"),
        "A should be kept (degree 2 == 2)"
    );
    assert!(
        stdout.contains("C\t0\t50"),
        "C should be kept (degree 2 == 2)"
    );
}

#[test]
fn command_paf_query_min_degree_skips_region() {
    // Only 2 distinct queries; --min-degree 3 skips the whole region
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M
";
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "--min-degree",
            "3",
            "-o",
            "bed",
        ])
        .stdin(paf)
        .run();
    assert!(
        stderr.contains("skipped") && stderr.contains("min-degree 3"),
        "missing degree-skip warning"
    );
    assert!(
        !stdout.contains("A\t0\t100") && !stdout.contains("C\t0\t50"),
        "no BED lines should be emitted when region is skipped"
    );
    assert!(
        stderr.contains("No results found"),
        "missing no-results notice"
    );
}

#[test]
fn command_paf_query_min_chain_length_filters_short() {
    // A: 100bp chain; C: 30bp chain. --min-chain-length 50 drops C
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t30\t+\tB\t100\t0\t30\t25\t30\t255\tcg:Z:30M
";
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "--min-chain-length",
            "50",
            "-o",
            "bed",
        ])
        .stdin(paf)
        .run();
    assert!(
        stdout.contains("A\t0\t100"),
        "A (100bp >= 50) should be kept"
    );
    assert!(
        !stdout.contains("C\t"),
        "C (30bp < 50) should be filtered out"
    );
}

#[test]
fn command_paf_query_min_chain_length_noop_when_zero() {
    // --min-chain-length 0 (default) keeps everything
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t30\t+\tB\t100\t0\t30\t25\t30\t255\tcg:Z:30M
";
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "--min-chain-length",
            "0",
            "-o",
            "bed",
        ])
        .stdin(paf)
        .run();
    assert!(
        stdout.contains("A\t0\t100"),
        "A should be kept (filter off)"
    );
    assert!(stdout.contains("C\t0\t30"), "C should be kept (filter off)");
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

// ── real-data validation ────────────────────────────────────────

#[test]
fn command_maf_to_paf_real_multiz_spar() {
    let (stdout, _) = PgrCmd::new()
        .args(&["maf", "to-paf", "tests/multiz/S288cvsSpar.maf"])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields[0], "Spar.gi_29362594");
    assert_eq!(fields[4], "-");
    assert_eq!(fields[5], "S288c.I");
}

#[test]
fn command_maf_to_paf_real_multiz_rm11() {
    let (stdout, _) = PgrCmd::new()
        .args(&["maf", "to-paf", "tests/multiz/S288cvsRM11_1a.maf"])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields[0], "RM11_1a.scaffold_17");
    assert_eq!(fields[4], "+");
    assert_eq!(fields[5], "S288c.I");
}

#[test]
fn command_paf_query_real_multiz_transitive() {
    use std::fs;
    use std::process::Command;
    let paf_path = "/tmp/pgr_real_test_merged.paf";
    let idx_path = "/tmp/pgr_real_test_merged.paf.idx";
    let pgr = std::env::current_dir().unwrap().join("target/debug/pgr");
    let spar_out = Command::new(&pgr)
        .args(["maf", "to-paf", "tests/multiz/S288cvsSpar.maf"])
        .output()
        .unwrap();
    let rm11_out = Command::new(&pgr)
        .args(["maf", "to-paf", "tests/multiz/S288cvsRM11_1a.maf"])
        .output()
        .unwrap();
    let mut merged = spar_out.stdout.clone();
    merged.extend_from_slice(&rm11_out.stdout);
    fs::write(paf_path, &merged).unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", paf_path, "-o", idx_path])
        .run();
    assert!(stderr.contains("saved to"));
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
