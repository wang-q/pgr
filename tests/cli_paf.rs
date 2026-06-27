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

// ── paf graph (V4a coarse GFA induction) ──────────────────────

fn write_temp_fasta(path: &str, records: &[(&str, &str)]) {
    use std::fs;
    let mut content = String::new();
    for (name, seq) in records {
        content.push_str(">");
        content.push_str(name);
        content.push('\n');
        content.push_str(seq);
        content.push('\n');
    }
    fs::write(path, content).unwrap();
}

#[test]
fn command_paf_graph_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "graph", "--help"]).run();
    assert!(stdout.contains("Induces a coarse GFA graph"));
    assert!(stdout.contains("--min-var-len"));
    assert!(stdout.contains("--fasta"));
}

#[test]
fn command_paf_graph_basic_forward() {
    // A and B share a 100bp alignment → one shared node + trailing novel segments.
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let fa = "/tmp/pgr_graph_basic.fa";
    write_temp_fasta(fa, &[("A", &"ACGT".repeat(25)), ("B", &"TGCA".repeat(25))]);
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa])
        .stdin(paf)
        .run();
    // At least one S line, one P line for each sequence.
    let s_count = stdout.lines().filter(|l| l.starts_with("S\t")).count();
    let p_count = stdout.lines().filter(|l| l.starts_with("P\t")).count();
    assert!(s_count >= 1, "expected >= 1 S line, got {s_count}");
    assert_eq!(p_count, 2, "expected 2 P lines (A, B), got {p_count}");
    assert!(stdout.contains("\nP\tA\t"), "missing P line for A");
    assert!(stdout.contains("\nP\tB\t"), "missing P line for B");
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_split_at_large_indel() {
    // 50M 200I 50M: 200I >= 100 → split. B has an insertion (novel node in B path).
    let paf = "A\t300\t0\t100\t+\tB\t300\t0\t300\t95\t300\t255\tcg:Z:50M200I50M\n";
    let fa = "/tmp/pgr_graph_split.fa";
    write_temp_fasta(fa, &[("A", &"A".repeat(300)), ("B", &"G".repeat(300))]);
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa, "--min-var-len", "100"])
        .stdin(paf)
        .run();
    // B path should have >= 3 steps (aligned, novel insertion, aligned).
    let b_line = stdout
        .lines()
        .find(|l| l.starts_with("P\tB\t"))
        .expect("missing P line for B");
    // P line format: P\tname\tpath\toverlaps — path is the 3rd field.
    let path_field: &str = b_line.split('\t').nth(2).unwrap();
    let step_count = path_field.split(',').count();
    assert!(
        step_count >= 3,
        "B path should have >= 3 steps (aligned, novel, aligned), got {step_count}: {path_field}"
    );
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_small_indel_no_split() {
    // 50M 30I 50M: 30I < 100 → no split. A and B share exactly one aligned node.
    let paf = "A\t200\t0\t130\t+\tB\t200\t0\t160\t95\t160\t255\tcg:Z:50M30I50M\n";
    let fa = "/tmp/pgr_graph_nosplit.fa";
    write_temp_fasta(fa, &[("A", &"A".repeat(200)), ("B", &"G".repeat(200))]);
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa, "--min-var-len", "100"])
        .stdin(paf)
        .run();
    // Find shared nodes between A and B paths.
    let a_line = stdout
        .lines()
        .find(|l| l.starts_with("P\tA\t"))
        .expect("missing P line for A");
    let b_line = stdout
        .lines()
        .find(|l| l.starts_with("P\tB\t"))
        .expect("missing P line for B");
    let a_steps: Vec<&str> = a_line.split('\t').nth(2).unwrap().split(',').collect();
    let b_steps: Vec<&str> = b_line.split('\t').nth(2).unwrap().split(',').collect();
    // Strip orientation suffix to compare node ids.
    let a_nodes: Vec<&str> = a_steps
        .iter()
        .map(|s| s.trim_end_matches(['+', '-']))
        .collect();
    let b_nodes: Vec<&str> = b_steps
        .iter()
        .map(|s| s.trim_end_matches(['+', '-']))
        .collect();
    let shared: Vec<&str> = a_nodes
        .iter()
        .filter(|n| b_nodes.contains(n))
        .copied()
        .collect();
    assert_eq!(
        shared.len(),
        1,
        "expected exactly 1 shared node (no split), got {shared:?}"
    );
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_reverse_strand() {
    // Reverse strand alignment: query coords flipped, but A and B still share a node.
    let paf = "A\t100\t0\t100\t-\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let fa = "/tmp/pgr_graph_rc.fa";
    write_temp_fasta(fa, &[("A", &"ACGT".repeat(25)), ("B", &"TGCA".repeat(25))]);
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa])
        .stdin(paf)
        .run();
    let a_line = stdout
        .lines()
        .find(|l| l.starts_with("P\tA\t"))
        .expect("missing P line for A");
    let b_line = stdout
        .lines()
        .find(|l| l.starts_with("P\tB\t"))
        .expect("missing P line for B");
    let a_nodes: Vec<&str> = a_line
        .split('\t')
        .nth(2)
        .unwrap()
        .split(',')
        .map(|s| s.trim_end_matches(['+', '-']))
        .collect();
    let b_nodes: Vec<&str> = b_line
        .split('\t')
        .nth(2)
        .unwrap()
        .split(',')
        .map(|s| s.trim_end_matches(['+', '-']))
        .collect();
    let shared: Vec<&str> = a_nodes
        .iter()
        .filter(|n| b_nodes.contains(n))
        .copied()
        .collect();
    assert!(
        !shared.is_empty(),
        "reverse-strand alignment should still produce a shared node"
    );
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_min_var_len_filter() {
    // 50M 150I 50M with --min-var-len 200: 150I < 200 → no split.
    // Same alignment with --min-var-len 100: 150I >= 100 → split.
    let paf = "A\t300\t0\t100\t+\tB\t300\t0\t250\t95\t250\t255\tcg:Z:50M150I50M\n";
    let fa = "/tmp/pgr_graph_filter.fa";
    write_temp_fasta(fa, &[("A", &"A".repeat(300)), ("B", &"G".repeat(300))]);

    // With threshold 200: no split, B path has 1 shared node + trailing novel.
    let (stdout_no_split, _) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa, "--min-var-len", "200"])
        .stdin(paf)
        .run();
    let b_line = stdout_no_split
        .lines()
        .find(|l| l.starts_with("P\tB\t"))
        .expect("missing P line for B");
    let steps_no_split = b_line.split('\t').nth(2).unwrap().split(',').count();

    // With threshold 100: split, B path has >= 3 steps.
    let (stdout_split, _) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa, "--min-var-len", "100"])
        .stdin(paf)
        .run();
    let b_line = stdout_split
        .lines()
        .find(|l| l.starts_with("P\tB\t"))
        .expect("missing P line for B");
    let steps_split = b_line.split('\t').nth(2).unwrap().split(',').count();

    assert!(
        steps_split > steps_no_split,
        "split path ({steps_split}) should have more steps than no-split ({steps_no_split})"
    );
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_missing_fasta_fails() {
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let (_stdout, stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", "/nonexistent/path.fa"])
        .stdin(paf)
        .run_fail();
    // Should fail with a friendly error, not panic.
    assert!(
        stderr.contains("could not open") || stderr.contains("No such file"),
        "expected file-not-found error, got: {stderr}"
    );
}
