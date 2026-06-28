#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── paf query bidirectional index (mirror in reverse_trees) ──────
// These tests verify the mirror index: when only A→B exists in the PAF,
// the index inserts a synthetic B→A entry into `reverse_trees[A]` so that
// BFS from A can still find B via the --transitive flag.

#[test]
fn command_paf_query_bidirectional_mirror_finds_target() {
    // Only A→B record (A is query, B is target). Without mirror index,
    // querying A would find nothing (trees[A] is empty). With mirror index,
    // reverse_trees[A] contains B, so BFS from A finds B.
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "A:0-100", "--transitive"])
        .stdin(paf)
        .run();
    assert!(
        stdout.contains("B\t0\t0\t100\t+\tA"),
        "B should be found via mirror index when querying from A"
    );
}

#[test]
fn command_paf_query_single_hop_does_not_use_mirror() {
    // Same PAF as above, but without --transitive. Single-hop query only
    // searches `trees[A]`, which is empty, so no results.
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "A:0-100"])
        .stdin(paf)
        .run();
    assert!(
        stderr.contains("No results found"),
        "single-hop should not use mirror index"
    );
}

#[test]
fn command_paf_query_bidirectional_multi_hop_via_mirror() {
    // A→B and C→B (both target B). Query from A should find B (via mirror
    // in reverse_trees[A]) and then C (via trees[B]). Without mirror index,
    // query from A would find nothing.
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t100\t+\tB\t100\t0\t100\t90\t100\t255\tcg:Z:100M
";
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "A:0-100", "--transitive"])
        .stdin(paf)
        .run();
    assert!(
        stdout.contains("B\t0\t0\t100\t+\tA"),
        "B should be found via mirror (1-hop)"
    );
    assert!(
        stdout.contains("C\t"),
        "C should be found via trees[B] (2-hop)"
    );
}

#[test]
fn command_paf_query_mirror_cigar_reversed() {
    // A has 200bp insertion relative to B. Original CIGAR: 50M200I50M.
    // Mirror entry reverses CIGAR and swaps I/D: 50M200D50M.
    let paf = "A\t300\t0\t300\t+\tB\t100\t0\t100\t95\t300\t255\tcg:Z:50M200I50M\n";
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "A:0-300", "--transitive"])
        .stdin(paf)
        .run();
    assert!(
        stdout.contains("cg:Z:50M200D50M"),
        "mirror entry should have reversed CIGAR with I/D swapped"
    );
}

#[test]
fn command_paf_query_reverse_strand_no_mirror() {
    // Minus-strand records do not get mirror entries (coordinate transform
    // is non-trivial). Query from A should find nothing via mirror.
    let paf = "A\t100\t0\t100\t-\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "A:0-100", "--transitive"])
        .stdin(paf)
        .run();
    assert!(
        stderr.contains("No results found"),
        "minus-strand should not have mirror entry"
    );
}

#[test]
fn command_paf_query_bidirectional_persists_across_save_load() {
    use std::fs;
    let paf_path = "/tmp/pgr_cli_test_bidir.paf";
    let idx_path = "/tmp/pgr_cli_test_bidir.paf.idx";
    fs::write(
        paf_path,
        "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n",
    )
    .unwrap();
    PgrCmd::new()
        .args(&["paf", "index", paf_path, "-o", idx_path])
        .run();
    // Query from A — only works if reverse_trees persisted across save/load.
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", idx_path, "A:0-100", "--transitive"])
        .run();
    assert!(
        stdout.contains("B\t0\t0\t100\t+\tA"),
        "bidirectional index should persist across save/load"
    );
    let _ = fs::remove_file(paf_path);
    let _ = fs::remove_file(idx_path);
}

// ── transitive closure invariants ────────────────────────────────
// Borrowed from impg test_transitive_integrity.rs: each test uses a tiny
// 2-4 line PAF to construct a precise graph topology and assert that the
// BFS transitive closure preserves the expected invariants.

#[test]
fn command_paf_transitive_non_overlapping_regions_stay_separate() {
    use std::collections::HashSet;
    // A:0-100 → B:0-100 and A:500-600 → C:0-100 (two non-overlapping A regions).
    // Query A:0-100 should find B but NOT C; query A:500-600 should find C but NOT B.
    let paf = "\
A\t1000\t0\t100\t+\tB\t1000\t0\t100\t100\t100\t60\tcg:Z:100=
A\t1000\t500\t600\t+\tC\t1000\t0\t100\t100\t100\t60\tcg:Z:100=";

    // Query A:0-100 — should find B only
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "to-bed",
            "stdin",
            "A:0-100",
            "--transitive",
            "--min-len",
            "0",
        ])
        .stdin(paf)
        .run();
    let names: HashSet<&str> = stdout
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert!(names.contains("B"), "query A:0-100 should find B");
    assert!(
        !names.contains("C"),
        "query A:0-100 should NOT find C (non-overlapping)"
    );

    // Query A:500-600 — should find C only
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "to-bed",
            "stdin",
            "A:500-600",
            "--transitive",
            "--min-len",
            "0",
        ])
        .stdin(paf)
        .run();
    let names: HashSet<&str> = stdout
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert!(names.contains("C"), "query A:500-600 should find C");
    assert!(
        !names.contains("B"),
        "query A:500-600 should NOT find B (non-overlapping)"
    );
}

#[test]
fn command_paf_transitive_coordinate_accuracy_subregion() {
    // A:0-100 → B:0-100 → C:0-100 (transitive chain).
    // Query A:25-75 should project to B:25-75 and C:25-75, not the full 0-100.
    let paf = "\
A\t1000\t0\t100\t+\tB\t1000\t0\t100\t100\t100\t60\tcg:Z:100=
B\t1000\t0\t100\t+\tC\t1000\t0\t100\t100\t100\t60\tcg:Z:100=";

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "to-bed",
            "stdin",
            "A:25-75",
            "--transitive",
            "--min-len",
            "0",
        ])
        .stdin(paf)
        .run();

    let mut found_b = false;
    let mut found_c = false;
    for line in stdout.lines().filter(|l| !l.is_empty()) {
        let fields: Vec<&str> = line.split('\t').collect();
        let name = fields[0];
        let start: i64 = fields[1].parse().unwrap();
        let end: i64 = fields[2].parse().unwrap();
        // Coordinates should be roughly 25-75, not 0-100
        assert!(
            (20..=30).contains(&start),
            "start on {name} should be ~25, got {start}"
        );
        assert!(
            (70..=80).contains(&end),
            "end on {name} should be ~75, got {end}"
        );
        if name == "B" {
            found_b = true;
        }
        if name == "C" {
            found_c = true;
        }
    }
    assert!(found_b, "should find B via transitive chain");
    assert!(found_c, "should find C via transitive chain");
}

#[test]
fn command_paf_transitive_distant_regions_no_collapse() {
    // A:0-100 → B:0-100 → D:0-100
    // A:1000-1100 → C:0-100 → D:500-600
    // Query A:0-100 should find D:0-100 (via B), NOT D:500-600.
    // Query A:1000-1100 should find D:500-600 (via C), NOT D:0-100.
    let paf = "\
A\t2000\t0\t100\t+\tB\t1000\t0\t100\t100\t100\t60\tcg:Z:100=
A\t2000\t1000\t1100\t+\tC\t1000\t0\t100\t100\t100\t60\tcg:Z:100=
B\t1000\t0\t100\t+\tD\t1000\t0\t100\t100\t100\t60\tcg:Z:100=
C\t1000\t0\t100\t+\tD\t1000\t500\t600\t100\t100\t60\tcg:Z:100=";

    // Query A:0-100 — D should be near 0, not 500
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "to-bed",
            "stdin",
            "A:0-100",
            "--transitive",
            "--max-depth",
            "3",
            "--min-len",
            "0",
        ])
        .stdin(paf)
        .run();
    let d_lines: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.is_empty() && l.starts_with("D\t"))
        .collect();
    assert!(!d_lines.is_empty(), "should find D via transitive path");
    for line in d_lines {
        let start: i64 = line.split('\t').nth(1).unwrap().parse().unwrap();
        assert!(
            start < 200,
            "D from A:0-100 path should be near 0, got {start}"
        );
    }

    // Query A:1000-1100 — D should be near 500, not 0
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "to-bed",
            "stdin",
            "A:1000-1100",
            "--transitive",
            "--max-depth",
            "3",
            "--min-len",
            "0",
        ])
        .stdin(paf)
        .run();
    let d_lines: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.is_empty() && l.starts_with("D\t"))
        .collect();
    assert!(!d_lines.is_empty(), "should find D via transitive path");
    for line in d_lines {
        let start: i64 = line.split('\t').nth(1).unwrap().parse().unwrap();
        assert!(
            start >= 400,
            "D from A:1000-1100 path should be near 500, got {start}"
        );
    }
}

#[test]
fn command_paf_transitive_multiple_alignments_to_same_target_stay_separate() {
    use std::collections::HashSet;
    // A:0-100 → B:0-100 and A:0-100 → B:500-600 (two alignments from same A region
    // to different B regions). Query A:0-100 should report TWO separate B results,
    // not one merged.
    let paf = "\
A\t1000\t0\t100\t+\tB\t1000\t0\t100\t100\t100\t60\tcg:Z:100=
A\t1000\t0\t100\t+\tB\t1000\t500\t600\t100\t100\t60\tcg:Z:100=";

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "to-bed",
            "stdin",
            "A:0-100",
            "--transitive",
            "--min-len",
            "0",
        ])
        .stdin(paf)
        .run();
    let b_lines: Vec<&str> = stdout
        .lines()
        .filter(|l| !l.is_empty() && l.starts_with("B\t"))
        .collect();
    assert_eq!(
        b_lines.len(),
        2,
        "should have 2 separate B results, got {}",
        b_lines.len()
    );
    // The two B results should be at different positions
    let starts: HashSet<i64> = b_lines
        .iter()
        .map(|l| l.split('\t').nth(1).unwrap().parse::<i64>().unwrap())
        .collect();
    assert_eq!(
        starts.len(),
        2,
        "the two B results should be at different positions"
    );
}
