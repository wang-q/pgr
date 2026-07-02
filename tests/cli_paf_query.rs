#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

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
fn command_paf_query_max_depth_short_m() {
    // 3-hop chain: D(target) <- C <- B <- A. `-m 1` == `--max-depth 1`
    // -> only 1-hop neighbor C reaches; B (2-hop) and A (3-hop) excluded.
    let paf = "C\t100\t0\t100\t+\tD\t100\t0\t100\t90\t100\t255\tcg:Z:100M\n\
B\t100\t0\t100\t+\tC\t100\t0\t100\t90\t100\t255\tcg:Z:100M\n\
A\t100\t0\t100\t+\tB\t100\t0\t100\t90\t100\t255\tcg:Z:100M\n";
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "D:0-100",
            "--transitive",
            "--max-depth",
            "1",
        ])
        .stdin(paf)
        .run();
    assert!(stdout.contains("C\t0\t0\t100\t+\tD"), "C (1-hop) not found");
    assert!(!stdout.contains("B\t"), "B (2-hop) should NOT appear: -m 1");
    assert!(!stdout.contains("A\t"), "A (3-hop) should NOT appear: -m 1");
}

#[test]
fn command_paf_query_max_depth_unlimited() {
    // 3-hop chain: D(target) <- C <- B <- A.
    // Default --max-depth 2 reaches C (1-hop) and B (2-hop) but NOT A (3-hop).
    // --max-depth 0 (unlimited) reaches A.
    let paf = "C\t100\t0\t100\t+\tD\t100\t0\t100\t90\t100\t255\tcg:Z:100M\n\
B\t100\t0\t100\t+\tC\t100\t0\t100\t90\t100\t255\tcg:Z:100M\n\
A\t100\t0\t100\t+\tB\t100\t0\t100\t90\t100\t255\tcg:Z:100M\n";
    // Default depth 2: A (3-hop) excluded.
    let (stdout_default, _) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "D:0-100", "--transitive"])
        .stdin(paf)
        .run();
    assert!(stdout_default.contains("C\t0\t0\t100\t+\tD"), "C (1-hop)");
    assert!(stdout_default.contains("B\t0\t0\t100\t+\tC"), "B (2-hop)");
    assert!(
        !stdout_default.contains("A\t"),
        "A (3-hop) should NOT appear at default max-depth=2"
    );
    // Unlimited: A (3-hop) included.
    let (stdout_unlim, _) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "D:0-100",
            "--transitive",
            "--max-depth",
            "0",
        ])
        .stdin(paf)
        .run();
    assert!(
        stdout_unlim.contains("A\t0\t0\t100\t+\tB"),
        "A (3-hop) at -m 0"
    );
}

#[test]
fn command_paf_query_syntenic_filter() {
    use std::fs;
    // PAF: A-B and C-B alignments, target B:0-100.
    // Chain: only B->A (covers A's 0-100). No B->C chain.
    // Without filter: A and C both present.
    // With filter: only A present (C dropped, no syntenic chain B->C).
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n\
C\t100\t0\t100\t+\tB\t100\t0\t100\t90\t100\t255\tcg:Z:100M\n";
    // chain header: chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    // then one block line "100" (size 100, no gap), then blank line.
    let chain = "chain 100 B 100 + 0 100 A 100 + 0 100 1\n100\n\n";
    let temp = tempfile::TempDir::new().unwrap();
    let chain_path = temp.path().join("syn.chain");
    fs::write(&chain_path, chain).unwrap();

    // Without filter: both A and C.
    let (stdout_no, _) = PgrCmd::new()
        .args(&["paf", "query", "stdin", "B:0-100"])
        .stdin(paf)
        .run();
    assert!(stdout_no.contains("A\t0\t0\t100\t+\tB"), "A without filter");
    assert!(stdout_no.contains("C\t0\t0\t100\t+\tB"), "C without filter");

    // With filter: only A (C has no B->C chain).
    let (stdout_f, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "--syntenic-filter",
            chain_path.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(
        stdout_f.contains("A\t0\t0\t100\t+\tB"),
        "A should be syntenic"
    );
    assert!(
        !stdout_f.contains("C\t"),
        "C should be dropped by syntenic-filter"
    );
    assert!(
        stderr.contains("syntenic-filter: dropped"),
        "should log dropped count"
    );
}

#[test]
fn command_paf_query_syntenic_filter_no_overlap() {
    use std::fs;
    // Chain B->A exists but query span (200-300) does NOT overlap A's query interval (0-100).
    // With filter: A also dropped (no chain covers its query interval).
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let chain = "chain 100 B 100 + 0 100 A 1000 + 200 300 1\n100\n\n";
    let temp = tempfile::TempDir::new().unwrap();
    let chain_path = temp.path().join("syn_nooverlap.chain");
    fs::write(&chain_path, chain).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "--syntenic-filter",
            chain_path.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(
        !stdout.contains("A\t0\t0\t100\t+\tB"),
        "A should be dropped: chain query span 200-300 does not overlap A's 0-100"
    );
}

#[test]
fn command_paf_query_subset_filter() {
    use std::fs;
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\nC\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M\n";
    let temp = tempfile::TempDir::new().unwrap();
    let list = temp.path().join("subset.txt");
    fs::write(&list, "A\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "--subset-sequence-list",
            list.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();
    assert!(stdout.contains("A"), "A should be included");
    assert!(!stdout.contains("C"), "C should be excluded");
}

// ── paf query -b (batch BED regions) ─────────────────────────────

#[test]
fn command_paf_query_batch_bed() {
    use std::fs;
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tD\t100\t50\t100\t45\t50\t255\tcg:Z:50M
";
    let temp = tempfile::TempDir::new().unwrap();
    let bed = temp.path().join("regions.bed");
    fs::write(&bed, "B\t0\t100\nD\t50\t100\n# comment line\n\n").unwrap();
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-bed", "stdin", "-b", bed.to_str().unwrap()])
        .stdin(paf)
        .run();
    assert!(stdout.contains("A\t0\t100"), "A (from region B) missing");
    assert!(stdout.contains("C\t0\t50"), "C (from region D) missing");
    assert!(stderr.contains("Total results: 2"), "total count missing");
}

#[test]
fn command_paf_query_batch_bed_skips_unknown_target() {
    use std::fs;
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let temp = tempfile::TempDir::new().unwrap();
    let bed = temp.path().join("unknown.bed");
    fs::write(&bed, "B\t0\t100\nZ\t0\t100\n").unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "to-bed", "stdin", "-b", bed.to_str().unwrap()])
        .stdin(paf)
        .run();
    assert!(
        stderr.contains("not found in index, skipping"),
        "missing skip warning for unknown target"
    );
}

#[test]
fn command_paf_query_region_and_bed_mutually_exclusive() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let bed = temp.path().join("mutex.bed");
    fs::write(&bed, "B\t0\t100\n").unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "-b",
            bed.to_str().unwrap(),
        ])
        .stdin("A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n")
        .run_fail();
    assert!(
        stderr.contains("mutually exclusive"),
        "missing mutual-exclusion error"
    );
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
        .args(&["paf", "to-bed", "stdin", "B:0-100", "--min-degree", "2"])
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
        .args(&["paf", "to-bed", "stdin", "B:0-100", "--min-degree", "3"])
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
            "to-bed",
            "stdin",
            "B:0-100",
            "--min-chain-length",
            "50",
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
            "to-bed",
            "stdin",
            "B:0-100",
            "--min-chain-length",
            "0",
        ])
        .stdin(paf)
        .run();
    assert!(
        stdout.contains("A\t0\t100"),
        "A should be kept (filter off)"
    );
    assert!(stdout.contains("C\t0\t30"), "C should be kept (filter off)");
}
