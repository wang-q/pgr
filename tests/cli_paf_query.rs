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
        .args(&["paf", "query", "stdin", "D:0-100", "-t", "-m", "1"])
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
        .args(&["paf", "query", "stdin", "D:0-100", "-t"])
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
        .args(&["paf", "query", "stdin", "D:0-100", "-t", "-m", "0"])
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
    let chain_path = "/tmp/pgr_syntenic.chain";
    fs::write(chain_path, chain).unwrap();

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
            chain_path,
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

    let _ = fs::remove_file(chain_path);
}

#[test]
fn command_paf_query_syntenic_filter_no_overlap() {
    use std::fs;
    // Chain B->A exists but query span (200-300) does NOT overlap A's query interval (0-100).
    // With filter: A also dropped (no chain covers its query interval).
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let chain = "chain 100 B 100 + 0 100 A 1000 + 200 300 1\n100\n\n";
    let chain_path = "/tmp/pgr_syntenic_nooverlap.chain";
    fs::write(chain_path, chain).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            "stdin",
            "B:0-100",
            "--syntenic-filter",
            chain_path,
        ])
        .stdin(paf)
        .run();
    assert!(
        !stdout.contains("A\t0\t0\t100\t+\tB"),
        "A should be dropped: chain query span 200-300 does not overlap A's 0-100"
    );

    let _ = fs::remove_file(chain_path);
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

// ── paf to-bed (BED3 output) ─────────────────────────────────────

#[test]
fn command_paf_to_bed_output() {
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M
";
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "to-bed", "stdin", "B:0-100"])
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
fn command_paf_to_bed_output_reverse_strand() {
    // Reverse-strand alignment: query coords should still be emitted as (min, max)
    let paf = "A\t100\t0\t100\t-\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "to-bed", "stdin", "B:0-100"])
        .stdin(paf)
        .run();
    assert!(
        stdout.contains("A\t0\t100"),
        "A BED3 line missing (reverse strand)"
    );
}

// ── indel coordinate accuracy at query layer ─────────────────────
// Borrowed from impg test_transitive_integrity.rs::test_indel_coordinate_accuracy
// (Test 6): with indels in the CIGAR, the query layer must project target
// sub-intervals onto the query without coordinate drift. pgr's indel
// coordinate tests were all at the to-maf layer; this covers the to-bed
// (query) layer.

#[test]
fn command_paf_to_bed_insertion_coordinate_accuracy() {
    // CIGAR: 50= 10I 50= → A:0-110 (query) → B:0-100 (target).
    //   - 50= : A:0-50  ↔ B:0-50
    //   - 10I : A:50-60 (insertion in A, no B consumption)
    //   - 50= : A:60-110 ↔ B:50-100
    // Query B:0-50 (before insertion) → A:0-50.
    // Query B:50-100 (after insertion) → A:60-110 (skip the 10bp insertion).
    let paf = "A\t110\t0\t110\t+\tB\t100\t0\t100\t100\t110\t60\tcg:Z:50=10I50=\n";

    // Query B:0-50 — should project to A:0-50 (before the insertion).
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "to-bed", "stdin", "B:0-50", "--min-len", "0"])
        .stdin(paf)
        .run();
    let a_line = stdout
        .lines()
        .find(|l| l.starts_with("A\t"))
        .unwrap_or_else(|| panic!("missing A line in BED output:\n{stdout}"));
    let fields: Vec<&str> = a_line.split('\t').collect();
    let start: i64 = fields[1].parse().unwrap();
    let end: i64 = fields[2].parse().unwrap();
    assert!(
        (0..=5).contains(&start) && (45..=55).contains(&end),
        "B:0-50 (before insertion) should map to A:~0-50, got A:{start}-{end}"
    );

    // Query B:50-100 — should project to A:60-110 (after the insertion).
    // Note: querying exactly at the insertion boundary (B:50) may include the
    // adjacent insertion bases (A:50-60) in the projected range; the end
    // coordinate (110) is the strong invariant. Query B:60-100 (well inside
    // the post-insertion region) for a clean start-coordinate check.
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "to-bed", "stdin", "B:50-100", "--min-len", "0"])
        .stdin(paf)
        .run();
    let a_line = stdout
        .lines()
        .find(|l| l.starts_with("A\t"))
        .unwrap_or_else(|| panic!("missing A line in BED output:\n{stdout}"));
    let fields: Vec<&str> = a_line.split('\t').collect();
    let end: i64 = fields[2].parse().unwrap();
    assert!(
        end >= 105,
        "B:50-100 (after insertion) should map end near 110, got A end={end}"
    );

    // Query B:60-100 (10bp inside the post-insertion region) — start should
    // be ~70 (60 + 10), cleanly after the insertion with no boundary effect.
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "to-bed", "stdin", "B:60-100", "--min-len", "0"])
        .stdin(paf)
        .run();
    let a_line = stdout
        .lines()
        .find(|l| l.starts_with("A\t"))
        .unwrap_or_else(|| panic!("missing A line in BED output:\n{stdout}"));
    let fields: Vec<&str> = a_line.split('\t').collect();
    let start: i64 = fields[1].parse().unwrap();
    let end: i64 = fields[2].parse().unwrap();
    assert!(
        (65..=75).contains(&start),
        "B:60-100 (inside post-insertion) should map start near 70, got A:{start}"
    );
    assert!(
        end >= 105,
        "B:60-100 (inside post-insertion) should map end near 110, got A:{end}"
    );
}

#[test]
fn command_paf_to_bed_deletion_coordinate_accuracy() {
    // CIGAR: 50= 10D 50= → A:0-100 (query) → B:0-110 (target).
    //   - 50= : A:0-50   ↔ B:0-50
    //   - 10D : B:50-60  (deletion in A, 10bp in B not in A)
    //   - 50= : A:50-100 ↔ B:60-110
    // Query B:0-50 (before deletion) → A:0-50.
    // Query B:60-110 (after deletion) → A:50-100.
    let paf = "A\t100\t0\t100\t+\tB\t110\t0\t110\t100\t100\t60\tcg:Z:50=10D50=\n";

    // Query B:0-50 — should project to A:0-50 (before the deletion).
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "to-bed", "stdin", "B:0-50", "--min-len", "0"])
        .stdin(paf)
        .run();
    let a_line = stdout
        .lines()
        .find(|l| l.starts_with("A\t"))
        .unwrap_or_else(|| panic!("missing A line in BED output:\n{stdout}"));
    let fields: Vec<&str> = a_line.split('\t').collect();
    let start: i64 = fields[1].parse().unwrap();
    let end: i64 = fields[2].parse().unwrap();
    assert!(
        (0..=5).contains(&start) && (45..=55).contains(&end),
        "B:0-50 (before deletion) should map to A:~0-50, got A:{start}-{end}"
    );

    // Query B:60-110 — should project to A:50-100 (after the deletion).
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "to-bed", "stdin", "B:60-110", "--min-len", "0"])
        .stdin(paf)
        .run();
    let a_line = stdout
        .lines()
        .find(|l| l.starts_with("A\t"))
        .unwrap_or_else(|| panic!("missing A line in BED output:\n{stdout}"));
    let fields: Vec<&str> = a_line.split('\t').collect();
    let start: i64 = fields[1].parse().unwrap();
    let end: i64 = fields[2].parse().unwrap();
    assert!(
        (45..=55).contains(&start) && end >= 95,
        "B:60-110 (after deletion) should map to A:~50-100, got A:{start}-{end}"
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
        .args(&["paf", "to-bed", "stdin", "-b", bed])
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
        .args(&["paf", "to-bed", "stdin", "-b", bed])
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

// ── paf query bidirectional index (mirror in reverse_trees) ──────

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

// ── paf to-maf (pairwise MAF from CIGAR) ─────────────────────────

/// Write `content` to a plain `path`, then BGZF-compress it via `pgr fa gz`
/// (which also creates the .gzi index required for random access).
fn write_bgzf_fa(path_no_gz: &str, content: &str) -> String {
    use std::fs;
    fs::write(path_no_gz, content).unwrap();
    let (out, _) = PgrCmd::new().args(&["fa", "gz", path_no_gz]).run();
    let _ = out;
    let gz_path = format!("{path_no_gz}.gz");
    // Sanity: compression produced the .gz file
    assert!(
        std::path::Path::new(&gz_path).exists(),
        "pgr fa gz failed to produce {gz_path}"
    );
    gz_path
}

#[test]
fn command_paf_to_maf_strict_name_validation() {
    use std::fs;
    // PAF references A and B; TSV only has A.
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_strict_A.fa", ">A\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_maf_strict.tsv";
    fs::write(tsv, format!("A\t{a_fa}\n")).unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:0-10", "-f", tsv])
        .stdin(paf)
        .run_fail();
    assert!(
        stderr.contains("FASTA TSV is missing") && stderr.contains("B"),
        "missing strict validation error for B"
    );
    let _ = fs::remove_file("/tmp/pgr_maf_strict_A.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
}

#[test]
fn command_paf_to_maf_perfect_match() {
    use std::fs;
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_perfect_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_perfect_B.fa", ">B\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_maf_perfect.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:0-10", "-f", tsv])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(stdout.contains("##maf version=1"), "missing MAF header");
    assert!(stdout.contains("a"), "missing alignment header");
    // target line first, query line second
    assert!(
        stdout.contains("s\tB\t0\t10\t+\t10\tACGTACGTAC"),
        "missing/incorrect target line"
    );
    assert!(
        stdout.contains("s\tA\t0\t10\t+\t10\tACGTACGTAC"),
        "missing/incorrect query line"
    );
    let _ = fs::remove_file("/tmp/pgr_maf_perfect_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_perfect_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

#[test]
fn command_paf_to_maf_with_insertion() {
    use std::fs;
    // CIGAR: 4= 3I 3= → target 0-7, query 0-10
    // target: ACGT---ACG  (4 match + 3 gaps + 3 match)
    // query:  ACGTACGTAC  (4 match + 3 bases + 3 match)
    // But query[7..10] should be "TAC" (query = ACGTACGTAC, idx 7,8,9 = T,A,C)
    // and query[4..7] = "ACG"
    // So query alignment = ACGT + ACG + TAC = ACGTACGTAC
    // target alignment = ACGT + --- + ACG = ACGT---ACG
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t7\t7\t10\t255\tcg:Z:4=3I3=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_ins_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_ins_B.fa", ">B\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_maf_ins.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:0-7", "-f", tsv])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // target has gaps where query inserted
    assert!(
        stdout.contains("ACGT---ACG"),
        "target alignment should contain gaps for insertion"
    );
    assert!(
        stdout.contains("ACGTACGTAC"),
        "query alignment should contain full sequence"
    );
    // sizes: target 7 non-gap, query 10 non-gap
    assert!(
        stdout.contains("s\tB\t0\t7\t+\t10\tACGT---ACG"),
        "target size should be 7"
    );
    assert!(
        stdout.contains("s\tA\t0\t10\t+\t10\tACGTACGTAC"),
        "query size should be 10"
    );
    let _ = fs::remove_file("/tmp/pgr_maf_ins_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_ins_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

#[test]
fn command_paf_to_maf_with_deletion() {
    use std::fs;
    // CIGAR: 4= 3D 3= → target 0-10, query 0-7
    // target: ACGTACGTAC (4 match + 3 bases + 3 match)
    // query:  ACGT---ACG (4 match + 3 gaps + 3 match)
    let paf = "A\t7\t0\t7\t+\tB\t10\t0\t10\t7\t10\t255\tcg:Z:4=3D3=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_del_A.fa", ">A\nACGTACG\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_del_B.fa", ">B\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_maf_del.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:0-10", "-f", tsv])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(
        stdout.contains("ACGTACGTAC"),
        "target alignment should contain full sequence"
    );
    assert!(
        stdout.contains("ACGT---ACG"),
        "query alignment should contain gaps for deletion"
    );
    let _ = fs::remove_file("/tmp/pgr_maf_del_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_del_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

#[test]
fn command_paf_to_maf_trimmed_subregion() {
    use std::fs;
    // Full alignment 10= over B:0-10. Query B:2-8 should trim CIGAR to 6=.
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_trim_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_trim_B.fa", ">B\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_maf_trim.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:2-8", "-f", tsv])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // B[2..8) = GTACGT, A[2..8) = GTACGT
    assert!(
        stdout.contains("s\tB\t2\t6\t+\t10\tGTACGT"),
        "target should be trimmed to B:2-8"
    );
    assert!(
        stdout.contains("s\tA\t2\t6\t+\t10\tGTACGT"),
        "query should be trimmed to A:2-8"
    );
    let _ = fs::remove_file("/tmp/pgr_maf_trim_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_trim_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

#[test]
fn command_paf_to_maf_reverse_strand_perfect_match() {
    use std::fs;
    // '-' strand perfect match: target B forward == RC(query A forward).
    // A forward = GTACGTACGT, RC = ACGTACGTAC = B forward.
    // CIGAR 10= describes 10 alignment columns of target vs RC(query).
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_rev_A.fa", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_rev_B.fa", ">B\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_maf_rev.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:0-10", "-f", tsv])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    assert!(stdout.contains("##maf version=1"), "missing MAF header");
    // Target line: forward strand, original sequence.
    assert!(
        stdout.contains("s\tB\t0\t10\t+\t10\tACGTACGTAC"),
        "missing/incorrect target line for '-' strand record"
    );
    // Query line: '-' strand, displayed sequence is RC of A forward.
    // q_start_maf = srcSize - qe = 10 - 10 = 0; q_size = 10.
    assert!(
        stdout.contains("s\tA\t0\t10\t-\t10\tACGTACGTAC"),
        "missing/incorrect query line for '-' strand record (RC not applied)"
    );
    let _ = fs::remove_file("/tmp/pgr_maf_rev_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_rev_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

#[test]
fn command_paf_to_maf_reverse_strand_with_insertion() {
    use std::fs;
    // '-' strand alignment with insertion: CIGAR 4=3I3= (7 target, 10 query cols).
    // A forward = GTACGTACGT, RC(A) = ACGTACGTAC.
    // target B = ACGT (RC(A)[0:4]) + TAC (RC(A)[7:10]) = ACGTTAC (7 bp).
    // Expected alignment columns:
    //   target: ACGT---TAC  (4 match + 3 gaps + 3 match)
    //   query:  ACGTACGTAC  (RC of A forward, walked left-to-right)
    let paf = "A\t10\t0\t10\t-\tB\t7\t0\t7\t7\t7\t255\tcg:Z:4=3I3=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_revins_A.fa", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_revins_B.fa", ">B\nACGTTAC\n");
    let tsv = "/tmp/pgr_maf_revins.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:0-7", "-f", tsv])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // Target has gaps where query inserted.
    assert!(
        stdout.contains("ACGT---TAC"),
        "target alignment should contain gaps for insertion on '-' strand"
    );
    // Query alignment is RC of A forward, walked left-to-right.
    assert!(
        stdout.contains("ACGTACGTAC"),
        "query alignment should be RC of A forward on '-' strand"
    );
    // Query line should be '-' strand with q_start = srcSize - qe = 0.
    assert!(
        stdout.contains("s\tA\t0\t10\t-\t10\tACGTACGTAC"),
        "missing/incorrect query s-line for '-' strand with insertion"
    );
    let _ = fs::remove_file("/tmp/pgr_maf_revins_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_revins_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

#[test]
fn command_paf_to_maf_reverse_strand_subinterval_first_half() {
    use std::fs;
    // '-' strand perfect match, sub-interval query on the first half.
    // A forward = GTACGTACGT (G0 T1 A2 C3 G4 T5 A6 C7 G8 T9), RC = ACGTACGTAC = B.
    // PAF CIGAR 10= aligns RC(A) vs B left-to-right.
    // Query B:0-5 → CIGAR first 5 query bases = RC(A)[0:5] = ACGTA, which
    // corresponds to forward A[5:10) = TACGT (RC = ACGTA).
    // Before the project() fix this returned forward A[0:5) = GTACG
    // (RC = CGTAC), which is wrong.
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_revsub1_A.fa", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_revsub1_B.fa", ">B\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_maf_revsub1.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:0-5", "-f", tsv])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // Target sub-interval B[0:5] = ACGTA, +strand, start 0, size 5, srcSize 10.
    assert!(
        stdout.contains("s\tB\t0\t5\t+\t10\tACGTA"),
        "missing/incorrect target line for '-' strand sub-interval (first half)"
    );
    // Query: forward A[5:10] RC'd = ACGTA. q_start_maf = srcSize - qe = 10-10 = 0.
    assert!(
        stdout.contains("s\tA\t0\t5\t-\t10\tACGTA"),
        "missing/incorrect query line for '-' strand sub-interval (first half); \
         this verifies project() maps RC offset back to forward [5,10)"
    );
    // Sanity: the buggy forward A[0:5] mapping must NOT appear.
    assert!(
        !stdout.contains("CGTAC"),
        "regression: query sequence looks like RC of forward A[0:5] — project() \
         did not convert RC offset to forward coordinates on '-' strand"
    );
    let _ = fs::remove_file("/tmp/pgr_maf_revsub1_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_revsub1_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

#[test]
fn command_paf_to_maf_reverse_strand_subinterval_second_half() {
    use std::fs;
    // Same setup as the first-half test, but query B:5-10.
    // CIGAR last 5 query bases = RC(A)[5:10] = CGTAC, corresponding to
    // forward A[0:5) = GTACG (RC = CGTAC).
    // q_start_maf = srcSize - qe = 10 - 5 = 5.
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_revsub2_A.fa", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_revsub2_B.fa", ">B\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_maf_revsub2.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:5-10", "-f", tsv])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // Target sub-interval B[5:10] = CGTAC.
    assert!(
        stdout.contains("s\tB\t5\t5\t+\t10\tCGTAC"),
        "missing/incorrect target line for '-' strand sub-interval (second half)"
    );
    // Query: forward A[0:5] RC'd = CGTAC. q_start_maf = srcSize - qe = 10-5 = 5.
    assert!(
        stdout.contains("s\tA\t5\t5\t-\t10\tCGTAC"),
        "missing/incorrect query line for '-' strand sub-interval (second half); \
         this verifies project() maps RC offset back to forward [0,5)"
    );
    let _ = fs::remove_file("/tmp/pgr_maf_revsub2_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_revsub2_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

#[test]
fn command_paf_to_maf_reverse_strand_subinterval_with_insertion() {
    use std::fs;
    // '-' strand with insertion, sub-interval query on the trailing target
    // segment (which spans op2 3I + op3 2=).
    // A forward = GTACGTACGT (G0 T1 A2 C3 G4 T5 A6 C7 G8 T9), RC = ACGTACGTAC.
    // CIGAR 5=3I2= aligns RC(A) vs B (7 bp): B[0:5]=ACGTA, B[5:7]=AC, with
    // RC(A)[5:8]=CGT inserted between them.
    // Query B:5-7 → project returns forward A[0,5) (union of op2 RC[5,8)→fwd
    // [2,5) and op3 RC[8,10)→fwd [0,2)).
    // q_seq = RC(A[0:5]) = RC(GTACG) = CGTAC. qs_eff = rec_qe - qe = 10-5 = 5.
    // Alignment columns: op2 3I at ct=5 emits q_seq[0..3]=CGT with target
    // gaps; op3 2= at ct=[5,7) emits q_seq[3..5]=AC paired with B[5:7]=AC.
    let paf = "A\t10\t0\t10\t-\tB\t7\t0\t7\t7\t7\t255\tcg:Z:5=3I2=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_revsubi_A.fa", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_revsubi_B.fa", ">B\nACGTAAC\n");
    let tsv = "/tmp/pgr_maf_revsubi.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:5-7", "-f", tsv])
        .stdin(paf)
        .run();
    assert!(stderr.contains("Total results: 1"), "expected 1 result");
    // Target sub-interval B[5:7] = AC, with 3 gap columns before it from the
    // insertion op (which sits at target position 5, inside [5,7)).
    assert!(
        stdout.contains("s\tB\t5\t2\t+\t7\t---AC"),
        "missing/incorrect target line for '-' strand sub-interval with insertion"
    );
    // Query: RC(A[0:5]) = CGTAC. q_start_maf = srcSize - qe = 10-5 = 5.
    assert!(
        stdout.contains("s\tA\t5\t5\t-\t10\tCGTAC"),
        "missing/incorrect query line for '-' strand sub-interval with insertion"
    );
    let _ = fs::remove_file("/tmp/pgr_maf_revsubi_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_revsubi_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

// ── paf to-maf --msa (multi-way MSA via POA) ─────────────────────

#[test]
fn command_paf_to_maf_msa_three_genomes_transitive() {
    use std::fs;
    // Three genomes A/B/C, all 10 bp, identical sequence ACGTACGTAC.
    // A-B and A-C alignments → query B:0-10 with --transitive gathers
    // {B(target), A, C} into one region. --msa merges them into a single
    // 3-sequence MAF block via POA. Since all sequences are identical, the
    // MSA columns should be gap-free and all three `s` lines equal.
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_msa_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_msa_B.fa", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_maf_msa_C.fa", ">C\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_maf_msa.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:0-10", "-t", "--msa", "-f", tsv])
        .stdin(paf)
        .run();
    assert!(
        stderr.contains("Total results:") && !stderr.contains("Total results: 0"),
        "expected non-zero results"
    );
    assert!(stdout.contains("##maf version=1"), "missing MAF header");
    // Exactly one `a` block (multi-way).
    let a_count = stdout.matches("\na\n").count() + if stdout.starts_with("a\n") { 1 } else { 0 };
    assert_eq!(a_count, 1, "expected exactly one MAF block, got {a_count}");
    // Three `s` lines (target B + queries A, C).
    let s_count = stdout.lines().filter(|l| l.starts_with("s\t")).count();
    assert_eq!(s_count, 3, "expected 3 s-lines, got {s_count}");
    // All identical → each s-line should end with ACGTACGTAC (no gaps).
    for line in stdout.lines().filter(|l| l.starts_with("s\t")) {
        assert!(
            line.ends_with("ACGTACGTAC"),
            "expected gap-free ACGTACGTAC in s-line: {line}"
        );
    }
    // Target B should appear first.
    let first_s = stdout.lines().find(|l| l.starts_with("s\t")).unwrap();
    assert!(
        first_s.starts_with("s\tB\t"),
        "target B should be first: {first_s}"
    );

    let _ = fs::remove_file("/tmp/pgr_maf_msa_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_msa_B.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_msa_C.fa");
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
fn command_paf_to_maf_msa_with_snp() {
    use std::fs;
    // B = ACGTACGTAC (target)
    // A = ACGTACGTAC (identical to B)
    // C = ACGTTCGTAC (SNP at position 4: A→T)
    // A-B and A-C alignments, query B:0-10 --transitive --msa.
    // POA should produce a 3-sequence MSA with one SNP column.
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_msa_snp_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_msa_snp_B.fa", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_maf_msa_snp_C.fa", ">C\nACGTTCGTAC\n");
    let tsv = "/tmp/pgr_maf_msa_snp.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:0-10", "-t", "--msa", "-f", tsv])
        .stdin(paf)
        .run();
    let s_count = stdout.lines().filter(|l| l.starts_with("s\t")).count();
    assert_eq!(s_count, 3, "expected 3 s-lines, got {s_count}");
    // All three s-lines should have length 10 (no gaps introduced for a SNP).
    for line in stdout.lines().filter(|l| l.starts_with("s\t")) {
        let aln = line.split('\t').next_back().unwrap();
        assert_eq!(aln.len(), 10, "expected 10-char alignment, got '{aln}'");
    }
    // C should differ from B at position 4 (0-indexed).
    let b_line = stdout.lines().find(|l| l.starts_with("s\tB\t")).unwrap();
    let c_line = stdout.lines().find(|l| l.starts_with("s\tC\t")).unwrap();
    let b_aln = b_line.split('\t').next_back().unwrap();
    let c_aln = c_line.split('\t').next_back().unwrap();
    let diffs: Vec<usize> = b_aln
        .chars()
        .zip(c_aln.chars())
        .enumerate()
        .filter_map(|(i, (a, b))| if a != b { Some(i) } else { None })
        .collect();
    assert_eq!(
        diffs,
        vec![4],
        "expected single SNP at pos 4, got {diffs:?}"
    );

    let _ = fs::remove_file("/tmp/pgr_maf_msa_snp_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_msa_snp_B.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_msa_snp_C.fa");
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
fn command_paf_to_maf_msa_reverse_strand_query() {
    use std::fs;
    // B = ACGTACGTAC (target, forward)
    // A = GTACGTACGT, aligned to B on '-' strand. RC(A) = ACGTACGTAC = B,
    // so after reverse-complementation A's aligned sequence equals B.
    // Query B:0-10 --transitive --msa: A is RC'd before POA, then both
    // sequences are identical → gap-free MSA.
    let paf = "A\t10\t0\t10\t-\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_maf_msa_rev_A.fa", ">A\nGTACGTACGT\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_maf_msa_rev_B.fa", ">B\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_maf_msa_rev.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-maf", "stdin", "B:0-10", "-t", "--msa", "-f", tsv])
        .stdin(paf)
        .run();
    let s_count = stdout.lines().filter(|l| l.starts_with("s\t")).count();
    assert_eq!(s_count, 2, "expected 2 s-lines (B + A), got {s_count}");
    // A should be emitted with strand '-'.
    let a_line = stdout.lines().find(|l| l.starts_with("s\tA\t")).unwrap();
    assert!(
        a_line.contains("\t-\t"),
        "A should be on '-' strand: {a_line}"
    );
    // A's aligned sequence should be RC(GTACGTACGT) = ACGTACGTAC, gap-free.
    let a_aln = a_line.split('\t').next_back().unwrap();
    assert_eq!(a_aln, "ACGTACGTAC", "expected RC(A) gap-free: {a_aln}");

    let _ = fs::remove_file("/tmp/pgr_maf_msa_rev_A.fa");
    let _ = fs::remove_file("/tmp/pgr_maf_msa_rev_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

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

#[test]
fn command_paf_to_gfa_identical() {
    use std::fs;
    // Two identical sequences -> unchopped to a single 10-bp segment, no
    // edges, 2 paths traversing that one segment.
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_gfa_id_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_gfa_id_B.fa", ">B\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_gfa_id.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-gfa", "stdin", "B:0-10", "-t", "-f", tsv])
        .stdin(paf)
        .run();

    // GFA v1.0 header present.
    assert!(
        stdout.lines().any(|l| l == "H\tVN:Z:1.0"),
        "missing GFA H header line"
    );

    let s_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("S\t")).collect();
    let l_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("L\t")).collect();
    let p_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("P\t")).collect();

    // Unchopping collapses the 10 identical bases into one segment.
    assert_eq!(s_lines.len(), 1, "expected 1 S line, got {}", s_lines.len());
    assert_eq!(
        l_lines.len(),
        0,
        "expected 0 L lines, got {}",
        l_lines.len()
    );
    assert_eq!(
        p_lines.len(),
        2,
        "expected 2 P lines, got {}",
        p_lines.len()
    );

    // The single segment carries the full 10-bp sequence + LN tag.
    let s_fields: Vec<&str> = s_lines[0].split('\t').collect();
    assert_eq!(s_fields[1], "1", "segment id should be 1");
    assert_eq!(s_fields[2], "ACGTACGTAC", "segment sequence mismatch");
    assert!(
        s_fields.contains(&"LN:i:10"),
        "missing LN:i:10 tag in S line: {}",
        s_lines[0]
    );

    // Each path visits exactly one node (segment 1), zero overlaps.
    for p in &p_lines {
        let fields: Vec<&str> = p.split('\t').collect();
        assert_eq!(fields[2], "1+", "path should visit only segment 1: {p}");
        assert!(fields[3].is_empty(), "path should have no overlaps: {p}");
    }

    let _ = fs::remove_file("/tmp/pgr_gfa_id_A.fa");
    let _ = fs::remove_file("/tmp/pgr_gfa_id_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

#[test]
fn command_paf_to_gfa_with_snp() {
    use std::fs;
    // B = ACGTACGTAC (target)
    // A = ACGTACGTAC (identical to B)
    // C = ACGTTCGTAC (SNP at pos 4: A->T)
    // After unchopping: 4 segments (ACGT, A, T, CGTAC), 4 edges, 3 paths.
    // The SNP forms a bubble: seg2(A) and seg3(T) share in={1}, out={4}.
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_gfa_snp_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_gfa_snp_B.fa", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_gfa_snp_C.fa", ">C\nACGTTCGTAC\n");
    let tsv = "/tmp/pgr_gfa_snp.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-gfa", "stdin", "B:0-10", "-t", "-f", tsv])
        .stdin(paf)
        .run();

    let s_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("S\t")).collect();
    let l_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("L\t")).collect();
    let p_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("P\t")).collect();

    // 4 segments: ACGT, A, T, CGTAC.
    assert_eq!(
        s_lines.len(),
        4,
        "expected 4 S lines, got {}",
        s_lines.len()
    );
    // 4 edges: 1->2, 1->3, 2->4, 3->4.
    assert_eq!(
        l_lines.len(),
        4,
        "expected 4 L lines, got {}",
        l_lines.len()
    );
    // 3 paths (B, A, C).
    assert_eq!(
        p_lines.len(),
        3,
        "expected 3 P lines, got {}",
        p_lines.len()
    );

    // Collect segment sequences (id -> seq).
    let mut seg_seq: std::collections::HashMap<&str, &str> = std::collections::HashMap::new();
    for s in &s_lines {
        let f: Vec<&str> = s.split('\t').collect();
        seg_seq.insert(f[1], f[2]);
    }
    assert_eq!(seg_seq.get("1"), Some(&"ACGT"));
    assert_eq!(seg_seq.get("2"), Some(&"A"));
    assert_eq!(seg_seq.get("3"), Some(&"T"));
    assert_eq!(seg_seq.get("4"), Some(&"CGTAC"));

    // B and C paths should differ at exactly one segment (the SNP), gap-free.
    let b_path: Vec<&str> = p_lines
        .iter()
        .find(|p| p.starts_with("P\tB\t"))
        .unwrap()
        .split('\t')
        .nth(2)
        .unwrap()
        .split(',')
        .collect();
    let c_path: Vec<&str> = p_lines
        .iter()
        .find(|p| p.starts_with("P\tC\t"))
        .unwrap()
        .split('\t')
        .nth(2)
        .unwrap()
        .split(',')
        .collect();

    assert_eq!(
        b_path.len(),
        c_path.len(),
        "B and C paths should have equal length (gap-free SNP)"
    );
    let diffs: usize = b_path
        .iter()
        .zip(c_path.iter())
        .filter(|(a, b)| a != b)
        .count();
    assert_eq!(
        diffs, 1,
        "B and C paths should differ at 1 segment (SNP), got {diffs}"
    );

    // B's path should go through the A allele (seg 2), C's through T (seg 3).
    assert!(
        b_path.iter().any(|s| s.starts_with("2+")),
        "B should traverse segment 2 (A allele): {b_path:?}"
    );
    assert!(
        c_path.iter().any(|s| s.starts_with("3+")),
        "C should traverse segment 3 (T allele): {c_path:?}"
    );

    let _ = fs::remove_file("/tmp/pgr_gfa_snp_A.fa");
    let _ = fs::remove_file("/tmp/pgr_gfa_snp_B.fa");
    let _ = fs::remove_file("/tmp/pgr_gfa_snp_C.fa");
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
fn command_paf_to_gfa_crush() {
    use std::fs;
    // Same setup as command_paf_to_gfa_with_snp, but with --crush.
    // The SNP bubble (seg2=A, seg3=T) collapses to one segment (A, the
    // higher-weight allele: B+A=2 vs C=1). Paths through T are rewritten
    // to A, losing base-level ALT info.
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_gfa_crush_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_gfa_crush_B.fa", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_gfa_crush_C.fa", ">C\nACGTTCGTAC\n");
    let tsv = "/tmp/pgr_gfa_crush.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf", "to-gfa", "stdin", "B:0-10", "-t", "-f", tsv, "--crush",
        ])
        .stdin(paf)
        .run();

    let s_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("S\t")).collect();
    let l_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("L\t")).collect();
    let p_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("P\t")).collect();

    // Crushed: 3 segments (ACGT, A, CGTAC), 2 edges, 3 identical paths.
    assert_eq!(
        s_lines.len(),
        3,
        "expected 3 S lines after crush, got {}",
        s_lines.len()
    );
    assert_eq!(
        l_lines.len(),
        2,
        "expected 2 L lines after crush, got {}",
        l_lines.len()
    );
    assert_eq!(
        p_lines.len(),
        3,
        "expected 3 P lines, got {}",
        p_lines.len()
    );

    // No 'T' segment should remain (the SNP ALT was crushed out).
    let has_t_seg = s_lines.iter().any(|s| s.split('\t').nth(2) == Some("T"));
    assert!(
        !has_t_seg,
        "T allele segment should be crushed out: {s_lines:?}"
    );

    // All three paths should be identical (the ALT path was rewritten to REF).
    let paths: Vec<&str> = p_lines
        .iter()
        .map(|p| p.split('\t').nth(2).unwrap())
        .collect();
    let first = paths[0];
    assert!(
        paths.iter().all(|p| *p == first),
        "all paths should be identical after crush: {paths:?}"
    );

    let _ = fs::remove_file("/tmp/pgr_gfa_crush_A.fa");
    let _ = fs::remove_file("/tmp/pgr_gfa_crush_B.fa");
    let _ = fs::remove_file("/tmp/pgr_gfa_crush_C.fa");
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

// ── GFA path spelling round-trip ─────────────────────────────────
// Borrowed from impg `test_graph_poa.rs::assert_gfa_paths_match_records`:
// parse GFA P lines, concatenate segment sequences along each path
// (respecting orientation), and verify the spelled sequence matches the
// input FASTA record byte-for-byte. This is the strongest correctness
// invariant for graph construction — every path must reconstruct its
// original sequence.

/// Parse GFA and spell each path's sequence by concatenating visited
/// segments (reverse-complementing for `-` steps). Returns
/// `Vec<(path_name, spelled_sequence)>`.
fn spell_gfa_paths(gfa: &str) -> Vec<(String, String)> {
    use std::collections::HashMap;
    let mut seg_seq: HashMap<String, String> = HashMap::new();
    let mut paths: Vec<(String, Vec<(String, char)>)> = Vec::new();

    for line in gfa.lines() {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields[0] == "S" && fields.len() >= 3 {
            seg_seq.insert(fields[1].to_string(), fields[2].to_string());
        } else if fields[0] == "P" && fields.len() >= 3 {
            let name = fields[1].to_string();
            let steps: Vec<(String, char)> = fields[2]
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    let orient = s.chars().last().unwrap_or('+');
                    let id = s.trim_end_matches(['+', '-']).to_string();
                    (id, orient)
                })
                .collect();
            paths.push((name, steps));
        }
    }

    paths
        .into_iter()
        .map(|(name, steps)| {
            let mut spelled = String::new();
            for (id, orient) in steps {
                let seq = seg_seq.get(&id).cloned().unwrap_or_default();
                if orient == '-' {
                    spelled.push_str(&revcomp(&seq));
                } else {
                    spelled.push_str(&seq);
                }
            }
            (name, spelled)
        })
        .collect()
}

/// Reverse-complement a DNA string (ACGTN, case-insensitive).
fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            other => other,
        })
        .collect()
}

#[test]
fn command_paf_to_gfa_roundtrip_identical() {
    use std::collections::BTreeMap;
    use std::fs;
    // Two identical sequences -> one segment, two paths. Each path must
    // spell back the original 10-bp sequence.
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_gfa_rt_id_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_gfa_rt_id_B.fa", ">B\nACGTACGTAC\n");
    let tsv = "/tmp/pgr_gfa_rt_id.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-gfa", "stdin", "B:0-10", "-t", "-f", tsv])
        .stdin(paf)
        .run();

    let spelled: BTreeMap<String, String> = spell_gfa_paths(&stdout).into_iter().collect();
    assert_eq!(
        spelled.get("A"),
        Some(&"ACGTACGTAC".to_string()),
        "path A should spell back the original sequence"
    );
    assert_eq!(
        spelled.get("B"),
        Some(&"ACGTACGTAC".to_string()),
        "path B should spell back the original sequence"
    );

    let _ = fs::remove_file("/tmp/pgr_gfa_rt_id_A.fa");
    let _ = fs::remove_file("/tmp/pgr_gfa_rt_id_B.fa");
    let _ = fs::remove_file(&a_fa);
    let _ = fs::remove_file(&b_fa);
    let _ = fs::remove_file(format!("{a_fa}.gzi"));
    let _ = fs::remove_file(format!("{b_fa}.gzi"));
    let _ = fs::remove_file(tsv);
    let _ = fs::remove_file(format!("{a_fa}.loc"));
    let _ = fs::remove_file(format!("{b_fa}.loc"));
}

#[test]
fn command_paf_to_gfa_roundtrip_snp_bubble() {
    use std::collections::BTreeMap;
    use std::fs;
    // B = ACGTACGTAC, A = ACGTACGTAC, C = ACGTTCGTAC (SNP at pos 4).
    // The SNP forms a bubble; each path must still spell its own sequence.
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_gfa_rt_snp_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_gfa_rt_snp_B.fa", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_gfa_rt_snp_C.fa", ">C\nACGTTCGTAC\n");
    let tsv = "/tmp/pgr_gfa_rt_snp.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-gfa", "stdin", "B:0-10", "-t", "-f", tsv])
        .stdin(paf)
        .run();

    let spelled: BTreeMap<String, String> = spell_gfa_paths(&stdout).into_iter().collect();
    assert_eq!(
        spelled.get("A"),
        Some(&"ACGTACGTAC".to_string()),
        "path A should spell ACGTACGTAC"
    );
    assert_eq!(
        spelled.get("B"),
        Some(&"ACGTACGTAC".to_string()),
        "path B should spell ACGTACGTAC"
    );
    assert_eq!(
        spelled.get("C"),
        Some(&"ACGTTCGTAC".to_string()),
        "path C should spell ACGTTCGTAC (SNP allele preserved)"
    );

    let _ = fs::remove_file("/tmp/pgr_gfa_rt_snp_A.fa");
    let _ = fs::remove_file("/tmp/pgr_gfa_rt_snp_B.fa");
    let _ = fs::remove_file("/tmp/pgr_gfa_rt_snp_C.fa");
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
fn command_paf_to_gfa_roundtrip_indel_bubble() {
    use std::collections::BTreeMap;
    use std::fs;
    // B = ACGTACGTAC (10bp)
    // A = ACGTACGTAC (10bp, identical to B)
    // C = ACGTGGGTAC (10bp, 2bp substitution) — keep simple, indels in POA
    //   are harder to predict; use a 2bp insertion instead:
    // C = ACGTACGGGTAC (12bp, 2bp insertion after pos 6)
    // C-A alignment: 6= 2I 4= (C has 2bp insertion relative to A)
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t12\t0\t12\t+\tA\t10\t0\t10\t10\t12\t255\tcg:Z:6=2I4=\n";
    let a_fa = write_bgzf_fa("/tmp/pgr_gfa_rt_indel_A.fa", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa("/tmp/pgr_gfa_rt_indel_B.fa", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa("/tmp/pgr_gfa_rt_indel_C.fa", ">C\nACGTACGGGTAC\n");
    let tsv = "/tmp/pgr_gfa_rt_indel.tsv";
    fs::write(tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "to-gfa", "stdin", "B:0-10", "-t", "-f", tsv])
        .stdin(paf)
        .run();

    let spelled: BTreeMap<String, String> = spell_gfa_paths(&stdout).into_iter().collect();
    assert_eq!(
        spelled.get("A"),
        Some(&"ACGTACGTAC".to_string()),
        "path A should spell ACGTACGTAC (no insertion)"
    );
    assert_eq!(
        spelled.get("B"),
        Some(&"ACGTACGTAC".to_string()),
        "path B should spell ACGTACGTAC (no insertion)"
    );
    assert_eq!(
        spelled.get("C"),
        Some(&"ACGTACGGGTAC".to_string()),
        "path C should spell ACGTACGGGTAC (2bp insertion preserved)"
    );

    let _ = fs::remove_file("/tmp/pgr_gfa_rt_indel_A.fa");
    let _ = fs::remove_file("/tmp/pgr_gfa_rt_indel_B.fa");
    let _ = fs::remove_file("/tmp/pgr_gfa_rt_indel_C.fa");
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
