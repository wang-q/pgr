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
