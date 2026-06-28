#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── paf to-bed (BED3 output from CIGAR) ──────────────────────────

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
