#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_maf_to_fas() {
    let (stdout, _) = PgrCmd::new()
        .args(&["maf", "to-fas", "tests/maf/example.maf"])
        .run();

    assert!(stdout.contains(">S288c.VIII(+):13377-13410"));
    assert!(stdout.contains("TTACTCGTCTTGCGGCCAAAACTCGAAGAAAAAC"));
    assert!(stdout.contains(">Spar.gi_29362578(-):72853-72885"));
    assert!(stdout.contains("TTACCCGTCTTGCGTCCAAAACTCGAA-AAAAAC"));
    assert_eq!(stdout.matches(">").count(), 8); // 2 blocks * 4 sequences
    assert_eq!(stdout.lines().count(), 18);
    assert!(stdout.contains("S288c.VIII"), "name list");
    assert!(stdout.contains(":42072-42168"), "coordinate transformed");
}

#[test]
fn command_maf_to_paf_basic() {
    let maf = "##maf version=1 scoring=multiz\n\
               a score=12345\n\
               s ref  100 10 + 1000 ACGTACGT--\n\
               s qry   50 10 +  500 ACGTACGTAA\n";

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["maf", "to-paf", "stdin"])
        .stdin(maf)
        .run();

    // 12 mandatory columns
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields.len(), 16, "PAF should have 12 cols + 4 tags = 16 fields");

    // Column checks
    assert_eq!(fields[0], "qry");           // query name
    assert_eq!(fields[1], "500");           // query length
    assert_eq!(fields[2], "50");            // query start
    assert_eq!(fields[3], "60");            // query end (50 + 10)
    assert_eq!(fields[4], "+");             // strand
    assert_eq!(fields[5], "ref");           // target name
    assert_eq!(fields[6], "1000");          // target length
    assert_eq!(fields[7], "100");           // target start
    assert_eq!(fields[8], "110");           // target end (100 + 10)
    assert_eq!(fields[9], "8");             // matches (8M)
    assert_eq!(fields[10], "10");           // block length
    assert_eq!(fields[11], "255");          // mapq

    // Custom tags
    assert!(stdout.contains("gi:f:"), "gi tag missing");
    assert!(stdout.contains("bi:f:"), "bi tag missing");
    assert!(stdout.contains("cg:Z:"), "cg tag missing");
    assert!(stdout.contains("ms:i:12345"), "score tag missing");
}

#[test]
fn command_maf_to_paf_multi_sequence_skipped() {
    // Multi-sequence MAF — should print warnings on stderr, nothing on stdout
    let maf = "##maf version=1\n\
               a score=100\n\
               s ref  100 10 + 1000 ACGTACGT--\n\
               s qry1  50 10 +  500 ACGTACGTAA\n\
               s qry2  30 10 +  300 ACGTACGTAA\n";

    let (stdout, stderr) = PgrCmd::new()
        .args(&["maf", "to-paf", "stdin"])
        .stdin(maf)
        .run();

    assert!(stdout.trim().is_empty(), "multi-seq should produce no output");
    assert!(stderr.contains("skipping"), "should warn about skipping");
    assert!(stderr.contains("3 sequences"), "should mention count");
}

#[test]
fn command_maf_to_paf_strand_minus() {
    // Reverse strand query
    let maf = "##maf version=1\n\
               a score=50\n\
               s ref  100 5 + 1000 ACGTA\n\
               s qry   50 5 -  500 ACGTA\n";

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["maf", "to-paf", "stdin"])
        .stdin(maf)
        .run();

    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields[0], "qry");
    assert_eq!(fields[4], "-", "reverse strand should be '-'");
    assert_eq!(fields[9], "5", "all 5 bases match");
    assert_eq!(fields[10], "5", "block length");
}

#[test]
fn command_maf_to_paf_with_gaps() {
    // Alignment with both insertion and deletion
    let maf = "##maf version=1\n\
               a score=200\n\
               s ref  100 6 + 1000 ACG-TG\n\
               s qry   50 6 +  500 ACGTT-\n";

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["maf", "to-paf", "stdin"])
        .stdin(maf)
        .run();

    assert!(stdout.contains("cg:Z:3M1I1M1D"), "CIGAR should be 3M1I1M1D");
    // gi: 4/(4+0+2) = 0.666..., bi: 4/(4+0+2) = 0.666...
    assert!(stdout.contains("gi:f:"), "gi tag present");
    assert!(stdout.contains("bi:f:"), "bi tag present");
    assert!(stdout.contains("ms:i:200"), "score tag present");
}
