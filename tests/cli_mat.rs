#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_mat_to_phylip() {
    let (stdout, _) = PgrCmd::new()
        .args(&["mat", "to-phylip", "tests/mat/IBPA.fa.tsv"])
        .run();

    assert_eq!(stdout.lines().count(), 11);
    assert!(stdout.contains("IBPA_ECOLI\t0\t0.0669"));
}

#[test]
fn command_mat_to_pair() {
    let (stdout, _) = PgrCmd::new()
        .args(&["mat", "to-pair", "tests/mat/IBPA.phy"])
        .run();

    assert_eq!(stdout.lines().count(), 55);
    assert!(stdout.contains("IBPA_ECOLI\tIBPA_ECOLI\t0\n"));
    assert!(stdout.contains("IBPA_ECOLI\tIBPA_ECOLI_GA\t0.058"));
}

#[test]
fn command_mat_format_full() {
    let (stdout, _) = PgrCmd::new()
        .args(&["mat", "format", "tests/mat/IBPA.phy"])
        .run();

    assert_eq!(stdout.lines().count(), 11);
    assert!(stdout.contains("IBPA_ECOLI\t0\t0.058394\t0.160584"));
    assert!(stdout.contains("IBPA_ECOLI_GA\t0.058394\t0\t0.10219"));
}

#[test]
fn command_mat_format_lower() {
    let (stdout, _) = PgrCmd::new()
        .args(&["mat", "format", "tests/mat/IBPA.phy", "--mode", "lower"])
        .run();

    assert_eq!(stdout.lines().count(), 11);
    assert!(stdout.contains("IBPA_ECOLI\n"));
    assert!(stdout.contains("IBPA_ECOLI_GA\t0.058394\n"));
}

#[test]
fn command_mat_format_strict() {
    let (stdout, _) = PgrCmd::new()
        .args(&["mat", "format", "tests/mat/IBPA.phy", "--mode", "strict"])
        .run();

    assert_eq!(stdout.lines().count(), 11);

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines[0].trim(), "10"); // Number of sequences line

    // Check format of the first sequence
    let first_seq = lines[1];
    assert!(first_seq.starts_with("IBPA_ECOLI"));
    assert_eq!(first_seq.chars().take(10).count(), 10); // Name length limit
    assert!(first_seq.contains(" 0.000000")); // Formatted distance value
}

#[test]
fn command_mat_subset() {
    let (stdout, _) = PgrCmd::new()
        .args(&["mat", "subset", "tests/mat/IBPA.phy", "tests/mat/IBPA.list"])
        .run();

    // Verify output
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines[0].trim(), "3"); // Number of sequences
    assert!(lines[1].starts_with("IBPA_ECOLI_GA\t0\t0.10219\t0.058394"));
    assert!(lines[3].starts_with("IBPA_ESCF3\t0.058394"));
}

#[test]
fn command_mat_compare() {
    // Test single method
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "mat",
            "compare",
            "tests/mat/IBPA.phy",
            "tests/mat/IBPA.71.phy",
            "--method",
            "pearson",
        ])
        .run();

    // Verify output format and approximate value
    assert!(stdout.contains("Method\tScore"));
    assert!(stdout.contains("pearson\t0.93"));

    // Test all methods
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "mat",
            "compare",
            "tests/mat/IBPA.phy",
            "tests/mat/IBPA.71.phy",
            "--method",
            "all",
        ])
        .run();

    // Verify matrix information in stderr
    assert!(stderr.contains("Sequences in matrices: 10 and 10"));
    assert!(stderr.contains("Common sequences: 10"));

    // Verify all methods are present with approximate values
    assert!(stdout.contains("pearson\t0.93"));
    assert!(stdout.contains("spearman\t0.91"));
    assert!(stdout.contains("mae\t0.11"));
    assert!(stdout.contains("cosine\t0.97"));
    assert!(stdout.contains("jaccard\t0.75"));
    assert!(stdout.contains("euclid\t1.22"));
}
