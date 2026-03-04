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

#[test]
fn command_mat_transform_linear() {
    // Input: A-B=0.1
    // Linear: x*2 + 1
    // Output: A-B=1.2
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "mat",
            "transform",
            "tests/mat/IBPA.phy",
            "--op",
            "linear",
            "--scale",
            "2.0",
            "--offset",
            "1.0",
        ])
        .run();

    // Original: IBPA_ECOLI vs IBPA_ECOLI_GA is 0.058394
    // Transformed: 0.058394 * 2 + 1 = 1.116788
    assert!(stdout.contains("1.116788"));
}

#[test]
fn command_mat_transform_inv_linear() {
    // Input: A-B=0.1
    // Inv-linear: 1.0 - x
    // Output: A-B=0.9
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "mat",
            "transform",
            "tests/mat/IBPA.phy",
            "--op",
            "inv-linear",
            "--max",
            "1.0",
        ])
        .run();

    // Original: IBPA_ECOLI vs IBPA_ECOLI_GA is 0.058394
    // Transformed: 1.0 - 0.058394 = 0.941606
    assert!(stdout.contains("0.941606"));
}

#[test]
fn command_mat_transform_log() {
    // Input: A-B=0.1
    // Log: -ln(x)
    // Output: -ln(0.1) = 2.302585
    let (stdout, _) = PgrCmd::new()
        .args(&["mat", "transform", "tests/mat/IBPA.phy", "--op", "log"])
        .run();

    // Original: IBPA_ECOLI vs IBPA_ECOLI_GA is 0.058394
    // Transformed: -ln(0.058394) = 2.8405
    assert!(stdout.contains("2.8405"));
}

#[test]
fn command_mat_transform_normalize() {
    // Create a dummy matrix with non-zero diagonals for testing normalization
    // 3
    // A 1.0 0.5 0.5
    // B 0.5 4.0 1.0
    // C 0.5 1.0 9.0
    //
    // Norm(A,B) = 0.5 / sqrt(1*4) = 0.5/2 = 0.25
    // Norm(B,C) = 1.0 / sqrt(4*9) = 1.0/6 = 0.166667
    // Norm(A,C) = 0.5 / sqrt(1*9) = 0.5/3 = 0.166667

    let input = "3\nA\t1.0\t0.5\t0.5\nB\t0.5\t4.0\t1.0\nC\t0.5\t1.0\t9.0\n";

    let (stdout, _) = PgrCmd::new()
        .args(&["mat", "transform", "stdin", "--normalize"])
        .stdin(input)
        .run();

    // Check normalized values
    assert!(stdout.contains("0.250000")); // A-B
    assert!(stdout.contains("0.166667")); // B-C or A-C
    assert!(stdout.contains("1.000000")); // Diagonals
}

#[test]
fn command_mat_transform_normalize_inv() {
    // Combine normalize and inv-linear (Sim -> Dist)
    // Input same as above.
    // Norm(A,B) = 0.25
    // Inv: 1.0 - 0.25 = 0.75

    let input = "3\nA\t1.0\t0.5\t0.5\nB\t0.5\t4.0\t1.0\nC\t0.5\t1.0\t9.0\n";

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "mat",
            "transform",
            "stdin",
            "--normalize",
            "--op",
            "inv-linear",
        ])
        .stdin(input)
        .run();

    assert!(stdout.contains("0.750000")); // 1 - 0.25
    assert!(stdout.contains("0.000000")); // Diagonals: 1 - 1.0
}
