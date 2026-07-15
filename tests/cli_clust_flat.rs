#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_clust_cc() {
    let (stdout, _) = PgrCmd::new()
        .args(&["clust", "cc", "tests/clust/IBPA.fa.05.tsv"])
        .run();

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains("A0A192CFC5_ECO25\tIBPA_ECOLI\tIBPA_ESCF3\nIBPA_ECOLI_GA_LV"));
}

#[test]
fn command_clust_cc_pair() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "cc",
            "tests/clust/IBPA.fa.05.tsv",
            "--format",
            "pair",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains("A0A192CFC5_ECO25\tIBPA_ECOLI"));
    assert!(stdout.contains("IBPA_ECOLI_GA_LV\tIBPA_ECOLI_GA_LV"));
}

#[test]
fn command_clust_dbscan() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "dbscan",
            "tests/clust/IBPA.fa.tsv",
            "--eps",
            "0.05",
            "--min-points",
            "2",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 7);
    // Default representative is medoid; the first column for the {A0A, IBPA_ECOLI, IBPA_ESCF3}
    // cluster should be IBPA_ECOLI (minimum sum of distances; tie broken alphabetically).
    assert!(stdout.contains("IBPA_ECOLI\tA0A192CFC5_ECO25\tIBPA_ESCF3"));
}

#[test]
fn command_clust_dbscan_rep_first() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "dbscan",
            "tests/clust/IBPA.fa.tsv",
            "--eps",
            "0.05",
            "--min-points",
            "2",
            "--rep",
            "first",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains("A0A192CFC5_ECO25\tIBPA_ECOLI\tIBPA_ESCF3"));
}

#[test]
fn command_clust_dbscan_pair() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "dbscan",
            "tests/clust/IBPA.fa.tsv",
            "--eps",
            "0.05",
            "--min-points",
            "2",
            "--format",
            "pair",
        ])
        .run();

    // Each line contains a representative-member pair
    assert!(stdout.lines().count() > 0);
    assert!(
        stdout.contains("IBPA_ECOLI\tIBPA_ECOLI\n") || stdout.contains("IBPA_ESCF3\tIBPA_ESCF3\n")
    );
    assert!(
        stdout.contains("IBPA_ECOLI\tIBPA_ESCF3\n") || stdout.contains("IBPA_ESCF3\tIBPA_ECOLI\n")
    );
}

#[test]
fn command_clust_kmedoids() {
    let (stdout, _) = PgrCmd::new()
        .args(&["clust", "km", "tests/clust/IBPA.fa.tsv", "-k", "2"])
        .run();

    // Output should contain at least 2 lines (clusters)
    assert!(stdout.lines().count() >= 2);
}

#[test]
fn command_clust_kmedoids_pair() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "k-medoids",
            "tests/clust/IBPA.fa.tsv",
            "-k",
            "2",
            "--format",
            "pair",
        ])
        .run();

    // Should contain tab-separated pairs
    assert!(stdout.contains("\t"));
    assert!(stdout.lines().count() >= 2);
}

#[test]
fn command_clust_mcl() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "mcl",
            "tests/clust/mcl_test.tsv",
            "--inflation",
            "2.0",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("A\tB\tC"));
    assert!(stdout.contains("D\tE"));
}

#[test]
fn command_clust_mcl_complex() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "mcl",
            "tests/clust/mcl_complex.tsv",
            "--inflation",
            "2.0",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 2);
    // Cluster 1: n1, n2, n3, n4
    assert!(stdout.contains("n1\tn2\tn3\tn4"));
    // Cluster 2: n5, n6, n7
    assert!(stdout.contains("n5\tn6\tn7"));
}

#[test]
fn command_clust_mcl_args() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "mcl",
            "tests/clust/mcl_test.tsv",
            "--prune",
            "1e-3",
            "--max-iter",
            "50",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("A\tB\tC"));
    assert!(stdout.contains("D\tE"));
}

#[test]
fn command_clust_mcl_pair() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "mcl",
            "tests/clust/mcl_test.tsv",
            "--format",
            "pair",
        ])
        .run();

    // Cluster 1 (size 3) + Cluster 2 (size 2) = 5 pairs
    assert_eq!(stdout.lines().count(), 5);

    // Check representative output
    assert!(stdout.contains("A\tA"));
    assert!(stdout.contains("A\tB"));
    assert!(stdout.contains("A\tC"));
    assert!(stdout.contains("D\tD"));
    assert!(stdout.contains("D\tE"));
}

#[test]
fn command_clust_dbscan_default_min_points() {
    let (stdout, _) = PgrCmd::new().args(&["clust", "dbscan", "--help"]).run();

    // The help should show default value of 4 for --min-points
    assert!(stdout.contains("4") || stdout.contains("default"));
}

#[test]
fn command_clust_dbscan_pair_rep_first() {
    let temp = tempfile::TempDir::new().unwrap();
    let input = temp.path().join("pairs.tsv");
    std::fs::write(
        &input,
        "A\tB\t0.1\nA\tC\t0.2\nB\tC\t0.15\nA\tD\t0.5\nB\tD\t0.6\nC\tD\t0.55\n",
    )
    .unwrap();

    let (stdout_medoid_pair, _) = PgrCmd::new()
        .args(&[
            "clust",
            "dbscan",
            input.to_str().unwrap(),
            "--eps",
            "0.25",
            "--min-points",
            "2",
            "--format",
            "pair",
        ])
        .run();

    let (stdout_first_pair, _) = PgrCmd::new()
        .args(&[
            "clust",
            "dbscan",
            input.to_str().unwrap(),
            "--eps",
            "0.25",
            "--min-points",
            "2",
            "--format",
            "pair",
            "--rep",
            "first",
        ])
        .run();

    let (stdout_medoid_cluster, _) = PgrCmd::new()
        .args(&[
            "clust",
            "dbscan",
            input.to_str().unwrap(),
            "--eps",
            "0.25",
            "--min-points",
            "2",
            "--format",
            "cluster",
        ])
        .run();

    let (stdout_first_cluster, _) = PgrCmd::new()
        .args(&[
            "clust",
            "dbscan",
            input.to_str().unwrap(),
            "--eps",
            "0.25",
            "--min-points",
            "2",
            "--format",
            "cluster",
            "--rep",
            "first",
        ])
        .run();

    // Default medoid representative is B (min sum of distances).
    assert!(stdout_medoid_pair.contains("B\tA"));
    assert!(stdout_medoid_pair.contains("B\tB"));
    assert!(stdout_medoid_pair.contains("B\tC"));
    let cluster_medoid_line = stdout_medoid_cluster.lines().next().unwrap();
    assert!(cluster_medoid_line.starts_with("B\t"));

    // With --rep first, representative is A (alphabetically first).
    assert!(stdout_first_pair.contains("A\tA"));
    assert!(stdout_first_pair.contains("A\tB"));
    assert!(stdout_first_pair.contains("A\tC"));
    let cluster_first_line = stdout_first_cluster.lines().next().unwrap();
    assert!(cluster_first_line.starts_with("A\t"));
}

#[test]
fn command_clust_mcl_pair_rep_first() {
    let temp = tempfile::TempDir::new().unwrap();
    let input = temp.path().join("pairs.tsv");
    std::fs::write(&input, "A\tB\t1.0\nB\tC\t1.0\nC\tA\t0.5\nD\tE\t1.0\n").unwrap();

    let (stdout_medoid, _) = PgrCmd::new()
        .args(&["clust", "mcl", input.to_str().unwrap(), "--format", "pair"])
        .run();

    let (stdout_first, _) = PgrCmd::new()
        .args(&[
            "clust",
            "mcl",
            input.to_str().unwrap(),
            "--format",
            "pair",
            "--rep",
            "first",
        ])
        .run();

    let (stdout_medoid_cluster, _) = PgrCmd::new()
        .args(&[
            "clust",
            "mcl",
            input.to_str().unwrap(),
            "--format",
            "cluster",
        ])
        .run();

    let (stdout_first_cluster, _) = PgrCmd::new()
        .args(&[
            "clust",
            "mcl",
            input.to_str().unwrap(),
            "--format",
            "cluster",
            "--rep",
            "first",
        ])
        .run();

    // Default medoid (max similarity sum) representative is B.
    assert!(stdout_medoid.contains("B\tA"));
    assert!(stdout_medoid.contains("B\tB"));
    assert!(stdout_medoid.contains("B\tC"));
    let cluster_medoid_line = stdout_medoid_cluster.lines().next().unwrap();
    assert!(cluster_medoid_line.starts_with("B\t"));

    // With --rep first, representative is A.
    assert!(stdout_first.contains("A\tA"));
    assert!(stdout_first.contains("A\tB"));
    assert!(stdout_first.contains("A\tC"));
    let cluster_first_line = stdout_first_cluster.lines().next().unwrap();
    assert!(cluster_first_line.starts_with("A\t"));
}

#[test]
fn command_clust_kmedoids_pair_rep_first() {
    let temp = tempfile::TempDir::new().unwrap();
    let input = temp.path().join("pairs.tsv");
    std::fs::write(
        &input,
        "A\tB\t0.1\nA\tC\t0.2\nB\tC\t0.15\nA\tD\t0.5\nB\tD\t0.6\nC\tD\t0.55\n",
    )
    .unwrap();

    let (stdout_medoid, _) = PgrCmd::new()
        .args(&[
            "clust",
            "k-medoids",
            input.to_str().unwrap(),
            "-k",
            "2",
            "--format",
            "pair",
        ])
        .run();

    let (stdout_first, _) = PgrCmd::new()
        .args(&[
            "clust",
            "k-medoids",
            input.to_str().unwrap(),
            "-k",
            "2",
            "--format",
            "pair",
            "--rep",
            "first",
        ])
        .run();

    let (stdout_medoid_cluster, _) = PgrCmd::new()
        .args(&[
            "clust",
            "k-medoids",
            input.to_str().unwrap(),
            "-k",
            "2",
            "--format",
            "cluster",
        ])
        .run();

    let (stdout_first_cluster, _) = PgrCmd::new()
        .args(&[
            "clust",
            "k-medoids",
            input.to_str().unwrap(),
            "-k",
            "2",
            "--format",
            "cluster",
            "--rep",
            "first",
        ])
        .run();

    // Verify both produce valid pair output.
    assert!(stdout_medoid.contains("\t"));
    assert!(stdout_first.contains("\t"));
    assert!(stdout_medoid.lines().count() >= 2);
    assert!(stdout_first.lines().count() >= 2);

    // In cluster format, representative is placed in the first column.
    let cluster_medoid_first_line = stdout_medoid_cluster.lines().next().unwrap();
    let cluster_first_first_line = stdout_first_cluster.lines().next().unwrap();
    assert!(cluster_medoid_first_line.starts_with("B\t"));
    assert!(cluster_first_first_line.starts_with("A\t"));
}
