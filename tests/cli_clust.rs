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
            "--min_points",
            "2",
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
            "--min_points",
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
            "--max_iter",
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
            "--inflation",
            "1.2",
            "--format",
            "pair",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 9);
    // Should contain self-loops for centers and connections
    assert!(stdout.contains("0\t5\n") || stdout.contains("5\t0\n"));
}

#[test]
fn command_clust_mcl_weighted() {
    // 0-1 (weight 2), 1-2 (weight 1)
    // inflation 1.2
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "mcl",
            "tests/clust/mcl_weighted.tsv",
            "--inflation",
            "1.2",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("0\t1\t2"));
}
