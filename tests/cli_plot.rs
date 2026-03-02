#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_plot_venn2() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "plot",
            "venn",
            "tests/plot/rocauc.result.tsv",
            "tests/plot/mcox.05.result.tsv",
        ])
        .run();

    assert!(stdout.contains("(-2.8, -1.8) { rocauc }"));
    assert!(stdout.contains("(-2,    0) { 669 }"));
}

#[test]
fn command_plot_venn3() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "plot",
            "venn",
            "tests/plot/rocauc.result.tsv",
            "tests/plot/mcox.05.result.tsv",
            "tests/plot/mcox.result.tsv",
        ])
        .run();

    assert!(stdout.contains("(-2.8, -1.8) { rocauc }"));
    assert!(stdout.contains("(-2,   -0.2) { 161 }"));
}

#[test]
fn command_plot_venn4() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "plot",
            "venn",
            "tests/plot/rocauc.result.tsv",
            "tests/plot/rocauc.result.tsv",
            "tests/plot/mcox.05.result.tsv",
            "tests/plot/mcox.result.tsv",
        ])
        .run();

    assert!(stdout.contains("(-2.2, -2.6) { rocauc }"));
    assert!(stdout.contains("(-2.2,  1.5) { 161 }"));
}

#[test]
fn command_plot_hh() {
    let (stdout, _) = PgrCmd::new()
        .args(&["plot", "hh", "tests/plot/hist.tsv"])
        .run();

    assert!(stdout.contains("31   0 0.0200"));
    assert!(stdout.contains("31   1 0.0000"));

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "plot",
            "hh",
            "tests/plot/hist.tsv",
            "-g",
            "2",
            "--bins",
            "20",
        ])
        .run();

    assert!(stdout.contains("11   0 0.0600"));
    assert!(stdout.contains("11   1 0.1600"));
}
