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
