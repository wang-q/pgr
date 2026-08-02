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
fn command_plot_dot() {
    let (stdout, _) = PgrCmd::new()
        .args(&["plot", "dot", "tests/plot/dot.paf", "--min-len", "400"])
        .run();

    assert!(stdout.starts_with("<?xml"));
    assert!(stdout.contains(r#"<svg xmlns="http://www.w3.org/2000/svg""#));
    // 4 records; two pass --min-len 400 (block 800 and 1250)
    let seg = stdout.split(r#"id="segments""#).nth(1).unwrap();
    assert_eq!(seg.matches("<line ").count(), 2);
    assert!(stdout.contains(r#"id="segments""#));
    // 2 identity-colored line strokes plus the border rect stroke
    assert!(stdout.matches("stroke=\"#").count() >= 3);
}

#[test]
fn command_plot_dot_filters_everything() {
    let (_, stderr) = PgrCmd::new()
        .args(&["plot", "dot", "tests/plot/dot.paf", "--min-len", "100000"])
        .run_fail();

    assert!(stderr.contains("no alignments pass filters"));
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

#[test]
fn command_plot_nrps() {
    let (stdout, _) = PgrCmd::new()
        .args(&["plot", "nrps", "tests/plot/srf.tsv"])
        .run();

    assert!(stdout.contains("(-0.4cm,0) -- (\\x1 + 0.2cm,0)"));
    assert!(!stdout.contains("\\textbf{M}ethyltransferase"));

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "plot",
            "nrps",
            "tests/plot/srf.tsv",
            "--legend",
            "--color",
            "black",
        ])
        .run();

    assert!(stdout.contains("     draw=black,"));
    assert!(stdout.contains("\\textbf{M}ethyltransferase"));
}
