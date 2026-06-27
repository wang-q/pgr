#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_paf_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "--help"]).run();
    assert!(stdout.contains("Manipulate PAF"));
    assert!(stdout.contains("index"));
}

#[test]
fn command_paf_index_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "index", "--help"]).run();
    assert!(stdout.contains("Build interval-tree index"));
    assert!(stdout.contains("infiles"));
}

#[test]
fn command_paf_index_single_file() {
    let paf = "\
q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255\tcg:Z:50M\tgi:f:0.9
q2\t300\t10\t60\t-\tt1\t200\t10\t60\t45\t50\t255\tcg:Z:50M
";
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", "stdin"])
        .stdin(paf)
        .run();
    assert!(stderr.contains("sequences: 3"));
    assert!(stderr.contains("targets:   1"));
}

#[test]
fn command_paf_index_no_cigar() {
    let paf = "\
q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255
q2\t300\t10\t60\t+\tt2\t400\t10\t60\t45\t50\t255
";
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", "stdin"])
        .stdin(paf)
        .run();
    assert!(stderr.contains("sequences: 4"));
    assert!(stderr.contains("targets:   2"));
}

#[test]
fn command_paf_index_empty() {
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", "stdin"])
        .stdin("")
        .run();
    assert!(stderr.contains("sequences: 0"));
    assert!(stderr.contains("targets:   0"));
}

#[test]
fn command_paf_index_comments_and_blanks() {
    let paf = "\
# header comment

q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255\tcg:Z:50M

# another comment
q2\t300\t10\t60\t-\tt1\t200\t10\t60\t45\t50\t255\tcg:Z:50M
";
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", "stdin"])
        .stdin(paf)
        .run();
    assert!(stderr.contains("sequences: 3"));
    assert!(stderr.contains("targets:   1"));
}

#[test]
fn command_paf_index_invalid() {
    PgrCmd::new()
        .args(&["paf", "index", "stdin"])
        .stdin("invalid line\n")
        .run_fail();
}
