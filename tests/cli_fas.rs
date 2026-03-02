#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use tempfile::TempDir;

#[test]
fn command_invalid() {
    let (_, stderr) = PgrCmd::new().args(&["fas", "foobar"]).run_fail();
    assert!(stderr.contains("recognized"));
}

#[test]
fn command_name() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "name", "tests/fas/example.fas", "-c"])
        .run();

    assert_eq!(stdout.lines().count(), 4);
    assert!(stdout.contains("S288c\t3"), "count");
    assert!(stdout.contains("S288c\t3\nYJM789\t3\nRM11"), "name order");
}

#[test]
fn command_cover() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "cover", "tests/fas/example.fas"])
        .run();

    assert_eq!(stdout.lines().count(), 16);
    assert!(stdout.contains("S288c"), "name list");
    assert!(stdout.contains("I"), "chr list");
    assert!(stdout.contains("13267-13287"), "runlist");

    // --name, --trim
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "cover",
            "tests/fas/example.fas",
            "--name",
            "S288c",
            "--trim",
            "10",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 3);
    assert!(!stdout.contains("S288c"), "name list");
    assert!(stdout.contains("I"), "chr list");
    assert!(stdout.contains("13277,184906"), "trimmed");
}

#[test]
fn command_concat() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "concat",
            "tests/fas/example.fas",
            "-r",
            "tests/fas/name.lst",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 4);
    assert_eq!(stdout.lines().next().unwrap().len(), 5); // >Spar
    assert_eq!(stdout.lines().last().unwrap().len(), 239);
    assert!(stdout.contains("Spar"), "name list");
    assert!(!stdout.contains("S288c"), "name list");
}

#[test]
fn command_concat_phylip() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "concat",
            "tests/fas/example.fas",
            "-r",
            "tests/fas/name.lst",
            "--phylip",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 3);
    assert_eq!(
        stdout.lines().last().unwrap().len(),
        "YJM789".to_string().len() + 1 + 239
    );
}

#[test]
fn command_subset() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "subset",
            "tests/fas/example.fas",
            "-r",
            "tests/fas/name.lst",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 15);
    assert!(stdout.lines().next().unwrap().contains("Spar")); // >Spar.

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "subset",
            "tests/fas/example.fas",
            "-r",
            "tests/fas/name.lst",
            "--strict",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 15);
    assert!(stdout.lines().next().unwrap().contains("Spar")); // >Spar.
}

#[test]
fn command_link() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "link", "tests/fas/example.fas"])
        .run();

    assert_eq!(stdout.lines().count(), 3);
    assert_eq!(stdout.lines().next().unwrap().split_whitespace().count(), 4);

    // --pair
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "link", "tests/fas/example.fas", "--pair"])
        .run();

    assert_eq!(stdout.lines().count(), 18);
    assert_eq!(stdout.lines().next().unwrap().split_whitespace().count(), 2);

    // --best
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "link", "tests/fas/example.fas", "--best"])
        .run();

    assert_eq!(stdout.lines().count(), 9);
    assert_eq!(stdout.lines().next().unwrap().split_whitespace().count(), 2);
}

#[test]
fn command_replace() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "replace",
            "tests/fas/example.fas",
            "-r",
            "tests/fas/replace.tsv",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 36);
    assert!(stdout.contains(">query.VIII(+)"));

    // fail
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "replace",
            "tests/fas/example.fas",
            "-r",
            "tests/fas/replace.fail.tsv",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 27);
    assert!(!stdout.contains("query"), "not replaced");
    assert!(
        stderr.contains("records") || stderr.contains("multiple records"),
        "error message"
    );

    // remove
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "replace",
            "tests/fas/example.fas",
            "-r",
            "tests/fas/replace.remove.tsv",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 18);
    assert!(!stdout.contains("13267-13287"), "block removed");
}

#[test]
fn command_check() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "check",
            "tests/fas/A_tha.pair.fas",
            "-r",
            "tests/fas/NC_000932.fa",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 3);
    assert!(stdout.lines().next().unwrap().contains("\tOK"));
    assert!(stdout.lines().last().unwrap().contains("\tFAILED"));

    // --name
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "check",
            "tests/fas/A_tha.pair.fas",
            "-r",
            "tests/fas/NC_000932.fa",
            "--name",
            "A_tha",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.lines().next().unwrap().contains("\tOK"));
    assert!(stdout.lines().last().unwrap().contains("\tOK"));
}

#[test]
fn command_create() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "create",
            "tests/fas/I.connect.tsv",
            "-r",
            "tests/fas/genome.fa",
            "--name",
            "S288c",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains("tgtgtgggtgtggtgtgg"), "revcom sequences");
    assert!(stdout.lines().next().unwrap().contains(">S288c."));
}

#[test]
fn command_separate() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "separate", "tests/fas/example.fas", "--rc"])
        .run();

    assert_eq!(stdout.lines().count(), 24);
    assert_eq!(
        stdout.lines().last().unwrap().len(),
        57,
        "length after remove dashes"
    );
    assert!(!stdout.contains("(-)"), "all strands are +");
    assert!(!stdout.contains("T-C"), "no dash, line 24");
}

#[test]
fn command_separate_to() {
    let tempdir = TempDir::new().unwrap();
    let tempdir_str = tempdir.path().to_str().unwrap();

    PgrCmd::new()
        .args(&[
            "fas",
            "separate",
            "tests/fas/example.fas",
            "--suffix",
            ".tmp",
            "-o",
            tempdir_str,
        ])
        .assert()
        .success()
        .stdout(predicates::str::is_empty());

    assert!(&tempdir.path().join("S288c.tmp").is_file());
    assert!(!&tempdir.path().join("YJM789.fasta").exists());

    tempdir.close().unwrap();
}

#[test]
fn command_split() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "split", "tests/fas/example.fas"])
        .run();

    assert_eq!(stdout.lines().count(), 27);

    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "split", "tests/fas/example.fas", "--simple"])
        .run();

    assert!(stdout.contains(">S288c\n"), "simple headers");
    assert!(!stdout.contains("I(+)"), "no positions");
}

#[test]
fn command_split_to() {
    let tempdir = TempDir::new().unwrap();
    let tempdir_str = tempdir.path().to_str().unwrap();

    PgrCmd::new()
        .args(&[
            "fas",
            "split",
            "tests/fas/example.fas",
            "--suffix",
            ".tmp",
            "--chr",
            "-o",
            tempdir_str,
        ])
        .assert()
        .success()
        .stdout(predicates::str::is_empty());

    assert!(&tempdir.path().join("S288c.I.tmp").is_file());
    assert!(!&tempdir.path().join("YJM789.fasta").exists());

    tempdir.close().unwrap();
}

#[test]
fn command_refine() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "refine", "tests/fas/example.fas", "--msa", "none"])
        .run();

    assert_eq!(stdout.lines().count(), 27);

    // --parallel 2
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "refine",
            "tests/fas/example.fas",
            "--msa",
            "none",
            "-p",
            "2",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 27);

    // --parallel 2
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "refine",
            "tests/fas/refine2.fas",
            "--msa",
            "none",
            "-p",
            "2",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 7);

    // --chop 10
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "refine",
            "tests/fas/example.fas",
            "--msa",
            "none",
            "--chop",
            "10",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 27);
    assert!(stdout.contains("185276-185332"), "new header"); // 185273-185334
    assert!(stdout.contains("156668-156724"), "new header"); // 156665-156726
    assert!(stdout.contains("3670-3727"), "new header"); // (-):3668-3730
    assert!(stdout.contains("2102-2159"), "new header"); // (-):2102-2161
}

#[test]
fn command_join() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "join",
            "tests/fas/S288cvsSpar.slice.fas",
            "--name",
            "Spar",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 5);
    assert!(
        stdout.lines().next().unwrap().contains(">Spar"),
        "Selected name first"
    );

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "join",
            "tests/fas/S288cvsRM11_1a.slice.fas",
            "tests/fas/S288cvsYJM789.slice.fas",
            "tests/fas/S288cvsSpar.slice.fas",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 9);
    assert!(
        stdout.lines().next().unwrap().contains(">S288c."),
        "First name first"
    );
}

#[test]
fn command_slice() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "slice",
            "tests/fas/slice.fas",
            "-r",
            "tests/fas/slice.json",
            "--name",
            "S288c",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains("13301-13400"), "sliced S288c");
    assert!(stdout.contains("2511-2636"), "sliced Spar");
    assert!(stdout.contains("\nTAGTCATCTCAG"), "sliced S288c seq");
}

#[test]
fn command_stat() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "stat", "tests/fas/example.fas"])
        .run();

    assert_eq!(stdout.lines().count(), 4);
    assert!(stdout.contains("0.192\t6\n"), "all together");

    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "stat", "tests/fas/example.fas", "--outgroup"])
        .run();

    assert_eq!(stdout.lines().count(), 4);
    assert!(stdout.contains("0.12\t3\n"), "exclude outgroup");
}

#[test]
fn command_filter() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "filter", "tests/fas/example.fas"])
        .run();

    assert_eq!(stdout.lines().count(), 27);

    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "filter", "tests/fas/example.fas", "--ge", "30"])
        .run();

    assert_eq!(stdout.lines().count(), 18);

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "filter",
            "tests/fas/example.fas",
            "--ge",
            "30",
            "--le",
            "100",
            "--name",
            "S288c",
            "--dash",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 9);
    assert!(stdout.contains("\nGCTAAAATATGAACG"), "no dash");
}
