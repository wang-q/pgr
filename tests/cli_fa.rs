#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use tempfile::TempDir;

#[test]
fn command_invalid() {
    let (_, stderr) = PgrCmd::new().args(&["fa", "foobar"]).run_fail();
    assert!(stderr.contains("recognized"));
}

#[test]
fn file_doesnt_provided() {
    let (_, stderr) = PgrCmd::new().args(&["fa", "size"]).run_fail();
    assert!(stderr.contains("not provided"));
}

#[test]
fn file_doesnt_exist() {
    let (_, stderr) = PgrCmd::new()
        .args(&["fa", "size", "tests/file/doesnt/exist"])
        .run_fail();
    assert!(stderr.contains("could not open"));
}

#[test]
fn command_fa_size() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("test.fa");

    fs::write(&input, ">seq1\nACGT\n>seq2\nACGTACGT\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "size", input.to_str().unwrap()])
        .run();

    assert!(stdout.contains("seq1\t4\n"));
    assert!(stdout.contains("seq2\t8\n"));
}

#[test]
fn command_fa_size_file() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "size", "tests/fasta/ufasta.fa"])
        .run();

    assert_eq!(stdout.lines().count(), 50);
    assert!(stdout.contains("read0\t359"), "read0");
    assert!(stdout.contains("read1\t106"), "read1");

    let mut sum = 0;
    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 2 {
            sum += fields[1].parse::<i32>().unwrap();
        }
    }
    assert_eq!(sum, 9317, "sum length");
}

#[test]
fn command_fa_size_gz() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "size",
            "tests/fasta/ufasta.fa",
            "tests/fasta/ufasta.fa.gz",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout.contains("read0\t359"), "read0");
    assert!(stdout.contains("read1\t106"), "read1");
}

#[test]
fn command_fa_size_no_ns() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("test_nons.fa");

    // seq1: 12 bases, 4 Ns (ACGT NNNN ACGT) -> 8 bases
    // seq2: 4 bases, 0 Ns -> 4 bases
    fs::write(&input, ">seq1\nACGTNNNNACGT\n>seq2\nACGT\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "size", input.to_str().unwrap(), "--no-ns"])
        .run();

    assert!(stdout.contains("seq1\t8\n"));
    assert!(stdout.contains("seq2\t4\n"));
}

#[test]
fn command_fa_some() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("test.fa");
    let list = temp.path().join("list.txt");
    let output = temp.path().join("out.fa");

    fs::write(&input, ">seq1\nACGT\n>seq2\nACGTACGT\n>seq3\nTTTT\n").unwrap();
    fs::write(&list, "seq1\nseq3\n").unwrap();

    PgrCmd::new()
        .args(&[
            "fa",
            "some",
            input.to_str().unwrap(),
            list.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .assert()
        .success();

    let content = fs::read_to_string(&output).unwrap();
    assert!(content.contains(">seq1"));
    assert!(content.contains(">seq3"));
    assert!(!content.contains(">seq2"));
}

#[test]
fn command_fa_some_invert() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("test.fa");
    let list = temp.path().join("list.txt");
    let output = temp.path().join("out.fa");

    fs::write(&input, ">seq1\nACGT\n>seq2\nACGTACGT\n>seq3\nTTTT\n").unwrap();
    fs::write(&list, "seq1\nseq3\n").unwrap();

    PgrCmd::new()
        .args(&[
            "fa",
            "some",
            input.to_str().unwrap(),
            list.to_str().unwrap(),
            "--invert",
            "-o",
            output.to_str().unwrap(),
        ])
        .assert()
        .success();

    let content = fs::read_to_string(&output).unwrap();
    assert!(!content.contains(">seq1"));
    assert!(!content.contains(">seq3"));
    assert!(content.contains(">seq2"));
}

#[test]
fn command_order() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "order",
            "tests/fasta/ufasta.fa",
            "tests/fasta/list.txt",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 4);
    assert!(stdout.contains("read12\n"), "read12");
    assert!(stdout.contains("read0\n"), "read0");
}

#[test]
fn command_one() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "one", "tests/fasta/ufasta.fa", "read12"])
        .run();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("read12\n"), "read12");
}

#[test]
fn command_masked() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "masked", "tests/fasta/ufasta.fa"])
        .run();

    assert!(stdout.contains("read46:3-4"), "read46");
}

#[test]
fn command_mask() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "mask",
            "tests/fasta/ufasta.fa",
            "tests/fasta/mask.json",
        ])
        .run();

    assert!(stdout.contains("read0\ntcgtttaacccaaatcaagg"), "read0");
    assert!(stdout.contains("read2\natagcaagct"), "read2");

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "mask",
            "--hard",
            "tests/fasta/ufasta.fa",
            "tests/fasta/mask.json",
        ])
        .run();

    assert!(stdout.contains("read0\nNNNNNNNNNNNNNNNNNNNN"), "read0");
    assert!(stdout.contains("read2\nNNNNNNNNNN"), "read2");
}

#[test]
fn command_rc() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "rc", "tests/fasta/ufasta.fa"])
        .run();

    assert!(stdout.contains("GgacTgcggCTagAA"), "read46");

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "rc", "tests/fasta/ufasta.fa", "tests/fasta/list.txt"])
        .run();

    assert!(stdout.contains(">RC_read12"), "read12");
    assert!(!stdout.contains(">RC_read46"), "read46");
    assert!(!stdout.contains("GgacTgcggCTagAA"), "read46");
}

#[test]
fn command_count() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "count", "tests/fasta/ufasta.fa"])
        .run();

    assert!(stdout.contains("read45\t0\t0"), "empty");
    assert!(stdout.contains("total\t9317\t2318"), "total");
}

#[test]
fn command_replace() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "replace",
            "tests/fasta/ufasta.fa",
            "tests/fasta/replace.tsv",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 95);
    assert!(stdout.contains(">359"), "read0");
    assert!(!stdout.contains(">read0"), "read0");

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "replace",
            "tests/fasta/ufasta.fa",
            "tests/fasta/replace.tsv",
            "--some",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains(">359"), "read0");
    assert!(!stdout.contains(">read0"), "read0");
}

#[test]
fn command_filter() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "filter",
            "tests/fasta/ufasta.fa",
            "-a",
            "10",
            "-z",
            "50",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 12);
    assert!(!stdout.contains(">read0"), "read0");
    assert!(stdout.contains(">read20"), "read20");

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "filter",
            "tests/fasta/ufasta.fa",
            "tests/fasta/ufasta.fa.gz",
            "--uniq",
            "-a",
            "1",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 90);
}

#[test]
fn command_filter_fmt() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "filter", "tests/fasta/filter.fa", "--iupac"])
        .run();

    assert!(!stdout.contains(">iupac\nAMRG"), "iupac");
    assert!(stdout.contains(">iupac\nANNG"), "iupac");
    assert!(stdout.contains(">dash\nA-NG"), "dash not changed");

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "filter", "tests/fasta/filter.fa", "--dash"])
        .run();

    assert!(!stdout.contains(">dash\nA-RG"), "dash");
    assert!(stdout.contains(">dash\nARG"), "dash");

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "filter", "tests/fasta/filter.fa", "--upper"])
        .run();

    assert!(!stdout.contains(">upper\nAtcG"), "upper");
    assert!(stdout.contains(">upper\nATCG"), "upper");

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "filter", "tests/fasta/filter.fa", "--simplify"])
        .run();

    assert!(!stdout.contains(">read.1 simplify\nAGGG"), "simplify");
    assert!(stdout.contains(">read\nAGGG"), "simplify");
}

#[test]
fn command_dedup() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "dedup", "tests/fasta/dedup.fa"])
        .run();

    assert_eq!(stdout.lines().count(), 8);
    assert!(!stdout.contains(">read0 some text"));

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "dedup", "tests/fasta/dedup.fa", "--desc"])
        .run();

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains(">read0 some text"));

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "dedup", "tests/fasta/dedup.fa", "--seq"])
        .run();

    assert_eq!(stdout.lines().count(), 6);
    assert!(!stdout.contains(">read1"));

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "dedup", "tests/fasta/dedup.fa", "--seq", "--case"])
        .run();

    assert_eq!(stdout.lines().count(), 4);
    assert!(!stdout.contains(">read2"));

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "dedup", "tests/fasta/dedup.fa", "--seq", "--both"])
        .run();

    assert_eq!(stdout.lines().count(), 2);
    assert!(!stdout.contains(">read3"));

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "dedup",
            "tests/fasta/dedup.fa",
            "--seq",
            "--both",
            "--file",
            "stdout",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains(">read0"));
    assert!(stdout.contains("read0\tread3"));
}

#[test]
fn command_split_name() {
    let tempdir = TempDir::new().unwrap();
    let tempdir_str = tempdir.path().to_str().unwrap();

    PgrCmd::new()
        .args(&[
            "fa",
            "split",
            "name",
            "tests/fasta/ufasta.fa",
            "-o",
            tempdir_str,
        ])
        .assert()
        .success()
        .stdout(predicates::str::is_empty());

    assert!(&tempdir.path().join("read0.fa").is_file());
    assert!(!&tempdir.path().join("000.fa").exists());

    tempdir.close().unwrap();
}

#[test]
fn command_split_about() {
    let tempdir = TempDir::new().unwrap();
    let tempdir_str = tempdir.path().to_str().unwrap();

    PgrCmd::new()
        .args(&[
            "fa",
            "split",
            "about",
            "tests/fasta/ufasta.fa",
            "-c",
            "2000",
            "-o",
            tempdir_str,
        ])
        .assert()
        .success()
        .stdout(predicates::str::is_empty());

    assert!(!&tempdir.path().join("read0.fa").is_file());
    assert!(&tempdir.path().join("000.fa").exists());
    assert!(&tempdir.path().join("004.fa").exists());
    assert!(!&tempdir.path().join("005.fa").exists());

    tempdir.close().unwrap();
}

#[test]
fn command_fa_n50() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("test.fa");

    fs::write(
        &input,
        ">seq1\nN\n>seq2\nN\n>seq3\nNN\n>seq4\nNN\n>seq5\nNNNN\n"
            .replace("N", "N".repeat(100).as_str()),
    )
    .unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "n50", input.to_str().unwrap()])
        .run();

    assert!(stdout.contains("N50\t200\n"));
}

#[test]
fn command_fa_n50_stats() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("test.fa");

    fs::write(
        &input,
        ">seq1\nN\n>seq2\nN\n>seq3\nNN\n>seq4\nNN\n>seq5\nNNNN\n"
            .replace("N", "N".repeat(100).as_str()),
    )
    .unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "n50", input.to_str().unwrap(), "-S", "-A", "-C", "-H"])
        .run();

    assert!(stdout.contains("200\n"));
    assert!(stdout.contains("1000\n"));
    assert!(stdout.contains("200.00\n"));
    assert!(stdout.contains("5\n"));
}

#[test]
fn command_fa_n50_comprehensive() {
    // display header
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "n50", "tests/fasta/ufasta.fa"])
        .run();

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("N50\t314"), "line 1");

    // doesn't display header
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "n50", "tests/fasta/ufasta.fa", "-H"])
        .run();

    assert_eq!(stdout.lines().count(), 1);
    assert!(!stdout.contains("N50\t314"), "line 1");
    assert!(stdout.contains("314"), "line 1");

    // set genome size (NG50)
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "n50", "tests/fasta/ufasta.fa", "-H", "-g", "10000"])
        .run();

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("297"), "line 1");

    // sum and average of size
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "n50", "tests/fasta/ufasta.fa", "-H", "-S", "-A"])
        .run();

    assert_eq!(stdout.lines().count(), 3);
    assert!(stdout.contains("314\n9317\n186.34"), "line 1,2,3");

    // N10, N90, E-size
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "n50",
            "tests/fasta/ufasta.fa",
            "-H",
            "-E",
            "-N",
            "10",
            "-N",
            "90",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 3);
    assert!(stdout.contains("516\n112\n314.70\n"), "line 1,2,3");

    // transposed
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "n50",
            "tests/fasta/ufasta.fa",
            "-E",
            "-N",
            "10",
            "-N",
            "90",
            "-t",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("N10\tN90\tE\n"), "line 1");
    assert!(stdout.contains("516\t112\t314.70\n"), "line 2");
}

#[test]
fn command_six_frame() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "six-frame", "tests/fasta/trans.fa"])
        .run();

    assert_eq!(stdout.lines().count(), 16);
    assert!(stdout.contains(">seq1(+):1-15|frame=0"));
    assert!(stdout.contains("MGMG*"));
    assert!(stdout.contains(">seq1(-):3-26|frame=2"));
    assert!(stdout.contains("TIYLYPIP"));

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "six-frame", "tests/fasta/trans.fa", "--len", "3"])
        .run();

    assert_eq!(stdout.lines().count(), 12);

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "six-frame",
            "tests/fasta/trans.fa",
            "--len",
            "3",
            "--end",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 4);

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "six-frame",
            "tests/fasta/trans.fa",
            "--len",
            "3",
            "--start",
            "--end",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 2);
}
