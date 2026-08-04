#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

/// Return the absolute path to a fixture in `tests/fasta/input`.
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta/input")
        .join(name)
}

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
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "size", fixture("basic.fa").to_str().unwrap()])
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
    // seq1: 12 bases, 4 Ns (ACGT NNNN ACGT) -> 8 bases
    // seq2: 4 bases, 0 Ns -> 4 bases
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "size",
            fixture("nons.fa").to_str().unwrap(),
            "--no-ns",
        ])
        .run();

    assert!(stdout.contains("seq1\t8\n"));
    assert!(stdout.contains("seq2\t4\n"));
}

#[test]
fn command_fa_some() {
    let temp = TempDir::new().unwrap();
    let list = temp.path().join("list.txt");
    let output = temp.path().join("out.fa");

    fs::write(&list, "seq1\nseq3\n").unwrap();

    PgrCmd::new()
        .args(&[
            "fa",
            "some",
            fixture("some.fa").to_str().unwrap(),
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
    let list = temp.path().join("list.txt");
    let output = temp.path().join("out.fa");

    fs::write(&list, "seq1\nseq3\n").unwrap();

    PgrCmd::new()
        .args(&[
            "fa",
            "some",
            fixture("some.fa").to_str().unwrap(),
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
fn command_fa_some_ignores_hash_comments() {
    let temp = TempDir::new().unwrap();
    let list = temp.path().join("list.txt");
    let output = temp.path().join("out.fa");

    // A `#` comment line must be ignored by the name list reader.
    fs::write(&list, "# comment\nseq1\nseq3\n").unwrap();

    PgrCmd::new()
        .args(&[
            "fa",
            "some",
            fixture("some.fa").to_str().unwrap(),
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
            "--runlist",
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
            "--runlist",
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
            "--replace-tsv",
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
            "--replace-tsv",
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
            "--min-len",
            "10",
            "--max-len",
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
            "--min-len",
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
    assert!(stdout.contains(">read simplify\nAGGG"), "simplify");
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
            "--dups-file",
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
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "n50", fixture("n50.fa").to_str().unwrap()])
        .run();

    assert!(stdout.contains("N50\t200\n"));
}

#[test]
fn command_fa_n50_stats() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "n50",
            fixture("n50.fa").to_str().unwrap(),
            "-S",
            "-A",
            "-C",
            "-H",
        ])
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
            "--transpose",
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
        .args(&["fa", "six-frame", "tests/fasta/trans.fa", "--min-len", "3"])
        .run();

    assert_eq!(stdout.lines().count(), 12);

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "six-frame",
            "tests/fasta/trans.fa",
            "--min-len",
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
            "--min-len",
            "3",
            "--start-met",
            "--end",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 2);
}

#[test]
fn command_six_frame_short_sequence_no_panic() {
    // A 1-base sequence must not panic (regression: reverse frames used to
    // underflow `dna_len - frame` and the forward frame-2 slice used to panic).
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("short.fa");
    fs::write(&input, ">s\nA\n").unwrap();

    PgrCmd::new()
        .args(&["fa", "six-frame", input.to_str().unwrap()])
        .assert()
        .success();
}

#[test]
fn command_fa_output_same_as_input_rejected() {
    // `-o` pointing at an input file must be rejected before the writer
    // truncates the input (regression: it used to silently empty the file).
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("in.fa");
    let original = ">seq\nACGTACGT\n";
    fs::write(&input, original).unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fa",
            "filter",
            input.to_str().unwrap(),
            "--min-len",
            "1",
            "-o",
            input.to_str().unwrap(),
        ])
        .run_fail();

    assert!(stderr.contains("is also an input file"));
    assert_eq!(fs::read_to_string(&input).unwrap(), original);
}

#[test]
fn command_six_frame_output_same_as_input_rejected() {
    // Same data-safety guarantee as the other fa subcommands: `-o` must not
    // truncate the input before it is read.
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("in.fa");
    let original = ">seq\nATGACGTAG\n";
    fs::write(&input, original).unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fa",
            "six-frame",
            input.to_str().unwrap(),
            "-o",
            input.to_str().unwrap(),
        ])
        .run_fail();

    assert!(stderr.contains("is also an input file"));
    assert_eq!(fs::read_to_string(&input).unwrap(), original);
}

#[test]
fn command_fa_split_output_not_overwrite_input() {
    // `split` writes into -o as a directory; an output file whose resolved
    // path collides with an input must be rejected before it is truncated.
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("chr.fa");
    let original = ">chr\nACGTACGTACGT\n";
    fs::write(&input, original).unwrap();

    // `name` mode: sequence name `chr` -> output `outdir/chr.fa`, which is the
    // input itself.
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fa",
            "split",
            "name",
            input.to_str().unwrap(),
            "-o",
            temp.path().to_str().unwrap(),
        ])
        .run_fail();

    assert!(stderr.contains("would overwrite input file"));
    assert_eq!(fs::read_to_string(&input).unwrap(), original);
}

#[test]
fn command_fa_one_not_found() {
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fa",
            "one",
            fixture("basic.fa").to_str().unwrap(),
            "nonexistent",
        ])
        .run_fail();

    assert!(stderr.contains("not found"));
}

#[test]
fn command_fa_one_success() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "one", fixture("basic.fa").to_str().unwrap(), "seq2"])
        .run();

    assert!(stdout.contains(">seq2"));
    assert!(stdout.contains("ACGTACGT"));
}

#[test]
fn command_fa_filter_uniq_ucsc_gold_names() {
    // faFilter test: duplicate basic.fa, then -uniq keeps the first of each
    // duplicated id. UCSC's gold output differs only by dropping header
    // descriptions, so compare record names.
    let dir = TempDir::new().unwrap();
    let dup = dir.path().join("dup.fa");
    let basic = fs::read_to_string(fixture("ucsc_basic.fa").to_str().unwrap()).unwrap();
    fs::write(&dup, format!("{}{}", basic, basic)).unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["fa", "filter", dup.to_str().unwrap(), "--uniq"])
        .run();

    let names: Vec<&str> = stdout
        .lines()
        .filter(|l| l.starts_with('>'))
        .map(|l| l.trim_start_matches('>').split_whitespace().next().unwrap())
        .collect();
    assert_eq!(names.len(), 10, "10 unique records after dedup");
    assert_eq!(names[0], "size9");
    assert_eq!(names[3], "foo9baz");
    assert_eq!(names[9], "1acc");
}

#[test]
fn command_fa_filter_max_len_ucsc_gold_names() {
    // faFilter -maxSize=20 keeps 7 of the 10 basic.fa records.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fa",
            "filter",
            fixture("ucsc_basic.fa").to_str().unwrap(),
            "--max-len",
            "20",
        ])
        .run();

    let names: Vec<&str> = stdout
        .lines()
        .filter(|l| l.starts_with('>'))
        .map(|l| l.trim_start_matches('>').split_whitespace().next().unwrap())
        .collect();
    assert_eq!(names.len(), 7);
    assert_eq!(names[6], "size20");
}
