//! CLI tests for `pgr rg` (migrated from `pgr runlist cover/coverage`).

#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::path::PathBuf;
use tempfile::TempDir;

fn fixture_dir() -> (TempDir, PathBuf) {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    std::fs::write(&rg, "chr1:1-10\nchr1:5-15\nchr2(+):100-200\nbad line\n").unwrap();
    let rg2 = dir.path().join("b.rg");
    std::fs::write(&rg2, "chr1:20-25\nchr2:150-160\n").unwrap();
    (dir, rg)
}

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/runlist")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn cmd(args: &[&str]) -> PgrCmd {
    let mut full = vec!["rg"];
    full.extend_from_slice(args);
    PgrCmd::new().args(&full)
}

#[test]
fn command_rg_cover() {
    let (dir, rg) = fixture_dir();
    let rg2 = dir.path().join("b.rg");
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", rg.to_str().unwrap(), rg2.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"chr1\": \"1-15,20-25\",\n  \"chr2\": \"100-200\"\n}\n"
    );
}

#[test]
fn command_rg_cover_stdin() {
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", "stdin"])
        .stdin("chr1:1-5\nchr1:3-8\n")
        .run();
    assert_eq!(stdout, "{\n  \"chr1\": \"1-8\"\n}\n");
}

#[test]
fn command_rg_cover_skips_reversed_ranges() {
    // `chr1:10-5` (start > end) must be skipped, not panic.
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", "stdin"])
        .stdin("chr1:10-5\nchr1:1-10\n")
        .run();
    assert_eq!(stdout, "{\n  \"chr1\": \"1-10\"\n}\n");
    // Coordinates above the representable maximum must be skipped too.
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", "stdin"])
        .stdin("chr1:2147483647-2147483647\n")
        .run();
    assert_eq!(stdout, "{}\n");
}

// Comment lines starting with `#` must be skipped by every rg subcommand
// (cover/coverage already did; count/span/sort/prop/runlist/merge used to
// treat `# chr1:1-10` as data).
#[test]
fn command_rg_comments_skipped() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("cmt.rg");
    std::fs::write(&rg, "# chr1:1-10\nchr1:5-15\n#chr2:1-5\nbad line\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", rg.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "{\n  \"chr1\": \"5-15\"\n}\n");

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "count", rg.to_str().unwrap(), rg.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "chr1:5-15\t1\n");

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "span", rg.to_str().unwrap(), "-n", "5"])
        .run();
    assert_eq!(stdout, "chr1:10\n");

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "sort", rg.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "chr1:5-15\nbad line\n");

    let json = dir.path().join("in.json");
    std::fs::write(&json, r#"{"chr1":"5-15"}"#).unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "rg",
            "runlist",
            json.to_str().unwrap(),
            rg.to_str().unwrap(),
        ])
        .run();
    assert_eq!(stdout, "chr1:5-15\n");

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "prop", json.to_str().unwrap(), rg.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "chr1:5-15\t1.0000\n");

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "merge", rg.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "");
}

// An oversized `end` coordinate (`chr1:5-99999999999`) used to be parsed as
// the point range `chr1:5` (overflow defaulted to start); it must be
// treated as an invalid line instead of silently corrupting coordinates.
#[test]
fn command_rg_overflow_end_skipped() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("ov.rg");
    std::fs::write(
        &rg,
        "chr1:5-99999999999\nchr1:1-10\nS288c.I(-):5-99999999999\n",
    )
    .unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", rg.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "{\n  \"chr1\": \"1-10\"\n}\n");

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "span", rg.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "chr1:1-10\n");

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "count", rg.to_str().unwrap(), rg.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "chr1:1-10\t1\n");

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "sort", rg.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "chr1:1-10\nchr1:5-99999999999\nS288c.I(-):5-99999999999\n"
    );
}

// Unicode word characters in contig names used to fail the ASCII-only
// scanner and were silently dropped by every rg subcommand; they must parse
// like the reference regex's Unicode-mode `\w` did.
#[test]
fn command_rg_unicode_contig_names_parsed() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("uni.rg");
    std::fs::write(&rg, "chr\u{3b1}:1-10\n\u{4e2d}(+):20-30\nbad line\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", rg.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"chr\u{3b1}\": \"1-10\",\n  \"\u{4e2d}\": \"20-30\"\n}\n"
    );

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "span", rg.to_str().unwrap(), "-n", "5"])
        .run();
    // `chrα:1-10` trimmed by 5 collapses to an invalid (0, 0) range whose
    // Display keeps the bare chromosome (vendored Range semantics), while
    // `中(+):20-30` becomes the point range `中(+):25`.
    assert_eq!(stdout, "chr\u{3b1}\n\u{4e2d}(+):25\n");

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "count", rg.to_str().unwrap(), rg.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "chr\u{3b1}:1-10\t1\n\u{4e2d}(+):20-30\t1\n");

    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "sort", rg.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "chr\u{3b1}:1-10\n\u{4e2d}(+):20-30\nbad line\n");
}

#[test]
fn command_rg_coverage() {
    let (dir, rg) = fixture_dir();
    let rg2 = dir.path().join("b.rg");
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "rg",
            "coverage",
            rg.to_str().unwrap(),
            rg2.to_str().unwrap(),
            "-m",
            "2",
        ])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"chr1\": \"5-10\",\n  \"chr2\": \"150-160\"\n}\n"
    );
}

#[test]
fn command_rg_coverage_detailed() {
    let (_dir, rg) = fixture_dir();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "coverage", rg.to_str().unwrap(), "-m", "1", "-d"])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"1\": {\n    \"chr1\": \"1-4,11-15\",\n    \"chr2\": \"100-200\"\n  },\n  \"2\": {\n    \"chr1\": \"5-10\"\n  }\n}\n"
    );
}

#[test]
fn command_rg_count() {
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("target.rg");
    let intervals = dir.path().join("iv.rg");
    std::fs::write(&target, "chr1:1-10\nchr1:5-15\nchr2:100-200\nbad line\n").unwrap();
    std::fs::write(&intervals, "chr1:1-10\nchr1:5-15\nchr2:150-160\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "rg",
            "count",
            target.to_str().unwrap(),
            intervals.to_str().unwrap(),
        ])
        .run();
    assert_eq!(stdout, "chr1:1-10\t2\nchr1:5-15\t2\nchr2:100-200\t1\n");
}

#[test]
fn command_rg_count_stdin_intervals() {
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("target.rg");
    std::fs::write(&target, "chr1:1-5\nchr1:10-20\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "count", target.to_str().unwrap(), "stdin"])
        .stdin("chr1:1-10\nchr1:15-16\n")
        .run();
    assert_eq!(stdout, "chr1:1-5\t1\nchr1:10-20\t2\n");
}

// Migrated from intspan `tests/cli_rgr.rs` `command_count` (.rg part; the
// TSV `-f` part tests functionality dropped by the rg-family input contract).
#[test]
fn command_rg_count_fixture() {
    let (stdout, _) = cmd(&["count", &fixture("S288c.rg"), &fixture("S288c.rg")]).run();
    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains("I:1-100\t2"), "got: {stdout}");
    assert!(stdout.contains("21294-22075\t1"), "got: {stdout}");
}

#[test]
fn command_rg_prop() {
    let dir = TempDir::new().unwrap();
    let json = dir.path().join("in.json");
    let rg = dir.path().join("a.rg");
    std::fs::write(&json, r#"{"chr1":"1-10,20-30"}"#).unwrap();
    std::fs::write(&rg, "chr1:5-25\nchr1:50-60\nchr2:1-10\nbad line\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "prop", json.to_str().unwrap(), rg.to_str().unwrap()])
        .run();
    // chr1:5-25 intersects 5-10 (6 bp) + 20-25 (6 bp) = 12 of 21 bp.
    assert_eq!(
        stdout,
        "chr1:5-25\t0.5714\nchr1:50-60\t0.0000\nchr2:1-10\t0.0000\n"
    );
}

#[test]
fn command_rg_prop_full() {
    let dir = TempDir::new().unwrap();
    let json = dir.path().join("in.json");
    let rg = dir.path().join("a.rg");
    std::fs::write(&json, r#"{"chr1":"1-10"}"#).unwrap();
    std::fs::write(&rg, "chr1:1-20\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "rg",
            "prop",
            json.to_str().unwrap(),
            rg.to_str().unwrap(),
            "--full",
        ])
        .run();
    assert_eq!(stdout, "chr1:1-20\t0.5000\t20\t10\n");
}

// Migrated from intspan `tests/cli_rgr.rs` `command_prop` (.rg part).
#[test]
fn command_rg_prop_fixture() {
    let (stdout, _) = cmd(&["prop", &fixture("intergenic.json"), &fixture("S288c.rg")]).run();
    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains("I:1-100\t0.0000"), "got: {stdout}");
    assert!(stdout.contains("II:21294-22075\t1.0000"), "got: {stdout}");
}

#[test]
fn command_rg_sort() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    std::fs::write(
        &rg,
        "chr2:100-200\nchr1:50-60\nchr1:100-110\nbad line\nchr1:20-30\n",
    )
    .unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "sort", rg.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "chr1:20-30\nchr1:50-60\nchr1:100-110\nchr2:100-200\nbad line\n"
    );
}

#[test]
fn command_rg_sort_strand_order() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    std::fs::write(&rg, "chr1(+):50-60\nchr1(-):50-60\nchr1(+):20-30\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "sort", rg.to_str().unwrap()])
        .run();
    // Key is (chr, start, strand): start 20 first, then start 50 with "+"
    // before "-" (ASCII '+' < '-').
    assert_eq!(stdout, "chr1(+):20-30\nchr1(+):50-60\nchr1(-):50-60\n");
}

// Migrated from intspan `tests/cli_rgr.rs` `command_sort` (.rg part).
#[test]
fn command_rg_sort_fixture() {
    let (stdout, _) = cmd(&["sort", &fixture("S288c.rg")]).run();
    assert_eq!(stdout.lines().count(), 6);
    assert_eq!(
        stdout.lines().next().unwrap().split('\t').count(),
        1,
        "field count"
    );
    assert!(
        stdout.contains("S288c.I(-):190-200\nS288c"),
        "got: {stdout}"
    );
}

// Migrated from intspan `tests/cli_rgr.rs` `command_runlist` (.rg parts).
#[test]
fn command_rg_runlist() {
    let rl = fixture("intergenic.json");
    let rg = fixture("S288c.rg");
    let (stdout, _) = cmd(&["runlist", &rl, &rg]).run();
    assert_eq!(stdout.lines().count(), 2);
    assert!(!stdout.contains("S288c"), "got: {stdout}");
    assert!(stdout.contains("21294-22075"), "got: {stdout}");

    let (stdout, _) = cmd(&["runlist", &rl, &rg, "--op", "non-overlap"]).run();
    assert_eq!(stdout.lines().count(), 4);
    assert!(stdout.contains("S288c"), "got: {stdout}");
    assert!(!stdout.contains("21294-22075"), "got: {stdout}");

    let (stdout, _) = cmd(&["runlist", &rl, &rg, "--op", "superset"]).run();
    assert_eq!(stdout.lines().count(), 2);
    assert!(!stdout.contains("S288c"), "got: {stdout}");
    assert!(stdout.contains("21294-22075"), "got: {stdout}");
}

// Migrated from intspan `tests/cli_rgr.rs` `command_runlist_invalid`.
#[test]
fn command_rg_runlist_invalid() {
    let (_, stderr) = cmd(&[
        "runlist",
        &fixture("intergenic.json"),
        &fixture("S288c.rg"),
        "--op",
        "invalid",
    ])
    .run_fail();
    assert!(stderr.contains("invalid value"), "got: {stderr}");
}

// `superset` keeps only ranges fully contained in the runlist. A range that
// merely overlaps or *contains* a runlist span is rejected — this locks in
// the "range inside runlist" reading (rgr `set.superset(range)`) and rules
// out the alternative "range contains a span" interpretation.
#[test]
fn command_rg_runlist_superset_partial_overlap() {
    let dir = TempDir::new().unwrap();
    let rl = dir.path().join("rl.json");
    std::fs::write(&rl, "{\"chr1\": \"10-20\"}").unwrap();
    let rg = dir.path().join("a.rg");
    std::fs::write(&rg, "chr1:10-20\nchr1:5-15\nchr1:5-25\nchr1:25-30\n").unwrap();
    let rl = rl.to_str().unwrap();
    let rg = rg.to_str().unwrap();

    // overlap: exact, partial and containing ranges all intersect the span.
    let (stdout, _) = cmd(&["runlist", rl, rg]).run();
    assert_eq!(stdout, "chr1:10-20\nchr1:5-15\nchr1:5-25\n");

    // superset: only the exact-containment line is kept.
    let (stdout, _) = cmd(&["runlist", rl, rg, "--op", "superset"]).run();
    assert_eq!(stdout, "chr1:10-20\n");

    // non-overlap: only the disjoint line is kept.
    let (stdout, _) = cmd(&["runlist", rl, rg, "--op", "non-overlap"]).run();
    assert_eq!(stdout, "chr1:25-30\n");
}

// Migrated from intspan `tests/cli_rgr.rs` `command_span` (.rg parts).
#[test]
fn command_rg_span() {
    let rg = fixture("S288c.rg");
    let (stdout, _) = cmd(&["span", &rg, "--op", "trim", "-n", "10"]).run();
    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains("I:11-90"), "got: {stdout}");
    assert!(stdout.contains("II:21304-22065"), "got: {stdout}");

    let (stdout, _) = cmd(&["span", &rg, "--op", "shift", "-m", "3p", "-n", "10"]).run();
    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains("I:11-110"), "got: {stdout}");
    assert!(stdout.contains("S288c.I(-):180-190"), "got: {stdout}");

    let (stdout, _) = cmd(&["span", &rg, "--op", "flank", "-m", "3p", "-n=-1", "-a"]).run();
    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains("I:1-100\tI:100"), "got: {stdout}");
    assert!(
        stdout.contains("S288c.I(-):190-200|Species=Yeast\tS288c.I(-):190"),
        "got: {stdout}"
    );

    let (stdout, _) = cmd(&["span", &rg, "--op", "excise", "-n", "20"]).run();
    assert_eq!(stdout.lines().count(), 6);
    assert_eq!(
        stdout.lines().filter(|e| e.is_empty()).count(),
        2,
        "empty lines"
    );
}

// Extreme `-n` values used to overflow the Range op arithmetic (debug panic,
// release wrap); they must not crash and must stay deterministic.
#[test]
fn command_rg_span_extreme_no_panic() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("mx.rg");
    std::fs::write(
        &rg,
        "chr1:2147483645-2147483645\nchr2:2000000000-2100000000\nchr3(-):1-10\n",
    )
    .unwrap();
    for (op, mode) in [
        ("trim", "both"),
        ("trim", "5p"),
        ("trim", "3p"),
        ("pad", "both"),
        ("shift", "5p"),
        ("shift", "3p"),
        ("flank", "5p"),
        ("flank", "3p"),
    ] {
        let (stdout, _) = PgrCmd::new()
            .args(&[
                "rg",
                "span",
                rg.to_str().unwrap(),
                "--op",
                op,
                "-m",
                mode,
                "-n",
                "2147483647",
            ])
            .run();
        assert_eq!(stdout.lines().count(), 3, "{op}/{mode}: {stdout}");
    }
    // `pad -n i32::MIN` used to panic on `-number` negation.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "rg",
            "span",
            rg.to_str().unwrap(),
            "--op",
            "pad",
            "-n=-2147483648",
        ])
        .run();
    assert_eq!(stdout.lines().count(), 3);
}

// The shift/flank + `--mode both` combination used to be validated only
// while processing lines, so an empty (or comment-only) input silently
// succeeded while a non-empty one failed. The validation must happen before
// any input is read.
#[test]
fn command_rg_span_invalid_mode_checked_before_input() {
    let dir = TempDir::new().unwrap();
    let empty = dir.path().join("empty.rg");
    std::fs::write(&empty, "").unwrap();
    let comments = dir.path().join("comments.rg");
    std::fs::write(&comments, "# chr1:1-10\nbad line\n").unwrap();
    for op in ["shift", "flank"] {
        for infile in [&empty, &comments] {
            let (_, stderr) = PgrCmd::new()
                .args(&[
                    "rg",
                    "span",
                    infile.to_str().unwrap(),
                    "--op",
                    op,
                    "-m",
                    "both",
                ])
                .run_fail();
            assert!(
                stderr.contains(&format!("invalid for {op}")),
                "{op}: got: {stderr}"
            );
        }
    }
}

// Writing the output over an input file used to truncate the input before it
// was read (span/prop/runlist streamed line by line; count built its index
// first but then truncated the target). The command must refuse instead.
#[test]
fn command_rg_output_same_as_input_rejected() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    let json = dir.path().join("in.json");
    std::fs::write(&rg, "chr1:1-10\nchr1:20-30\n").unwrap();
    std::fs::write(&json, r#"{"chr1":"1-5"}"#).unwrap();

    let cases: Vec<Vec<&str>> = vec![
        vec![
            "rg",
            "span",
            rg.to_str().unwrap(),
            "-o",
            rg.to_str().unwrap(),
        ],
        vec![
            "rg",
            "prop",
            json.to_str().unwrap(),
            rg.to_str().unwrap(),
            "-o",
            rg.to_str().unwrap(),
        ],
        vec![
            "rg",
            "runlist",
            json.to_str().unwrap(),
            rg.to_str().unwrap(),
            "-o",
            rg.to_str().unwrap(),
        ],
        vec![
            "rg",
            "count",
            rg.to_str().unwrap(),
            rg.to_str().unwrap(),
            "-o",
            rg.to_str().unwrap(),
        ],
    ];
    for args in cases {
        let (_, stderr) = PgrCmd::new().args(&args).run_fail();
        assert!(
            stderr.contains("also an input file"),
            "{args:?}: got: {stderr}"
        );
    }
    // The input must be untouched.
    assert_eq!(
        std::fs::read_to_string(&rg).unwrap(),
        "chr1:1-10\nchr1:20-30\n"
    );
}

// A file literally named `stdout` must not be mistaken for the screen
// sentinel: writing to the screen cannot overwrite it, so the command must
// read it normally instead of rejecting the run (`same_path` compares the
// literal `stdout` output value against the input path).
#[test]
fn command_rg_stdout_named_input_allowed() {
    let dir = TempDir::new().unwrap();
    std::fs::write(dir.path().join("in.json"), r#"{"chr1":"1-10"}"#).unwrap();
    std::fs::write(dir.path().join("stdout"), "chr1:1-10\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .current_dir(dir.path())
        .args(&["rg", "runlist", "in.json", "stdout"])
        .run();
    assert_eq!(stdout, "chr1:1-10\n");

    let (stdout, _) = PgrCmd::new()
        .current_dir(dir.path())
        .args(&["rg", "runlist", "in.json", "stdout", "-o", "stdout"])
        .run();
    assert_eq!(stdout, "chr1:1-10\n");
}

// The `stdin` input sentinel denotes the standard-input stream, never a
// file; an output path literally named `stdin` (a fresh file) cannot
// overwrite it, so the same-path guard must not reject the run (mirror of
// the `stdout` exemption).
#[test]
fn command_rg_stdin_sentinel_output_allowed() {
    let dir = TempDir::new().unwrap();
    std::fs::write(dir.path().join("in.json"), r#"{"chr1":"1-10"}"#).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .current_dir(dir.path())
        .stdin("chr1:1-10\n")
        .args(&["rg", "runlist", "in.json", "stdin", "-o", "stdin"])
        .run();
    assert_eq!(stdout, "", "stdout: {stdout}");
    assert_eq!(stderr, "", "stderr: {stderr}");
    assert_eq!(
        std::fs::read_to_string(dir.path().join("stdin")).unwrap(),
        "chr1:1-10\n"
    );
}

// Streaming commands used to create the output before opening the later
// inputs; a missing input then truncated an existing output file and left
// partial results. All inputs must be validated before the output is touched.
#[test]
fn command_rg_output_preserved_on_missing_input() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    let json = dir.path().join("in.json");
    std::fs::write(&rg, "chr1:1-10\n").unwrap();
    std::fs::write(&json, r#"{"chr1":"1-5"}"#).unwrap();
    let missing = dir.path().join("missing.rg");

    let cases: Vec<Vec<&str>> = vec![
        vec![
            "rg",
            "runlist",
            json.to_str().unwrap(),
            rg.to_str().unwrap(),
            missing.to_str().unwrap(),
        ],
        vec![
            "rg",
            "prop",
            json.to_str().unwrap(),
            rg.to_str().unwrap(),
            missing.to_str().unwrap(),
        ],
        vec![
            "rg",
            "span",
            rg.to_str().unwrap(),
            missing.to_str().unwrap(),
        ],
        vec![
            "rg",
            "count",
            missing.to_str().unwrap(),
            rg.to_str().unwrap(),
        ],
    ];
    for args in cases {
        let out = dir.path().join("out.txt");
        std::fs::write(&out, "PRECIOUS\n").unwrap();
        let mut full = args.clone();
        full.extend_from_slice(&["-o", out.to_str().unwrap()]);
        let (_, stderr) = PgrCmd::new().args(&full).run_fail();
        assert!(stderr.contains("could not open"), "{full:?}: got: {stderr}");
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            "PRECIOUS\n",
            "{full:?}: output was truncated"
        );
    }
}

// `File::open` succeeds on directories on Unix, so the pre-output open probe
// used to pass for a directory input and the "Is a directory" error only
// surfaced on the first read — after the output file had been truncated.
// The reader must reject directories up front, before any output is touched.
#[test]
fn command_rg_output_preserved_on_directory_input() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    let json = dir.path().join("in.json");
    std::fs::write(&rg, "chr1:1-10\n").unwrap();
    std::fs::write(&json, r#"{"chr1":"1-5"}"#).unwrap();
    let adir = dir.path().join("a_dir");
    std::fs::create_dir(&adir).unwrap();

    let cases: Vec<Vec<&str>> = vec![
        vec![
            "rg",
            "runlist",
            json.to_str().unwrap(),
            rg.to_str().unwrap(),
            adir.to_str().unwrap(),
        ],
        vec![
            "rg",
            "prop",
            json.to_str().unwrap(),
            rg.to_str().unwrap(),
            adir.to_str().unwrap(),
        ],
        vec!["rg", "span", rg.to_str().unwrap(), adir.to_str().unwrap()],
        vec!["rg", "count", adir.to_str().unwrap(), rg.to_str().unwrap()],
    ];
    for args in cases {
        let out = dir.path().join("out.txt");
        std::fs::write(&out, "PRECIOUS\n").unwrap();
        let mut full = args.clone();
        full.extend_from_slice(&["-o", out.to_str().unwrap()]);
        let (_, stderr) = PgrCmd::new().args(&full).run_fail();
        assert!(stderr.contains("is a directory"), "{full:?}: got: {stderr}");
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            "PRECIOUS\n",
            "{full:?}: output was truncated"
        );
    }
}

// `-o` pointing at the input through a symlink (or a hard link) used to
// bypass the lexical `same_path` check and truncate the input while it was
// still being read. Both aliases must be refused.
#[test]
fn command_rg_output_alias_of_input_rejected() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    std::fs::write(&rg, "chr1:1-10\nchr1:20-30\n").unwrap();
    let symlink = dir.path().join("sym.rg");
    std::fs::hard_link(&rg, &symlink).unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "rg",
            "span",
            rg.to_str().unwrap(),
            "-o",
            symlink.to_str().unwrap(),
        ])
        .run_fail();
    assert!(stderr.contains("also an input file"), "got: {stderr}");

    #[cfg(unix)]
    {
        use std::os::unix::fs::symlink;
        let link = dir.path().join("ln.rg");
        symlink(&rg, &link).unwrap();
        let (_, stderr) = PgrCmd::new()
            .args(&[
                "rg",
                "span",
                rg.to_str().unwrap(),
                "-o",
                link.to_str().unwrap(),
            ])
            .run_fail();
        assert!(stderr.contains("also an input file"), "got: {stderr}");
    }
    // Neither alias must have truncated the input.
    assert_eq!(
        std::fs::read_to_string(&rg).unwrap(),
        "chr1:1-10\nchr1:20-30\n"
    );
}

// New test: rgr's `command_merge` used a multi-part TSV fixture
// (II.links.tsv), which decision A excludes; this covers the .rg adaptation.
#[test]
fn command_rg_merge() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    std::fs::write(&rg, "chr1:100-200\nchr1:105-205\nchr1:1000-2000\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "merge", rg.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "chr1:100-200\tchr1(+):100-205\nchr1:105-205\tchr1(+):100-205\n"
    );
    // Looser threshold still leaves the disjoint range unmerged.
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "merge", rg.to_str().unwrap(), "-c", "0.5"])
        .run();
    assert_eq!(stdout.lines().count(), 2);
}

// Identical lines are deduplicated per chromosome before clustering, so two
// copies of the same range form a single part (no self-cluster) while a
// third, overlapping range still joins it.
#[test]
fn command_rg_merge_dedups_identical_lines() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("dup.rg");
    std::fs::write(&rg, "chr1:100-200\nchr1:100-200\nchr1:105-205\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "merge", rg.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "chr1:100-200\tchr1(+):100-205\nchr1:105-205\tchr1(+):100-205\n"
    );
}

// rgr parity: both ratios are computed in f32, so a pair whose reciprocal
// overlap equals the f32-rounded `--coverage` exactly (40/50 with 0.8) must
// still be joined. Comparing in f64 used to round the threshold slightly
// above the exact ratio and split such pairs, dropping 208-247 here.
#[test]
fn command_rg_merge_exact_threshold_parity() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("edge.rg");
    std::fs::write(
        &rg,
        "chr1(+):189-255\nchr1(+):192-246\nchr1(+):200-249\nchr1(+):208-247\n",
    )
    .unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "merge", rg.to_str().unwrap(), "--coverage", "0.8"])
        .run();
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 3);
    assert!(lines.contains(&"chr1(+):208-247\tchr1(+):189-255"));
    assert!(lines.contains(&"chr1(+):192-246\tchr1(+):189-255"));
    assert!(lines.contains(&"chr1(+):200-249\tchr1(+):189-255"));
}

#[test]
fn command_cover() {
    let (stdout, _) = cmd(&["cover", &fixture("S288c.rg")]).run();
    let lines = stdout.lines().count();
    assert!(lines == 3 || lines == 4, "line count {lines}");
    assert!(!stdout.contains("S288c"), "species name: {stdout}");
    assert!(!stdout.contains("1-100"), "merged: {stdout}");
    assert!(stdout.contains("1-150"), "covered: {stdout}");

    let (stdout, _) = cmd(&["cover", &fixture("dazzname.rg")]).run();
    let lines = stdout.lines().count();
    assert!(lines == 2 || lines == 3, "line count {lines}");
    assert!(stdout.contains("infile_0/1/0_514"), "chr name: {stdout}");
    assert!(stdout.contains("19-499"), "covered: {stdout}");
}

#[test]
fn command_coverage() {
    let (stdout, _) = cmd(&["coverage", &fixture("S288c.rg"), "-m", "2"]).run();
    let lines = stdout.lines().count();
    assert!(lines == 3 || lines == 4, "line count {lines}");
    assert!(!stdout.contains("S288c"), "species name: {stdout}");
    assert!(!stdout.contains("1-150"), "coverage 1: {stdout}");
    assert!(stdout.contains("90-100"), "coverage 2: {stdout}");
}

#[test]
fn command_coverage_detailed() {
    let (stdout, _) = cmd(&["coverage", &fixture("S288c.rg"), "-m", "1", "-d"]).run();
    let lines = stdout.lines().count();
    assert!(lines == 9 || lines == 10, "line count {lines}");
    assert!(!stdout.contains("S288c"), "species name: {stdout}");
    assert!(stdout.contains("1-89"), "coverage 1: {stdout}");
    assert!(stdout.contains("90-100"), "coverage 2: {stdout}");
    assert!(stdout.contains("190-200"), "coverage 2: {stdout}");
}
