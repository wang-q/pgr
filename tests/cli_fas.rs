#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use std::io::Write;
use tempfile::TempDir;

#[test]
fn command_invalid() {
    let (_, stderr) = PgrCmd::new().args(&["fas", "foobar"]).run_fail();
    assert!(stderr.contains("recognized"));
}

#[test]
fn command_name() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "name", "tests/fas/example.fas", "-C"])
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
            "-R",
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
            "-R",
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
            "-R",
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
            "-R",
            "tests/fas/name.lst",
            "--strict",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 15);
    assert!(stdout.lines().next().unwrap().contains("Spar")); // >Spar.
}

#[test]
fn command_fas_subset_duplicate_required_no_dup() {
    // A species listed twice in --required must not be emitted as duplicate
    // entries in the output blocks (mirrors the `concat` dedup).
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("dedup.fas");
    fs::write(
        &fas_file,
        ">speciesA.chr1:1-5\nACGTA\n>speciesB.chr1:1-5\nACGTG\n\n",
    )
    .unwrap();

    let name_lst = temp.path().join("names.lst");
    fs::write(&name_lst, "speciesA\nspeciesB\nspeciesA\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "subset",
            fas_file.to_str().unwrap(),
            "-R",
            name_lst.to_str().unwrap(),
        ])
        .run();

    assert_eq!(stdout.matches(">speciesA").count(), 1, "got: {stdout}");
    assert_eq!(stdout.matches(">speciesB").count(), 1, "got: {stdout}");
    assert_eq!(stdout.matches("ACGTA").count(), 1, "got: {stdout}");
}

#[test]
fn command_fas_subset_duplicate_species_in_block_keeps_first() {
    // A block containing the same species name twice keeps the first
    // occurrence (matching `concat`), so the duplicate sequence is not
    // silently dropped in favor of the last one.
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("dup_species.fas");
    fs::write(
        &fas_file,
        ">speciesA.chr1(+):1-5\nACGTA\n>speciesA.chr1(+):6-10\nTTTTT\n>speciesB.chr1(+):1-5\nACGTG\n\n",
    )
    .unwrap();
    let name_lst = temp.path().join("names.lst");
    fs::write(&name_lst, "speciesA\nspeciesB\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "subset",
            fas_file.to_str().unwrap(),
            "-R",
            name_lst.to_str().unwrap(),
        ])
        .run();

    assert_eq!(stdout.matches(">speciesA").count(), 1, "got: {stdout}");
    assert!(stdout.contains("ACGTA"), "got: {stdout}");
    assert!(!stdout.contains("TTTTT"), "got: {stdout}");
    assert_eq!(stdout.matches(">speciesB").count(), 1, "got: {stdout}");
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
            "--replace-tsv",
            "tests/fas/replace.tsv",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 32);
    assert!(stdout.contains(">query.VIII(+)"));

    // fail
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "replace",
            "tests/fas/example.fas",
            "--replace-tsv",
            "tests/fas/replace.fail.tsv",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 24);
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
            "--replace-tsv",
            "tests/fas/replace.remove.tsv",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 16);
    assert!(!stdout.contains("13267-13287"), "block removed");
}

#[test]
fn command_check() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "check",
            "tests/fas/A_tha.pair.fas",
            "-g",
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
            "-g",
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
fn command_fas_check_out_of_range_is_failed_not_abort() {
    // A block FA coordinate beyond the reference chromosome length must be
    // reported as FAILED for that line, not abort the whole check command.
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("oob.fas");
    fs::write(&fas_file, ">NC_000932:1-999999\nATGGGCGAAC\n\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "check",
            fas_file.to_str().unwrap(),
            "-g",
            "tests/fas/NC_000932.fa",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 1);
    assert!(
        stdout.lines().next().unwrap().contains("\tFAILED"),
        "out-of-range coordinate must be FAILED, got {}",
        stdout
    );
}

#[test]
fn command_create() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "create",
            "tests/fas/I.connect.tsv",
            "-g",
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
fn command_separate_rc_preserves_non_iupac() {
    // `separate --rc` reverse-complements minus-strand sequences. Non-IUPAC
    // bytes such as `*` must be preserved as-is, not mangled to the 255
    // sentinel (which previously rendered as `ÿ`).
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("rc.fas");
    fs::write(&fas_file, ">sp.chr1(-):1-5\nACG*T\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "separate", fas_file.to_str().unwrap(), "--rc"])
        .run();

    // Reverse complement of "ACG*T" is "A*CGT" (dashes removed in both).
    assert!(stdout.contains("A*CGT"), "got: {stdout}");
    assert!(!stdout.contains('ÿ'), "got: {stdout}");
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
            "-s",
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
            "-s",
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
fn command_split_to_simple() {
    let tempdir = TempDir::new().unwrap();
    let tempdir_str = tempdir.path().to_str().unwrap();

    PgrCmd::new()
        .args(&[
            "fas",
            "split",
            "tests/fas/example.fas",
            "--chr",
            "--simple",
            "-o",
            tempdir_str,
        ])
        .assert()
        .success()
        .stdout(predicates::str::is_empty());

    let content = fs::read_to_string(tempdir.path().join("S288c.I.fas")).unwrap();
    assert!(content.contains(">S288c\n"), "simple header in file output");
    assert!(!content.contains("I(+)"), "no positions in file output");

    tempdir.close().unwrap();
}

#[test]
fn command_refine() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "refine", "tests/fas/example.fas", "--engine", "none"])
        .run();

    assert_eq!(stdout.lines().count(), 27);

    // --parallel 2
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "refine",
            "tests/fas/example.fas",
            "--engine",
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
            "--engine",
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
            "--engine",
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
            "--runlist",
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

/// A reference species whose own sequence has a gap covers fewer genomic
/// positions than its range length. A runlist spanning the full range used
/// to abort the whole slice immediately; it must instead skip the subspan
/// and keep processing.
#[test]
fn command_slice_gapped_reference_no_abort() {
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("gap_ref.fas");
    fs::write(
        &fas_file,
        ">Ref.chr1(+):1-5\nA-CGT\n>Oth.chr2(+):1-5\nATGCA\n",
    )
    .unwrap();
    let runlist = temp.path().join("all.json");
    fs::write(&runlist, "{\"chr1\": \"1-5\"}").unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "slice",
            fas_file.to_str().unwrap(),
            "--runlist",
            runlist.to_str().unwrap(),
            "--name",
            "Ref",
        ])
        .run();

    // The command must not abort; the out-of-range subspan is skipped.
    assert!(
        stderr.contains("skipping slice subspan"),
        "stderr: {stderr}"
    );
    assert!(
        stdout.trim().is_empty(),
        "stdout should be empty, got: {stdout}"
    );
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
fn command_fas_stat_unequal_length_no_panic() {
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("malformed.fas");
    fs::write(
        &fas_file,
        ">a.chr1(+):1-10\nAAAATTTTGG\n>b.chr2(+):1-10\nAAAATTTTAG\n>c.chr3(+):1-8\nAAAATTTT\n",
    )
    .unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&["fas", "stat", fas_file.to_str().unwrap()])
        .run();
    assert!(
        stderr.contains("unequal lengths"),
        "expected a friendly error, not a panic; stderr={}",
        stderr
    );
    assert!(stdout.trim().is_empty() || stdout.starts_with("target"));
}

#[test]
fn command_fas_refine_unequal_length_no_panic() {
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("malformed.fas");
    fs::write(
        &fas_file,
        ">a.chr1(+):1-10\nAAAATTTTGG\n>b.chr2(+):1-10\nAAAATTTTAG\n>c.chr3(+):1-8\nAAAATTTT\n",
    )
    .unwrap();

    // `refine` realigns sequences of unequal length (its whole purpose), so it
    // must succeed and emit the (shortest-anchored) alignment rather than
    // error or panic.
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "refine",
            fas_file.to_str().unwrap(),
            "--engine",
            "none",
            "--chop",
            "2",
        ])
        .run();
    assert!(
        stderr.is_empty(),
        "expected success, not an error or panic; stderr={}",
        stderr
    );
    assert!(
        stdout.contains(">a.chr1") && stdout.contains(">b.chr2") && stdout.contains(">c.chr3"),
        "expected all three species in the refined output; stdout={}",
        stdout
    );
}

#[test]
fn command_filter() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "filter", "tests/fas/example.fas"])
        .run();

    assert_eq!(stdout.lines().count(), 27);

    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "filter", "tests/fas/example.fas", "--min-len", "30"])
        .run();

    assert_eq!(stdout.lines().count(), 18);

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "filter",
            "tests/fas/example.fas",
            "--min-len",
            "30",
            "--max-len",
            "100",
            "--name",
            "S288c",
            "--dash",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 9);
    assert!(stdout.contains("\nGCTAAAATATGAACG"), "no dash");
}

#[test]
fn command_fas_concat_phylip_unequal_lengths() {
    let temp = TempDir::new().unwrap();
    // Create a fas file where species have different total lengths
    let fas_file = temp.path().join("unequal.fas");
    fs::write(
        &fas_file,
        ">speciesA.chr1:1-10\nACGTACGTAC\n>speciesB.chr1:1-8\nACGTACGT\n",
    )
    .unwrap();

    let name_lst = temp.path().join("names.lst");
    fs::write(&name_lst, "speciesA\nspeciesB\n").unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "concat",
            "--required",
            name_lst.to_str().unwrap(),
            fas_file.to_str().unwrap(),
            "--phylip",
        ])
        .run_fail();

    assert!(stderr.contains("PHYLIP requires equal-length sequences"));
}

#[test]
fn command_fas_stat_outgroup_single_entry() {
    let temp = TempDir::new().unwrap();
    // Single-entry block
    let fas_file = temp.path().join("single.fas");
    fs::write(&fas_file, ">chr1.speciesA 1-10\nACGTACGTAC\n").unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&["fas", "stat", fas_file.to_str().unwrap(), "--outgroup"])
        .run_fail();

    assert!(stderr.contains("cannot apply --outgroup"));
}

#[test]
fn command_fas_stat_outgroup_two_entries() {
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("two.fas");
    fs::write(&fas_file, ">target.chr1:1-5\nACGTA\n>out.chr1:1-5\nACGTC\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "stat", fas_file.to_str().unwrap(), "--outgroup"])
        .run();

    assert_eq!(stdout.lines().count(), 2, "header plus one block");
    let data_line = stdout.lines().nth(1).unwrap();
    let cols: Vec<&str> = data_line.split('\t').collect();
    assert_eq!(cols[1], "5", "length");
    assert_eq!(cols[6], "0", "D after excluding outgroup");
}

#[test]
fn command_fas_concat_required_order() {
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("order.fas");
    fs::write(
        &fas_file,
        ">speciesA.chr1:1-5\nACGTA\n>speciesB.chr1:1-5\nACGTG\n\n",
    )
    .unwrap();

    let name_lst = temp.path().join("names.lst");
    fs::write(&name_lst, "speciesB\nspeciesA\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "concat",
            fas_file.to_str().unwrap(),
            "-R",
            name_lst.to_str().unwrap(),
        ])
        .run();

    let first_header = stdout.lines().next().unwrap();
    assert!(
        first_header.starts_with(">speciesB"),
        "output should follow --required order, got {}",
        first_header
    );
}

#[test]
fn command_fas_concat_duplicate_required_no_dup() {
    // A species listed twice in --required must not be concatenated twice nor
    // emitted as duplicate output lines.
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("dedup.fas");
    fs::write(
        &fas_file,
        ">speciesA.chr1:1-5\nACGTA\n>speciesB.chr1:1-5\nACGTG\n\n",
    )
    .unwrap();

    let name_lst = temp.path().join("names.lst");
    fs::write(&name_lst, "speciesA\nspeciesB\nspeciesA\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "concat",
            fas_file.to_str().unwrap(),
            "-R",
            name_lst.to_str().unwrap(),
        ])
        .run();

    assert_eq!(stdout.matches(">speciesA").count(), 1, "got: {stdout}");
    assert_eq!(stdout.matches(">speciesB").count(), 1, "got: {stdout}");
    // No duplicated sequence: speciesA appears exactly once.
    assert_eq!(stdout.matches("ACGTA").count(), 1, "got: {stdout}");
}

#[test]
fn command_fas_replace_duplicate_header() {
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("dup.fas");
    fs::write(
        &fas_file,
        ">target.chr1:1-5\nACGTA\n>target.chr1:1-5\nACGTC\n\n",
    )
    .unwrap();

    let tsv = temp.path().join("dup.tsv");
    fs::write(&tsv, "target.chr1:1-5\n").unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "replace",
            fas_file.to_str().unwrap(),
            "--replace-tsv",
            tsv.to_str().unwrap(),
        ])
        .run();

    assert_eq!(
        stdout.matches(">target.").count(),
        2,
        "duplicate header block should be kept unchanged"
    );
    assert!(
        stderr.contains("appears") || stderr.contains("keeping block unchanged"),
        "expected warning about duplicate header, got {}",
        stderr
    );
}

#[test]
fn command_fas_concat_empty_required() {
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("empty_required.fas");
    fs::write(
        &fas_file,
        ">speciesA.chr1:1-5\nACGTA\n>speciesB.chr1:1-5\nACGTG\n",
    )
    .unwrap();

    let name_lst = temp.path().join("empty.lst");
    fs::write(&name_lst, "").unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "concat",
            "--required",
            name_lst.to_str().unwrap(),
            fas_file.to_str().unwrap(),
        ])
        .run_fail();

    assert!(
        stderr.contains("required file is empty"),
        "expected empty --required error, got {}",
        stderr
    );
}

#[test]
fn command_consensus() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "consensus",
            "tests/fas/example.fas",
            "--engine",
            "builtin",
        ])
        .run();

    assert!(stdout.contains(">consensus"), "consensus header");
    assert!(stdout.lines().count() > 2, "has header and sequence");
}

#[test]
fn command_variation() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "variation", "tests/fas/example.fas"])
        .run();

    assert!(stdout.contains("#target"), "header line");
    assert!(stdout.lines().count() > 1, "has data rows");
}

#[test]
fn command_to_vcf() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "to-vcf",
            "--sizes",
            "tests/fas_vcf/S288c.chr.sizes",
            "tests/fas_vcf/YDL184C.fas",
        ])
        .run();

    assert!(stdout.starts_with("##fileformat=VCFv4.2"), "vcf header");
    assert!(stdout.contains("##contig=<ID="), "contig header");
    assert!(stdout.contains("\nIV\t"), "data line");
}

#[test]
fn command_to_xlsx() {
    let temp = TempDir::new().unwrap();
    let out = temp.path().join("variations.xlsx");

    PgrCmd::new()
        .args(&[
            "fas",
            "to-xlsx",
            "tests/fas/example.fas",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert!(out.is_file());
    assert!(fs::metadata(&out).unwrap().len() > 0);
}

#[test]
fn command_multiz() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "multiz",
            "-r",
            "S288c",
            "tests/fas/S288cvsRM11_1a.slice.fas",
            "tests/fas/S288cvsSpar.slice.fas",
        ])
        .run();

    assert!(stdout.contains(">S288c."), "has S288c entry");
    assert!(stdout.lines().count() > 2, "has output blocks");
}

#[test]
fn command_multiz_gzip() {
    let temp = TempDir::new().unwrap();

    let rm11 = fs::read("tests/fas/S288cvsRM11_1a.slice.fas").unwrap();
    let spar = fs::read("tests/fas/S288cvsSpar.slice.fas").unwrap();

    let rm11_gz = temp.path().join("rm11.fas.gz");
    let spar_gz = temp.path().join("spar.fas.gz");

    {
        let mut encoder = flate2::write::GzEncoder::new(
            fs::File::create(&rm11_gz).unwrap(),
            flate2::Compression::default(),
        );
        encoder.write_all(&rm11).unwrap();
        encoder.finish().unwrap();
    }
    {
        let mut encoder = flate2::write::GzEncoder::new(
            fs::File::create(&spar_gz).unwrap(),
            flate2::Compression::default(),
        );
        encoder.write_all(&spar).unwrap();
        encoder.finish().unwrap();
    }

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "multiz",
            "-r",
            "S288c",
            rm11_gz.to_str().unwrap(),
            spar_gz.to_str().unwrap(),
        ])
        .run();

    assert!(stdout.contains(">S288c."), "has S288c entry");
    assert!(stdout.lines().count() > 2, "has output blocks");
}

#[test]
fn command_fas_stat_outgroup_length_consistent() {
    let (stdout_no, _) = PgrCmd::new()
        .args(&["fas", "stat", "tests/fas/example.fas"])
        .run();
    let (stdout_og, _) = PgrCmd::new()
        .args(&["fas", "stat", "tests/fas/example.fas", "--outgroup"])
        .run();

    let lines_no: Vec<&str> = stdout_no.lines().collect();
    let lines_og: Vec<&str> = stdout_og.lines().collect();
    assert_eq!(lines_no.len(), lines_og.len());

    for (no, og) in lines_no.iter().skip(1).zip(lines_og.iter().skip(1)) {
        let cols_no: Vec<&str> = no.split('\t').collect();
        let cols_og: Vec<&str> = og.split('\t').collect();
        assert_eq!(
            cols_no[1], cols_og[1],
            "length should be consistent with and without --outgroup"
        );
    }
}

#[test]
fn command_refine_skips_malformed_block() {
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("malformed.fas");
    fs::write(
        &fas_file,
        ">target.chr1:1-5\nACGTA\nACGT\n\n>target.chr1:1-5\nACGTA\n>out.chr1:1-5\nACGTC\n",
    )
    .unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "refine",
            fas_file.to_str().unwrap(),
            "--engine",
            "none",
        ])
        .run();

    assert_eq!(
        stdout.lines().count(),
        5,
        "only the valid block should be output"
    );
    assert!(
        stderr.contains("skipping malformed fas block"),
        "expected warning about malformed block, got {}",
        stderr
    );
}

#[test]
fn command_filter_upper() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "filter", "tests/fas/example.fas", "--upper"])
        .run();

    assert_eq!(stdout.lines().count(), 27);
    let seq_lower = stdout
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.chars())
        .filter(|c| c.is_ascii_lowercase())
        .count();
    assert_eq!(seq_lower, 0, "all sequence characters should be uppercase");
}

#[test]
fn command_slice_default_name() {
    let (stdout_with_name, _) = PgrCmd::new()
        .args(&[
            "fas",
            "slice",
            "tests/fas/slice.fas",
            "--runlist",
            "tests/fas/slice.json",
            "--name",
            "S288c",
        ])
        .run();

    let (stdout_default, _) = PgrCmd::new()
        .args(&[
            "fas",
            "slice",
            "tests/fas/slice.fas",
            "--runlist",
            "tests/fas/slice.json",
        ])
        .run();

    assert_eq!(stdout_with_name, stdout_default);
}

#[test]
fn command_replace_three_fields() {
    let temp = TempDir::new().unwrap();
    let fas_file = temp.path().join("replace3.fas");
    fs::write(
        &fas_file,
        ">target.chr1:1-5\nACGTA\n>query.chr1:1-5\nACGTC\n",
    )
    .unwrap();

    let tsv = temp.path().join("replace3.tsv");
    fs::write(&tsv, "target.chr1:1-5\tnewA\tnewB\n").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fas",
            "replace",
            fas_file.to_str().unwrap(),
            "--replace-tsv",
            tsv.to_str().unwrap(),
        ])
        .run();

    assert_eq!(
        stdout.matches(">newA").count(),
        1,
        "first replacement name should appear once"
    );
    assert_eq!(
        stdout.matches(">newB").count(),
        1,
        "second replacement name should appear once"
    );
    assert_eq!(
        stdout.matches(">target.").count(),
        0,
        "original name should be replaced"
    );
}

#[test]
fn command_create_skips_invalid_range() {
    let temp = TempDir::new().unwrap();
    let connect = temp.path().join("connect.tsv");
    fs::write(&connect, "S288c.I:1-10\tinvalid_range\nS288c.I:11-20\n").unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "create",
            connect.to_str().unwrap(),
            "-g",
            "tests/fas/genome.fa",
            "--name",
            "S288c",
        ])
        .run();

    assert!(
        stdout.contains(">S288c.I:11-20"),
        "valid range should produce output"
    );
    assert!(
        !stdout.contains("invalid_range"),
        "invalid range should be skipped"
    );
    assert!(
        stderr.contains("skipping invalid range"),
        "expected warning about invalid range, got {}",
        stderr
    );
}

#[test]
fn command_create_output_not_overwrite_loc_index() {
    // `create` opens the output writer (truncating) before reading the
    // reference genome's `.loc` sidecar index. If `-o` names the `.loc` file,
    // the index is truncated before `open_indexed` loads it, silently dropping
    // every link. It must be rejected.
    let temp = TempDir::new().unwrap();
    let genome = temp.path().join("genome.fa");
    fs::write(&genome, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    let loc = format!("{}.loc", genome.display());

    let connect = temp.path().join("connect.tsv");
    fs::write(&connect, "S288c.chr1:1-5\n").unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "create",
            connect.to_str().unwrap(),
            "-g",
            genome.to_str().unwrap(),
            "-o",
            &loc,
        ])
        .run_fail();

    assert!(
        stderr.contains("is also an input file"),
        "expected rejection of -o matching .loc, got {}",
        stderr
    );
}

#[test]
fn command_create_skips_out_of_range_link() {
    // A link coordinate beyond the reference chromosome length must be
    // skipped (with a warning) rather than aborting the whole create run.
    let temp = TempDir::new().unwrap();
    let genome = temp.path().join("genome.fa");
    fs::write(&genome, ">chr1\nACGTACGTACGTACGT\n").unwrap();

    let connect = temp.path().join("connect.tsv");
    fs::write(
        &connect,
        "A.chr1:1-4\tB.chr1:1-4\nA.chr1:100-110\tB.chr1:1-4\n",
    )
    .unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "fas",
            "create",
            connect.to_str().unwrap(),
            "-g",
            genome.to_str().unwrap(),
        ])
        .run();

    // The in-range block is emitted; the out-of-range range is skipped.
    assert!(stdout.starts_with(">A.chr1:1-4\n"), "got: {stdout}");
    assert!(
        !stdout.contains(">A.chr1:100-110"),
        "out-of-range range must be skipped, got: {stdout}"
    );
    // Both valid ranges (line 1 and line 2 `B.chr1:1-4`) are fetched; a run
    // that aborts on the out-of-range range would emit only one.
    assert_eq!(
        stdout.matches(">B.chr1:1-4\n").count(),
        2,
        "expected both valid ranges to be emitted, got: {stdout}"
    );
    assert!(
        stderr.contains("out-of-range"),
        "expected a warning about the out-of-range range, got: {stderr}"
    );
}

#[test]
fn command_split_simple_stdout() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fas", "split", "tests/fas/example.fas", "--simple"])
        .run();

    assert!(
        stdout.contains(">S288c\n"),
        "simple headers should use species names only"
    );
    assert!(
        !stdout.contains("I(+)"),
        "simple headers should not contain coordinates"
    );
}
