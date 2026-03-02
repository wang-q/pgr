#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use tempfile::TempDir;

#[test]
fn command_axt_help() {
    let (stdout, _) = PgrCmd::new().args(&["axt", "--help"]).run();
    assert!(stdout.contains("Usage: pgr") || stdout.contains(" axt"));
}

#[test]
fn command_axt_sort_help() {
    let (stdout, _) = PgrCmd::new().args(&["axt", "sort", "--help"]).run();
    assert!(stdout.contains("Sort axt files"));
}

#[test]
fn command_axt_to_maf_help() {
    let (stdout, _) = PgrCmd::new().args(&["axt", "to-maf", "--help"]).run();
    assert!(stdout.contains("Convert from axt to maf format"));
}

#[test]
fn command_axt_sort_basic() {
    let dir = TempDir::new().unwrap();
    let input_path = dir.path().join("input.axt");
    let output_path = dir.path().join("output.axt");

    let input_content = "\
0 chr1 11 21 chr2 11 21 - 100
ACTG
ACTG

1 chr1 6 16 chr2 31 41 + 50
AAAA
AAAA

2 chr1 31 41 chr2 6 16 + 200
TTTT
TTTT
";
    fs::write(&input_path, input_content).unwrap();

    PgrCmd::new()
        .args(&[
            "axt",
            "sort",
            input_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success();

    let output = fs::read_to_string(&output_path).unwrap();
    let lines: Vec<&str> = output.lines().filter(|l| !l.is_empty()).collect();

    // Expected order: 1 (start 6), 0 (start 11), 2 (start 31)
    assert!(lines[0].contains("chr1 6 16"));
    assert!(lines[3].contains("chr1 11 21"));
    assert!(lines[6].contains("chr1 31 41"));
}

#[test]
fn command_axt_to_maf_zero_score() {
    let dir = TempDir::new().unwrap();
    let input_path = dir.path().join("zeroScore.axt");
    let t_sizes_path = dir.path().join("rn3.sizes");
    let q_sizes_path = dir.path().join("hg15.sizes");
    let output_path = dir.path().join("zeroScore.maf");

    // input/zeroScore.axt
    let input_content = "\
7579 chr8 16067581 16067602 chr19 54649198 54649224 - 0
ATTCGCTGGTATACAT---------GAGGTC
----GCTGATGTATCTACAAAAGAAGAGGTC
";
    // input/rn3.sizes (subset)
    let t_sizes_content = "\
chr8	129061546
";
    // input/hg15.sizes (subset)
    let q_sizes_content = "\
chr19	63790860
";

    fs::write(&input_path, input_content).unwrap();
    fs::write(&t_sizes_path, t_sizes_content).unwrap();
    fs::write(&q_sizes_path, q_sizes_content).unwrap();

    PgrCmd::new()
        .args(&[
            "axt",
            "to-maf",
            input_path.to_str().unwrap(),
            "--t-sizes",
            t_sizes_path.to_str().unwrap(),
            "--q-sizes",
            q_sizes_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success();

    let output = fs::read_to_string(&output_path).unwrap();

    // Expected output check
    // s chr8  16067580 22 + 129061546 ATTCGCTGGTATACAT---------GAGGTC
    // s chr19 54649197 27 -  63790860 ----GCTGATGTATCTACAAAAGAAGAGGTC

    // Check score=0.0 (pgr output format)
    assert!(output.contains("score=0.0"));

    // Check target line components
    assert!(output.contains("chr8"));
    assert!(output.contains("16067580"));
    assert!(output.contains("22"));
    assert!(output.contains("+"));
    assert!(output.contains("129061546"));
    assert!(output.contains("ATTCGCTGGTATACAT---------GAGGTC"));

    // Check query line components
    assert!(output.contains("chr19"));
    assert!(output.contains("54649197"));
    assert!(output.contains("27"));
    assert!(output.contains("-"));
    assert!(output.contains("63790860"));
    assert!(output.contains("----GCTGATGTATCTACAAAAGAAGAGGTC"));
}

#[test]
fn command_axt_sort_by_score() {
    let dir = TempDir::new().unwrap();
    let input_path = dir.path().join("input.axt");
    let output_path = dir.path().join("output.axt");

    let input_content = "\
0 chr1 11 21 chr2 11 21 - 100
ACTG
ACTG

1 chr1 6 16 chr2 31 41 + 50
AAAA
AAAA

2 chr1 31 41 chr2 6 16 + 200
TTTT
TTTT
";
    fs::write(&input_path, input_content).unwrap();

    PgrCmd::new()
        .args(&[
            "axt",
            "sort",
            "--by-score",
            input_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success();

    let output = fs::read_to_string(&output_path).unwrap();
    let lines: Vec<&str> = output.lines().filter(|l| !l.is_empty()).collect();

    // Expected order: 2 (score 200), 0 (score 100), 1 (score 50)
    assert!(lines[0].contains("2 chr1 31 41"));
    assert!(lines[3].contains("0 chr1 11 21"));
    assert!(lines[6].contains("1 chr1 6 16"));
}

#[test]
fn command_axt_to_maf_basic() {
    let dir = TempDir::new().unwrap();
    let axt_path = dir.path().join("input.axt");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let output_path = dir.path().join("output.maf");

    let axt_content = "\
0 chr1 11 14 chr2 11 14 - 100
ACTG
ACTG
";
    fs::write(&axt_path, axt_content).unwrap();
    fs::write(&t_sizes_path, b"chr1 1000\n").unwrap();
    fs::write(&q_sizes_path, b"chr2 2000\n").unwrap();

    PgrCmd::new()
        .args(&[
            "axt",
            "to-maf",
            axt_path.to_str().unwrap(),
            "-t",
            t_sizes_path.to_str().unwrap(),
            "-q",
            q_sizes_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success();

    let output = fs::read_to_string(&output_path).unwrap();
    assert!(output.contains("scoring=blastz"));
    assert!(output.contains("s chr1"));
    assert!(output.contains("s chr2"));

    // AXT: chr1 11 14 (1-based, inclusive). Length 4.
    // MAF: start 10 (0-based), size 4.
    assert!(output.contains("chr1                         10          4"));
}

#[test]
fn command_axt_to_maf_split() {
    let dir = TempDir::new().unwrap();
    let axt_path = dir.path().join("input.axt");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let output_dir = dir.path().join("output_split");

    let axt_content = "\
0 chr1 11 14 chr2 11 14 - 100
ACTG
ACTG

1 chr2 21 24 chr1 21 24 + 100
ACTG
ACTG
";
    fs::write(&axt_path, axt_content).unwrap();
    fs::write(&t_sizes_path, b"chr1 1000\nchr2 1000\n").unwrap();
    fs::write(&q_sizes_path, b"chr1 2000\nchr2 2000\n").unwrap();

    PgrCmd::new()
        .args(&[
            "axt",
            "to-maf",
            "--t-split",
            axt_path.to_str().unwrap(),
            "-t",
            t_sizes_path.to_str().unwrap(),
            "-q",
            q_sizes_path.to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert!(output_dir.exists());
    assert!(output_dir.join("chr1.maf").exists());
    assert!(output_dir.join("chr2.maf").exists());

    let output_chr1 = fs::read_to_string(output_dir.join("chr1.maf")).unwrap();
    assert!(output_chr1.contains("s chr1                         10"));
    assert!(!output_chr1.contains("s chr2                         20")); // Should not contain the second record

    let output_chr2 = fs::read_to_string(output_dir.join("chr2.maf")).unwrap();
    assert!(output_chr2.contains("s chr2                         20"));
    assert!(!output_chr2.contains("s chr1                         10")); // Should not contain the first record
}

#[test]
fn command_axt_to_psl_blockbug() {
    let dir = TempDir::new().unwrap();
    let input_path = dir.path().join("blockBug.axt");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let output_path = dir.path().join("output.psl");

    // BlockBug content
    let input_content = "\
183 MmUn_161829_35 2396 3008 MmUn_71944_35 586 1136 + 3228
GTCAAAGTGCTAAACTGTGGTGTATCAGTACTGCAATATTCTGTTACAGTCCTCCTCTGTCGTGTCACAGTTCTCAAAACTGCTGTCCCACTGTTCCAGTTTTGTATCACTGTGCAGCATTTTTCTGTAACAGTCATCTACTGCGGTGTCAAGGTGCTCACCTGTGGTGCCACAGTGCT-CCACTGTGGTGTCACATTGCTCCACTGTGGTGTCAGAGTGATCCgctgtggtgtcacagtgctccactgtgctgtcacagtgctgcactctgttgtcactgtgctccaatgtggtgtcacagtgctc---cactgtgctgtcacagtcctccactgtggtgtCACATTGCTCCACTGTTGTGTCAGAGTGATCCCCTGTGGTTTCACAGTGCTCCACTGT-GATGTCGCAGTGCTCCATGTTGTGTCACAGTTCT-CCACTGTGCTGTCACGGTGCACCCTGTGGAGTCACAGTGCTCCCCTGTATTGTCAAAGTGCTCCACTTTCATGTCAAAGTGCTCCACTGTGCTATCAAAGTGCTCTATTGTGATGTCACAGTGTTCAAGTTTGTTTTAACAGTGCTCCACTGTGGTGTGACAGTTTTCCACTTTGTTGTGACA
GTCACAGTTCACCAATGTGGTGTCACAGTGCTCCACTGTGGTGTCACAGTGCTCCACTGTGATGTCACAATGCTCCACTGTGATGTCAGA--------GTGCTCCAATAATGTG-------------TCAAAG-----------AGCTCCATTTTGTTGTCATGGTACTCCACTGTGCTGTCACAGTGCTCCTATATGGTGTCATTGTGCTCCACTTTTAAGTCACAGTGCTCCTCTGTGGCTTCACAATGCTGCACTCTTATTTCACAGTGCTCCAACTGTGCTGTCACATTGCTTCACTCTGCTTTCACAGTGCACCACTATGATGTCATAATGCTCCACTGTGATGTCACAGTGCTCTGCT---GTGGTGTCACAGTGCTCCACTG-----------T----------TGAGTCACAGTGCTCCAGTGTGGTGTCACAGTGC--TCAATGTG------GTG---TCAAAGTGCTCCACTTTTGTTTCACAG-ACTCCACTGTGGAGACCGCATTCTCTATTGTACTGTCACAGTGCACCACTGTGATGACATAGTTCTCCCCTGTGATGTCAGAGTCTTCCAGCTAGATGTTACAGTGTTCCATTGTGCTGTAACA
";

    // blockBug.sizes:
    // MmUn_71944_35	1218
    // MmUn_161829_35	7971
    let sizes_content = "\
MmUn_71944_35\t1218
MmUn_161829_35\t7971
";

    fs::write(&input_path, input_content).unwrap();
    fs::write(&t_sizes_path, sizes_content).unwrap();
    fs::write(&q_sizes_path, sizes_content).unwrap();

    PgrCmd::new()
        .args(&[
            "axt",
            "to-psl",
            input_path.to_str().unwrap(),
            "-t",
            t_sizes_path.to_str().unwrap(),
            "-q",
            q_sizes_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success();

    let output = fs::read_to_string(&output_path).unwrap();

    // Check against expected
    let expected = "261\t231\t53\t0\t4\t6\t10\t68\t+\tMmUn_71944_35\t1218\t585\t1136\tMmUn_161829_35\t7971\t2395\t3008\t13\t90,16,6,35,127,54,22,24,9,8,3,27,124,\t585,675,691,697,733,863,917,940,965,974,982,985,1012,\t2395,2493,2522,2539,2574,2701,2758,2801,2825,2836,2850,2856,2884,";

    assert!(output.trim().contains(expected.trim()));
}

#[test]
fn command_axt_to_fas_basic() {
    let dir = TempDir::new().unwrap();
    let input_path = dir.path().join("input.axt");
    let sizes_path = dir.path().join("q.sizes");
    let output_path = dir.path().join("output.fas");

    // input.axt
    // 0 chr1 11 14 chr2 21 24 + 100
    // ACGT
    // ACGT
    let input_content = "\
0 chr1 11 14 chr2 21 24 + 100
ACGT
ACGT
";
    // q.sizes
    let sizes_content = "chr2\t100\n";

    fs::write(&input_path, input_content).unwrap();
    fs::write(&sizes_path, sizes_content).unwrap();

    PgrCmd::new()
        .args(&[
            "axt",
            "to-fas",
            sizes_path.to_str().unwrap(),
            input_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success();

    let output = fs::read_to_string(&output_path).unwrap();

    // Check for expected FASTA headers and sequences
    assert!(output.contains(">target.chr1(+):11-14"));
    assert!(output.contains(">query.chr2(+):21-24"));
}

#[test]
fn command_axt_to_fas_example() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "axt",
            "to-fas",
            "tests/axt/RM11_1a.sizes",
            "tests/axt/example.axt",
            "--qname",
            "RM11_1a",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains("target.I(+)"), "name list");
    assert!(stdout.contains("RM11_1a.scaffold_14"), "name list");
    assert!(stdout.contains("(+):3634-3714"), "positive strand");
    assert!(stdout.contains("(-):22732-22852"), "coordinate transformed");
}
