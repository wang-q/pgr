use assert_cmd::Command;
use predicates::prelude::*;
use std::fs::File;
use std::io::Write;
use tempfile::tempdir;

#[test]
fn command_axt_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt").arg("--help");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Axt tools"));
    Ok(())
}

#[test]
fn command_axt_sort_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt").arg("sort").arg("--help");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Sort axt files"));
    Ok(())
}

#[test]
fn command_axt_tomaf_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt").arg("tomaf").arg("--help");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Convert from axt to maf format"));
    Ok(())
}

#[test]
fn command_axt_sort_basic() -> anyhow::Result<()> {
    let dir = tempdir()?;
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
    {
        let mut f = File::create(&input_path)?;
        f.write_all(input_content.as_bytes())?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt")
        .arg("sort")
        .arg(&input_path)
        .arg("-o")
        .arg(&output_path);

    cmd.assert().success();

    let output = std::fs::read_to_string(&output_path)?;
    let lines: Vec<&str> = output.lines().filter(|l| !l.is_empty()).collect();

    // Expected order: 1 (start 6), 0 (start 11), 2 (start 31)
    assert!(lines[0].contains("chr1 6 16"));
    assert!(lines[3].contains("chr1 11 21"));
    assert!(lines[6].contains("chr1 31 41"));

    Ok(())
}

#[test]
fn command_axt_tomaf_zero_score() -> anyhow::Result<()> {
    let dir = tempdir()?;
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

    {
        let mut f = File::create(&input_path)?;
        f.write_all(input_content.as_bytes())?;
    }
    {
        let mut f = File::create(&t_sizes_path)?;
        f.write_all(t_sizes_content.as_bytes())?;
    }
    {
        let mut f = File::create(&q_sizes_path)?;
        f.write_all(q_sizes_content.as_bytes())?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt")
        .arg("tomaf")
        .arg(&input_path)
        .arg("--t-sizes")
        .arg(&t_sizes_path)
        .arg("--q-sizes")
        .arg(&q_sizes_path)
        .arg("-o")
        .arg(&output_path);

    cmd.assert().success();

    let output = std::fs::read_to_string(&output_path)?;

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

    Ok(())
}

#[test]
fn command_axt_sort_by_score() -> anyhow::Result<()> {
    let dir = tempdir()?;
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
    {
        let mut f = File::create(&input_path)?;
        f.write_all(input_content.as_bytes())?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt")
        .arg("sort")
        .arg("--by-score")
        .arg(&input_path)
        .arg("-o")
        .arg(&output_path);

    cmd.assert().success();

    let output = std::fs::read_to_string(&output_path)?;
    let lines: Vec<&str> = output.lines().filter(|l| !l.is_empty()).collect();

    // Expected order: 2 (score 200), 0 (score 100), 1 (score 50)
    assert!(lines[0].contains("2 chr1 31 41"));
    assert!(lines[3].contains("0 chr1 11 21"));
    assert!(lines[6].contains("1 chr1 6 16"));

    Ok(())
}

#[test]
fn command_axt_tomaf_basic() -> anyhow::Result<()> {
    let dir = tempdir()?;
    let axt_path = dir.path().join("input.axt");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let output_path = dir.path().join("output.maf");

    let axt_content = "\
0 chr1 11 14 chr2 11 14 - 100
ACTG
ACTG
";
    {
        let mut f = File::create(&axt_path)?;
        f.write_all(axt_content.as_bytes())?;
    }
    {
        let mut f = File::create(&t_sizes_path)?;
        f.write_all(b"chr1 1000\n")?;
    }
    {
        let mut f = File::create(&q_sizes_path)?;
        f.write_all(b"chr2 2000\n")?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt")
        .arg("tomaf")
        .arg(&axt_path)
        .arg("-t")
        .arg(&t_sizes_path)
        .arg("-q")
        .arg(&q_sizes_path)
        .arg("-o")
        .arg(&output_path);

    cmd.assert().success();

    let output = std::fs::read_to_string(&output_path)?;
    assert!(output.contains("scoring=blastz"));
    assert!(output.contains("s chr1"));
    assert!(output.contains("s chr2"));

    // AXT: chr1 11 14 (1-based, inclusive). Length 4.
    // MAF: start 10 (0-based), size 4.
    assert!(output.contains("chr1                         10          4"));

    Ok(())
}

#[test]
fn command_axt_tomaf_split() -> anyhow::Result<()> {
    let dir = tempdir()?;
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
    {
        let mut f = File::create(&axt_path)?;
        f.write_all(axt_content.as_bytes())?;
    }
    {
        let mut f = File::create(&t_sizes_path)?;
        f.write_all(b"chr1 1000\nchr2 1000\n")?;
    }
    {
        let mut f = File::create(&q_sizes_path)?;
        f.write_all(b"chr1 2000\nchr2 2000\n")?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt")
        .arg("tomaf")
        .arg("--t-split")
        .arg(&axt_path)
        .arg("-t")
        .arg(&t_sizes_path)
        .arg("-q")
        .arg(&q_sizes_path)
        .arg("-o")
        .arg(&output_dir);

    cmd.assert().success();

    assert!(output_dir.exists());
    assert!(output_dir.join("chr1.maf").exists());
    assert!(output_dir.join("chr2.maf").exists());

    let output_chr1 = std::fs::read_to_string(output_dir.join("chr1.maf"))?;
    assert!(output_chr1.contains("s chr1                         10"));
    assert!(!output_chr1.contains("s chr2                         20")); // Should not contain the second record

    let output_chr2 = std::fs::read_to_string(output_dir.join("chr2.maf"))?;
    assert!(output_chr2.contains("s chr2                         20"));
    assert!(!output_chr2.contains("s chr1                         10")); // Should not contain the first record

    Ok(())
}

#[test]
fn command_axt_topsl_blockbug() -> anyhow::Result<()> {
    let dir = tempdir()?;
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

    {
        let mut f = File::create(&input_path)?;
        f.write_all(input_content.as_bytes())?;
    }
    {
        let mut f = File::create(&t_sizes_path)?;
        f.write_all(sizes_content.as_bytes())?;
    }
    {
        let mut f = File::create(&q_sizes_path)?;
        f.write_all(sizes_content.as_bytes())?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt")
        .arg("topsl")
        .arg(&input_path)
        .arg("-t")
        .arg(&t_sizes_path)
        .arg("-q")
        .arg(&q_sizes_path)
        .arg("-o")
        .arg(&output_path);

    cmd.assert().success();

    let output = std::fs::read_to_string(&output_path)?;

    // Check against expected
    let expected = "261\t231\t53\t0\t4\t6\t10\t68\t+\tMmUn_71944_35\t1218\t585\t1136\tMmUn_161829_35\t7971\t2395\t3008\t13\t90,16,6,35,127,54,22,24,9,8,3,27,124,\t585,675,691,697,733,863,917,940,965,974,982,985,1012,\t2395,2493,2522,2539,2574,2701,2758,2801,2825,2836,2850,2856,2884,";

    assert!(output.trim().contains(expected.trim()));

    Ok(())
}

#[test]
fn command_axt_tofas_basic() -> anyhow::Result<()> {
    let dir = tempdir()?;
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

    {
        let mut f = File::create(&input_path)?;
        f.write_all(input_content.as_bytes())?;
    }
    {
        let mut f = File::create(&sizes_path)?;
        f.write_all(sizes_content.as_bytes())?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt")
        .arg("tofas")
        .arg(&sizes_path)
        .arg(&input_path)
        .arg("-o")
        .arg(&output_path);

    cmd.assert().success();

    let output = std::fs::read_to_string(&output_path)?;

    // Check for expected FASTA headers and sequences
    assert!(output.contains(">target.chr1(+):11-14"));
    assert!(output.contains(">query.chr2(+):21-24"));

    Ok(())
}

#[test]
fn command_axt_tofas_example() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("axt")
        .arg("tofas")
        .arg("tests/fasr/RM11_1a.chr.sizes")
        .arg("tests/fasr/example.axt")
        .arg("--qname")
        .arg("RM11_1a")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains("target.I(+)"), "name list");
    assert!(stdout.contains("RM11_1a.scaffold_14"), "name list");
    assert!(stdout.contains("(+):3634-3714"), "positive strand");
    assert!(stdout.contains("(-):22732-22852"), "coordinate transformed");

    Ok(())
}
