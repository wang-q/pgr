use assert_cmd::prelude::*;
use std::process::Command;

fn run_vcf(args: &[&str]) -> anyhow::Result<String> {
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fas").arg("to-vcf");
    for a in args {
        cmd.arg(a);
    }
    let output = cmd.output()?;
    Ok(String::from_utf8(output.stdout)?)
}

#[test]
fn command_vcf_basic() -> anyhow::Result<()> {
    let stdout = run_vcf(&["tests/fas/example.fas"])?;

    assert!(stdout.starts_with("##fileformat=VCF"), "vcf header");
    assert!(
        stdout.contains(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS288c\tYJM789\tRM11\tSpar"
        ),
        "samples in header"
    );

    let first_row = stdout
        .lines()
        .find(|l| !l.starts_with('#'))
        .unwrap()
        .to_string();
    assert!(first_row.starts_with("I\t"), "chrom is target chr");
    assert_eq!(
        first_row.split('\t').count(),
        13,
        "columns including 4 samples"
    );

    Ok(())
}

#[test]
fn command_vcf_example_fields() -> anyhow::Result<()> {
    let stdout = run_vcf(&["tests/fas/example.fas"])?;

    let mut rows: std::collections::HashMap<i32, Vec<String>> = std::collections::HashMap::new();
    for line in stdout.lines() {
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        rows.insert(0, cols);
    }

    let r = rows.get(&0).unwrap();
    assert_eq!(r[2], ".", "ID is '.'");
    assert_eq!(r[5], ".", "QUAL is '.'");
    assert_eq!(r[6], ".", "FILTER is '.'");
    assert_eq!(r[7], ".", "INFO is '.'");
    assert_ne!(r[4], ".", "ALT is not '.'");

    Ok(())
}

#[test]
fn command_vcf_ydl_basic() -> anyhow::Result<()> {
    let stdout = run_vcf(&["tests/fas_vcf/YDL184C.fas"])?;

    assert!(stdout.starts_with("##fileformat=VCF"), "vcf header");
    assert!(
        stdout.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS288c\twild006\twine003\tbeer050\tbeer007\tSpar"),
        "samples in header"
    );

    let first_row = stdout
        .lines()
        .find(|l| !l.starts_with('#'))
        .unwrap()
        .to_string();
    assert!(first_row.starts_with("IV\t"), "chrom is target chr IV");
    assert_eq!(
        first_row.split('\t').count(),
        15,
        "columns including 6 samples"
    );

    Ok(())
}

#[test]
fn command_vcf_ydl_fields() -> anyhow::Result<()> {
    let stdout = run_vcf(&["tests/fas_vcf/YDL184C.fas"])?;

    let mut rows: std::collections::HashMap<i32, Vec<String>> = std::collections::HashMap::new();
    for line in stdout.lines() {
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        rows.insert(0, cols);
    }

    let r = rows.get(&0).unwrap();
    assert_eq!(r[2], ".", "ID is '.'");
    assert_eq!(r[5], ".", "QUAL is '.'");
    assert_eq!(r[6], ".", "FILTER is '.'");
    assert_eq!(r[7], ".", "INFO is '.'");
    assert_ne!(r[4], ".", "ALT is not '.'");

    Ok(())
}

#[test]
fn command_vcf_sizes_contig_real() -> anyhow::Result<()> {
    let stdout = run_vcf(&[
        "--sizes",
        "tests/fas_vcf/S288c.chr.sizes",
        "tests/fas_vcf/YDL184C.fas",
    ])?;

    assert!(
        stdout.contains("##contig=<ID=IV,length=1531933>"),
        "contig header contains IV from sizes"
    );

    Ok(())
}

#[test]
fn command_vcf_ydl_expected_rows() -> anyhow::Result<()> {
    let stdout = run_vcf(&["tests/fas_vcf/YDL184C.fas"])?;

    let mut rows: std::collections::HashMap<i32, Vec<String>> = std::collections::HashMap::new();
    for line in stdout.lines() {
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        let pos = cols[1].parse::<i32>().unwrap();
        rows.insert(pos, cols);
    }

    let r1 = rows.get(&130401).expect("row at 130401");
    assert_eq!(r1[0], "IV");
    assert_eq!(r1[3], "A");
    assert!(r1[4].contains('G'));
    let gt1 = &r1[9..];
    assert_eq!(gt1, ["0", "0", "0", "0", "0", "1"]);

    let r2 = rows.get(&130402).expect("row at 130402");
    assert_eq!(r2[0], "IV");
    assert_eq!(r2[3], "T");
    assert!(r2[4].contains('C'));
    let gt2 = &r2[9..];
    assert_eq!(gt2, ["0", "0", "0", "0", "0", "1"]);

    let r3 = rows.get(&130495).expect("row at 130495");
    assert_eq!(r3[0], "IV");
    assert_eq!(r3[3], "T");
    assert!(r3[4].contains('G'));
    let gt3 = &r3[9..];
    assert_eq!(gt3, ["0", "0", "0", "0", "1", "0"]);

    Ok(())
}
