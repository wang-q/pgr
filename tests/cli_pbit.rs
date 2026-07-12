#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

/// Return the absolute path to a fixture in `tests/pbit/input`.
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/pbit/input")
        .join(name)
}

/// Read a FASTA file and return the concatenated sequence (all non-header
/// lines joined, uppercased).
fn read_fasta_seq(path: &std::path::Path) -> String {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect::<String>()
}

/// Read a FASTA file and return its records as (name, sequence) pairs.
fn read_fasta_records(path: &std::path::Path) -> Vec<(String, String)> {
    let content = fs::read_to_string(path).unwrap();
    let mut records = Vec::new();
    let mut name = String::new();
    let mut seq = String::new();
    for line in content.lines() {
        if let Some(stripped) = line.strip_prefix('>') {
            if !name.is_empty() {
                records.push((name.clone(), std::mem::take(&mut seq)));
            }
            name = stripped.to_string();
        } else {
            seq.push_str(line);
        }
    }
    if !name.is_empty() {
        records.push((name, seq));
    }
    records
}

#[test]
fn test_pbit_create_basic() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    assert!(out_pbit.exists());
    assert!(fs::metadata(&out_pbit).unwrap().len() > 0);
}

#[test]
fn test_pbit_stat_overview() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap()])
        .run();

    assert!(stdout.contains("Version:"));
    assert!(stdout.contains("Segment size: 4096"));
    assert!(stdout.contains("K-mer length: 15"));
    assert!(stdout.contains("Reference groups: 1"));
    assert!(stdout.contains("Samples: 1"));
}

#[test]
fn test_pbit_stat_samples() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let name_tsv = temp.path().join("names.tsv");
    let sample = fixture("sample_2000_identical.fa");
    fs::write(
        &name_tsv,
        format!("s1\t{}\ns2\t{}\n", sample.display(), sample.display()),
    )
    .unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "--name",
            name_tsv.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--samples"])
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines.contains(&"s1"));
    assert!(lines.contains(&"s2"));
    assert_eq!(lines.len(), 2);
}

#[test]
fn test_pbit_stat_refs() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_5000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_5000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--refs"])
        .run();

    // 5000 bp / 4096 segment_size -> 2 segments for chr1.
    assert_eq!(stdout.trim(), "chr1\t2");
}

#[test]
fn test_pbit_stat_contigs() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_1000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_1000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--contigs"])
        .run();

    let sample_path = fixture("sample_1000_identical.fa");
    let stem = sample_path.file_stem().unwrap().to_str().unwrap();
    assert_eq!(stdout.trim(), format!("{}\tchr1", stem));
}

#[test]
fn test_pbit_to_fa_roundtrip() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");
    let sample = fixture("sample_2000_snps100.fa");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            sample.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    let stem = sample.file_stem().unwrap().to_str().unwrap();
    let out_fa = out_dir.join(format!("{}.fa", stem));
    assert!(out_fa.exists());

    let expected = read_fasta_seq(&sample).to_ascii_uppercase();
    let extracted = read_fasta_seq(&out_fa);
    assert_eq!(extracted, expected);
}

#[test]
fn test_pbit_range_full_contig() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "range", out_pbit.to_str().unwrap(), "chr1"])
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines[0].starts_with('>') && lines[0].contains(" chr1"));
    let seq: String = lines[1..].concat();
    assert_eq!(seq.len(), 2000);
}

#[test]
fn test_pbit_range_slice() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "range", out_pbit.to_str().unwrap(), "chr1:1-100"])
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines[0].starts_with('>') && lines[0].contains("chr1:1-100(+)"));
    let seq: String = lines[1..].concat();
    assert_eq!(seq.len(), 100);
    let expected = read_fasta_seq(&fixture("ref_2000.fa"))[..100].to_uppercase();
    assert_eq!(seq, expected);
}

#[test]
fn test_pbit_range_neg_strand() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "range", out_pbit.to_str().unwrap(), "chr1(-):1-100"])
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines[0].starts_with('>') && lines[0].contains("chr1:1-100(-)"));
    let seq: String = lines[1..].concat();
    assert_eq!(seq.len(), 100);

    let fwd: Vec<u8> = read_fasta_seq(&fixture("ref_2000.fa"))[..100]
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let expected: Vec<u8> = pgr::libs::nt::rev_comp(&fwd).collect();
    assert_eq!(seq.as_bytes(), expected);
}

#[test]
fn test_pbit_range_multi_ranges() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "pbit",
            "range",
            out_pbit.to_str().unwrap(),
            "chr1:1-10",
            "chr1:11-20",
        ])
        .run();

    // Two ranges -> two FASTA entries.
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 2);
}

#[test]
fn test_pbit_some_basic() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let list_file = temp.path().join("list.txt");
    let out_fa = temp.path().join("out.fa");

    fs::write(&list_file, "chr1\n").unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000_2contig_identical.fa").to_str().unwrap(),
            "-i",
            fixture("ref_2000_2contig_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "some",
            out_pbit.to_str().unwrap(),
            list_file.to_str().unwrap(),
            "-o",
            out_fa.to_str().unwrap(),
        ])
        .run();

    // Read the output file (not stdout -- output goes to the file).
    let content = fs::read_to_string(&out_fa).unwrap();
    let headers: Vec<&str> = content.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 1);
    assert!(headers[0].contains("chr1"));
}

#[test]
fn test_pbit_some_invert() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let list_file = temp.path().join("list.txt");

    fs::write(&list_file, "chr1\n").unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000_2contig_identical.fa").to_str().unwrap(),
            "-i",
            fixture("ref_2000_2contig_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "pbit",
            "some",
            out_pbit.to_str().unwrap(),
            list_file.to_str().unwrap(),
            "-i",
            "-o",
            "stdout",
        ])
        .run();

    // Invert: should only contain chr2 (not chr1).
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 1);
    assert!(headers[0].contains("chr2"));
}

#[test]
fn test_pbit_create_custom_params() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
            "-s",
            "8192",
            "-k",
            "10",
            "-l",
            "15",
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap()])
        .run();

    assert!(stdout.contains("Segment size: 8192"));
    assert!(stdout.contains("K-mer length: 10"));
    assert!(stdout.contains("Min match length: 15"));
}

#[test]
fn test_pbit_custom_min_match_len_roundtrip() {
    // Regression test: custom -l must be stored in the header and used
    // correctly by the decompressor. Previously min_match_len was inferred
    // as kmer_len + 3, which mismatched the compressor's value for non-default
    // -l, causing decode errors.
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");
    let sample = fixture("sample_2000_snps100.fa");

    // Use non-default -k and -l where kmer_len + 3 != min_match_len.
    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            sample.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
            "-k",
            "10",
            "-l",
            "15",
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    let stem = sample.file_stem().unwrap().to_str().unwrap();
    let out_fa = out_dir.join(format!("{}.fa", stem));
    let expected = read_fasta_seq(&sample).to_ascii_uppercase();
    let extracted = read_fasta_seq(&out_fa);
    assert_eq!(extracted, expected);
}

#[test]
fn test_pbit_no_match_contig() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_1000.fa").to_str().unwrap(),
            "-i",
            fixture("unknown_1000_seed100.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // The sample should be registered (sample_count = 1) but with no contigs.
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap()])
        .run();
    assert!(stdout.contains("Samples: 1"));

    // stat --contigs should output nothing (no contigs for this sample).
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--contigs"])
        .run();
    assert!(stdout.trim().is_empty());
}

#[test]
fn test_pbit_multi_contig_reference() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("multi_500.fa").to_str().unwrap(),
            "-i",
            fixture("multi_500.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--refs"])
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 3);
    assert!(stdout.contains("chr1\t1"));
    assert!(stdout.contains("chr2\t1"));
    assert!(stdout.contains("chr3\t1"));
}

#[test]
fn test_pbit_multi_segment_contig() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");
    let sample = fixture("sample_5000_snp4090g.fa");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_5000.fa").to_str().unwrap(),
            "-i",
            sample.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    let stem = sample.file_stem().unwrap().to_str().unwrap();
    let out_fa = out_dir.join(format!("{}.fa", stem));
    let content = fs::read_to_string(&out_fa).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    let seq: String = lines[1..].concat();
    assert_eq!(seq.len(), 5000);
    let expected = read_fasta_seq(&sample).to_ascii_uppercase();
    assert_eq!(seq, expected);
}

#[test]
fn test_pbit_identical_samples_dedup() {
    let temp = TempDir::new().unwrap();
    let name_tsv = temp.path().join("names.tsv");
    let out_pbit = temp.path().join("out.pbit");
    let sample = fixture("sample_1000_identical.fa");

    // Use --name TSV to give two different sample names to the same file.
    fs::write(
        &name_tsv,
        format!("s1\t{}\ns2\t{}\n", sample.display(), sample.display()),
    )
    .unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_1000.fa").to_str().unwrap(),
            "--name",
            name_tsv.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // Both samples should be registered.
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--samples"])
        .run();
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2);
    assert!(lines.contains(&"s1"));
    assert!(lines.contains(&"s2"));

    // to-fa should produce two files.
    let out_dir = temp.path().join("outdir");
    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    let f1 = out_dir.join("s1.fa");
    let f2 = out_dir.join("s2.fa");
    assert!(f1.exists());
    assert!(f2.exists());

    // Both should have identical content (identical input -> identical output).
    let c1 = fs::read_to_string(&f1).unwrap();
    let c2 = fs::read_to_string(&f2).unwrap();
    assert_eq!(c1, c2);
}

#[test]
fn test_pbit_with_snp_sample() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");
    let sample = fixture("sample_2000_snps100.fa");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            sample.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    let stem = sample.file_stem().unwrap().to_str().unwrap();
    let out_fa = out_dir.join(format!("{}.fa", stem));
    let expected = read_fasta_seq(&sample).to_ascii_uppercase();
    let seq = read_fasta_seq(&out_fa);
    assert_eq!(seq, expected);
}

#[test]
fn test_pbit_create_with_name_tsv() {
    let temp = TempDir::new().unwrap();
    let name_tsv = temp.path().join("names.tsv");
    let out_pbit = temp.path().join("out.pbit");
    let sample = fixture("sample_1000_identical.fa");

    fs::write(&name_tsv, format!("custom_name\t{}\n", sample.display())).unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_1000.fa").to_str().unwrap(),
            "--name",
            name_tsv.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--samples"])
        .run();
    assert_eq!(stdout.trim(), "custom_name");
}

#[test]
fn test_pbit_to_fa_single_sample() {
    let temp = TempDir::new().unwrap();
    let name_tsv = temp.path().join("names.tsv");
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");
    let sample = fixture("sample_1000_identical.fa");

    fs::write(
        &name_tsv,
        format!("s1\t{}\ns2\t{}\n", sample.display(), sample.display()),
    )
    .unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_1000.fa").to_str().unwrap(),
            "--name",
            name_tsv.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "-s",
            "s2",
        ])
        .run();

    // Only s2.fa should exist in the output directory.
    assert!(out_dir.join("s2.fa").exists());
    assert!(!out_dir.join("s1.fa").exists());
}

#[test]
fn test_pbit_pseudocat_roundtrip() {
    // Smoke test with a real-world FASTA file.
    let temp = TempDir::new().unwrap();
    let ref_fa = std::path::PathBuf::from("tests/pgr/pseudocat.fa");
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
            "-i",
            ref_fa.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // Verify stat output.
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap()])
        .run();
    assert!(stdout.contains("Reference contigs: 1"));
    assert!(stdout.contains("Samples: 1"));
    // pseudocat.fa is ~18803 bp -> 5 segments (4*4096 + remainder).
    assert!(stdout.contains("Reference groups: 5"));

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    let out_fa = out_dir.join("pseudocat.fa");
    assert!(out_fa.exists());

    // Compare extracted sequence length with original.
    let original = fs::read_to_string(&ref_fa).unwrap();
    let orig_seq: String = original.lines().filter(|l| !l.starts_with('>')).collect();
    let extracted = fs::read_to_string(&out_fa).unwrap();
    let ext_seq: String = extracted.lines().filter(|l| !l.starts_with('>')).collect();
    assert_eq!(ext_seq.len(), orig_seq.len());
}

#[test]
fn test_pbit_append() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let name_tsv = temp.path().join("names.tsv");
    let s2_fa = temp.path().join("s2.fa");
    let sample = fixture("sample_2000_identical.fa");

    fs::write(&name_tsv, format!("s1\t{}\n", sample.display())).unwrap();
    fs::copy(&sample, &s2_fa).unwrap();

    // Create with 1 sample named s1.
    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "--name",
            name_tsv.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // Append a second sample named s2.
    PgrCmd::new()
        .args(&[
            "pbit",
            "append",
            out_pbit.to_str().unwrap(),
            "-i",
            s2_fa.to_str().unwrap(),
        ])
        .run();

    // Both samples should be registered.
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--samples"])
        .run();
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2);
    assert!(lines.contains(&"s1"));
    assert!(lines.contains(&"s2"));
}

#[test]
fn test_pbit_append_overwrite() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let new_pbit = temp.path().join("new.pbit");
    let name_tsv = temp.path().join("names.tsv");
    let s2_fa = temp.path().join("s2.fa");
    let sample = fixture("sample_2000_identical.fa");

    fs::write(&name_tsv, format!("s1\t{}\n", sample.display())).unwrap();
    fs::copy(&sample, &s2_fa).unwrap();

    // Create with 1 sample named s1.
    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "--name",
            name_tsv.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // Append to a new archive via -o.
    PgrCmd::new()
        .args(&[
            "pbit",
            "append",
            out_pbit.to_str().unwrap(),
            "-i",
            s2_fa.to_str().unwrap(),
            "-o",
            new_pbit.to_str().unwrap(),
        ])
        .run();

    // Original archive unchanged (still 1 sample).
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--samples"])
        .run();
    assert_eq!(stdout.lines().count(), 1);

    // New archive has 2 samples.
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", new_pbit.to_str().unwrap(), "--samples"])
        .run();
    assert_eq!(stdout.lines().count(), 2);
}

#[test]
fn test_pbit_append_in_place() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let name_tsv = temp.path().join("names.tsv");
    let s2_fa = temp.path().join("s2.fa");

    fs::write(
        &name_tsv,
        format!("s1\t{}\n", fixture("sample_2000_snps100.fa").display()),
    )
    .unwrap();
    fs::copy(fixture("sample_2000_snps200.fa"), &s2_fa).unwrap();

    // Create with 1 sample named s1.
    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "--name",
            name_tsv.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // Append in-place (no -o).
    PgrCmd::new()
        .args(&[
            "pbit",
            "append",
            out_pbit.to_str().unwrap(),
            "-i",
            s2_fa.to_str().unwrap(),
        ])
        .run();

    // Both samples extractable with correct content.
    let out_dir = temp.path().join("outdir");
    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    // s1 content preserved.
    let c1 = fs::read_to_string(out_dir.join("s1.fa")).unwrap();
    let seq1: String = c1.lines().filter(|l| !l.starts_with('>')).collect();
    let expected1 = read_fasta_seq(&fixture("sample_2000_snps100.fa")).to_ascii_uppercase();
    assert_eq!(seq1, expected1);

    // s2 content correct.
    let c2 = fs::read_to_string(out_dir.join("s2.fa")).unwrap();
    let seq2: String = c2.lines().filter(|l| !l.starts_with('>')).collect();
    let expected2 = read_fasta_seq(&fixture("sample_2000_snps200.fa")).to_ascii_uppercase();
    assert_eq!(seq2, expected2);
}

#[test]
fn test_pbit_range_multicontig() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("multi_500.fa").to_str().unwrap(),
            "-i",
            fixture("multi_500.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // Extract multiple contig ranges in one call.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "pbit",
            "range",
            out_pbit.to_str().unwrap(),
            "chr1:1-50",
            "chr2:1-50",
            "chr3",
        ])
        .run();

    // 3 ranges -> 3 FASTA entries.
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 3);
    assert!(headers.iter().any(|h| h.contains("chr1:1-50")));
    assert!(headers.iter().any(|h| h.contains("chr2:1-50")));
    assert!(headers.iter().any(|h| h.contains("chr3")));
}

#[test]
fn test_pbit_empty_contig() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_empty_1000.fa").to_str().unwrap(),
            "-i",
            fixture("ref_empty_1000.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // stat should show 1 sample without panicking.
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap()])
        .run();
    assert!(stdout.contains("Samples: 1"));

    // stat --contigs should list both contigs.
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--contigs"])
        .run();
    assert!(stdout.contains("chr1"));
    assert!(stdout.contains("chr_empty"));
}

#[test]
fn test_pbit_mask_roundtrip() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("mask_ref.fa").to_str().unwrap(),
            "-i",
            fixture("mask_ref.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    // Extracted sequence is uppercase (mask not applied), content correct.
    let mask_path = fixture("mask_ref.fa");
    let stem = mask_path.file_stem().unwrap().to_str().unwrap();
    let content = fs::read_to_string(out_dir.join(format!("{}.fa", stem))).unwrap();
    let seq: String = content.lines().filter(|l| !l.starts_with('>')).collect();
    let expected = read_fasta_seq(&fixture("mask_ref.fa")).to_ascii_uppercase();
    assert_eq!(seq, expected);
}

#[test]
fn test_pbit_n_roundtrip() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");
    let sample = fixture("n_sample.fa");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("n_ref.fa").to_str().unwrap(),
            "-i",
            sample.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    // Extracted sequence preserves N (uppercase).
    let stem = sample.file_stem().unwrap().to_str().unwrap();
    let content = fs::read_to_string(out_dir.join(format!("{}.fa", stem))).unwrap();
    let seq: String = content.lines().filter(|l| !l.starts_with('>')).collect();
    let expected = read_fasta_seq(&sample).to_ascii_uppercase();
    assert_eq!(seq, expected);
}

#[test]
fn test_pbit_random_roundtrip() {
    // Property test: random reference (multi-contig, with N) + random SNP
    // samples -> create -> to-fa -> compare. Repeated 5 times.
    for seed in 0..5u64 {
        let temp = TempDir::new().unwrap();
        let out_pbit = temp.path().join("out.pbit");
        let out_dir = temp.path().join("outdir");
        let ref_fa = fixture(&format!("prop_ref_{}.fa", seed));
        let sample_fa = fixture(&format!("prop_sample_{}.fa", seed));

        PgrCmd::new()
            .args(&[
                "pbit",
                "create",
                "-r",
                ref_fa.to_str().unwrap(),
                "-i",
                sample_fa.to_str().unwrap(),
                "-o",
                out_pbit.to_str().unwrap(),
            ])
            .run();

        PgrCmd::new()
            .args(&[
                "pbit",
                "to-fa",
                out_pbit.to_str().unwrap(),
                "-o",
                out_dir.to_str().unwrap(),
            ])
            .run();

        // Verify each sample contig matches (uppercase).
        let sample_records = read_fasta_records(&sample_fa);
        let out_fa = out_dir.join(format!("prop_sample_{}.fa", seed));
        let got_records = read_fasta_records(&out_fa);
        assert_eq!(got_records.len(), sample_records.len());
        for (i, (name, seq)) in got_records.iter().enumerate() {
            assert_eq!(*name, sample_records[i].0);
            let expected: String = sample_records[i]
                .1
                .chars()
                .map(|c| c.to_ascii_uppercase())
                .collect();
            assert_eq!(seq.len(), expected.len());
            assert_eq!(*seq, expected);
        }
    }
}

#[test]
fn test_pbit_large_contig_segment_boundary() {
    // Reference contig length = segment_size * 3 + 1 (spans 4 segments).
    // Default segment_size = 4096 -> contig length = 12289.
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let sample = fixture("sample_12289_snps100.fa");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_12289.fa").to_str().unwrap(),
            "-i",
            sample.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // Extract ranges at segment boundaries (1-based inclusive).
    // Segment 0: [1, 4096], Segment 1: [4097, 8192], etc.
    let boundary_ranges = ["chr1:4095-4098", "chr1:8191-8194", "chr1:12288-12289"];

    let mut v = vec!["pbit", "range", out_pbit.to_str().unwrap()];
    v.extend_from_slice(&boundary_ranges);
    let (stdout, _) = PgrCmd::new().args(&v).run();

    // 3 ranges -> 3 FASTA entries.
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 3);

    // Verify each slice matches the original sample.
    let sample_seq = read_fasta_seq(&sample);
    let lines: Vec<&str> = stdout.lines().collect();
    let mut seq_idx = 0;
    for range in &boundary_ranges {
        // Parse start-end from range string.
        let coords: &str = range.split(':').nth(1).unwrap();
        let mut parts = coords.split('-');
        let start: usize = parts.next().unwrap().parse().unwrap();
        let end: usize = parts.next().unwrap().parse().unwrap();

        // Skip header line, then collect sequence lines until next '>'.
        seq_idx += 1;
        let mut seq = String::new();
        while seq_idx < lines.len() && !lines[seq_idx].starts_with('>') {
            seq.push_str(lines[seq_idx]);
            seq_idx += 1;
        }
        // 1-based inclusive -> 0-based half-open.
        let expected: String = sample_seq
            .chars()
            .skip(start - 1)
            .take(end - start + 1)
            .map(|c| c.to_ascii_uppercase())
            .collect();
        assert_eq!(seq.len(), end - start + 1);
        assert_eq!(seq, expected);
    }
}

#[test]
fn test_pbit_range_invalid_range_warns() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let assert = PgrCmd::new()
        .args(&["pbit", "range", out_pbit.to_str().unwrap(), "chr1:abc-def"])
        .assert()
        .success();
    let output = assert.get_output();
    let stdout = String::from_utf8(output.stdout.clone()).unwrap();
    let stderr = String::from_utf8(output.stderr.clone()).unwrap();

    assert!(stdout.trim().is_empty());
    assert!(stderr.contains("invalid range format"));
}

#[test]
fn test_pbit_create_invalid_name_tsv_line_number() {
    let temp = TempDir::new().unwrap();
    let name_tsv = temp.path().join("names.tsv");
    let out_pbit = temp.path().join("out.pbit");

    // Line 1 is a comment (skipped); line 2 is missing the FASTA path.
    fs::write(&name_tsv, "# sample list\nsample1\n").unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_1000.fa").to_str().unwrap(),
            "--name",
            name_tsv.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run_fail();

    let combined = format!("{}{}", stdout, stderr);
    assert!(combined.contains("line 2"));
    assert!(combined.contains("missing FASTA path"));
}

#[test]
fn test_pbit_range_nonexistent_contig_warns() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let assert = PgrCmd::new()
        .args(&[
            "pbit",
            "range",
            out_pbit.to_str().unwrap(),
            "missing_contig",
        ])
        .assert()
        .success();
    let output = assert.get_output();
    let stdout = String::from_utf8(output.stdout.clone()).unwrap();
    let stderr = String::from_utf8(output.stderr.clone()).unwrap();

    assert!(stdout.trim().is_empty());
    assert!(
        stderr.contains("missing_contig") || stderr.contains("not found"),
        "expected warning about missing contig, got stderr: {}",
        stderr
    );
}

#[test]
fn test_pbit_empty_sample_fasta() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_1000.fa").to_str().unwrap(),
            "-i",
            fixture("empty.fa").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap()])
        .run();
    assert!(stdout.contains("Samples: 1"));

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--contigs"])
        .run();
    assert!(stdout.trim().is_empty());

    let out_dir = temp.path().join("outdir");
    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    // The command must not panic. The output file may be absent or empty.
    let out_file = out_dir.join("empty.fa");
    if out_file.exists() {
        let content = fs::read_to_string(&out_file).unwrap();
        assert!(
            content.trim().is_empty(),
            "expected empty output file for sample with no contigs, got: {}",
            content
        );
    }
}
