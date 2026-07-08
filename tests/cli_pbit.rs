#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use std::io::Write;
use tempfile::TempDir;

/// Write a FASTA file with the given (name, seq) records.
fn write_fasta(path: &std::path::Path, records: &[(&str, &str)]) {
    let mut f = fs::File::create(path).unwrap();
    for (name, seq) in records {
        writeln!(f, ">{}", name).unwrap();
        writeln!(f, "{}", seq).unwrap();
    }
}

/// Generate a deterministic random DNA sequence of the given length.
fn random_dna(len: usize, seed: u64) -> String {
    use rand::rngs::StdRng;
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = StdRng::seed_from_u64(seed);
    (0..len)
        .map(|_| match rng.random_range(0u8..4) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        })
        .collect()
}

/// Introduce SNPs at every 100th position.
fn introduce_snps(seq: &str, seed: u64) -> String {
    use rand::rngs::StdRng;
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = StdRng::seed_from_u64(seed);
    let mut out: Vec<char> = seq.chars().collect();
    for i in (0..out.len()).step_by(100) {
        let orig = out[i];
        let new = match orig {
            'A' => {
                let r = rng.random_range(0u8..3);
                if r == 0 {
                    'C'
                } else if r == 1 {
                    'G'
                } else {
                    'T'
                }
            }
            'C' => {
                let r = rng.random_range(0u8..3);
                if r == 0 {
                    'A'
                } else if r == 1 {
                    'G'
                } else {
                    'T'
                }
            }
            'G' => {
                let r = rng.random_range(0u8..3);
                if r == 0 {
                    'A'
                } else if r == 1 {
                    'C'
                } else {
                    'T'
                }
            }
            _ => {
                let r = rng.random_range(0u8..3);
                if r == 0 {
                    'A'
                } else if r == 1 {
                    'C'
                } else {
                    'G'
                }
            }
        };
        out[i] = new;
    }
    out.into_iter().collect()
}

#[test]
fn test_pbit_create_basic() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

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

    assert!(out_pbit.exists());
    assert!(fs::metadata(&out_pbit).unwrap().len() > 0);
}

#[test]
fn test_pbit_stat_overview() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

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
    let ref_fa = temp.path().join("ref.fa");
    let s1_fa = temp.path().join("s1.fa");
    let s2_fa = temp.path().join("s2.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&s1_fa, &[("chr1", &ref_seq)]);
    write_fasta(&s2_fa, &[("chr1", &ref_seq)]);

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
            "-i",
            s1_fa.to_str().unwrap(),
            "-i",
            s2_fa.to_str().unwrap(),
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
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(5000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

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

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--refs"])
        .run();

    // 5000 bp / 4096 segment_size → 2 segments for chr1.
    assert_eq!(stdout.trim(), "chr1\t2");
}

#[test]
fn test_pbit_stat_contigs() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(1000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

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

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap(), "--contigs"])
        .run();

    assert_eq!(stdout.trim(), "sample\tchr1");
}

#[test]
fn test_pbit_to_fa_roundtrip() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    let ref_seq = random_dna(2000, 42);
    let sample_seq = introduce_snps(&ref_seq, 100);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &sample_seq)]);

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

    let out_fa = out_dir.join("sample.fa");
    assert!(out_fa.exists());

    // Read the output FASTA and verify the sequence matches (uppercase).
    let content = fs::read_to_string(&out_fa).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines[0].starts_with(">chr1"));
    let extracted: String = lines[1..].concat();
    let expected: String = sample_seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    assert_eq!(extracted, expected);
}

#[test]
fn test_pbit_range_full_contig() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

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

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "range", out_pbit.to_str().unwrap(), "chr1"])
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines[0].starts_with(">sample chr1"));
    let seq: String = lines[1..].concat();
    assert_eq!(seq.len(), 2000);
}

#[test]
fn test_pbit_range_slice() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

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

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "range", out_pbit.to_str().unwrap(), "chr1:1-100"])
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines[0].starts_with(">sample chr1:1-100(+)"));
    let seq: String = lines[1..].concat();
    assert_eq!(seq.len(), 100);
    let expected: String = ref_seq[..100].to_uppercase();
    assert_eq!(seq, expected);
}

#[test]
fn test_pbit_range_neg_strand() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

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

    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "range", out_pbit.to_str().unwrap(), "chr1(-):1-100"])
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines[0].starts_with(">sample chr1:1-100(-)"));
    let seq: String = lines[1..].concat();
    assert_eq!(seq.len(), 100);

    // Compute the expected reverse complement of ref_seq[..100].
    let fwd: Vec<u8> = ref_seq[..100]
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let expected: Vec<u8> = pgr::libs::nt::rev_comp(&fwd).collect();
    assert_eq!(seq.as_bytes(), expected);
}

#[test]
fn test_pbit_range_multi_ranges() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

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

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "pbit",
            "range",
            out_pbit.to_str().unwrap(),
            "chr1:1-10",
            "chr1:11-20",
        ])
        .run();

    // Two ranges → two FASTA entries.
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 2);
}

#[test]
fn test_pbit_some_basic() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");
    let list_file = temp.path().join("list.txt");
    let out_fa = temp.path().join("out.fa");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq), ("chr2", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq), ("chr2", &ref_seq)]);

    fs::write(&list_file, "chr1\n").unwrap();

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
            "some",
            out_pbit.to_str().unwrap(),
            list_file.to_str().unwrap(),
            "-o",
            out_fa.to_str().unwrap(),
        ])
        .run();

    // Read the output file (not stdout — output goes to the file).
    let content = fs::read_to_string(&out_fa).unwrap();
    let headers: Vec<&str> = content.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 1);
    assert!(headers[0].contains("chr1"));
}

#[test]
fn test_pbit_some_invert() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");
    let list_file = temp.path().join("list.txt");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq), ("chr2", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq), ("chr2", &ref_seq)]);

    fs::write(&list_file, "chr1\n").unwrap();

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
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

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
}

#[test]
fn test_pbit_no_match_contig() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(1000, 42);
    let sample_seq = random_dna(1000, 100);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    // Sample has a contig name that does not match any reference contig.
    write_fasta(&sample_fa, &[("unknown_contig", &sample_seq)]);

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
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let seq1 = random_dna(500, 1);
    let seq2 = random_dna(500, 2);
    let seq3 = random_dna(500, 3);
    write_fasta(
        &ref_fa,
        &[("chr1", &seq1), ("chr2", &seq2), ("chr3", &seq3)],
    );
    write_fasta(
        &sample_fa,
        &[("chr1", &seq1), ("chr2", &seq2), ("chr3", &seq3)],
    );

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
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    // 5000 bp → 2 segments of 4096 + 904.
    let ref_seq = random_dna(5000, 42);
    let mut sample_seq = ref_seq.clone();
    // Introduce a SNP near the segment boundary.
    sample_seq.replace_range(4090..4091, "G");
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &sample_seq)]);

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

    let out_fa = out_dir.join("sample.fa");
    let content = fs::read_to_string(&out_fa).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    let seq: String = lines[1..].concat();
    assert_eq!(seq.len(), 5000);
    let expected: String = sample_seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    assert_eq!(seq, expected);
}

#[test]
fn test_pbit_identical_samples_dedup() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let name_tsv = temp.path().join("names.tsv");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(1000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

    // Use --name TSV to give two different sample names to the same file.
    fs::write(
        &name_tsv,
        format!("s1\t{}\ns2\t{}", sample_fa.display(), sample_fa.display()),
    )
    .unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
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

    // Both should have identical content (identical input → identical output).
    let c1 = fs::read_to_string(&f1).unwrap();
    let c2 = fs::read_to_string(&f2).unwrap();
    assert_eq!(c1, c2);
}

#[test]
fn test_pbit_with_snp_sample() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    let ref_seq = random_dna(2000, 42);
    let sample_seq = introduce_snps(&ref_seq, 100);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &sample_seq)]);

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

    let out_fa = out_dir.join("sample.fa");
    let content = fs::read_to_string(&out_fa).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    let seq: String = lines[1..].concat();
    let expected: String = sample_seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    assert_eq!(seq, expected);
}

#[test]
fn test_pbit_create_with_name_tsv() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let name_tsv = temp.path().join("names.tsv");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(1000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq)]);

    fs::write(&name_tsv, format!("custom_name\t{}", sample_fa.display())).unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
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
    let ref_fa = temp.path().join("ref.fa");
    let s1_fa = temp.path().join("s1.fa");
    let s2_fa = temp.path().join("s2.fa");
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    let ref_seq = random_dna(1000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&s1_fa, &[("chr1", &ref_seq)]);
    write_fasta(&s2_fa, &[("chr1", &ref_seq)]);

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
            "-i",
            s1_fa.to_str().unwrap(),
            "-i",
            s2_fa.to_str().unwrap(),
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
    // pseudocat.fa is ~18803 bp → 5 segments (4*4096 + remainder).
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
    let ref_fa = temp.path().join("ref.fa");
    let s1_fa = temp.path().join("s1.fa");
    let s2_fa = temp.path().join("s2.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&s1_fa, &[("chr1", &ref_seq)]);
    write_fasta(&s2_fa, &[("chr1", &ref_seq)]);

    // Create with 1 sample.
    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
            "-i",
            s1_fa.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // Append a second sample.
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
    let ref_fa = temp.path().join("ref.fa");
    let s1_fa = temp.path().join("s1.fa");
    let s2_fa = temp.path().join("s2.fa");
    let out_pbit = temp.path().join("out.pbit");
    let new_pbit = temp.path().join("new.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&s1_fa, &[("chr1", &ref_seq)]);
    write_fasta(&s2_fa, &[("chr1", &ref_seq)]);

    // Create with 1 sample.
    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
            "-i",
            s1_fa.to_str().unwrap(),
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
    let ref_fa = temp.path().join("ref.fa");
    let s1_fa = temp.path().join("s1.fa");
    let s2_fa = temp.path().join("s2.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    let s1_seq = introduce_snps(&ref_seq, 100);
    let s2_seq = introduce_snps(&ref_seq, 200);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&s1_fa, &[("chr1", &s1_seq)]);
    write_fasta(&s2_fa, &[("chr1", &s2_seq)]);

    // Create with 1 sample.
    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
            "-i",
            s1_fa.to_str().unwrap(),
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
    let expected1: String = s1_seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    assert_eq!(seq1, expected1);

    // s2 content correct.
    let c2 = fs::read_to_string(out_dir.join("s2.fa")).unwrap();
    let seq2: String = c2.lines().filter(|l| !l.starts_with('>')).collect();
    let expected2: String = s2_seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    assert_eq!(seq2, expected2);
}

#[test]
fn test_pbit_range_multicontig() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let seq1 = random_dna(500, 1);
    let seq2 = random_dna(500, 2);
    let seq3 = random_dna(500, 3);
    write_fasta(
        &ref_fa,
        &[("chr1", &seq1), ("chr2", &seq2), ("chr3", &seq3)],
    );
    write_fasta(
        &sample_fa,
        &[("chr1", &seq1), ("chr2", &seq2), ("chr3", &seq3)],
    );

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

    // 3 ranges → 3 FASTA entries.
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 3);
    assert!(headers.iter().any(|h| h.contains("chr1:1-50")));
    assert!(headers.iter().any(|h| h.contains("chr2:1-50")));
    assert!(headers.iter().any(|h| h.contains("chr3")));
}

#[test]
fn test_pbit_empty_contig() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(1000, 42);
    // Reference has chr1 (normal) + chr2 (empty).
    write_fasta(&ref_fa, &[("chr1", &ref_seq), ("chr2", "")]);
    write_fasta(&sample_fa, &[("chr1", &ref_seq), ("chr2", "")]);

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
    assert!(stdout.contains("chr2"));
}

#[test]
fn test_pbit_mask_roundtrip() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    // Mixed-case reference: lowercase regions are soft-masked.
    let ref_seq = "acgtACGTacgtACGTacgtACGTacgtACGT";
    write_fasta(&ref_fa, &[("chr1", ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", ref_seq)]);

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

    // Extracted sequence is uppercase (mask not applied), content correct.
    let content = fs::read_to_string(out_dir.join("sample.fa")).unwrap();
    let seq: String = content.lines().filter(|l| !l.starts_with('>')).collect();
    let expected: String = ref_seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    assert_eq!(seq, expected);
}

#[test]
fn test_pbit_n_roundtrip() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    // Reference and sample contain N runs.
    let ref_seq = "ACGTNNNNACGTACGTACGTNNNNACGTACGT";
    let sample_seq = "ACGTNNNNACGTACGTACGTNNNnACGTACGT";
    write_fasta(&ref_fa, &[("chr1", ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", sample_seq)]);

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

    // Extracted sequence preserves N (uppercase).
    let content = fs::read_to_string(out_dir.join("sample.fa")).unwrap();
    let seq: String = content.lines().filter(|l| !l.starts_with('>')).collect();
    let expected: String = sample_seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    assert_eq!(seq, expected);
}

/// Generate random DNA with occasional N runs.
fn random_dna_with_n(len: usize, seed: u64, n_freq: f64) -> String {
    use rand::rngs::StdRng;
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = StdRng::seed_from_u64(seed);
    (0..len)
        .map(|_| {
            if rng.random::<f64>() < n_freq {
                'N'
            } else {
                match rng.random_range(0u8..4) {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    _ => 'T',
                }
            }
        })
        .collect()
}

#[test]
fn test_pbit_random_roundtrip() {
    // Property test: random reference (multi-contig, with N) + random SNP
    // samples → create → to-fa → compare. Repeated 5 times.
    for seed in 0..5u64 {
        let temp = TempDir::new().unwrap();
        let ref_fa = temp.path().join("ref.fa");
        let sample_fa = temp.path().join("sample.fa");
        let out_pbit = temp.path().join("out.pbit");
        let out_dir = temp.path().join("outdir");

        // Random reference: 2-3 contigs, 500-1000 bp each, ~5% N.
        let n_contigs = 2 + (seed % 2) as usize;
        let mut ref_records: Vec<(String, String)> = Vec::new();
        for c in 0..n_contigs {
            let len = 500 + (seed as usize + c * 100) % 500;
            let seq = random_dna_with_n(len, seed * 100 + c as u64, 0.05);
            ref_records.push((format!("chr{}", c + 1), seq));
        }

        // Build sample FASTA with SNPs introduced per contig.
        let mut sample_records: Vec<(String, String)> = Vec::new();
        for (name, seq) in &ref_records {
            let mutated = introduce_snps(seq, seed * 7 + name.len() as u64);
            sample_records.push((name.clone(), mutated));
        }

        // Write reference and sample FASTA.
        let ref_pairs: Vec<(&str, &str)> = ref_records
            .iter()
            .map(|(n, s)| (n.as_str(), s.as_str()))
            .collect();
        write_fasta(&ref_fa, &ref_pairs);
        let sample_pairs: Vec<(&str, &str)> = sample_records
            .iter()
            .map(|(n, s)| (n.as_str(), s.as_str()))
            .collect();
        write_fasta(&sample_fa, &sample_pairs);

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
        let content = fs::read_to_string(out_dir.join("sample.fa")).unwrap();
        let mut current_name = String::new();
        let mut current_seq = String::new();
        let mut entries: Vec<(String, String)> = Vec::new();
        for line in content.lines() {
            if line.starts_with('>') {
                if !current_name.is_empty() {
                    entries.push((current_name.clone(), current_seq.clone()));
                    current_seq.clear();
                }
                current_name = line[1..].to_string();
            } else {
                current_seq.push_str(line);
            }
        }
        if !current_name.is_empty() {
            entries.push((current_name, current_seq));
        }

        assert_eq!(entries.len(), sample_records.len());
        for (i, (name, seq)) in entries.iter().enumerate() {
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
    // Default segment_size = 4096 → contig length = 12289.
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(12289, 42);
    let sample_seq = introduce_snps(&ref_seq, 100);
    write_fasta(&ref_fa, &[("chr1", &ref_seq)]);
    write_fasta(&sample_fa, &[("chr1", &sample_seq)]);

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

    // Extract ranges at segment boundaries (1-based inclusive).
    // Segment 0: [1, 4096], Segment 1: [4097, 8192], etc.
    let boundary_ranges = ["chr1:4095-4098", "chr1:8191-8194", "chr1:12288-12289"];

    let (stdout, _) = PgrCmd::new()
        .args({
            let mut v = vec!["pbit", "range", out_pbit.to_str().unwrap()];
            v.extend_from_slice(&boundary_ranges);
            v
        })
        .run();

    // 3 ranges → 3 FASTA entries.
    let headers: Vec<&str> = stdout.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers.len(), 3);

    // Verify each slice matches the original sample.
    let lines: Vec<&str> = stdout.lines().collect();
    let mut seq_idx = 1; // skip first header
    for range in &boundary_ranges {
        // Parse start-end from range string.
        let coords: &str = range.split(':').nth(1).unwrap();
        let mut parts = coords.split('-');
        let start: usize = parts.next().unwrap().parse().unwrap();
        let end: usize = parts.next().unwrap().parse().unwrap();

        // Skip header line.
        seq_idx += 1;
        let mut seq = String::new();
        while seq_idx < lines.len() && !lines[seq_idx].starts_with('>') {
            seq.push_str(lines[seq_idx]);
            seq_idx += 1;
        }
        // 1-based inclusive → 0-based half-open.
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
