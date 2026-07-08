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
