#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

/// Parses a FASTA into (header, sequence) pairs.
fn parse_fa(data: &[u8]) -> Vec<(String, String)> {
    let mut recs = Vec::new();
    let mut cur: Option<(String, String)> = None;
    for line in std::str::from_utf8(data).unwrap().lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(c) = cur.take() {
                recs.push(c);
            }
            cur = Some((rest.to_string(), String::new()));
        } else if let Some(c) = cur.as_mut() {
            c.1.push_str(line);
        }
    }
    if let Some(c) = cur {
        recs.push(c);
    }
    recs
}

/// A linear genome compresses into a single maximal unitig.
#[test]
fn command_asm_unitig_linear() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    // Random 60 bp (all 30 k-mers unique -> linear, not cyclic).
    let seq = "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGA";
    // 4 identical reads so every k-mer is solid (count 4 >= seed 3).
    let fa = format!(">r1\n{seq}\n>r2\n{seq}\n>r3\n{seq}\n>r4\n{seq}\n");
    std::fs::write(&infile, fa).unwrap();
    PgrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    assert_eq!(recs.len(), 1, "expected one unitig, got {}", recs.len());
    let rc: String = seq
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => c,
        })
        .collect();
    assert!(recs[0].1 == seq || recs[0].1 == rc, "got {}", recs[0].1);
}

/// A bubble (two parallel paths) stays split: each branch is its own unitig
/// instead of being merged into one representative path.
#[test]
fn command_asm_unitig_keeps_branches() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    let prefix = "ACGT".repeat(15); // 60 bp shared prefix
    let window_a = "ACGT".repeat(8); // 32 bp path A
    let window_b = "ACGA".repeat(8); // 32 bp path B (variant)
    let suffix = "TGCA".repeat(15); // 60 bp shared suffix
    let path_a = format!("{prefix}{window_a}{suffix}");
    let path_b = format!("{prefix}{window_b}{suffix}");
    // 10 reads per path: every k-mer has count 10 (solid at seed 3).
    let mut fa = String::new();
    for i in 0..10 {
        fa.push_str(&format!(">a{i}\n{path_a}\n>b{i}\n{path_b}\n"));
    }
    std::fs::write(&infile, fa).unwrap();
    PgrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    // Prefix, both variant windows, and suffix form separate unitigs; a
    // bubble-popping assembler would merge them into fewer, longer contigs.
    assert!(recs.len() >= 4, "expected >= 4 unitigs, got {}", recs.len());
    assert!(
        recs.iter().all(|(_, s)| s.len() < path_a.len()),
        "a unitig spans the whole bubble"
    );
}

/// Output is deterministic and non-empty on the Lambda dataset.
#[test]
fn command_asm_unitig_deterministic() {
    let out_dir = tempfile::tempdir().unwrap();
    let out1 = out_dir.path().join("u1.fa");
    let out2 = out_dir.path().join("u2.fa");
    for out in [&out1, &out2] {
        PgrCmd::new()
            .args(&[
                "asm",
                "unitig",
                "tests/bbtools/Lambda/R1.2k.fq.gz",
                "tests/bbtools/Lambda/R2.2k.fq.gz",
                "-o",
                out.to_str().unwrap(),
                "--kmer",
                "31",
            ])
            .assert()
            .success();
    }
    assert_eq!(std::fs::read(&out1).unwrap(), std::fs::read(&out2).unwrap());
    let recs = parse_fa(&std::fs::read(&out1).unwrap());
    assert!(!recs.is_empty());
}

/// Raising the solid threshold drops low-count k-mers (bcalm `-abundance-min`).
#[test]
fn command_asm_unitig_min_count_seed() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let seq = "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGA";
    // 2 identical reads: every k-mer has count 2 (not solid at the default
    // threshold of 3, solid at --min-count-seed 2).
    let fa = format!(">r1\n{seq}\n>r2\n{seq}\n");
    std::fs::write(&infile, fa).unwrap();
    let default_out = out_dir.path().join("default.fa");
    PgrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            default_out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    assert!(parse_fa(&std::fs::read(&default_out).unwrap()).is_empty());
    let strict_out = out_dir.path().join("strict.fa");
    PgrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            strict_out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
            "--min-count-seed",
            "2",
        ])
        .assert()
        .success();
    assert_eq!(parse_fa(&std::fs::read(&strict_out).unwrap()).len(), 1);
}

/// A branching graph emits GFA segments and (k-1)-overlap links.
#[test]
fn command_asm_unitig_gfa() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.gfa");
    let prefix = "ACGT".repeat(15); // 60 bp shared prefix
    let window_a = "ACGT".repeat(8); // 32 bp path A
    let window_b = "ACGA".repeat(8); // 32 bp path B (variant)
    let suffix = "TGCA".repeat(15); // 60 bp shared suffix
    let path_a = format!("{prefix}{window_a}{suffix}");
    let path_b = format!("{prefix}{window_b}{suffix}");
    let mut fa = String::new();
    for i in 0..10 {
        fa.push_str(&format!(">a{i}\n{path_a}\n>b{i}\n{path_b}\n"));
    }
    std::fs::write(&infile, fa).unwrap();
    PgrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
            "--gfa",
        ])
        .assert()
        .success();
    let gfa = std::fs::read_to_string(&out).unwrap();
    let mut headers = 0;
    let mut segments = 0;
    let mut links = 0;
    for line in gfa.lines() {
        match line.as_bytes().first() {
            Some(b'H') => headers += 1,
            Some(b'S') => segments += 1,
            Some(b'L') => {
                links += 1;
                assert!(line.ends_with("\t30M"), "overlap: {line}");
            }
            _ => {}
        }
    }
    assert_eq!(headers, 1);
    assert!(segments >= 4, "segments: {segments}");
    assert!(links >= 4, "links: {links}");
}

/// `--links` appends BCALM-style `L:` entries to FASTA headers.
#[test]
fn command_asm_unitig_links_header() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    let prefix = "ACGT".repeat(15);
    let window_a = "ACGT".repeat(8);
    let window_b = "ACGA".repeat(8);
    let suffix = "TGCA".repeat(15);
    let path_a = format!("{prefix}{window_a}{suffix}");
    let path_b = format!("{prefix}{window_b}{suffix}");
    let mut fa = String::new();
    for i in 0..10 {
        fa.push_str(&format!(">a{i}\n{path_a}\n>b{i}\n{path_b}\n"));
    }
    std::fs::write(&infile, fa).unwrap();
    PgrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
            "--links",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    assert!(recs.len() >= 4);
    let any_link = recs.iter().any(|(h, _)| {
        h.split_whitespace()
            .any(|f| f.starts_with("L:+:") || f.starts_with("L:-:"))
    });
    assert!(any_link, "no L: entries in headers");
}

/// A k-mer above the 128-base key limit must fail cleanly instead of
/// panicking in `Kmer::new().expect()` (zero-panic policy).
#[test]
fn command_asm_unitig_rejects_kmer_above_limit() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fa");
    PgrCmd::new()
        .args(&[
            "asm",
            "unitig",
            "tests/bbtools/Lambda/R1.2k.fq.gz",
            "tests/bbtools/Lambda/R2.2k.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "129",
        ])
        .assert()
        .failure();
}

/// `-o` must not overwrite an input file (the writer is opened before the
/// reads are consumed).
#[test]
fn command_asm_unitig_outfile_not_input() {
    let infile = "tests/bbtools/Lambda/R1.2k.fq.gz";
    PgrCmd::new()
        .args(&["asm", "unitig", infile, "-o", infile])
        .assert()
        .failure()
        .stderr(predicates::str::contains("is also an input file"));
}

/// A cyclic genome (periodic sequence) assembles into a circular unitig
/// flagged in the FASTA header.
#[test]
fn command_asm_unitig_circular() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    // Periodic 80 bp "genome": with k=31 the k-mer graph is a 4-kmer cycle.
    let genome = "ACGT".repeat(20);
    let mut fa = String::new();
    for i in 0..10 {
        let start = (i * 2) % 40;
        fa.push_str(&format!(">r{i}\n{}\n", &genome[start..start + 40]));
    }
    std::fs::write(&infile, fa).unwrap();
    PgrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    assert_eq!(
        recs.len(),
        1,
        "expected one circular unitig, got {}",
        recs.len()
    );
    assert!(recs[0].0.contains("circular"), "header: {}", recs[0].0);
}
