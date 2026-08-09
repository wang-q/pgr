#[macro_use]
#[path = "common/mod.rs"]
mod common;

use std::io::Write;
use std::path::Path;

const HIST_FILE_LEN: u64 = 28 + 32767 * 8;

fn write_fa(path: &Path) {
    std::fs::write(
        path,
        ">chr1\nACGTACGTACGTACGTACGTACGTACGT\n>chr2\nTTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAA\n",
    )
    .unwrap();
}

fn write_fq(path: &Path) {
    std::fs::write(
        path,
        "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@r2\nTTTTGGGGCCCCAAAA\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
}

#[test]
fn command_kmer_help() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd.arg("kmer").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.contains("Analyzes k-mer counts, profiles"));
    assert!(stdout.contains("table"));
    assert!(stdout.contains("profile"));
    assert!(stdout.contains("hist"));
    assert!(stdout.contains("gc"));
    assert!(stdout.contains("qhist"));
    assert!(stdout.contains("qcheck"));
    assert!(stdout.contains("gsize"));
    Ok(())
}

#[test]
fn command_kmer_table_end_to_end() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    let pkt = temp.path().join("t.pkt");
    write_fa(&fa);

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("14 unique 8-mers"), "stderr: {stderr}");
    assert!(pkt.exists());

    // The table can be reused: histogram from -t matches the sequence path.
    let h1 = temp.path().join("h1.hist");
    let h2 = temp.path().join("h2.hist");
    common::PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            h1.to_str().unwrap(),
        ])
        .run();
    common::PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            "-t",
            pkt.to_str().unwrap(),
            "-o",
            h2.to_str().unwrap(),
        ])
        .run();
    assert_eq!(h1.metadata()?.len(), HIST_FILE_LEN);
    assert_eq!(
        std::fs::read(&h1)?,
        std::fs::read(&h2)?,
        "hist from table must match hist from sequences"
    );
    Ok(())
}

#[test]
fn command_kmer_table_reads_gzip() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    write_fa(&fa);
    let plain = temp.path().join("plain.pkt");
    common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            plain.to_str().unwrap(),
        ])
        .run();

    // Gzip the same input: the resulting table must be identical.
    let gz = temp.path().join("in.fa.gz");
    let f = std::fs::File::create(&gz)?;
    let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
    enc.write_all(&std::fs::read(&fa)?)?;
    enc.finish()?;
    let gz_pkt = temp.path().join("gz.pkt");
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            gz.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            gz_pkt.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("14 unique 8-mers"), "stderr: {stderr}");
    assert_eq!(
        std::fs::read(&plain)?,
        std::fs::read(&gz_pkt)?,
        "gzipped input must produce the same table"
    );
    Ok(())
}

#[test]
fn command_kmer_profile_self_and_relative() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    let pkt = temp.path().join("t.pkt");
    let self_pkp = temp.path().join("self.pkp");
    let rel_pkp = temp.path().join("rel.pkp");
    write_fa(&fa);

    common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .run();

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "profile",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            self_pkp.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("2 profiles"), "stderr: {stderr}");
    assert_eq!(&std::fs::read(&self_pkp)?[0..4], b"PKPP");

    // Relative profile reuses the table; k is read from the table when the
    // command line omits --kmer.
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "profile",
            fa.to_str().unwrap(),
            "-t",
            pkt.to_str().unwrap(),
            "-o",
            rel_pkp.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("2 profiles"), "stderr: {stderr}");
    assert_eq!(&std::fs::read(&rel_pkp)?[0..4], b"PKPP");
    Ok(())
}

#[test]
fn command_kmer_reads_fastq_and_stdin() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fq = temp.path().join("in.fq");
    let pkt = temp.path().join("fq.pkt");
    write_fq(&fq);

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fq.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("8 unique 8-mers"), "stderr: {stderr}");

    let fa = temp.path().join("in.fa");
    write_fa(&fa);
    let stdin_pkt = temp.path().join("stdin.pkt");
    let input = std::fs::read_to_string(&fa)?;
    let (_, stderr) = common::PgrCmd::new()
        .stdin(input)
        .args(&[
            "kmer",
            "table",
            "stdin",
            "-k",
            "8",
            "-o",
            stdin_pkt.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("14 unique 8-mers"), "stderr: {stderr}");
    Ok(())
}

#[test]
fn command_kmer_argument_validation() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    let pkt = temp.path().join("t.pkt");
    write_fa(&fa);
    common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .run();

    // --kmer mismatching the table must fail.
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            "-t",
            pkt.to_str().unwrap(),
            "-k",
            "10",
            "-o",
            temp.path().join("x.hist").to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("does not match table k"),
        "stderr: {stderr}"
    );

    // No --kmer and no --table must fail.
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            fa.to_str().unwrap(),
            "-o",
            temp.path().join("y.hist").to_str().unwrap(),
        ])
        .run_fail();
    assert!(stderr.contains("--kmer is required"), "stderr: {stderr}");
    Ok(())
}

#[test]
fn command_kmer_gc_end_to_end() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    let pkt = temp.path().join("t.pkt");
    let kgc = temp.path().join("out.kgc");
    write_fa(&fa);

    // Sequence path: peak count 2 -> xmax = 2.1*2 = 4 -> 8 rows x 4 bins.
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "gc",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            kgc.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("peak 2"), "stderr: {stderr}");
    let text = std::fs::read_to_string(&kgc)?;
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "GCP\tKF\tCount");
    assert_eq!(lines.len(), 1 + 8 * 4);

    // Table path must produce the same matrix.
    common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .run();
    let kgc2 = temp.path().join("out2.kgc");
    common::PgrCmd::new()
        .args(&[
            "kmer",
            "gc",
            "-t",
            pkt.to_str().unwrap(),
            "-o",
            kgc2.to_str().unwrap(),
        ])
        .run();
    assert_eq!(
        std::fs::read(&kgc)?,
        std::fs::read(&kgc2)?,
        "gc from table must match gc from sequences"
    );
    Ok(())
}

#[test]
fn command_kmer_gc_tex_renders_heatmap() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    write_fa(&fa);
    let tex = temp.path().join("out.tex");
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "gc",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "--tex",
            "-o",
            tex.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("peak 2"), "stderr: {stderr}");
    let text = std::fs::read_to_string(&tex)?;
    assert!(text.contains("k-mer coverage"));
    assert!(text.contains("GC content"));
    assert!(text.contains("\\addplot"));
    assert!(text.contains("x y  C"));
    Ok(())
}

#[test]
fn command_kmer_qhist_end_to_end() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fq = temp.path().join("in.fq");
    write_fq(&fq);

    // Default threshold (detected Phred+33 + 5 = 38) and the explicit
    // threshold must produce the same quorum-format histogram.
    let out1 = temp.path().join("a.qhist");
    let out2 = temp.path().join("b.qhist");
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "qhist",
            fq.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            out1.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("threshold 38"), "stderr: {stderr}");
    common::PgrCmd::new()
        .args(&[
            "kmer",
            "qhist",
            fq.to_str().unwrap(),
            "-k",
            "8",
            "-q",
            "38",
            "-o",
            out2.to_str().unwrap(),
        ])
        .run();
    let text = std::fs::read_to_string(&out1)?;
    let lines: Vec<&str> = text.lines().collect();
    assert!(lines.iter().all(|l| {
        let cols: Vec<&str> = l.split_whitespace().collect();
        cols.len() == 3 && cols[0].parse::<u64>().is_ok()
    }));
    assert_eq!(
        std::fs::read(&out1)?,
        std::fs::read(&out2)?,
        "default and explicit threshold must match"
    );
    Ok(())
}

#[test]
fn command_kmer_qhist_rejects_fasta() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    write_fa(&fa);
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "qhist",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            temp.path().join("x.qhist").to_str().unwrap(),
        ])
        .run_fail();
    assert!(stderr.contains("requires FASTQ"), "stderr: {stderr}");
    Ok(())
}

#[test]
fn command_kmer_qcheck_end_to_end() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fq = temp.path().join("in.fq");
    let seq = "ACGTACGTACGTACGTACGTACGTACGT";
    let qual = "I".repeat(28);
    let mut fastq = String::new();
    for i in 1..=3 {
        fastq.push_str(&format!("@r{i}\n{seq}\n+\n{qual}\n"));
    }
    // One read with a single-base substitution at position 10.
    let mut bad = seq.to_string();
    bad.replace_range(10..11, "C");
    fastq.push_str(&format!("@bad\n{bad}\n+\n{qual}\n"));
    std::fs::write(&fq, fastq)?;

    let kept = temp.path().join("kept.fq");
    let discarded = temp.path().join("discarded.fq");
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "qcheck",
            fq.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            kept.to_str().unwrap(),
            "--discard-file",
            discarded.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("Kept 3 reads, flagged 1"),
        "stderr: {stderr}"
    );
    let kept_text = std::fs::read_to_string(&kept)?;
    let discarded_text = std::fs::read_to_string(&discarded)?;
    assert_eq!(kept_text.matches("@r").count(), 3);
    assert!(!kept_text.contains("@bad"));
    assert!(discarded_text.contains("@bad"));
    assert_eq!(discarded_text.matches('@').count(), 1);
    Ok(())
}

#[test]
fn command_kmer_gsize_estimates_coverage_and_size() -> anyhow::Result<()> {
    // Synthetic 30x coverage: 300 reads of 100 bp from a 1 kb genome.
    let temp = tempfile::TempDir::new()?;
    let fq = temp.path().join("reads.fq");
    let genome = random_dna(1000, 42);
    let mut fastq = String::new();
    for i in 1..=300 {
        let start = ((i as u64 * 2654435761) % 901) as usize;
        let read: String = genome[start..start + 100]
            .iter()
            .map(|&b| b as char)
            .collect();
        fastq.push_str(&format!("@r{i}\n{read}\n+\n{}\n", "I".repeat(100)));
    }
    std::fs::write(&fq, fastq)?;

    let (stdout, _) = common::PgrCmd::new()
        .args(&["kmer", "gsize", fq.to_str().unwrap(), "-k", "17"])
        .run();
    let mut peak = 0u64;
    let mut size = 0f64;
    for line in stdout.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() == 2 {
            match cols[0] {
                "peak_coverage" => peak = cols[1].parse()?,
                "genome_size" => size = cols[1].parse()?,
                _ => {}
            }
        }
    }
    assert!(
        (8..=60).contains(&peak),
        "peak coverage {peak} far from 30x"
    );
    assert!(
        (500.0..=2000.0).contains(&size),
        "genome size {size} far from 1000 bp"
    );
    Ok(())
}

#[test]
fn command_kmer_gsize_model_fit() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fq = temp.path().join("reads.fq");
    let genome = random_dna(1000, 42);
    let mut fastq = String::new();
    for i in 1..=600 {
        let start = ((i as u64 * 2654435761) % 901) as usize;
        let read: String = genome[start..start + 100]
            .iter()
            .map(|&b| b as char)
            .collect();
        fastq.push_str(&format!("@r{i}\n{read}\n+\n{}\n", "I".repeat(100)));
    }
    std::fs::write(&fq, fastq)?;

    let (stdout, _) = common::PgrCmd::new()
        .args(&["kmer", "gsize", fq.to_str().unwrap(), "-k", "17", "--model"])
        .run();
    let mut kmercov = 0f64;
    let mut size = 0f64;
    for line in stdout.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() == 2 {
            match cols[0] {
                "kmercov" => kmercov = cols[1].parse()?,
                "genome_size" => size = cols[1].parse()?,
                _ => {}
            }
        }
    }
    assert!(
        (30.0..=90.0).contains(&kmercov),
        "kmercov {kmercov} far from 60x"
    );
    // The two-component Poisson fit overestimates the genome size on real
    // reads (coverage heterogeneity from read ends); assert order of
    // magnitude only. Exact recovery is covered by the noiseless synthetic
    // unit test in libs/kmer/hist.rs.
    assert!(
        (300.0..=3000.0).contains(&size),
        "genome size {size} far from 1000 bp"
    );
    Ok(())
}

#[test]
fn command_kmer_gsize_model_writes_genomescope_outputs() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fq = temp.path().join("reads.fq");
    let genome = random_dna(1000, 42);
    let mut fastq = String::new();
    for i in 1..=600 {
        let start = ((i as u64 * 2654435761) % 901) as usize;
        let read: String = genome[start..start + 100]
            .iter()
            .map(|&b| b as char)
            .collect();
        fastq.push_str(&format!("@r{i}\n{read}\n+\n{}\n", "I".repeat(100)));
    }
    std::fs::write(&fq, fastq)?;

    let outdir = temp.path().join("gs");
    let (stdout, _) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "gsize",
            fq.to_str().unwrap(),
            "-k",
            "17",
            "--model",
            "-o",
            outdir.to_str().unwrap(),
        ])
        .run();
    let summary = outdir.join("summary.txt");
    let model = outdir.join("model.txt");
    assert!(summary.is_file(), "summary.txt must be written");
    assert!(model.is_file(), "model.txt must be written");
    let summary_text = std::fs::read_to_string(&summary)?;
    // anchr 2_fastk deletes the first six lines (`sed '1,6 d'`), so the
    // header must be the GenomeScope 5-liner plus a blank line.
    let head: Vec<&str> = summary_text.lines().take(6).collect();
    assert_eq!(head[0], "GenomeScope version 2.0");
    assert!(head[1].starts_with("input file = "));
    assert!(head[2].starts_with("output directory = "));
    assert_eq!(head[3], "p = 1");
    assert_eq!(head[4], "k = 17");
    assert_eq!(head[5], "");
    assert!(summary_text.contains("Genome Haploid Length"));
    assert!(summary_text.contains("Read Error Rate"));
    let model_text = std::fs::read_to_string(&model)?;
    // anchr 2_fastk parses `grep '^kmercov' model.txt | cut -f 2`.
    let kcov_line = model_text
        .lines()
        .find(|l| l.starts_with("kmercov"))
        .expect("model.txt must have a kmercov line");
    assert!(kcov_line.split_whitespace().nth(1).is_some());
    assert!(stdout.contains("kmercov\t"));
    assert!(stdout.contains("converged\ttrue"));
    Ok(())
}

/// Deterministic pseudo-random DNA (LCG).
fn random_dna(len: usize, seed: u64) -> Vec<u8> {
    let bases = *b"ACGT";
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[(x >> 33) as usize & 3]
        })
        .collect()
}
