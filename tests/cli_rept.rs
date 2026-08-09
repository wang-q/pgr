#[macro_use]
#[path = "common/mod.rs"]
mod common;

#[test]
fn command_rept_e_kmer_help() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd.arg("rept").arg("e-kmer").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Identifies repeats against an external library"));
    Ok(())
}

#[test]
fn command_rept_s_kmer_help() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd.arg("rept").arg("s-kmer").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Identifies repetitive regions in a genome"));
    Ok(())
}

#[test]
fn command_rept_trf_help() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd.arg("rept").arg("trf").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Identifies tandem repeats in a genome"));
    Ok(())
}

#[test]
fn command_rept_e_align_help() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd.arg("rept").arg("e-align").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Identifies repeats against an external library (alignment)"));
    Ok(())
}

/// Deterministic pseudo-random DNA (same LCG as cli_align_pgi).
fn random_seq(len: usize, seed: u64) -> String {
    let bases = *b"ACGT";
    let mut x = seed;
    let mut s = String::with_capacity(len);
    for _ in 0..len {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        s.push(bases[(x >> 33) as usize & 3] as char);
    }
    s
}

/// e-align on synthetic data: a 400 bp "repeat" inserted twice in a random
/// genome should be reported as covered intervals.
#[test]
fn command_rept_e_align_end_to_end() -> anyhow::Result<()> {
    let repeat = random_seq(400, 7);
    let genome = format!(
        "{}{}{}{}{}",
        random_seq(400, 1),
        repeat,
        random_seq(200, 2),
        repeat,
        random_seq(400, 3)
    );

    let temp = tempfile::TempDir::new()?;
    let lib = temp.path().join("lib.fa");
    let genome_fa = temp.path().join("genome.fa");
    let out = temp.path().join("out.json");
    std::fs::write(&lib, format!(">rep1\n{}\n", repeat))?;
    // Dotted contig name: the runlist parser truncates these to the last
    // '.' segment, and e-align must restore the full name in the runlist.
    std::fs::write(&genome_fa, format!(">chr1.1\n{}\n", genome))?;

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "e-align",
            lib.to_str().unwrap(),
            genome_fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("==> Outputs"), "pipeline failed: {stderr}");

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let spans = json["chr1.1"].as_str().expect("chr1.1 spans missing");
    let covered: usize = spans
        .split(',')
        .filter_map(|s| {
            let mut it = s.split('-');
            let start: usize = it.next()?.parse().ok()?;
            let end: usize = it.next()?.parse().ok()?;
            Some(end - start + 1)
        })
        .sum();
    // Both 400 bp copies overlap the reported intervals; allow partial
    // boundary trimming, so require at least one full copy worth of coverage.
    assert!(covered >= 400, "expected repeat coverage, got: {spans}");
    Ok(())
}

/// rm on synthetic data: a 400 bp "repeat" inserted twice in a random genome
/// should be reported as covered intervals (needs makeblastdb, rmblastn and
/// trf in $PATH).
#[test]
fn command_rept_masker_end_to_end() -> anyhow::Result<()> {
    for tool in ["makeblastdb", "rmblastn", "trf"] {
        if which::which(tool).is_err() {
            eprintln!("skipping: {tool} not found");
            return Ok(());
        }
    }

    let repeat = random_seq(400, 7);
    let genome = format!(
        "{}{}{}{}{}",
        random_seq(400, 1),
        repeat,
        random_seq(200, 2),
        repeat,
        random_seq(400, 3)
    );

    let temp = tempfile::TempDir::new()?;
    let lib = temp.path().join("lib.fa");
    let genome_fa = temp.path().join("genome.fa");
    let out = temp.path().join("out.json");
    std::fs::write(&lib, format!(">rep1\n{}\n", repeat))?;
    // Dotted contig name: the runlist parser truncates these to the last
    // '.' segment, and rmblast must restore the full name in the runlist.
    std::fs::write(&genome_fa, format!(">chr1.1\n{}\n", genome))?;

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "masker",
            lib.to_str().unwrap(),
            genome_fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("==> Outputs"), "pipeline failed: {stderr}");

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let spans = json["chr1.1"].as_str().expect("chr1.1 spans missing");
    let covered: usize = spans
        .split(',')
        .filter_map(|s| {
            let mut it = s.split('-');
            let start: usize = it.next()?.parse().ok()?;
            let end: usize = it.next()?.parse().ok()?;
            Some(end - start + 1)
        })
        .sum();
    // Both 400 bp copies overlap the reported intervals; allow partial
    // boundary trimming, so require at least one full copy worth of coverage.
    assert!(covered >= 400, "expected repeat coverage, got: {spans}");
    Ok(())
}

/// s-kmer on a genome with a dotted contig name must report the full name
/// in the runlist (regression for runlist-name truncation). Runs with no
/// external tools; the native k-mer pipeline must work out of the box.
#[test]
fn command_rept_s_kmer_dotted_name() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let genome_fa = temp.path().join("genome.fa");
    let out = temp.path().join("out.json");
    let seq = random_seq(2000, 11);
    std::fs::write(
        &genome_fa,
        // One duplicated block guarantees repetitive k-mers for Profex.
        format!(">NC_000913.1\n{}{}\n", seq, seq),
    )?;

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "s-kmer",
            genome_fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("==> Outputs"), "pipeline failed: {stderr}");

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let spans = json["NC_000913.1"]
        .as_str()
        .expect("NC_000913.1 key missing (runlist truncated the name?)");
    // A 2000 bp block duplicated head-to-tail: every k-mer appears twice, and
    // the native extractor closes the tail run (the old Profex wrapper
    // dropped it), so the merged runlist spans the full 4000 bp.
    assert_eq!(spans, "1-4000", "unexpected s-kmer spans: {spans}");
    Ok(())
}

/// e-kmer on a perfect tandem duplication must report the full duplicated
/// interval, including the chromosome-tail run (no external tools needed).
#[test]
fn command_rept_e_kmer_tandem_coordinates() -> anyhow::Result<()> {
    let seq = random_seq(2000, 13);
    let temp = tempfile::TempDir::new()?;
    let lib = temp.path().join("lib.fa");
    let genome_fa = temp.path().join("genome.fa");
    let out = temp.path().join("out.json");
    std::fs::write(&lib, format!(">rep\n{}\n", seq))?;
    std::fs::write(&genome_fa, format!(">NC_000913.1\n{}{}\n", seq, seq))?;

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "e-kmer",
            lib.to_str().unwrap(),
            genome_fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("==> Outputs"), "pipeline failed: {stderr}");

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let spans = json["NC_000913.1"]
        .as_str()
        .expect("NC_000913.1 key missing");
    assert_eq!(spans, "1-4000", "e-kmer missed the tandem copy: {spans}");
    Ok(())
}

/// e-kmer end-to-end on the repo's small fixtures. Guards the native k-mer
/// pipeline and the dotted-name mapping.
#[test]
fn command_rept_e_kmer_end_to_end() -> anyhow::Result<()> {
    let lib = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/pgr/tncentral.fa.gz");
    let genome = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz");
    let temp = tempfile::TempDir::new()?;
    let out = temp.path().join("out.json");

    let (_, stderr) = common::PgrCmd::new()
        .args(&["rept", "e-kmer", lib, genome, "-o", out.to_str().unwrap()])
        .run();
    assert!(stderr.contains("==> Outputs"), "pipeline failed: {stderr}");
    assert!(
        !stderr.contains("FastK") && !stderr.contains("Profex"),
        "native pipeline must not invoke external tools: {stderr}"
    );

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let spans = json["NC_000913"].as_str().expect("NC_000913 key missing");
    assert!(spans.contains('-'), "expected intervals, got: {spans}");
    Ok(())
}

/// e-kmer `--keep-index` writes a `.pkt` table next to the library and
/// reuses it on the next run instead of rebuilding.
#[test]
fn command_rept_e_kmer_keep_index_cache() -> anyhow::Result<()> {
    let seq = random_seq(2000, 13);
    let temp = tempfile::TempDir::new()?;
    let lib = temp.path().join("lib.fa");
    let genome_fa = temp.path().join("genome.fa");
    let out = temp.path().join("out.json");
    std::fs::write(&lib, format!(">rep\n{}\n", seq))?;
    std::fs::write(&genome_fa, format!(">chr1\n{}{}\n", seq, seq))?;

    let cache = temp.path().join("lib.pkt");
    assert!(!cache.exists());

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "e-kmer",
            lib.to_str().unwrap(),
            genome_fa.to_str().unwrap(),
            "--keep-index",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("==> Outputs"), "pipeline failed: {stderr}");
    assert!(cache.exists(), ".pkt cache must be written");

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "e-kmer",
            lib.to_str().unwrap(),
            genome_fa.to_str().unwrap(),
            "--keep-index",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("reusing repeat table"),
        "second run must reuse the cache: {stderr}"
    );
    Ok(())
}

/// trf end-to-end on MG1655 (needs `trf` in $PATH).
#[test]
fn command_rept_trf_end_to_end() -> anyhow::Result<()> {
    for tool in ["trf"] {
        if which::which(tool).is_err() {
            eprintln!("skipping: {tool} not found");
            return Ok(());
        }
    }

    let genome = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz");
    let temp = tempfile::TempDir::new()?;
    let out = temp.path().join("out.json");

    let (_, stderr) = common::PgrCmd::new()
        .args(&["rept", "trf", genome, "-o", out.to_str().unwrap()])
        .run();
    assert!(stderr.contains("==> Outputs"), "pipeline failed: {stderr}");

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let spans = json["NC_000913"].as_str().expect("NC_000913 key missing");
    assert!(spans.contains('-'), "expected intervals, got: {spans}");
    Ok(())
}

#[test]
fn command_rept_s_align_help() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd.arg("rept").arg("s-align").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Identifies repetitive regions by self alignment"));
    Ok(())
}

/// s-align end-to-end on MG1655 (needs `lastz` in $PATH).
#[test]
fn command_rept_s_align_end_to_end() -> anyhow::Result<()> {
    for tool in ["lastz"] {
        if which::which(tool).is_err() {
            eprintln!("skipping: {tool} not found");
            return Ok(());
        }
    }

    // A 300 kb NC_000913 fragment keeps the end-to-end lastz guard while
    // cutting the self-alignment from ~10 s to well under a second.
    let genome = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/genome/mg1655.300k.fa.gz"
    );
    let temp = tempfile::TempDir::new()?;
    let out = temp.path().join("out.json");

    let (_, stderr) = common::PgrCmd::new()
        .args(&["rept", "s-align", genome, "-o", out.to_str().unwrap()])
        .run();
    assert!(stderr.contains("==> Coverage"), "pipeline failed: {stderr}");

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let spans = json["NC_000913"].as_str().expect("NC_000913 key missing");
    assert!(spans.contains('-'), "expected intervals, got: {spans}");
    Ok(())
}

/// s-align must restore dotted contig names in the runlist (regression for
/// runlist-name truncation; needs lastz in $PATH).
#[test]
fn command_rept_s_align_dotted_name() -> anyhow::Result<()> {
    for tool in ["lastz"] {
        if which::which(tool).is_err() {
            eprintln!("skipping: {tool} not found");
            return Ok(());
        }
    }

    let temp = tempfile::TempDir::new()?;
    let genome_fa = temp.path().join("genome.fa");
    let out = temp.path().join("out.json");
    let dup = random_seq(500, 31);
    let seq = format!(
        "{}{}{}{}{}",
        random_seq(300, 32),
        dup,
        random_seq(200, 33),
        dup,
        random_seq(200, 34)
    );
    std::fs::write(&genome_fa, format!(">NC_000913.1\n{}\n", seq))?;

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "s-align",
            genome_fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("==> Coverage"), "pipeline failed: {stderr}");

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let spans = json["NC_000913.1"]
        .as_str()
        .expect("NC_000913.1 key missing (runlist truncated the name?)");
    // Regression: the lift used to misparse `NC_000913.1:1-200` as contig
    // "1" and silently skip, so coverage came from unlifted window
    // coordinates. With the lift fixed, the reported intervals are the real
    // duplicated regions (0-based 300-800 and 1000-1500 -> 301-800 and
    // 1001-1500).
    assert_eq!(
        spans, "301-800,1001-1500",
        "unexpected s-align spans: {spans}"
    );
    Ok(())
}

/// s-align on a contig whose name contains ':' (`chr1:alt`) must keep the
/// full name in the runlist key (regression: `split_once(':')` / the Range
/// parser truncated it to "alt").
#[test]
fn command_rept_s_align_colon_name() -> anyhow::Result<()> {
    if which::which("lastz").is_err() {
        eprintln!("skipping: lastz not found");
        return Ok(());
    }

    let temp = tempfile::TempDir::new()?;
    let genome_fa = temp.path().join("genome.fa");
    let out = temp.path().join("out.json");
    let dup = random_seq(400, 41);
    let seq = format!(
        "{}{}{}{}{}",
        random_seq(300, 42),
        dup,
        random_seq(200, 43),
        dup,
        random_seq(200, 44)
    );
    std::fs::write(&genome_fa, format!(">chr1:alt\n{}\n", seq))?;

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "s-align",
            genome_fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("==> Coverage"), "pipeline failed: {stderr}");

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let spans = json["chr1:alt"]
        .as_str()
        .expect("chr1:alt key missing (name truncated?)");
    // Copies at 0-based 300-700 and 900-1300 (400 bp dup).
    assert_eq!(
        spans, "301-700,901-1300",
        "unexpected s-align spans: {spans}"
    );
    Ok(())
}

/// trf must resolve `fa split` names with special characters (sanitized).
#[test]
fn command_rept_trf_special_chars() -> anyhow::Result<()> {
    if which::which("trf").is_err() {
        eprintln!("skipping: trf not found");
        return Ok(());
    }

    let temp = tempfile::TempDir::new()?;
    let genome_fa = temp.path().join("genome.fa");
    let out = temp.path().join("out.json");
    let seq = format!("{}{}", random_seq(500, 21), random_seq(500, 21));
    std::fs::write(&genome_fa, format!(">chr(1):x\n{}\n", seq))?;

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "trf",
            genome_fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("==> Outputs"), "pipeline failed: {stderr}");

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let _ = json["chr(1):x"]
        .as_str()
        .expect("sanitized name not restored");
    Ok(())
}

/// `-o` pointing at an input must be rejected before the output is written;
/// the genome/repeat FASTA would otherwise be silently overwritten with the
/// runlist JSON.
#[test]
fn command_rept_output_same_as_input_rejected() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("genome.fa");
    std::fs::write(&fa, format!(">chr\n{}\n", random_seq(1200, 71)))?;
    let before = std::fs::read(&fa)?;

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "s-kmer",
            fa.to_str().unwrap(),
            "-o",
            fa.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("also an input file"),
        "output-as-input must error, got: {stderr}"
    );
    assert_eq!(std::fs::read(&fa)?, before, "input must stay intact");
    Ok(())
}

/// e-align rejects out-of-range `--min-identity`.
#[test]
fn command_rept_e_align_invalid_identity() -> anyhow::Result<()> {
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "rept",
            "e-align",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/pgr/tncentral.fa.gz"),
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz"),
            "--min-identity",
            "1.5",
            "-o",
            "/tmp/never.json",
        ])
        .run_fail();
    assert!(
        stderr.contains("must be in (0, 1]"),
        "expected range error, got: {stderr}"
    );
    Ok(())
}

#[test]
fn command_rept_without_subcommand_errors_not_panics() {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd.arg("rept").output().unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(
        stderr.contains("requires a subcommand"),
        "expected a missing-subcommand error, got: {stderr}"
    );
    assert!(!stderr.contains("panicked"), "got {stderr}");
}
