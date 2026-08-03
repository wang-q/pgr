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
    let bases = [b'A', b'C', b'G', b'T'];
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
/// genome should be reported as covered intervals (needs `spanr` in $PATH).
#[test]
fn command_rept_e_align_end_to_end() -> anyhow::Result<()> {
    if which::which("spanr").is_err() {
        eprintln!("skipping: spanr not found");
        return Ok(());
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
    // Dotted contig name: `spanr cover` truncates these to the last '.'
    // segment, and e-align must restore the full name in the runlist.
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

/// s-kmer on a genome with a dotted contig name must report the full name
/// in the runlist (regression for `spanr cover` truncation; needs FastK,
/// Profex and spanr in $PATH).
#[test]
fn command_rept_s_kmer_dotted_name() -> anyhow::Result<()> {
    for tool in ["FastK", "Profex", "spanr"] {
        if which::which(tool).is_err() {
            eprintln!("skipping: {tool} not found");
            return Ok(());
        }
    }

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
    let _ = json["NC_000913.1"]
        .as_str()
        .expect("NC_000913.1 key missing (spanr truncated the name?)");
    Ok(())
}

/// e-kmer end-to-end on the repo's small fixtures (needs FastK, Profex and
/// spanr in $PATH). Guards the FastK pipeline and the dotted-name mapping.
#[test]
fn command_rept_e_kmer_end_to_end() -> anyhow::Result<()> {
    for tool in ["FastK", "Profex", "spanr"] {
        if which::which(tool).is_err() {
            eprintln!("skipping: {tool} not found");
            return Ok(());
        }
    }

    let lib = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/pgr/tncentral.fa.gz");
    let genome = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/genome/mg1655.fa.gz");
    let temp = tempfile::TempDir::new()?;
    let out = temp.path().join("out.json");

    let (_, stderr) = common::PgrCmd::new()
        .args(&["rept", "e-kmer", lib, genome, "-o", out.to_str().unwrap()])
        .run();
    assert!(stderr.contains("==> Outputs"), "pipeline failed: {stderr}");

    let json: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out)?)?;
    let spans = json["NC_000913"].as_str().expect("NC_000913 key missing");
    assert!(spans.contains('-'), "expected intervals, got: {spans}");
    Ok(())
}

/// trf end-to-end on MG1655 (needs `trf` and `spanr` in $PATH).
#[test]
fn command_rept_trf_end_to_end() -> anyhow::Result<()> {
    for tool in ["trf", "spanr"] {
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
