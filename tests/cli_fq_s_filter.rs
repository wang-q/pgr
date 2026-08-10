#[macro_use]
#[path = "common/mod.rs"]
mod common;

#[test]
fn command_fq_s_filter_real_lambda() -> anyhow::Result<()> {
    // Real Lambda reads carry sequencing errors; quorum-style self-checking
    // must flag a few percent but keep the overwhelming majority.
    let temp = tempfile::TempDir::new()?;
    let kept = temp.path().join("kept.fq");
    let discarded = temp.path().join("discarded.fq");
    common::PgrCmd::new()
        .args(&[
            "fq",
            "s-filter",
            "tests/bbtools/Lambda/golden/filter.fq.gz",
            "-k",
            "31",
            "-o",
            kept.to_str().unwrap(),
            "--discard-file",
            discarded.to_str().unwrap(),
        ])
        .run();
    let n_kept = std::fs::read_to_string(&kept)?.lines().count() / 4;
    let n_disc = std::fs::read_to_string(&discarded)?.lines().count() / 4;
    assert_eq!(n_kept + n_disc, 36384, "all reads must be classified");
    let frac = n_disc as f64 / (n_kept + n_disc) as f64;
    assert!(
        (0.02..=0.05).contains(&frac),
        "flagged fraction {frac:.3} ({n_disc}/{}) must be a few percent",
        n_kept + n_disc
    );
    Ok(())
}

#[test]
fn command_fq_s_filter_lambda_flag_count_matches_quorum() -> anyhow::Result<()> {
    // On the Lambda filter golden (k=24, default skip/good/anchor-count),
    // pgr flags 1267 reads; quorum error_correct_reads with the same
    // parameters flags 1264 (3 borderline reads differ). This pins the
    // backward-extension fix (an off-by-one that previously flagged
    // ~27k reads when skip/good were non-default).
    let temp = tempfile::TempDir::new()?;
    let kept = temp.path().join("kept.fq");
    let discarded = temp.path().join("discarded.fq");
    common::PgrCmd::new()
        .args(&[
            "fq",
            "s-filter",
            "tests/bbtools/Lambda/golden/filter.fq.gz",
            "-k",
            "24",
            "-o",
            kept.to_str().unwrap(),
            "--discard-file",
            discarded.to_str().unwrap(),
        ])
        .run();
    let flagged = std::fs::read_to_string(&discarded)?.lines().count() / 4;
    assert!(
        (1240..=1300).contains(&flagged),
        "flagged count {flagged} far from the quorum reference 1264"
    );
    Ok(())
}

#[test]
fn command_fq_s_filter_end_to_end() -> anyhow::Result<()> {
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
            "fq",
            "s-filter",
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
