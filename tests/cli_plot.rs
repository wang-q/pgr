#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_plot_venn2() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "plot",
            "venn",
            "tests/plot/rocauc.result.tsv",
            "tests/plot/mcox.05.result.tsv",
        ])
        .run();

    assert!(stdout.contains("(-2.8, -1.8) { rocauc }"));
    assert!(stdout.contains("(-2,    0) { 669 }"));
}

#[test]
fn command_plot_venn3() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "plot",
            "venn",
            "tests/plot/rocauc.result.tsv",
            "tests/plot/mcox.05.result.tsv",
            "tests/plot/mcox.result.tsv",
        ])
        .run();

    assert!(stdout.contains("(-2.8, -1.8) { rocauc }"));
    assert!(stdout.contains("(-2,   -0.2) { 161 }"));
}

#[test]
fn command_plot_venn4() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "plot",
            "venn",
            "tests/plot/rocauc.result.tsv",
            "tests/plot/rocauc.result.tsv",
            "tests/plot/mcox.05.result.tsv",
            "tests/plot/mcox.result.tsv",
        ])
        .run();

    assert!(stdout.contains("(-2.2, -2.6) { rocauc }"));
    assert!(stdout.contains("(-2.2,  1.5) { 161 }"));
}

#[test]
fn command_plot_dot() {
    let (stdout, _) = PgrCmd::new()
        .args(&["plot", "dot", "tests/plot/dot.paf", "--min-len", "400"])
        .run();

    assert!(stdout.starts_with("<?xml"));
    assert!(stdout.contains(r#"<svg xmlns="http://www.w3.org/2000/svg""#));
    // 4 records; two pass --min-len 400 (block 800 and 1250)
    let seg = stdout.split(r#"id="segments""#).nth(1).unwrap();
    assert_eq!(seg.matches("<line ").count(), 2);
    assert!(stdout.contains(r#"id="segments""#));
    // 2 identity-colored line strokes plus the border rect stroke
    assert!(stdout.matches("stroke=\"#").count() >= 3);
}

#[test]
fn command_plot_dot_filters_everything() {
    let (_, stderr) = PgrCmd::new()
        .args(&["plot", "dot", "tests/plot/dot.paf", "--min-len", "100000"])
        .run_fail();

    assert!(stderr.contains("no alignments pass filters"));
}

#[test]
fn command_plot_hh() {
    let (stdout, _) = PgrCmd::new()
        .args(&["plot", "hh", "tests/plot/hist.tsv"])
        .run();

    assert!(stdout.contains("31   0 0.0200"));
    assert!(stdout.contains("31   1 0.0000"));

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "plot",
            "hh",
            "tests/plot/hist.tsv",
            "-g",
            "2",
            "--bins",
            "20",
        ])
        .run();

    assert!(stdout.contains("11   0 0.0600"));
    assert!(stdout.contains("11   1 0.1600"));
}

#[test]
fn command_plot_nrps() {
    let (stdout, _) = PgrCmd::new()
        .args(&["plot", "nrps", "tests/plot/srf.tsv"])
        .run();

    assert!(stdout.contains("(-0.4cm,0) -- (\\x1 + 0.2cm,0)"));
    assert!(!stdout.contains("\\textbf{M}ethyltransferase"));

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "plot",
            "nrps",
            "tests/plot/srf.tsv",
            "--legend",
            "--color",
            "black",
        ])
        .run();

    assert!(stdout.contains("     draw=black,"));
    assert!(stdout.contains("\\textbf{M}ethyltransferase"));
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

#[test]
fn command_plot_heat_end_to_end() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    std::fs::write(
        &fa,
        ">g\nGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCC\n>a\nAAAATTTTAAAATTTTAAAATTTTAAAATTTTAAAATTTTAAAATTTT\n",
    )?;
    let kgc = temp.path().join("g.kgc");
    PgrCmd::new()
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
    let tex = temp.path().join("h.tex");
    PgrCmd::new()
        .args(&[
            "plot",
            "heat",
            kgc.to_str().unwrap(),
            "-o",
            tex.to_str().unwrap(),
        ])
        .run();
    let text = std::fs::read_to_string(&tex)?;
    assert!(text.contains("addplot"));
    assert!(text.contains("k-mer coverage"));
    Ok(())
}

#[test]
fn command_plot_spectra_end_to_end() -> anyhow::Result<()> {
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

    let pkt = temp.path().join("t.pkt");
    PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fq.to_str().unwrap(),
            "-k",
            "17",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .run();
    let hist = temp.path().join("t.hist");
    PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            "-t",
            pkt.to_str().unwrap(),
            "-o",
            hist.to_str().unwrap(),
        ])
        .run();
    let gs_out = temp.path().join("gs");
    PgrCmd::new()
        .args(&[
            "kmer",
            "gsize",
            "-t",
            pkt.to_str().unwrap(),
            "--model",
            "-o",
            gs_out.to_str().unwrap(),
        ])
        .run();
    let tex = temp.path().join("s.tex");
    PgrCmd::new()
        .args(&[
            "plot",
            "spectra",
            hist.to_str().unwrap(),
            gs_out.join("model.txt").to_str().unwrap(),
            "-o",
            tex.to_str().unwrap(),
        ])
        .run();
    let text = std::fs::read_to_string(&tex)?;
    assert!(text.contains("addplot"));
    assert!(text.contains("observed"));
    assert!(text.contains("kcov"));
    // Single standard view with correct percent escaping.
    assert_eq!(text.matches("\\begin{tikzpicture}").count(), 1);
    assert!(text.contains("\\%"), "percent must be escaped in the title");
    assert!(text.contains("full model"));
    Ok(())
}

#[test]
fn command_plot_heat_rejects_empty_matrix() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let kgc = temp.path().join("empty.kgc");
    std::fs::write(&kgc, "GCP\tKF\tCount\n")?;
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "plot",
            "heat",
            kgc.to_str().unwrap(),
            "-o",
            temp.path().join("h.tex").to_str().unwrap(),
        ])
        .run_fail();
    assert!(stderr.contains("no matrix rows"), "stderr: {stderr}");
    Ok(())
}

#[test]
fn command_plot_spectra_rejects_model_without_kmercov() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    std::fs::write(&fa, ">s\nACGTACGTACGTACGTACGT\n")?;
    let hist = temp.path().join("t.hist");
    PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            hist.to_str().unwrap(),
        ])
        .run();
    let bad_model = temp.path().join("bad_model.txt");
    std::fs::write(&bad_model, "d 0.1\nbias 0.5\n")?;
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "plot",
            "spectra",
            hist.to_str().unwrap(),
            bad_model.to_str().unwrap(),
            "-o",
            temp.path().join("s.tex").to_str().unwrap(),
        ])
        .run_fail();
    assert!(stderr.contains("kmercov"), "stderr: {stderr}");
    Ok(())
}
