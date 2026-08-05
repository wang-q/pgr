#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;

/// Deterministic pseudo-random DNA (LCG, no ACGT periodicity).
fn random_seq(len: usize, seed: u64) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[(x >> 33) as usize & 3] as char
        })
        .collect()
}

fn write_fa(dir: &std::path::Path, name: &str, text: &str) -> String {
    let path = dir.join(format!("{name}.fa"));
    fs::write(&path, text).unwrap();
    path.to_string_lossy().to_string()
}

/// A tandem repeat (two identical 1200 bp copies) flanked by unique sequence.
fn tandem_genome() -> String {
    let dup = random_seq(1200, 11);
    format!(">genome\n{dup}\n{dup}\n{}\n", random_seq(400, 12))
}

#[test]
fn command_sd_search_pgi_engine() {
    // The native pgi engine needs no external tools and must find the tandem
    // repeat as an SD hit.
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "genome", &tandem_genome());

    let out = temp.path().join("hits.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "sd",
            "search",
            &fa,
            "--engine",
            "pgi",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "pgi search failed: {stderr}");
    let hits = fs::read_to_string(&out).unwrap();
    assert!(!hits.trim().is_empty(), "pgi engine found no SD hits");

    // The T2T-CHM13 filters apply to the pgi engine too: a huge min-len
    // drops everything.
    let filtered = temp.path().join("hits2.psl");
    let _ = PgrCmd::new()
        .args(&[
            "sd",
            "search",
            &fa,
            "--engine",
            "pgi",
            "--min-len",
            "100000",
            "-o",
            filtered.to_str().unwrap(),
        ])
        .run();
    assert!(
        fs::read_to_string(&filtered).unwrap().trim().is_empty(),
        "min-len filter must drop short hits"
    );
}

/// An inverted repeat flanked by non-homologous sequence must be detected by
/// the pgi engine (regression: greedy chaining merged the two reciprocal
/// chains - they share one diagonal - into a chimeric low-identity block
/// that the SD filter dropped).
#[test]
fn command_sd_search_pgi_inverted_repeat() {
    let temp = tempfile::TempDir::new().unwrap();
    let dup = random_seq(1200, 21);
    let rc: String = dup
        .chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect();
    let genome = format!(
        "{}{}{}{}{}",
        random_seq(2000, 22),
        dup,
        random_seq(1800, 23),
        rc,
        random_seq(1500, 24)
    );
    let fa = write_fa(temp.path(), "genome", &format!(">chr\n{genome}\n"));

    let out = temp.path().join("hits.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "sd",
            "search",
            &fa,
            "--engine",
            "pgi",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "pgi search failed: {stderr}");

    let hits = fs::read_to_string(&out).unwrap();
    let blocks: Vec<Vec<&str>> = hits
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(blocks.len(), 2, "expected two reciprocal blocks: {hits}");
    for b in &blocks {
        // Minus-strand, >= 1000 bp, no mismatches.
        assert_eq!(b[8], "-", "inverted copies align on '-': {hits}");
        let qlen: i32 = b[12].parse::<i32>().unwrap() - b[11].parse::<i32>().unwrap();
        assert!(qlen >= 1000, "block too short: {hits}");
        assert_eq!(b[1], "0", "no mismatches expected: {hits}");
    }
}

/// Inverted copies closer than `max_gap` (1000): the two reciprocal chains
/// sit on the same diagonal and their seeds are within the greedy chaining
/// gap, so the loop used to bridge them into one chimeric chain whose diluted
/// identity failed the SD filter - losing the pair entirely. The greedy loop
/// must split at the non-homologous intervening span.
#[test]
fn command_sd_search_pgi_close_inverted_repeat() {
    let temp = tempfile::TempDir::new().unwrap();
    let dup = random_seq(1200, 61);
    let rc: String = dup
        .chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect();
    // 800 bp gap: the reciprocal chains' seeds stay within max_gap 1000.
    let genome = format!(
        "{}{}{}{}{}",
        random_seq(5000, 62),
        dup,
        random_seq(800, 63),
        rc,
        random_seq(5000, 64)
    );
    let fa = write_fa(temp.path(), "genome", &format!(">chr\n{genome}\n"));

    let out = temp.path().join("hits.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "sd",
            "search",
            &fa,
            "--engine",
            "pgi",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "pgi search failed: {stderr}");

    let hits = fs::read_to_string(&out).unwrap();
    let blocks: Vec<Vec<&str>> = hits
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(blocks.len(), 2, "expected two reciprocal blocks: {hits}");
    for b in &blocks {
        assert_eq!(b[8], "-", "inverted copies align on '-': {hits}");
        let qlen: i32 = b[12].parse::<i32>().unwrap() - b[11].parse::<i32>().unwrap();
        assert!(qlen >= 1000, "block too short: {hits}");
        assert_eq!(b[1], "0", "no mismatches expected: {hits}");
    }
}

/// Four copies of one repeat (two forward, two reverse) whose copy pairs sit
/// on nearby diagonals with gaps inside the merge window. Regression: the
/// adjacent-chain merge stitched two independent copy pairs into a chimeric
/// chain whose diluted identity was dropped by the SD filter, losing the
/// real hits (the elementary SDs of one copy were left non-core).
#[test]
fn command_sd_search_pgi_multi_copy_close_diagonals() {
    let temp = tempfile::TempDir::new().unwrap();
    let dup = random_seq(2000, 31);
    let rc: String = dup
        .chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect();
    let genome = format!(
        "{}{}{}{}{}{}{}{}{}",
        random_seq(2787, 32),
        dup,
        random_seq(1646, 33),
        dup,
        random_seq(2351, 34),
        rc,
        random_seq(3590, 35),
        rc,
        random_seq(2000, 36)
    );
    let fa = write_fa(temp.path(), "genome", &format!(">chr\n{genome}\n"));

    let out = temp.path().join("hits.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "sd",
            "search",
            &fa,
            "--engine",
            "pgi",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "pgi search failed: {stderr}");

    let hits = fs::read_to_string(&out).unwrap();
    let blocks: Vec<Vec<&str>> = hits
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split('\t').collect())
        .collect();
    // All six copy pairs, both directions: every copy must participate in a
    // minus-strand hit with each reverse copy and a plus-strand hit with each
    // forward copy. The chimeric merge used to swallow the cross pairs
    // (copy 1 <-> copy 4 and copy 2 <-> copy 3).
    assert_eq!(blocks.len(), 12, "expected all 12 reciprocal hits: {hits}");
    let minus: Vec<(u32, u32)> = blocks
        .iter()
        .filter(|b| b[8] == "-")
        .map(|b| (b[11].parse().unwrap(), b[12].parse().unwrap()))
        .collect();
    assert_eq!(
        minus.len(),
        8,
        "expected 8 minus-strand reciprocal hits: {hits}"
    );
    // The four copies (0-based intervals from the layout above) must each
    // appear as the query of a minus-strand hit: the cross pairs
    // (copy 1 <-> copy 4, copy 2 <-> copy 3) used to be swallowed by the
    // chimeric merge, leaving the reverse copies without minus hits.
    let copies = [
        (2787u32, 4787u32),
        (6433, 8433),
        (10784, 12784),
        (16374, 18374),
    ];
    for (s, e) in copies {
        let covered = minus.iter().any(|&(qs, qe)| qs < e && s < qe);
        assert!(covered, "copy {s}-{e} must appear in a minus hit: {hits}");
    }
}

/// A `.pgi` index as the genome input must be rejected with a friendly error
/// instead of silently returning an empty hit list (blocks align without
/// extension sequences and score 0).
#[test]
fn command_sd_search_rejects_pgi_input() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "genome", &tandem_genome());
    let idx = temp.path().join("genome.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&["pgi", "build", &fa, "-o", idx.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"), "pgi build failed: {stderr}");

    let out = temp.path().join("hits.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "sd",
            "search",
            idx.to_str().unwrap(),
            "--engine",
            "pgi",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("needs genome FASTA"),
        "pgi input must error, got: {stderr}"
    );
}

/// `-o` pointing at an input file must be rejected before the output is
/// written (the FASTA/PAF/BED inputs would otherwise be silently overwritten
/// with the transformed output).
#[test]
fn command_sd_output_same_as_input_rejected() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "genome", &tandem_genome());

    for args in [
        vec!["sd", "search", &fa, "-o", &fa],
        vec!["sd", "decompose", &fa, "-o", &fa],
    ] {
        let before = fs::read(&fa).unwrap();
        let (_, stderr) = PgrCmd::new().args(&args).run();
        assert!(
            stderr.contains("also an input file"),
            "output-as-input must error, got: {stderr}"
        );
        assert_eq!(
            fs::read(&fa).unwrap(),
            before,
            "input must stay intact after rejected run"
        );
    }
}

#[test]
fn command_sd_search_lastz_engine() {
    if which::which("lastz").is_err() {
        eprintln!("Skipping command_sd_search_lastz_engine: lastz not installed");
        return;
    }
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "genome", &tandem_genome());

    let out = temp.path().join("hits.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "sd",
            "search",
            &fa,
            "--engine",
            "lastz",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(!stderr.contains("error:"), "lastz search failed: {stderr}");
    assert!(
        !fs::read_to_string(&out).unwrap().trim().is_empty(),
        "lastz engine found no SD hits"
    );
}

/// Regression: `sd search --engine lastz` without `--preset` must use the
/// documented set01 default (the CLI used to omit preset parameters entirely,
/// running lastz with its own built-in defaults).
#[test]
fn command_sd_search_lastz_default_preset() {
    if which::which("lastz").is_err() {
        eprintln!("Skipping command_sd_search_lastz_default_preset: lastz not installed");
        return;
    }
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "genome", &tandem_genome());
    let out = temp.path().join("hits.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "sd",
            "search",
            &fa,
            "--engine",
            "lastz",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("K=3000"),
        "set01 preset params must be applied by default, got: {stderr}"
    );
}

#[test]
fn command_sd_cross_pgi_engine() {
    // Two genomes sharing a duplicated region: the pgi engine maps the
    // homology across genomes without lastz.
    let temp = tempfile::TempDir::new().unwrap();
    let dup = random_seq(3000, 11);
    let fa_a = write_fa(
        temp.path(),
        "a",
        &format!(">a\n{}\n{dup}\n", random_seq(500, 12)),
    );
    let fa_b = write_fa(
        temp.path(),
        "b",
        &format!(">b\n{}\n{dup}\n", random_seq(500, 13)),
    );

    let out = temp.path().join("cross.paf");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "sd",
            "cross",
            &fa_a,
            &fa_b,
            "--engine",
            "pgi",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(!stderr.contains("error:"), "pgi cross failed: {stderr}");
    assert!(
        !fs::read_to_string(&out).unwrap().trim().is_empty(),
        "pgi cross found no cross-genome homology"
    );
}

#[test]
fn command_sd_cross_pgi_engine_relative_paths() {
    // Regression: the query path must resolve before the pipeline chdirs
    // into its tempdir; relative query paths used to break pgi cross.
    let temp = tempfile::TempDir::new().unwrap();
    let dup = random_seq(3000, 11);
    fs::write(
        temp.path().join("a.fa"),
        format!(">a\n{}\n{dup}\n", random_seq(500, 12)),
    )
    .unwrap();
    fs::write(
        temp.path().join("b.fa"),
        format!(">b\n{}\n{dup}\n", random_seq(500, 13)),
    )
    .unwrap();

    let (_, stderr) = PgrCmd::new()
        .current_dir(temp.path())
        .args(&[
            "sd",
            "cross",
            "a.fa",
            "b.fa",
            "--engine",
            "pgi",
            "-o",
            "cross.paf",
        ])
        .run();
    assert!(
        !stderr.contains("error:"),
        "relative-path cross failed: {stderr}"
    );
    assert!(
        !fs::read_to_string(temp.path().join("cross.paf"))
            .unwrap()
            .trim()
            .is_empty(),
        "pgi cross found no cross-genome homology"
    );
}

/// `sd run` end-to-end on a synthetic duplicated genome: the full pipeline
/// (search -> align -> cluster -> decompose -> cover) must produce a
/// CORE-annotated elementary BED.
#[test]
fn command_sd_run_end_to_end() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "genome", &tandem_genome());
    let outdir = temp.path().join("sd_out");

    let (_, stderr) = PgrCmd::new()
        .args(&["sd", "run", &fa, "-o", outdir.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"), "sd run failed: {stderr}");

    let bed = fs::read_to_string(outdir.join("out.elem.bed")).unwrap();
    assert!(bed.contains("CORE"), "expected CORE rows, got: {bed}");
    assert!(
        bed.lines().count() >= 2,
        "expected 2+ fragments, got: {bed}"
    );
}

/// `sd run` on a plain-gzip (non-BGZF) genome must work end to end: the
/// `.loc`-based cluster step used to fail on plain gzip (it only accepted
/// plain or BGZF files), even though search/align already accept it.
#[test]
fn command_sd_run_gzipped_genome() {
    use std::io::Write;

    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "genome", &tandem_genome());
    let gz = temp.path().join("genome.fa.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        fs::File::create(&gz).unwrap(),
        flate2::Compression::default(),
    );
    encoder
        .write_all(fs::read_to_string(&fa).unwrap().as_bytes())
        .unwrap();
    encoder.finish().unwrap();

    let outdir = temp.path().join("sd_out");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "sd",
            "run",
            gz.to_str().unwrap(),
            "-o",
            outdir.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "sd run on gz failed: {stderr}");
    let bed = fs::read_to_string(outdir.join("out.elem.bed")).unwrap();
    assert!(bed.contains("CORE"), "expected CORE rows, got: {bed}");
    assert!(
        !temp.path().join("genome.fa.gz.loc").exists(),
        "no stray .loc may be written next to a gzip input"
    );
}

/// `sd run` must produce byte-identical output across runs: the cluster
/// numbering (and the global set_id renumbering) must not depend on the
/// HashMap iteration order (which is randomized per process). Regression:
/// two distinct repeat families used to swap set_id/row order between runs.
#[test]
fn command_sd_run_output_deterministic_across_runs() {
    let temp = tempfile::TempDir::new().unwrap();
    let dup1 = random_seq(1200, 81);
    let dup2 = random_seq(1100, 82);
    let genome = format!(
        "{}{}{}{}{}{}{}{}{}",
        random_seq(8000, 83),
        dup1,
        random_seq(3000, 84),
        dup1,
        random_seq(2000, 85),
        dup2,
        random_seq(1500, 86),
        dup2,
        random_seq(8000, 87)
    );
    let fa = write_fa(temp.path(), "genome", &format!(">chr\n{genome}\n"));

    let out1 = temp.path().join("out1");
    let out2 = temp.path().join("out2");
    let (_, stderr) = PgrCmd::new()
        .args(&["sd", "run", &fa, "-o", out1.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"), "sd run failed: {stderr}");
    let (_, stderr) = PgrCmd::new()
        .args(&["sd", "run", &fa, "-o", out2.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"), "sd run failed: {stderr}");

    let a = fs::read_to_string(out1.join("out.elem.bed")).unwrap();
    let b = fs::read_to_string(out2.join("out.elem.bed")).unwrap();
    assert_eq!(a, b, "sd run output must be deterministic");
}

/// Regression: `sd run --engine lastz --preset <p>` used to pass
/// `"--preset set01"` as a single argv element to the inner `sd search`,
/// which clap rejects ("unexpected argument"). The preset must expand to
/// two separate arguments. When lastz is installed the full pipeline runs;
/// otherwise the inner search still fails on the missing binary rather than
/// on argument parsing.
#[test]
fn command_sd_run_lastz_preset_parses() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "genome", &tandem_genome());
    let outdir = temp.path().join("sd_out");

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "sd",
            "run",
            &fa,
            "--engine",
            "lastz",
            "--preset",
            "set01",
            "-o",
            outdir.to_str().unwrap(),
        ])
        .run();
    assert!(
        !stderr.contains("unexpected argument"),
        "preset must expand to separate args, got: {stderr}"
    );
    if which::which("lastz").is_ok() {
        assert!(
            stderr.contains("wrote"),
            "lastz present: the full run must succeed, got: {stderr}"
        );
    }
}

#[test]
fn command_sd_without_subcommand_errors_not_panics() {
    let (_, stderr) = PgrCmd::new().args(&["sd"]).run_fail();
    assert!(
        stderr.contains("requires a subcommand"),
        "expected a missing-subcommand error, got: {stderr}"
    );
    assert!(!stderr.contains("panicked"), "got {stderr}");
}
