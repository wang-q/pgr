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
