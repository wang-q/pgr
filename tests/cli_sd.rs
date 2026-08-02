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
