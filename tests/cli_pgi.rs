#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;

/// Write a FASTA with one sequence `seq` named `name` into `dir/name.fa`.
fn write_fa(dir: &std::path::Path, name: &str, seq: &str) -> String {
    let path = dir.join(format!("{name}.fa"));
    fs::write(&path, format!(">{name}\n{seq}\n")).unwrap();
    path.to_string_lossy().to_string()
}

/// Parse the unique k-mer count from a `pgr pgi build` log line.
fn unique_from_build(stderr: &str) -> u64 {
    stderr
        .split(" unique k-mers")
        .next()
        .unwrap()
        .split_whitespace()
        .last()
        .unwrap()
        .parse()
        .unwrap()
}

#[test]
fn command_pgi_build_mask_fasta_2bit_equivalent() {
    // --mask must skip soft-masked regions identically whether the mask comes
    // from lowercase FASTA or from 2bit mask blocks.
    let temp = tempfile::TempDir::new().unwrap();
    let upper: String = (0..200u32)
        .map(|i| b"ACGT"[(i % 4) as usize] as char)
        .collect();
    let lower: String = (0..200u32)
        .map(|i| {
            let b = b"TGCA"[(i % 4) as usize];
            (b as char).to_ascii_lowercase()
        })
        .collect();
    let fa = write_fa(temp.path(), "g", &format!("{upper}{lower}"));

    // Plain build keeps every k-mer; both masked builds must agree and be smaller.
    let out_plain = temp.path().join("plain.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&["pgi", "build", &fa, "-o", out_plain.to_str().unwrap()])
        .run();
    let plain = unique_from_build(&stderr);

    let out_fa = temp.path().join("masked_fa.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            &fa,
            "--mask",
            "-o",
            out_fa.to_str().unwrap(),
        ])
        .run();
    let masked_fa = unique_from_build(&stderr);

    // The lowercase FASTA becomes 2bit mask blocks (fa to-2bit keeps masking).
    let tb = temp.path().join("g.2bit");
    let (_, stderr) = PgrCmd::new()
        .args(&["fa", "to-2bit", &fa, "-o", tb.to_str().unwrap()])
        .run();
    assert!(!stderr.contains("error:"), "to-2bit failed: {stderr}");
    let out_tb = temp.path().join("masked_tb.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            tb.to_str().unwrap(),
            "--mask",
            "-o",
            out_tb.to_str().unwrap(),
        ])
        .run();
    let masked_tb = unique_from_build(&stderr);

    assert!(masked_fa < plain, "mask must drop k-mers");
    assert_eq!(
        masked_fa, masked_tb,
        "FASTA lowercase and 2bit mask blocks must be equivalent"
    );
}

#[test]
fn command_pgi_build_stat() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "g", &"ACGT".repeat(100));
    let out = temp.path().join("g.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&["pgi", "build", &fa, "-o", out.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"), "missing build log: {stderr}");
    assert!(out.exists());

    let (stdout, _) = PgrCmd::new()
        .args(&["pgi", "stat", out.to_str().unwrap()])
        .run();
    assert!(stdout.contains("K-mer size: 40"), "got {stdout}");
    assert!(stdout.contains("Syncmer: 8/5"), "got {stdout}");
    assert!(stdout.contains("Contigs: 1"), "got {stdout}");
    assert!(stdout.contains("Unique k-mers:"), "got {stdout}");
}

#[test]
fn command_pgi_build_from_2bit() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "g", &"ACGT".repeat(100));
    let tb = temp.path().join("g.2bit");
    let (_, stderr) = PgrCmd::new()
        .args(&["fa", "to-2bit", &fa, "-o", tb.to_str().unwrap()])
        .run();
    assert!(!stderr.contains("error:"), "2bit build failed: {stderr}");

    let out = temp.path().join("g.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            tb.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("wrote"),
        "pgi build from 2bit failed: {stderr}"
    );
    assert!(out.exists());
}

#[test]
fn command_pgi_dist_identical_and_disjoint() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa1 = write_fa(temp.path(), "a", &"ACGT".repeat(100));
    let fa2 = write_fa(temp.path(), "a2", &"ACGT".repeat(100));
    let fa3 = write_fa(temp.path(), "b", &"TTTT".repeat(100));
    let idx1 = temp.path().join("a.pgi");
    let idx2 = temp.path().join("a2.pgi");
    let idx3 = temp.path().join("b.pgi");
    for (fa, out) in [(&fa1, &idx1), (&fa2, &idx2), (&fa3, &idx3)] {
        let (_, stderr) = PgrCmd::new()
            .args(&["pgi", "build", fa, "-o", out.to_str().unwrap()])
            .run();
        assert!(stderr.contains("wrote"), "build failed: {stderr}");
    }

    // Identical sequences -> Jaccard 1, Mash 0.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "pgi",
            idx1.to_str().unwrap(),
            idx2.to_str().unwrap(),
        ])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields.len(), 9, "unexpected output: {stdout}");
    assert_eq!(
        fields[4], fields[5],
        "identical indexes should have inter == union"
    );
    assert_eq!(fields[7], "1.0000", "jaccard should be 1: {stdout}");

    // Disjoint sequences -> Jaccard 0.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "pgi",
            idx1.to_str().unwrap(),
            idx3.to_str().unwrap(),
        ])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields[4], "0", "inter should be 0: {stdout}");
    assert_eq!(fields[7], "0.0000", "jaccard should be 0: {stdout}");
}

#[test]
fn command_pgi_dist_param_mismatch_fails() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa1 = write_fa(temp.path(), "a", &"ACGT".repeat(100));
    let fa2 = write_fa(temp.path(), "b", &"ACGT".repeat(100));
    let idx1 = temp.path().join("a.pgi");
    let idx2 = temp.path().join("b.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&["pgi", "build", &fa1, "-o", idx1.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"));
    // Different k-mer size.
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            &fa2,
            "-o",
            idx2.to_str().unwrap(),
            "--kmer",
            "20",
        ])
        .run();
    assert!(stderr.contains("wrote"));

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "dist",
            "pgi",
            idx1.to_str().unwrap(),
            idx2.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("k-mer size mismatch"),
        "expected mismatch error: {stderr}"
    );
}

#[test]
fn command_pgi_dist_empty_indexes() {
    // A sequence shorter than k (default k=40) yields an index with 0 unique
    // k-mers. Two empty indexes are identical -> jaccard 1, mash 0 (not NaN
    // and not a spurious distance of 1).
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "tiny", "ACGT");
    let idx1 = temp.path().join("a.pgi");
    let idx2 = temp.path().join("b.pgi");
    for out in [&idx1, &idx2] {
        let (_, stderr) = PgrCmd::new()
            .args(&["pgi", "build", &fa, "-o", out.to_str().unwrap()])
            .run();
        assert!(stderr.contains("wrote"), "build failed: {stderr}");
    }

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "pgi",
            idx1.to_str().unwrap(),
            idx2.to_str().unwrap(),
        ])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields.len(), 9, "unexpected output: {stdout}");
    assert_eq!(fields[2], "0", "total1 should be 0: {stdout}");
    assert_eq!(
        fields[6], "0.0000",
        "empty==empty mash should be 0: {stdout}"
    );
    assert_eq!(
        fields[7], "1.0000",
        "empty==empty jaccard should be 1: {stdout}"
    );
    assert_eq!(
        fields[8], "0.0000",
        "empty first-set containment should be 0: {stdout}"
    );
    assert!(!stdout.contains("NaN"), "no NaN allowed: {stdout}");
}

#[test]
fn command_pgi_dist_containment_directional() {
    // Containment uses the FIRST index as denominator (directional), like the
    // rest of the dist family. A smaller first set fully contained in a larger
    // second set gives containment 1; swapping the arguments must change it.
    let temp = tempfile::TempDir::new().unwrap();
    // fa_short is a 60 bp sequence; fa_long shares it and adds extra sequence.
    let shared = "ACGTGCAATGGCTTAGCGTACCGAT".repeat(3); // 78 bp
    let fa_short = write_fa(temp.path(), "short", &shared);
    let fa_long = write_fa(
        temp.path(),
        "long",
        &format!("{}{}", shared, "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTA".repeat(5)),
    );
    let short_idx = temp.path().join("short.pgi");
    let long_idx = temp.path().join("long.pgi");
    for (fa, out) in [(&fa_short, &short_idx), (&fa_long, &long_idx)] {
        let (_, stderr) = PgrCmd::new()
            .args(&["pgi", "build", fa, "-o", out.to_str().unwrap()])
            .run();
        assert!(stderr.contains("wrote"), "build failed: {stderr}");
    }

    // short in long: containment (short/long-first) ~ 1 (short's k-mers all
    // present in long since long is a superset, modulo window-boundary
    // syncmer positions).
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "pgi",
            short_idx.to_str().unwrap(),
            long_idx.to_str().unwrap(),
        ])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    let c_short_in_long = fields[8].parse::<f64>().unwrap();
    // long in short: containment is long's k-mers / short's k-mers -> smaller.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "pgi",
            long_idx.to_str().unwrap(),
            short_idx.to_str().unwrap(),
        ])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    let c_long_in_short = fields[8].parse::<f64>().unwrap();

    // Directional: the two containment values differ, and the short-in-long
    // view is the larger one.
    assert!(
        c_short_in_long > c_long_in_short,
        "containment should be directional: short-in-long={c_short_in_long} long-in-short={c_long_in_short}"
    );
    assert!(
        c_short_in_long > 0.9,
        "containment of a subset-first view should be near 1: {c_short_in_long}"
    );
}

#[test]
fn command_pgi_to_hv() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "g", &"ACGT".repeat(100));
    let idx = temp.path().join("g.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&["pgi", "build", &fa, "-o", idx.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"));

    let hv = temp.path().join("g.hv");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "to-hv",
            idx.to_str().unwrap(),
            "-o",
            hv.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("hypervector"), "to-hv failed: {stderr}");
    let bytes = fs::read(&hv).unwrap();
    assert_eq!(&bytes[0..4], b"PGV1", "bad hv magic");
    // 4 magic + 4 ver + 4 k + 4 dim + 4 sparse + 8 n_kmer + 4 name_len
    // + name(1) + 4096*4 hv bytes (default dim)
    assert_eq!(bytes.len(), 4 + 4 + 4 + 4 + 4 + 8 + 4 + 1 + 4096 * 4);
}

#[test]
fn command_pgi_build_invalid_k_fails() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "g", &"ACGT".repeat(100));
    let out = temp.path().join("g.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            &fa,
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "65",
        ])
        .run_fail();
    assert!(stderr.contains("k must be in 1..=64"), "got {stderr}");
}

#[test]
fn command_pgi_without_subcommand_errors_not_panics() {
    let (_, stderr) = PgrCmd::new().args(&["pgi"]).run_fail();
    assert!(
        stderr.contains("requires a subcommand"),
        "expected a missing-subcommand error, got: {stderr}"
    );
    assert!(!stderr.contains("panicked"), "got {stderr}");
}
