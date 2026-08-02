#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;

/// Deterministic pseudo-random DNA of length `len` (LCG, no ACGT periodicity).
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

fn write_fa(dir: &std::path::Path, name: &str, seq: &str) -> String {
    let path = dir.join(format!("{name}.fa"));
    fs::write(&path, format!(">{name}\n{seq}\n")).unwrap();
    path.to_string_lossy().to_string()
}

fn build_pgi(dir: &std::path::Path, name: &str) -> (String, String) {
    let fa = write_fa(dir, name, &random_seq(400, 42));
    let out = dir.join(format!("{name}.pgi"));
    let (_, stderr) = PgrCmd::new()
        .args(&["pgi", "build", &fa, "-o", out.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"), "build failed: {stderr}");
    (fa, out.to_string_lossy().to_string())
}

/// Parse PSL stdout into (strand, q_start, q_end, t_start, t_end, q_size).
fn parse_psl(stdout: &str) -> Vec<(String, u32, u32, u32, u32, u32)> {
    stdout
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            assert!(f.len() >= 18, "malformed PSL line: {l}");
            (
                f[8].to_string(),
                f[11].parse().unwrap(),
                f[12].parse().unwrap(),
                f[15].parse().unwrap(),
                f[16].parse().unwrap(),
                f[10].parse().unwrap(),
            )
        })
        .collect()
}

fn q_covered(records: &[(String, u32, u32, u32, u32, u32)]) -> u32 {
    records.iter().map(|r| r.2 - r.1).sum()
}

#[test]
fn command_align_pgi_identical() {
    let temp = tempfile::TempDir::new().unwrap();
    let (_, ref_idx) = build_pgi(temp.path(), "ref");
    let (_, query_idx) = build_pgi(temp.path(), "query");

    let out = temp.path().join("out.psl");
    let _ = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_idx,
            &query_idx,
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    let records = parse_psl(&fs::read_to_string(&out).unwrap());
    assert!(!records.is_empty(), "expected PSL blocks");
    assert!(
        records.iter().all(|r| r.0 == "+"),
        "identical sequences must be plus strand"
    );
    assert!(records.iter().all(|r| r.1 < r.2 && r.2 <= r.5));
    assert!(q_covered(&records) >= 200, "expected >50% query coverage");
}

#[test]
fn command_align_pgi_rc_query() {
    let temp = tempfile::TempDir::new().unwrap();
    let (_, ref_idx) = build_pgi(temp.path(), "ref");
    let rc: String = random_seq(400, 42)
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => unreachable!(),
        })
        .collect();
    let fa = write_fa(temp.path(), "query", &rc);
    let query_idx = temp.path().join("query.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&["pgi", "build", &fa, "-o", query_idx.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"));

    let out = temp.path().join("out.psl");
    let _ = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_idx,
            query_idx.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    let records = parse_psl(&fs::read_to_string(&out).unwrap());
    assert!(!records.is_empty(), "expected PSL blocks");
    assert!(
        records.iter().all(|r| r.0 == "-"),
        "RC query must be minus strand"
    );
    assert!(records.iter().all(|r| r.1 < r.2 && r.2 <= r.5));
    assert!(q_covered(&records) >= 200, "expected >50% query coverage");
}

#[test]
fn command_align_pgi_param_mismatch_fails() {
    let temp = tempfile::TempDir::new().unwrap();
    let (_, ref_idx) = build_pgi(temp.path(), "ref");
    let fa = write_fa(temp.path(), "query", &random_seq(400, 7));
    let query_idx = temp.path().join("query.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            &fa,
            "-o",
            query_idx.to_str().unwrap(),
            "--kmer",
            "20",
        ])
        .run();
    assert!(stderr.contains("wrote"));

    let out = temp.path().join("out.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_idx,
            query_idx.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("k-mer size mismatch"),
        "expected mismatch error: {stderr}"
    );
}

#[test]
fn command_align_pgi_with_sequences() {
    let temp = tempfile::TempDir::new().unwrap();
    let (ref_fa, ref_idx) = build_pgi(temp.path(), "ref");
    let (query_fa, query_idx) = build_pgi(temp.path(), "query");

    let out = temp.path().join("out.psl");
    let _ = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_idx,
            &query_idx,
            "--ref-seq",
            &ref_fa,
            "--query-seq",
            &query_fa,
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    let records = parse_psl(&fs::read_to_string(&out).unwrap());
    assert!(!records.is_empty(), "expected PSL blocks");
    // Extended blocks must carry match counts (field 0 > 0).
    let text = fs::read_to_string(&out).unwrap();
    assert!(
        text.lines()
            .map(|l| l.split_whitespace().next().unwrap().parse::<u32>().unwrap())
            .any(|m| m > 0),
        "expected a scored alignment: {text}"
    );
    assert!(q_covered(&records) >= 200, "expected >50% query coverage");
}

#[test]
fn command_align_pgi_sequences_direct() {
    // Genome inputs are indexed automatically; no `pgr pgi build` needed.
    let temp = tempfile::TempDir::new().unwrap();
    let ref_fa = write_fa(temp.path(), "ref", &random_seq(400, 42));
    let query_fa = write_fa(temp.path(), "query", &random_seq(400, 42));

    let out = temp.path().join("out.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_fa,
            &query_fa,
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("built reference index") && stderr.contains("built query index"),
        "expected automatic indexing: {stderr}"
    );
    let records = parse_psl(&fs::read_to_string(&out).unwrap());
    assert!(!records.is_empty(), "expected PSL blocks");
    assert!(q_covered(&records) >= 200, "expected >50% query coverage");
}

#[test]
fn command_align_pgi_reuses_sibling_index() {
    // A same-named .pgi next to the genome input is reused, not rebuilt.
    let temp = tempfile::TempDir::new().unwrap();
    let (ref_fa, _) = build_pgi(temp.path(), "ref");
    let (query_fa, _) = build_pgi(temp.path(), "query");

    let out = temp.path().join("out.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_fa,
            &query_fa,
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("reusing reference index") && stderr.contains("reusing query index"),
        "expected sibling index reuse: {stderr}"
    );
    assert!(
        !stderr.contains("built reference index"),
        "must not rebuild: {stderr}"
    );
    let records = parse_psl(&fs::read_to_string(&out).unwrap());
    assert!(!records.is_empty(), "expected PSL blocks");
}

#[test]
fn command_align_pgi_mixed_inputs() {
    // A genome input on one side and a .pgi on the other work together.
    let temp = tempfile::TempDir::new().unwrap();
    let ref_fa = write_fa(temp.path(), "ref", &random_seq(400, 42));
    let (query_fa, query_idx) = build_pgi(temp.path(), "query");

    let out = temp.path().join("out.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_fa,
            &query_idx,
            "--query-seq",
            &query_fa,
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("built reference index"),
        "expected one automatic build: {stderr}"
    );
    let records = parse_psl(&fs::read_to_string(&out).unwrap());
    assert!(!records.is_empty(), "expected PSL blocks");
    assert!(q_covered(&records) >= 200, "expected >50% query coverage");
}

#[test]
fn command_align_pgi_seq_flag_conflict() {
    // --ref-seq/--query-seq only apply to .pgi inputs.
    let temp = tempfile::TempDir::new().unwrap();
    let ref_fa = write_fa(temp.path(), "ref", &random_seq(400, 42));
    let query_fa = write_fa(temp.path(), "query", &random_seq(400, 42));

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_fa,
            &query_fa,
            "--ref-seq",
            &ref_fa,
            "-o",
            temp.path().join("out.psl").to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("applies only to .pgi inputs"),
        "expected conflict error: {stderr}"
    );
}

#[test]
fn command_align_pgi_sequence_validation() {
    // A .pgi input with a mismatched --ref-seq must be rejected.
    let temp = tempfile::TempDir::new().unwrap();
    let (_, ref_idx) = build_pgi(temp.path(), "ref");
    let (query_fa, query_idx) = build_pgi(temp.path(), "query");
    // Overwrite ref.fa with a same-named, shorter sequence: the index says
    // 400 bp, so validation must fail on the length.
    let _ = write_fa(temp.path(), "ref", &random_seq(300, 7));

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_idx,
            &query_idx,
            "--ref-seq",
            &ref_idx.replace(".pgi", ".fa"),
            "--query-seq",
            &query_fa,
            "-o",
            temp.path().join("out.psl").to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("length mismatch"),
        "expected contig validation error: {stderr}"
    );
}

#[test]
fn command_align_pgi_keep_index() {
    let temp = tempfile::TempDir::new().unwrap();
    let ref_fa = write_fa(temp.path(), "ref", &random_seq(400, 42));
    let query_fa = write_fa(temp.path(), "query", &random_seq(400, 42));

    let _ = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_fa,
            &query_fa,
            "--keep-index",
            "-o",
            temp.path().join("out.psl").to_str().unwrap(),
        ])
        .run();
    assert!(
        temp.path().join("ref.pgi").exists() && temp.path().join("query.pgi").exists(),
        "--keep-index must leave indexes next to the inputs"
    );
}

#[test]
fn command_align_pgi_cached_param_conflict() {
    // Explicit --kmer must agree with a reused sibling index.
    let temp = tempfile::TempDir::new().unwrap();
    let ref_fa = write_fa(temp.path(), "ref", &random_seq(400, 42));
    let query_fa = write_fa(temp.path(), "query", &random_seq(400, 42));
    let _ = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            &query_fa,
            "-o",
            temp.path().join("query.pgi").to_str().unwrap(),
            "--kmer",
            "20",
        ])
        .run();

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_fa,
            &query_fa,
            "--kmer",
            "40",
            "-o",
            temp.path().join("out.psl").to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("conflicts with the cached index"),
        "expected cached-index conflict error: {stderr}"
    );
}

#[test]
fn command_align_pgi_self() {
    // A tandem repeat (two copies of the same 400 bp sequence): self-alignment
    // must find the copy pair and must not emit an exact self-identity block.
    let temp = tempfile::TempDir::new().unwrap();
    let seq = random_seq(400, 42);
    let genome = format!(">genome\n{seq}\n{seq}\n");
    let fa = temp.path().join("genome.fa");
    fs::write(&fa, &genome).unwrap();

    let out = temp.path().join("out.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "self-align failed: {stderr}");
    let records = parse_psl(&fs::read_to_string(&out).unwrap());
    assert!(!records.is_empty(), "expected repeat blocks");
    assert!(
        records
            .iter()
            .all(|r| !(r.0 == "+" && r.1 == r.3 && r.2 == r.4)),
        "an exact self-identity block was emitted"
    );
}
