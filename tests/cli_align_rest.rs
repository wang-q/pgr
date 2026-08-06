#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;

/// Deterministic pseudo-random DNA of length `len` (LCG, no ACGT periodicity).
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

/// Mutate a sequence in place at `rate` divergence (substitutions only).
fn mutate(seq: &str, rate: f64) -> String {
    let mut x = 0x9E3779B97F4A7C15u64;
    seq.chars()
        .map(|b| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let r = (x >> 33) as f64 / (1u64 << 31) as f64;
            if r < rate {
                let bases = *b"ACGT";
                bases[(x >> 33) as usize & 3] as char
            } else {
                b
            }
        })
        .collect()
}

/// Build a genome pair with three conserved regions (2000 bp each) separated by
/// 3000 bp of random (non-homologous) filler. Returns (target, query).
struct Genomes {
    target: String,
    query: String,
    /// Target start (0-based) of r1 / r2 / r3.
    r: [usize; 3],
}

fn build_genomes(seed: u64) -> Genomes {
    let region = 2000;
    let gap = 3000;
    let r1 = random_seq(region, seed);
    let r2 = random_seq(region, seed + 1);
    let r3 = random_seq(region, seed + 2);
    let target = format!(
        "{r1}{}{r2}{}{r3}",
        random_seq(gap, seed + 3),
        random_seq(gap, seed + 4)
    );
    let q1 = mutate(&r1, 0.02);
    let q2 = mutate(&r2, 0.02);
    let q3 = mutate(&r3, 0.02);
    let query = format!(
        "{q1}{}{q2}{}{q3}",
        random_seq(gap, seed + 5),
        random_seq(gap, seed + 6)
    );
    Genomes {
        target,
        query,
        r: [0, region + gap, region + gap + region + gap],
    }
}

fn write_fa(dir: &std::path::Path, name: &str, seq: &str) -> String {
    let path = dir.join(format!("{name}.fa"));
    fs::write(&path, format!(">{name}\n{seq}\n")).unwrap();
    path.to_string_lossy().to_string()
}

/// Write a multi-record FASTA (records joined with newlines).
fn write_fa_multi(dir: &std::path::Path, name: &str, seqs: &[(&str, &str)]) -> String {
    let path = dir.join(format!("{name}.fa"));
    let text = seqs
        .iter()
        .map(|(n, s)| format!(">{n}\n{s}"))
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(&path, format!("{text}\n")).unwrap();
    path.to_string_lossy().to_string()
}

/// Parse PSL into (strand, q_start, q_end, t_start, t_end).
fn parse_psl(text: &str) -> Vec<(String, u32, u32, u32, u32)> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            assert!(f.len() >= 18, "malformed PSL line: {l}");
            (
                f[8].to_string(),
                f[11].parse().unwrap(),
                f[12].parse().unwrap(),
                f[15].parse().unwrap(),
                f[16].parse().unwrap(),
            )
        })
        .collect()
}

/// Does any record cover the half-open target interval [s, e)?
fn covers(records: &[(String, u32, u32, u32, u32)], s: u32, e: u32) -> bool {
    records.iter().any(|r| r.3 < e && r.4 > s && r.1 < r.2)
}

/// Run `align rest` on the given pgi PSL (or None to auto-run pgi).
/// Returns the output PSL text.
fn run_rest(
    temp: &tempfile::TempDir,
    target: &str,
    query: &str,
    pgi_psl: Option<&str>,
    extra: &[&str],
) -> (String, String) {
    let out = temp.path().join("out.psl");
    let mut args: Vec<String> = vec!["align".into(), "rest".into(), target.into(), query.into()];
    if let Some(p) = pgi_psl {
        args.push("--avail-psl".into());
        args.push(p.into());
    }
    for e in extra {
        args.push(e.to_string());
    }
    args.push("-o".into());
    args.push(out.to_str().unwrap().to_string());
    let arg_refs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
    let (_, stderr) = PgrCmd::new().args(&arg_refs).run();
    let text = fs::read_to_string(&out).unwrap();
    (text, stderr)
}

/// Gate the whole module on lastz being installed.
fn lastz_missing() -> bool {
    which::which("lastz").is_err()
}

/// Keep only the pgi blocks covering a target region (sorted by tStart).
fn pgi_blocks_on_target(
    temp: &tempfile::TempDir,
    target: &str,
    query: &str,
    keep: &[usize],
) -> String {
    let pgi_out = temp.path().join("pgi.psl");
    PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            target,
            query,
            "-o",
            pgi_out.to_str().unwrap(),
        ])
        .run();
    let text = fs::read_to_string(&pgi_out).unwrap();
    let mut lines: Vec<&str> = text.lines().collect();
    lines.sort_by(|a, b| {
        let fa: Vec<&str> = a.split_whitespace().collect();
        let fb: Vec<&str> = b.split_whitespace().collect();
        fa[16]
            .parse::<u32>()
            .unwrap()
            .cmp(&fb[16].parse::<u32>().unwrap())
    });
    let kept: Vec<&str> = keep.iter().map(|&i| lines[i]).collect();
    format!("{}\n", kept.join("\n"))
}

#[test]
fn command_align_rest_fills_complement_beyond_anchor_gaps() {
    if lastz_missing() {
        eprintln!("lastz not installed; skipping");
        return;
    }
    // Keep only the r1 anchor: rest must fill everything else (the r1-r2
    // gap, r2, the r2-r3 gap and r3) via the whole-genome complement.
    let temp = tempfile::TempDir::new().unwrap();
    let g = build_genomes(71);
    let t = write_fa(temp.path(), "target", &g.target);
    let q = write_fa(temp.path(), "query", &g.query);

    let reduced = temp.path().join("pgi_reduced.psl");
    fs::write(&reduced, pgi_blocks_on_target(&temp, &t, &q, &[0])).unwrap();

    let (rest, _) = run_rest(&temp, &t, &q, Some(reduced.to_str().unwrap()), &[]);
    let records = parse_psl(&rest);
    assert!(
        records.len() >= 3,
        "expected anchors plus LASTZ fills: {rest}"
    );
    for &s in &g.r {
        assert!(
            covers(&records, s as u32, (s + 2000) as u32),
            "region at {s} not covered by rest: {rest}"
        );
    }
}

#[test]
fn command_align_rest_multicontig_query_uses_all_query_holes() {
    if lastz_missing() {
        eprintln!("lastz not installed; skipping");
        return;
    }
    // Query is one file with two contigs: qa carries the r1/r2/r3 homologs
    // (plus fillers), qb is unrelated. The query-side holes span both
    // contigs; LASTZ must align against the merged multi-sequence FASTA and
    // lift coordinates per contig.
    let temp = tempfile::TempDir::new().unwrap();
    let g = build_genomes(72);
    let t = write_fa(temp.path(), "target", &g.target);
    let q = write_fa_multi(
        temp.path(),
        "query",
        &[("qa", &g.query), ("qb", &random_seq(4000, 73))],
    );

    let reduced = temp.path().join("pgi_reduced.psl");
    fs::write(&reduced, pgi_blocks_on_target(&temp, &t, &q, &[0])).unwrap();

    let (rest, _) = run_rest(&temp, &t, &q, Some(reduced.to_str().unwrap()), &[]);
    let records = parse_psl(&rest);
    assert!(
        records.len() >= 2,
        "expected anchors plus a LASTZ fill: {rest}"
    );
    for &s in &g.r {
        assert!(
            covers(&records, s as u32, (s + 2000) as u32),
            "region at {s} not covered: {rest}"
        );
    }
}

#[test]
fn command_align_rest_excises_small_anchors() {
    if lastz_missing() {
        eprintln!("lastz not installed; skipping");
        return;
    }
    // With a huge --min-anchor every anchor is excised, so the complement is
    // the whole target; LASTZ must still recover the conserved regions.
    let temp = tempfile::TempDir::new().unwrap();
    let g = build_genomes(74);
    let t = write_fa(temp.path(), "target", &g.target);
    let q = write_fa(temp.path(), "query", &g.query);

    let (rest, _) = run_rest(&temp, &t, &q, None, &["--min-anchor", "1000000"]);
    let records = parse_psl(&rest);
    for &s in &g.r {
        assert!(
            covers(&records, s as u32, (s + 2000) as u32),
            "region at {s} not covered after full excise: {rest}"
        );
    }
}

#[test]
fn command_align_rest_missing_subcommand_errors_not_panics() {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd.arg("align").arg("rest").output().unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(
        stderr.contains("required"),
        "expected missing-argument error, got: {stderr}"
    );
}
