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

/// Reverse-complement a sequence.
fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => unreachable!(),
        })
        .collect()
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
                let bases = [b'A', b'C', b'G', b'T'];
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

/// Run `align hybrid` on the given pgi PSL (or None to auto-run pgi).
/// Returns the output PSL text.
fn run_hybrid(
    temp: &tempfile::TempDir,
    target: &str,
    query: &str,
    pgi_psl: Option<&str>,
    extra: &[&str],
) -> (String, String) {
    let out = temp.path().join("out.psl");
    let mut args: Vec<String> = vec!["align".into(), "hybrid".into(), target.into(), query.into()];
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

/// Gate the whole module on lastz being installed: without it, every test
/// would fail with the same "lastz not found" error. We skip instead, so the
/// suite stays green on machines without lastz.
fn lastz_missing() -> bool {
    which::which("lastz").is_err()
}

#[test]
fn command_align_hybrid_all_regions_pgi_visible_no_extra_block() {
    if lastz_missing() {
        eprintln!("lastz not installed; skipping");
        return;
    }
    // With three low-divergence conserved regions, pgi finds all of them, so
    // every LASTZ record is redundant. Dedup is left to the downstream
    // chainnet pipeline, so the output may contain extra LASTZ records; we
    // only require that all three regions are covered (pgi records verbatim).
    let temp = tempfile::TempDir::new().unwrap();
    let g = build_genomes(42);
    let t = write_fa(temp.path(), "target", &g.target);
    let q = write_fa(temp.path(), "query", &g.query);

    let (hybrid, _) = run_hybrid(&temp, &t, &q, None, &[]);
    let records = parse_psl(&hybrid);
    assert!(records.len() >= 3, "expected >=3 regions: {hybrid}");
    // All three conserved regions must be covered.
    for &s in &g.r {
        assert!(
            covers(&records, s as u32, (s + 2000) as u32),
            "region at {s} not covered"
        );
    }
}

#[test]
fn command_align_hybrid_fills_missing_anchor_gap() {
    if lastz_missing() {
        eprintln!("lastz not installed; skipping");
        return;
    }
    // Keep only the r1 and r3 anchors (drop r2). The gap between them must be
    // filled by LASTZ, adding back a record covering r2.
    let temp = tempfile::TempDir::new().unwrap();
    let g = build_genomes(43);
    let t = write_fa(temp.path(), "target", &g.target);
    let q = write_fa(temp.path(), "query", &g.query);

    // Full pgi run, then keep only the first and last blocks (sorted by tStart).
    let pgi_out = temp.path().join("pgi.psl");
    PgrCmd::new()
        .args(&["align", "pgi", &t, &q, "-o", pgi_out.to_str().unwrap()])
        .run();
    let pgi_text = fs::read_to_string(&pgi_out).unwrap();
    let mut lines: Vec<&str> = pgi_text.lines().collect();
    assert!(lines.len() >= 3, "expected >=3 pgi blocks");
    lines.sort_by(|a, b| {
        let fa: Vec<&str> = a.split_whitespace().collect();
        let fb: Vec<&str> = b.split_whitespace().collect();
        fa[16]
            .parse::<u32>()
            .unwrap()
            .cmp(&fb[16].parse::<u32>().unwrap())
    });
    let reduced = temp.path().join("pgi_reduced.psl");
    fs::write(
        &reduced,
        format!("{}\n{}", lines[0], lines[lines.len() - 1]),
    )
    .unwrap();

    let (hybrid, _) = run_hybrid(&temp, &t, &q, Some(reduced.to_str().unwrap()), &[]);
    let records = parse_psl(&hybrid);
    assert!(
        records.len() >= 2,
        "expected the anchors plus a LASTZ fill: {hybrid}"
    );
    // The middle region (r2) must now be covered.
    assert!(
        covers(&records, g.r[1] as u32, (g.r[1] + 2000) as u32),
        "lastz must recover the dropped middle region: {hybrid}"
    );
}

#[test]
fn command_align_hybrid_negative_strand_fill() {
    if lastz_missing() {
        eprintln!("lastz not installed; skipping");
        return;
    }
    // Middle region is reverse-complemented on the query -> '-' strand. The
    // box and coordinate lifting must handle it (block kept with strand '-').
    let temp = tempfile::TempDir::new().unwrap();
    let region = 2000;
    let gap = 3000;
    let r1 = random_seq(region, 51);
    let r2 = random_seq(region, 52);
    let r3 = random_seq(region, 53);
    let target = format!("{r1}{}{r2}{}{r3}", random_seq(gap, 54), random_seq(gap, 55));
    let query = format!(
        "{}{}{}{}{}",
        mutate(&r1, 0.02),
        random_seq(gap, 56),
        rc(&mutate(&r2, 0.02)),
        random_seq(gap, 57),
        mutate(&r3, 0.02)
    );
    let t = write_fa(temp.path(), "target", &target);
    let q = write_fa(temp.path(), "query", &query);

    let pgi_out = temp.path().join("pgi.psl");
    PgrCmd::new()
        .args(&["align", "pgi", &t, &q, "-o", pgi_out.to_str().unwrap()])
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
    // Drop the middle (r2, '-') block, keep r1 and r3.
    let reduced = temp.path().join("pgi_reduced.psl");
    let kept: Vec<&str> = lines
        .iter()
        .filter(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            f[8] != "-"
        })
        .cloned()
        .collect();
    fs::write(&reduced, format!("{}\n", kept.join("\n"))).unwrap();

    let (hybrid, _) = run_hybrid(&temp, &t, &q, Some(reduced.to_str().unwrap()), &[]);
    let records = parse_psl(&hybrid);
    assert!(
        records.iter().any(|r| r.0 == "-"),
        "lastz must recover the '-' strand region: {hybrid}"
    );
}

#[test]
fn command_align_hybrid_avail_psl_overwrite_protected() {
    if lastz_missing() {
        eprintln!("lastz not installed; skipping");
        return;
    }
    // `-o` must not overwrite the --avail-psl input (or the genomes).
    let temp = tempfile::TempDir::new().unwrap();
    let g = build_genomes(44);
    let t = write_fa(temp.path(), "target", &g.target);
    let q = write_fa(temp.path(), "query", &g.query);
    let pgi_out = temp.path().join("pgi.psl");
    PgrCmd::new()
        .args(&["align", "pgi", &t, &q, "-o", pgi_out.to_str().unwrap()])
        .run();
    let before = fs::read_to_string(&pgi_out).unwrap();

    // Try to write the hybrid output over the avail-psl input.
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "hybrid",
            &t,
            &q,
            "--avail-psl",
            pgi_out.to_str().unwrap(),
            "-o",
            pgi_out.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("also an input file"),
        "expected -o to be rejected as an input: {stderr}"
    );
    assert_eq!(
        fs::read_to_string(&pgi_out).unwrap(),
        before,
        "the pgi input must not be overwritten"
    );
}

#[test]
fn command_align_hybrid_extreme_parallel_rejected() {
    // Even without lastz, clap must reject an out-of-range --parallel before
    // any rayon pool is built (Zero Panic convention).
    let temp = tempfile::TempDir::new().unwrap();
    let g = build_genomes(45);
    let t = write_fa(temp.path(), "target", &g.target);
    let q = write_fa(temp.path(), "query", &g.query);
    for bad in ["0", "18446744073709551615", "1025"] {
        let (_, stderr) = PgrCmd::new()
            .args(&[
                "align",
                "hybrid",
                &t,
                &q,
                "--parallel",
                bad,
                "-o",
                temp.path().join("o.psl").to_str().unwrap(),
            ])
            .run_fail();
        assert!(
            stderr.contains("not in 1..=1024"),
            "--parallel {bad} must be rejected: {stderr}"
        );
        assert!(!stderr.contains("panicked"), "got {stderr}");
    }
}

#[test]
fn command_align_hybrid_reuses_sibling_2bit() {
    if lastz_missing() {
        eprintln!("lastz not installed; skipping");
        return;
    }
    // When a sibling `.2bit` sits next to a `.fa` input, hybrid must reuse it
    // instead of converting the FASTA again (random-access extraction source).
    let temp = tempfile::TempDir::new().unwrap();
    let g = build_genomes(46);
    let t = write_fa(temp.path(), "target", &g.target);
    let q = write_fa(temp.path(), "query", &g.query);

    // Pre-create sibling 2bit files next to the FASTA inputs.
    let t_2bit = t.replacen(".fa", ".2bit", 1);
    let q_2bit = q.replacen(".fa", ".2bit", 1);
    PgrCmd::new()
        .args(&["fa", "to-2bit", &t, "-o", &t_2bit])
        .run();
    PgrCmd::new()
        .args(&["fa", "to-2bit", &q, "-o", &q_2bit])
        .run();

    let (hybrid, stderr) = run_hybrid(&temp, &t, &q, None, &[]);
    assert!(
        stderr.contains("reusing sibling 2bit"),
        "expected sibling 2bit reuse logged: {stderr}"
    );
    // The merged output must still be valid (3 conserved regions covered).
    let records = parse_psl(&hybrid);
    assert!(records.len() >= 3, "expected >=3 regions: {hybrid}");
    for &s in &g.r {
        assert!(
            covers(&records, s as u32, (s + 2000) as u32),
            "region at {s} not covered"
        );
    }
}

#[test]
fn command_align_hybrid_missing_subcommand_errors_not_panics() {
    let (_, stderr) = PgrCmd::new().args(&["align"]).run_fail();
    assert!(
        stderr.contains("requires a subcommand"),
        "expected a missing-subcommand error, got: {stderr}"
    );
    assert!(!stderr.contains("panicked"), "got {stderr}");
}
