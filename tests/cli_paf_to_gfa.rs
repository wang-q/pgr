#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::collections::BTreeMap;
use std::path::PathBuf;

/// Return the absolute path to a fixture in `tests/paf/input`.
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/paf/input")
        .join(name)
}

// ── paf to-gfa (V4b local graph from POA MSA) ────────────────────

#[test]
fn command_paf_to_gfa_identical() {
    // Two identical sequences -> unchopped to a single 10-bp segment, no
    // edges, 2 paths traversing that one segment.
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("AB.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    assert!(
        stdout.lines().any(|l| l == "H\tVN:Z:1.0"),
        "missing GFA H header line"
    );

    let s_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("S\t")).collect();
    let l_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("L\t")).collect();
    let p_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("P\t")).collect();

    assert_eq!(s_lines.len(), 1, "expected 1 S line, got {}", s_lines.len());
    assert_eq!(
        l_lines.len(),
        0,
        "expected 0 L lines, got {}",
        l_lines.len()
    );
    assert_eq!(
        p_lines.len(),
        2,
        "expected 2 P lines, got {}",
        p_lines.len()
    );

    let s_fields: Vec<&str> = s_lines[0].split('\t').collect();
    assert_eq!(s_fields[1], "1", "segment id should be 1");
    assert_eq!(s_fields[2], "ACGTACGTAC", "segment sequence mismatch");
    assert!(
        s_fields.contains(&"LN:i:10"),
        "missing LN:i:10 tag in S line: {}",
        s_lines[0]
    );

    for p in &p_lines {
        let fields: Vec<&str> = p.split('\t').collect();
        assert_eq!(fields[2], "1+", "path should visit only segment 1: {p}");
        assert!(fields[3].is_empty(), "path should have no overlaps: {p}");
    }
}

#[test]
fn command_paf_to_gfa_with_snp() {
    // B = ACGTACGTAC, A = ACGTACGTAC, C = ACGTTCGTAC (SNP at pos 4).
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("ABC_snp.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let s_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("S\t")).collect();
    let l_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("L\t")).collect();
    let p_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("P\t")).collect();

    assert_eq!(
        s_lines.len(),
        4,
        "expected 4 S lines, got {}",
        s_lines.len()
    );
    assert_eq!(
        l_lines.len(),
        4,
        "expected 4 L lines, got {}",
        l_lines.len()
    );
    assert_eq!(
        p_lines.len(),
        3,
        "expected 3 P lines, got {}",
        p_lines.len()
    );

    let mut seg_seq: std::collections::HashMap<&str, &str> = std::collections::HashMap::new();
    for s in &s_lines {
        let f: Vec<&str> = s.split('\t').collect();
        seg_seq.insert(f[1], f[2]);
    }
    assert_eq!(seg_seq.get("1"), Some(&"ACGT"));
    assert_eq!(seg_seq.get("2"), Some(&"A"));
    assert_eq!(seg_seq.get("3"), Some(&"T"));
    assert_eq!(seg_seq.get("4"), Some(&"CGTAC"));

    let b_path: Vec<&str> = p_lines
        .iter()
        .find(|p| p.starts_with("P\tB\t"))
        .unwrap()
        .split('\t')
        .nth(2)
        .unwrap()
        .split(',')
        .collect();
    let c_path: Vec<&str> = p_lines
        .iter()
        .find(|p| p.starts_with("P\tC\t"))
        .unwrap()
        .split('\t')
        .nth(2)
        .unwrap()
        .split(',')
        .collect();

    assert_eq!(
        b_path.len(),
        c_path.len(),
        "B and C paths should have equal length (gap-free SNP)"
    );
    let diffs: usize = b_path
        .iter()
        .zip(c_path.iter())
        .filter(|(a, b)| a != b)
        .count();
    assert_eq!(
        diffs, 1,
        "B and C paths should differ at 1 segment (SNP), got {diffs}"
    );

    assert!(
        b_path.iter().any(|s| s.starts_with("2+")),
        "B should traverse segment 2 (A allele): {b_path:?}"
    );
    assert!(
        c_path.iter().any(|s| s.starts_with("3+")),
        "C should traverse segment 3 (T allele): {c_path:?}"
    );
}

#[test]
fn command_paf_to_gfa_crush() {
    // Same setup as with_snp, but with --crush.
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("ABC_snp.tsv").to_str().unwrap(),
            "--crush",
        ])
        .stdin(paf)
        .run();

    let s_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("S\t")).collect();
    let l_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("L\t")).collect();
    let p_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("P\t")).collect();

    assert_eq!(
        s_lines.len(),
        3,
        "expected 3 S lines after crush, got {}",
        s_lines.len()
    );
    assert_eq!(
        l_lines.len(),
        2,
        "expected 2 L lines after crush, got {}",
        l_lines.len()
    );
    assert_eq!(
        p_lines.len(),
        3,
        "expected 3 P lines, got {}",
        p_lines.len()
    );

    let has_t_seg = s_lines.iter().any(|s| s.split('\t').nth(2) == Some("T"));
    assert!(
        !has_t_seg,
        "T allele segment should be crushed out: {s_lines:?}"
    );

    let paths: Vec<&str> = p_lines
        .iter()
        .map(|p| p.split('\t').nth(2).unwrap())
        .collect();
    let first = paths[0];
    assert!(
        paths.iter().all(|p| *p == first),
        "all paths should be identical after crush: {paths:?}"
    );
}

// ── GFA path spelling round-trip ─────────────────────────────────

/// Parse GFA and spell each path's sequence by concatenating visited
/// segments (reverse-complementing for `-` steps).
fn spell_gfa_paths(gfa: &str) -> Vec<(String, String)> {
    use std::collections::HashMap;
    let mut seg_seq: HashMap<String, String> = HashMap::new();
    let mut paths: Vec<(String, Vec<(String, char)>)> = Vec::new();

    for line in gfa.lines() {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields[0] == "S" && fields.len() >= 3 {
            seg_seq.insert(fields[1].to_string(), fields[2].to_string());
        } else if fields[0] == "P" && fields.len() >= 3 {
            let name = fields[1].to_string();
            let steps: Vec<(String, char)> = fields[2]
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    let orient = s.chars().last().unwrap_or('+');
                    let id = s.trim_end_matches(['+', '-']).to_string();
                    (id, orient)
                })
                .collect();
            paths.push((name, steps));
        }
    }

    paths
        .into_iter()
        .map(|(name, steps)| {
            let mut spelled = String::new();
            for (id, orient) in steps {
                let seq = seg_seq.get(&id).cloned().unwrap_or_default();
                if orient == '-' {
                    spelled.push_str(&revcomp(&seq));
                } else {
                    spelled.push_str(&seq);
                }
            }
            (name, spelled)
        })
        .collect()
}

/// Reverse-complement a DNA string (ACGTN, case-insensitive).
fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            other => other,
        })
        .collect()
}

#[test]
fn command_paf_to_gfa_roundtrip_identical() {
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("AB.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let spelled: BTreeMap<String, String> = spell_gfa_paths(&stdout).into_iter().collect();
    assert_eq!(
        spelled.get("A"),
        Some(&"ACGTACGTAC".to_string()),
        "path A should spell back the original sequence"
    );
    assert_eq!(
        spelled.get("B"),
        Some(&"ACGTACGTAC".to_string()),
        "path B should spell back the original sequence"
    );
}

#[test]
fn command_paf_to_gfa_roundtrip_snp_bubble() {
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("ABC_snp.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let spelled: BTreeMap<String, String> = spell_gfa_paths(&stdout).into_iter().collect();
    assert_eq!(
        spelled.get("A"),
        Some(&"ACGTACGTAC".to_string()),
        "path A should spell ACGTACGTAC"
    );
    assert_eq!(
        spelled.get("B"),
        Some(&"ACGTACGTAC".to_string()),
        "path B should spell ACGTACGTAC"
    );
    assert_eq!(
        spelled.get("C"),
        Some(&"ACGTTCGTAC".to_string()),
        "path C should spell ACGTTCGTAC (SNP allele preserved)"
    );
}

#[test]
fn command_paf_to_gfa_roundtrip_indel_bubble() {
    // B = ACGTACGTAC (10bp), A = ACGTACGTAC, C = ACGTACGGGTAC (12bp insertion).
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t12\t0\t12\t+\tA\t10\t0\t10\t10\t12\t255\tcg:Z:6=2I4=\n";
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("ABC_ins2.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let spelled: BTreeMap<String, String> = spell_gfa_paths(&stdout).into_iter().collect();
    assert_eq!(
        spelled.get("A"),
        Some(&"ACGTACGTAC".to_string()),
        "path A should spell ACGTACGTAC (no insertion)"
    );
    assert_eq!(
        spelled.get("B"),
        Some(&"ACGTACGTAC".to_string()),
        "path B should spell ACGTACGTAC (no insertion)"
    );
    assert_eq!(
        spelled.get("C"),
        Some(&"ACGTACGGGTAC".to_string()),
        "path C should spell ACGTACGGGTAC (2bp insertion preserved)"
    );
}

#[test]
fn command_paf_to_gfa_lowercase_roundtrip() {
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "--transitive",
            "-f",
            fixture("ABC_lower.tsv").to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let spelled: BTreeMap<String, String> = spell_gfa_paths(&stdout).into_iter().collect();
    assert_eq!(
        spelled.get("A"),
        Some(&"acgtacgtac".to_string()),
        "path A should spell back the lowercase input sequence"
    );
    assert_eq!(
        spelled.get("B"),
        Some(&"acgtacgtac".to_string()),
        "path B should spell back the lowercase input sequence"
    );
    assert_eq!(
        spelled.get("C"),
        Some(&"acgtttcgtac".to_string()),
        "path C should spell back the lowercase input sequence (SNP allele preserved)"
    );
}
