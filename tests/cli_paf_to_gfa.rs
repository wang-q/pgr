#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── helper ───────────────────────────────────────────────────────

/// Write `content` to `<dir>/<name>.fa`, then BGZF-compress it via `pgr fa gz`
/// (which also creates the .gzi index required for random access). Returns
/// the path to the produced `.fa.gz` file.
fn write_bgzf_fa(dir: &std::path::Path, name: &str, content: &str) -> String {
    use std::fs;
    let fa_path = dir.join(format!("{name}.fa"));
    fs::write(&fa_path, content).unwrap();
    let fa_str = fa_path.to_string_lossy().into_owned();
    let (out, _) = PgrCmd::new().args(&["fa", "gz", &fa_str]).run();
    let _ = out;
    let gz_path = format!("{fa_str}.gz");
    assert!(
        std::path::Path::new(&gz_path).exists(),
        "pgr fa gz failed to produce {gz_path}"
    );
    gz_path
}

// ── paf to-gfa (V4b local graph from POA MSA) ────────────────────

#[test]
fn command_paf_to_gfa_identical() {
    use std::fs;
    // Two identical sequences -> unchopped to a single 10-bp segment, no
    // edges, 2 paths traversing that one segment.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    // GFA v1.0 header present.
    assert!(
        stdout.lines().any(|l| l == "H\tVN:Z:1.0"),
        "missing GFA H header line"
    );

    let s_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("S\t")).collect();
    let l_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("L\t")).collect();
    let p_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("P\t")).collect();

    // Unchopping collapses the 10 identical bases into one segment.
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

    // The single segment carries the full 10-bp sequence + LN tag.
    let s_fields: Vec<&str> = s_lines[0].split('\t').collect();
    assert_eq!(s_fields[1], "1", "segment id should be 1");
    assert_eq!(s_fields[2], "ACGTACGTAC", "segment sequence mismatch");
    assert!(
        s_fields.contains(&"LN:i:10"),
        "missing LN:i:10 tag in S line: {}",
        s_lines[0]
    );

    // Each path visits exactly one node (segment 1), zero overlaps.
    for p in &p_lines {
        let fields: Vec<&str> = p.split('\t').collect();
        assert_eq!(fields[2], "1+", "path should visit only segment 1: {p}");
        assert!(fields[3].is_empty(), "path should have no overlaps: {p}");
    }
}

#[test]
fn command_paf_to_gfa_with_snp() {
    use std::fs;
    // B = ACGTACGTAC (target)
    // A = ACGTACGTAC (identical to B)
    // C = ACGTTCGTAC (SNP at pos 4: A->T)
    // After unchopping: 4 segments (ACGT, A, T, CGTAC), 4 edges, 3 paths.
    // The SNP forms a bubble: seg2(A) and seg3(T) share in={1}, out={4}.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nACGTTCGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
        ])
        .stdin(paf)
        .run();

    let s_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("S\t")).collect();
    let l_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("L\t")).collect();
    let p_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("P\t")).collect();

    // 4 segments: ACGT, A, T, CGTAC.
    assert_eq!(
        s_lines.len(),
        4,
        "expected 4 S lines, got {}",
        s_lines.len()
    );
    // 4 edges: 1->2, 1->3, 2->4, 3->4.
    assert_eq!(
        l_lines.len(),
        4,
        "expected 4 L lines, got {}",
        l_lines.len()
    );
    // 3 paths (B, A, C).
    assert_eq!(
        p_lines.len(),
        3,
        "expected 3 P lines, got {}",
        p_lines.len()
    );

    // Collect segment sequences (id -> seq).
    let mut seg_seq: std::collections::HashMap<&str, &str> = std::collections::HashMap::new();
    for s in &s_lines {
        let f: Vec<&str> = s.split('\t').collect();
        seg_seq.insert(f[1], f[2]);
    }
    assert_eq!(seg_seq.get("1"), Some(&"ACGT"));
    assert_eq!(seg_seq.get("2"), Some(&"A"));
    assert_eq!(seg_seq.get("3"), Some(&"T"));
    assert_eq!(seg_seq.get("4"), Some(&"CGTAC"));

    // B and C paths should differ at exactly one segment (the SNP), gap-free.
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

    // B's path should go through the A allele (seg 2), C's through T (seg 3).
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
    use std::fs;
    // Same setup as command_paf_to_gfa_with_snp, but with --crush.
    // The SNP bubble (seg2=A, seg3=T) collapses to one segment (A, the
    // higher-weight allele: B+A=2 vs C=1). Paths through T are rewritten
    // to A, losing base-level ALT info.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nACGTTCGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
            "--crush",
        ])
        .stdin(paf)
        .run();

    let s_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("S\t")).collect();
    let l_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("L\t")).collect();
    let p_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("P\t")).collect();

    // Crushed: 3 segments (ACGT, A, CGTAC), 2 edges, 3 identical paths.
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

    // No 'T' segment should remain (the SNP ALT was crushed out).
    let has_t_seg = s_lines.iter().any(|s| s.split('\t').nth(2) == Some("T"));
    assert!(
        !has_t_seg,
        "T allele segment should be crushed out: {s_lines:?}"
    );

    // All three paths should be identical (the ALT path was rewritten to REF).
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
// Borrowed from impg `test_graph_poa.rs::assert_gfa_paths_match_records`:
// parse GFA P lines, concatenate segment sequences along each path
// (respecting orientation), and verify the spelled sequence matches the
// input FASTA record byte-for-byte. This is the strongest correctness
// invariant for graph construction — every path must reconstruct its
// original sequence.

/// Parse GFA and spell each path's sequence by concatenating visited
/// segments (reverse-complementing for `-` steps). Returns
/// `Vec<(path_name, spelled_sequence)>`.
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
    use std::collections::BTreeMap;
    use std::fs;
    // Two identical sequences -> one segment, two paths. Each path must
    // spell back the original 10-bp sequence.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::collections::BTreeMap;
    use std::fs;
    // B = ACGTACGTAC, A = ACGTACGTAC, C = ACGTTCGTAC (SNP at pos 4).
    // The SNP forms a bubble; each path must still spell its own sequence.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nACGTTCGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::collections::BTreeMap;
    use std::fs;
    // B = ACGTACGTAC (10bp)
    // A = ACGTACGTAC (10bp, identical to B)
    // C = ACGTACGGGTAC (12bp, 2bp insertion after pos 6)
    // C-A alignment: 6= 2I 4= (C has 2bp insertion relative to A)
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t12\t0\t12\t+\tA\t10\t0\t10\t10\t12\t255\tcg:Z:6=2I4=\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nACGTACGTAC\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nACGTACGTAC\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nACGTACGGGTAC\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
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
    use std::collections::BTreeMap;
    use std::fs;
    // Inspired by seqwish `test/HLA/01_seqwish.t` masked-sequence test.
    // Unlike seqwish, pgr's POA does NOT normalize case (poa/align.rs
    // compares bases by strict equality), so an all-lowercase input yields
    // a different topology (more segments: lowercase a/c/g/t are distinct
    // nodes from each other just like uppercase A/C/G/T, but the SNP
    // bubble shape differs). We verify the round-trip invariant that still
    // holds: each path spells back its original (lowercase) sequence.
    let temp = tempfile::TempDir::new().unwrap();
    let dir = temp.path();
    let paf = "\
A\t10\t0\t10\t+\tB\t10\t0\t10\t10\t10\t255\tcg:Z:10=\n\
C\t10\t0\t10\t+\tA\t10\t0\t10\t9\t10\t255\tcg:Z:10M\n";
    let a_fa = write_bgzf_fa(dir, "A", ">A\nacgtacgtac\n");
    let b_fa = write_bgzf_fa(dir, "B", ">B\nacgtacgtac\n");
    let c_fa = write_bgzf_fa(dir, "C", ">C\nacgtttcgtac\n");
    let tsv = dir.join("in.tsv");
    fs::write(&tsv, format!("A\t{a_fa}\nB\t{b_fa}\nC\t{c_fa}\n")).unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "to-gfa",
            "stdin",
            "B:0-10",
            "-t",
            "-f",
            tsv.to_str().unwrap(),
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
