#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── paf graph-report (V6 graph quality) ───────────────────────

fn write_temp_fasta(dir: &std::path::Path, name: &str, records: &[(&str, &str)]) -> String {
    use std::fs;
    let mut content = String::new();
    for (rec_name, seq) in records {
        content.push('>');
        content.push_str(rec_name);
        content.push('\n');
        content.push_str(seq);
        content.push('\n');
    }
    let path = dir.join(format!("{name}.fa"));
    fs::write(&path, content).unwrap();
    path.to_string_lossy().into_owned()
}

#[test]
fn command_paf_graph_report_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "graph-report", "--help"]).run();
    assert!(stdout.contains("Reports coarse GFA topology metrics"));
    assert!(stdout.contains("--min-var-len"));
    assert!(stdout.contains("--fasta"));
}

#[test]
fn command_paf_graph_report_basic_forward() {
    // A and B share a 100bp alignment → one shared node + trailing novel segments.
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_temp_fasta(
        temp.path(),
        "basic",
        &[("A", &"ACGT".repeat(25)), ("B", &"TGCA".repeat(25))],
    );
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph-report", "stdin", "-f", &fa])
        .stdin(paf)
        .run();

    // Parse TSV output into a map.
    let metrics: std::collections::HashMap<&str, String> = stdout
        .lines()
        .filter_map(|l| {
            let mut it = l.split('\t');
            let k = it.next()?;
            let v = it.next()?;
            Some((k, v.to_string()))
        })
        .collect();

    assert_eq!(metrics["paths"], "2", "expected 2 paths (A, B)");
    // At least one shared node between A and B → reused_nodes_cross_path >= 1.
    let reused: usize = metrics["reused_nodes_cross_path"].parse().unwrap();
    assert!(
        reused >= 1,
        "expected >= 1 cross-path reused node, got {reused}"
    );
    // One connected component (A and B share a node).
    assert_eq!(metrics["components"], "1", "expected 1 component");
    // path_steps >= 2 (each path has at least the shared node + novel tail).
    let path_steps: usize = metrics["path_steps"].parse().unwrap();
    assert!(
        path_steps >= 2,
        "expected >= 2 path_steps, got {path_steps}"
    );
}

#[test]
fn command_paf_graph_report_split_at_large_indel() {
    // 50M 200I 50M: 200I >= 100 → split. B has an insertion (novel node in B path).
    let paf = "A\t300\t0\t100\t+\tB\t300\t0\t300\t95\t300\t255\tcg:Z:50M200I50M\n";
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_temp_fasta(
        temp.path(),
        "split",
        &[("A", &"A".repeat(300)), ("B", &"G".repeat(300))],
    );
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "graph-report",
            "stdin",
            "-f",
            &fa,
            "--min-var-len",
            "100",
        ])
        .stdin(paf)
        .run();

    let metrics: std::collections::HashMap<&str, String> = stdout
        .lines()
        .filter_map(|l| {
            let mut it = l.split('\t');
            let k = it.next()?;
            let v = it.next()?;
            Some((k, v.to_string()))
        })
        .collect();

    // B has a novel insertion node → singleton_nodes >= 1 (only B visits it).
    let singletons: usize = metrics["singleton_nodes"].parse().unwrap();
    assert!(
        singletons >= 1,
        "expected >= 1 singleton (novel insertion), got {singletons}"
    );
    // segments >= 3: shared-left + novel-insertion + shared-right (plus trailing novel).
    let segments: usize = metrics["segments"].parse().unwrap();
    assert!(
        segments >= 3,
        "expected >= 3 segments after split, got {segments}"
    );
}

#[test]
fn command_paf_graph_report_small_indel_no_split() {
    // 50M 30I 50M: 30I < 100 → no split. A and B share exactly one aligned node.
    let paf = "A\t200\t0\t130\t+\tB\t200\t0\t160\t95\t160\t255\tcg:Z:50M30I50M\n";
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_temp_fasta(
        temp.path(),
        "nosplit",
        &[("A", &"A".repeat(200)), ("B", &"G".repeat(200))],
    );
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "graph-report",
            "stdin",
            "-f",
            &fa,
            "--min-var-len",
            "100",
        ])
        .stdin(paf)
        .run();

    let metrics: std::collections::HashMap<&str, String> = stdout
        .lines()
        .filter_map(|l| {
            let mut it = l.split('\t');
            let k = it.next()?;
            let v = it.next()?;
            Some((k, v.to_string()))
        })
        .collect();

    // No split → exactly 1 cross-path reused node (the single aligned segment).
    let reused: usize = metrics["reused_nodes_cross_path"].parse().unwrap();
    assert_eq!(
        reused, 1,
        "expected exactly 1 cross-path reused node (no split), got {reused}"
    );
}

#[test]
fn command_paf_graph_report_threshold_comparison() {
    // Same alignment, different thresholds: stricter threshold yields fewer splits.
    // 50M 200I 50M
    let paf = "A\t300\t0\t100\t+\tB\t300\t0\t300\t95\t300\t255\tcg:Z:50M200I50M\n";
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_temp_fasta(
        temp.path(),
        "thr",
        &[("A", &"A".repeat(300)), ("B", &"G".repeat(300))],
    );

    // Threshold 100: 200I >= 100 → split.
    let (out_strict, _) = PgrCmd::new()
        .args(&[
            "paf",
            "graph-report",
            "stdin",
            "-f",
            &fa,
            "--min-var-len",
            "100",
        ])
        .stdin(paf)
        .run();
    // Threshold 500: 200I < 500 → no split.
    let (out_loose, _) = PgrCmd::new()
        .args(&[
            "paf",
            "graph-report",
            "stdin",
            "-f",
            &fa,
            "--min-var-len",
            "500",
        ])
        .stdin(paf)
        .run();

    let parse_seg = |s: &str| -> usize {
        s.lines()
            .find_map(|l| {
                let mut it = l.split('\t');
                if it.next() == Some("segments") {
                    it.next().and_then(|v| v.parse().ok())
                } else {
                    None
                }
            })
            .unwrap_or(0)
    };
    let seg_strict = parse_seg(&out_strict);
    let seg_loose = parse_seg(&out_loose);
    assert!(
        seg_strict > seg_loose,
        "stricter threshold should yield more segments: strict={seg_strict} loose={seg_loose}"
    );
}

#[test]
fn command_paf_graph_report_no_alignment() {
    // Empty PAF → each sequence becomes one isolated novel node.
    let paf = "";
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_temp_fasta(
        temp.path(),
        "empty",
        &[("A", &"ACGT".repeat(10)), ("B", &"TGCA".repeat(10))],
    );
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph-report", "stdin", "-f", &fa])
        .stdin(paf)
        .run();

    let metrics: std::collections::HashMap<&str, String> = stdout
        .lines()
        .filter_map(|l| {
            let mut it = l.split('\t');
            let k = it.next()?;
            let v = it.next()?;
            Some((k, v.to_string()))
        })
        .collect();

    assert_eq!(metrics["segments"], "2", "expected 2 isolated novel nodes");
    assert_eq!(metrics["paths"], "2");
    assert_eq!(metrics["links"], "0", "no alignments → no edges");
    assert_eq!(
        metrics["components"], "2",
        "2 isolated nodes → 2 components"
    );
    assert_eq!(metrics["isolated_nodes"], "2");
    assert_eq!(metrics["tips"], "0");
}
