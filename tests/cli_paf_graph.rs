#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── paf graph (V4a coarse GFA induction) ──────────────────────

fn write_temp_fasta(path: &str, records: &[(&str, &str)]) {
    use std::fs;
    let mut content = String::new();
    for (name, seq) in records {
        content.push('>');
        content.push_str(name);
        content.push('\n');
        content.push_str(seq);
        content.push('\n');
    }
    fs::write(path, content).unwrap();
}

#[test]
fn command_paf_graph_help() {
    let (stdout, _) = PgrCmd::new().args(&["paf", "graph", "--help"]).run();
    assert!(stdout.contains("Induces a coarse GFA graph"));
    assert!(stdout.contains("--min-var-len"));
    assert!(stdout.contains("--fasta"));
}

#[test]
fn command_paf_graph_basic_forward() {
    // A and B share a 100bp alignment → one shared node + trailing novel segments.
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let fa = "/tmp/pgr_graph_basic.fa";
    write_temp_fasta(fa, &[("A", &"ACGT".repeat(25)), ("B", &"TGCA".repeat(25))]);
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa])
        .stdin(paf)
        .run();
    // At least one S line, one P line for each sequence.
    let s_count = stdout.lines().filter(|l| l.starts_with("S\t")).count();
    let p_count = stdout.lines().filter(|l| l.starts_with("P\t")).count();
    assert!(s_count >= 1, "expected >= 1 S line, got {s_count}");
    assert_eq!(p_count, 2, "expected 2 P lines (A, B), got {p_count}");
    assert!(stdout.contains("\nP\tA\t"), "missing P line for A");
    assert!(stdout.contains("\nP\tB\t"), "missing P line for B");
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_split_at_large_indel() {
    // 50M 200I 50M: 200I >= 100 → split. B has an insertion (novel node in B path).
    let paf = "A\t300\t0\t100\t+\tB\t300\t0\t300\t95\t300\t255\tcg:Z:50M200I50M\n";
    let fa = "/tmp/pgr_graph_split.fa";
    write_temp_fasta(fa, &[("A", &"A".repeat(300)), ("B", &"G".repeat(300))]);
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa, "--min-var-len", "100"])
        .stdin(paf)
        .run();
    // B path should have >= 3 steps (aligned, novel insertion, aligned).
    let b_line = stdout
        .lines()
        .find(|l| l.starts_with("P\tB\t"))
        .expect("missing P line for B");
    // P line format: P\tname\tpath\toverlaps — path is the 3rd field.
    let path_field: &str = b_line.split('\t').nth(2).unwrap();
    let step_count = path_field.split(',').count();
    assert!(
        step_count >= 3,
        "B path should have >= 3 steps (aligned, novel, aligned), got {step_count}: {path_field}"
    );
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_small_indel_no_split() {
    // 50M 30I 50M: 30I < 100 → no split. A and B share exactly one aligned node.
    let paf = "A\t200\t0\t130\t+\tB\t200\t0\t160\t95\t160\t255\tcg:Z:50M30I50M\n";
    let fa = "/tmp/pgr_graph_nosplit.fa";
    write_temp_fasta(fa, &[("A", &"A".repeat(200)), ("B", &"G".repeat(200))]);
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa, "--min-var-len", "100"])
        .stdin(paf)
        .run();
    // Find shared nodes between A and B paths.
    let a_line = stdout
        .lines()
        .find(|l| l.starts_with("P\tA\t"))
        .expect("missing P line for A");
    let b_line = stdout
        .lines()
        .find(|l| l.starts_with("P\tB\t"))
        .expect("missing P line for B");
    let a_steps: Vec<&str> = a_line.split('\t').nth(2).unwrap().split(',').collect();
    let b_steps: Vec<&str> = b_line.split('\t').nth(2).unwrap().split(',').collect();
    // Strip orientation suffix to compare node ids.
    let a_nodes: Vec<&str> = a_steps
        .iter()
        .map(|s| s.trim_end_matches(['+', '-']))
        .collect();
    let b_nodes: Vec<&str> = b_steps
        .iter()
        .map(|s| s.trim_end_matches(['+', '-']))
        .collect();
    let shared: Vec<&str> = a_nodes
        .iter()
        .filter(|n| b_nodes.contains(n))
        .copied()
        .collect();
    assert_eq!(
        shared.len(),
        1,
        "expected exactly 1 shared node (no split), got {shared:?}"
    );
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_reverse_strand() {
    // Reverse strand alignment: query coords flipped, but A and B still share a node.
    let paf = "A\t100\t0\t100\t-\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let fa = "/tmp/pgr_graph_rc.fa";
    write_temp_fasta(fa, &[("A", &"ACGT".repeat(25)), ("B", &"TGCA".repeat(25))]);
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa])
        .stdin(paf)
        .run();
    let a_line = stdout
        .lines()
        .find(|l| l.starts_with("P\tA\t"))
        .expect("missing P line for A");
    let b_line = stdout
        .lines()
        .find(|l| l.starts_with("P\tB\t"))
        .expect("missing P line for B");
    let a_nodes: Vec<&str> = a_line
        .split('\t')
        .nth(2)
        .unwrap()
        .split(',')
        .map(|s| s.trim_end_matches(['+', '-']))
        .collect();
    let b_nodes: Vec<&str> = b_line
        .split('\t')
        .nth(2)
        .unwrap()
        .split(',')
        .map(|s| s.trim_end_matches(['+', '-']))
        .collect();
    let shared: Vec<&str> = a_nodes
        .iter()
        .filter(|n| b_nodes.contains(n))
        .copied()
        .collect();
    assert!(
        !shared.is_empty(),
        "reverse-strand alignment should still produce a shared node"
    );
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_min_var_len_filter() {
    // 50M 150I 50M with --min-var-len 200: 150I < 200 → no split.
    // Same alignment with --min-var-len 100: 150I >= 100 → split.
    let paf = "A\t300\t0\t100\t+\tB\t300\t0\t250\t95\t250\t255\tcg:Z:50M150I50M\n";
    let fa = "/tmp/pgr_graph_filter.fa";
    write_temp_fasta(fa, &[("A", &"A".repeat(300)), ("B", &"G".repeat(300))]);

    // With threshold 200: no split, B path has 1 shared node + trailing novel.
    let (stdout_no_split, _) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa, "--min-var-len", "200"])
        .stdin(paf)
        .run();
    let b_line = stdout_no_split
        .lines()
        .find(|l| l.starts_with("P\tB\t"))
        .expect("missing P line for B");
    let steps_no_split = b_line.split('\t').nth(2).unwrap().split(',').count();

    // With threshold 100: split, B path has >= 3 steps.
    let (stdout_split, _) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa, "--min-var-len", "100"])
        .stdin(paf)
        .run();
    let b_line = stdout_split
        .lines()
        .find(|l| l.starts_with("P\tB\t"))
        .expect("missing P line for B");
    let steps_split = b_line.split('\t').nth(2).unwrap().split(',').count();

    assert!(
        steps_split > steps_no_split,
        "split path ({steps_split}) should have more steps than no-split ({steps_no_split})"
    );
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_missing_fasta_fails() {
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let (_stdout, stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", "/nonexistent/path.fa"])
        .stdin(paf)
        .run_fail();
    // Should fail with a friendly error, not panic.
    assert!(
        stderr.contains("could not open") || stderr.contains("No such file"),
        "expected file-not-found error, got: {stderr}"
    );
}

#[test]
fn command_paf_graph_rgfa_tags() {
    // PAF: query=A, target=B. Target is registered first → B has seq_id 0.
    // Shared aligned node originates from B (target) at offset 0.
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let fa = "/tmp/pgr_graph_rgfa.fa";
    write_temp_fasta(fa, &[("A", &"ACGT".repeat(25)), ("B", &"TGCA".repeat(25))]);
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa])
        .stdin(paf)
        .run();

    // Every S line should carry SN:Z, SO:i, SR:i tags.
    let s_lines: Vec<&str> = stdout.lines().filter(|l| l.starts_with("S\t")).collect();
    assert!(!s_lines.is_empty(), "expected at least one S line");
    for line in &s_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        // S\tid\tseq\tSN:Z:name\tSO:i:pos\tSR:i:0
        assert!(fields.len() >= 6, "S line missing rGFA tags: {line}");
        let tags = &fields[3..];
        assert!(
            tags.iter().any(|t| t.starts_with("SN:Z:")),
            "missing SN:Z tag in S line: {line}"
        );
        assert!(
            tags.iter().any(|t| t.starts_with("SO:i:")),
            "missing SO:i tag in S line: {line}"
        );
        assert!(
            tags.iter().any(|t| *t == "SR:i:0"),
            "missing SR:i:0 tag in S line: {line}"
        );
    }

    // The shared aligned node originates from B (target, seq_id 0) at offset 0.
    let shared_line = s_lines
        .iter()
        .find(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f.iter().any(|t| *t == "SN:Z:B") && f.iter().any(|t| *t == "SO:i:0")
        })
        .expect("missing shared node with SN:Z:B and SO:i:0");
    let _ = shared_line;
    let _ = std::fs::remove_file(fa);
}

#[test]
fn command_paf_graph_rgfa_novel_node_origin() {
    // PAF: query=A, target=B. CIGAR 50M200I50M → A (query) has 200bp insertion.
    // The novel insertion node in A's path spans A:50-250, origin SN:Z:A SO:i:50.
    let paf = "A\t300\t0\t100\t+\tB\t300\t0\t300\t95\t300\t255\tcg:Z:50M200I50M\n";
    let fa = "/tmp/pgr_graph_rgfa_novel.fa";
    write_temp_fasta(fa, &[("A", &"A".repeat(300)), ("B", &"G".repeat(300))]);
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["paf", "graph", "stdin", "-f", fa, "--min-var-len", "100"])
        .stdin(paf)
        .run();

    // Novel insertion node in A's path: SN:Z:A and SO:i:50 (gap 50-250 in A).
    let has_novel_a = stdout
        .lines()
        .any(|l| l.starts_with("S\t") && l.contains("SN:Z:A") && l.contains("SO:i:50"));
    assert!(
        has_novel_a,
        "expected a novel node with SN:Z:A and SO:i:50 (A's 200bp insertion at offset 50)"
    );
    let _ = std::fs::remove_file(fa);
}
