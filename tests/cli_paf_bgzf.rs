#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── Gzipped PAF input ───────────────────────────────────────────

/// Write `content` to `path` as a plain gzip file using flate2.
fn write_gz(path: &std::path::Path, content: &str) {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::fs::File;
    let f = File::create(path).unwrap();
    let mut enc = GzEncoder::new(f, Compression::default());
    std::io::Write::write_all(&mut enc, content.as_bytes()).unwrap();
    enc.finish().unwrap();
}

/// Write `content` to `path` as a BGZF file using noodles-bgzf.
fn write_bgzf(path: &std::path::Path, content: &str) {
    use noodles_bgzf as bgzf;
    let f = std::fs::File::create(path).unwrap();
    let mut writer = bgzf::io::Writer::new(f);
    std::io::Write::write_all(&mut writer, content.as_bytes()).unwrap();
    writer.finish().unwrap();
}

const GZ_PAF: &str = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t50\t+\tB\t100\t50\t100\t45\t50\t255\tcg:Z:50M
";

#[test]
fn command_paf_index_gzipped() {
    let temp = tempfile::TempDir::new().unwrap();
    let path = temp.path().join("gz_index.paf.gz");
    write_gz(&path, GZ_PAF);
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", path.to_str().unwrap()])
        .run();
    assert!(stderr.contains("sequences: 3"));
    assert!(stderr.contains("targets:   1"));
}

#[test]
fn command_paf_index_bgzf() {
    let temp = tempfile::TempDir::new().unwrap();
    let path = temp.path().join("bgzf_index.paf.gz");
    write_bgzf(&path, GZ_PAF);
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", path.to_str().unwrap()])
        .run();
    assert!(stderr.contains("sequences: 3"));
    assert!(stderr.contains("targets:   1"));
}

#[test]
fn command_paf_query_gzipped() {
    let temp = tempfile::TempDir::new().unwrap();
    let path = temp.path().join("gz_query.paf.gz");
    write_gz(&path, GZ_PAF);
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "query", path.to_str().unwrap(), "B:0-100"])
        .run();
    assert!(stderr.contains("Building index"));
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
    assert!(stdout.contains("C\t0\t0\t50\t+\tB"), "C not found");
}

#[test]
fn command_paf_query_bgzf() {
    let temp = tempfile::TempDir::new().unwrap();
    let path = temp.path().join("bgzf_query.paf.gz");
    write_bgzf(&path, GZ_PAF);
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", path.to_str().unwrap(), "B:0-100"])
        .run();
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
}

#[test]
fn command_paf_graph_gzipped() {
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("gz_graph.paf.gz");
    let fa_path = temp.path().join("gz_graph.fa");
    write_gz(
        &paf_path,
        "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n",
    );
    fs::write(&fa_path, ">A\nACGTACGTAC\n>B\nACGTACGTAC\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "paf",
            "graph",
            paf_path.to_str().unwrap(),
            "-f",
            fa_path.to_str().unwrap(),
        ])
        .run();
    assert!(stdout.contains("S\t1\t"), "missing S line");
    assert!(stdout.contains("P\tA\t"), "missing P line for A");
}

#[test]
fn command_paf_index_gz_save_and_query() {
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("gz_persist.paf.gz");
    let idx_path = temp.path().join("gz_persist.paf.idx");
    write_gz(&paf_path, GZ_PAF);
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "index",
            paf_path.to_str().unwrap(),
            "-o",
            idx_path.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("saved to"));
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", idx_path.to_str().unwrap(), "B:0-100"])
        .run();
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
}

// ── Lazy CIGAR loading (BGZF virtual-position) ──────────────────
//
// These tests verify the IndexedReader-based lazy CIGAR path:
// BGZF PAF → build_from_path records vpos per line → CIGAR fetched
// on-demand during query. Compared with plain gzip (in-memory CIGAR)
// and direct vs persisted index, results must be identical.

#[test]
fn command_paf_index_bgzf_lazy_mode_reported() {
    // `paf index` on a BGZF file should report "lazy" mode in stderr.
    let temp = tempfile::TempDir::new().unwrap();
    let path = temp.path().join("bgzf_lazy_mode.paf.gz");
    write_bgzf(&path, GZ_PAF);
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", path.to_str().unwrap()])
        .run();
    assert!(
        stderr.contains("mode:      lazy"),
        "expected lazy mode banner, got: {stderr}"
    );
}

#[test]
fn command_paf_index_bgzf_plain_gzip_not_lazy() {
    // Plain gzip must NOT trigger lazy mode (no virtual positions).
    let temp = tempfile::TempDir::new().unwrap();
    let path = temp.path().join("gz_not_lazy.paf.gz");
    write_gz(&path, GZ_PAF);
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", path.to_str().unwrap()])
        .run();
    assert!(
        !stderr.contains("mode:      lazy"),
        "plain gzip should not report lazy mode, got: {stderr}"
    );
}

#[test]
fn command_paf_query_bgzf_lazy_cigar_resolved() {
    // Querying a BGZF file directly must fetch CIGAR on-demand and
    // produce the same output as the in-memory path.
    let temp = tempfile::TempDir::new().unwrap();
    let path = temp.path().join("bgzf_lazy_query.paf.gz");
    write_bgzf(&path, GZ_PAF);
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "query", path.to_str().unwrap(), "B:0-100"])
        .run();
    assert!(
        stderr.contains("mode: lazy"),
        "expected lazy mode banner on query, got: {stderr}"
    );
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
    assert!(stdout.contains("C\t0\t0\t50\t+\tB"), "C not found");
    // CIGAR tag must be present (proves lazy fetch worked).
    assert!(stdout.contains("cg:Z:100M"), "CIGAR not resolved");
}

#[test]
fn command_paf_bgzf_lazy_persist_and_reload() {
    // End-to-end: BGZF → index → save .paf.idx → load → query.
    // The persisted index must reopen the BGZF file and resolve CIGAR.
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("bgzf_lazy_persist.paf.gz");
    let idx_path = temp.path().join("bgzf_lazy_persist.paf.idx");
    write_bgzf(&paf_path, GZ_PAF);

    // Build + save.
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "index",
            paf_path.to_str().unwrap(),
            "-o",
            idx_path.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("mode:      lazy"), "lazy mode at build");
    assert!(stderr.contains("saved to"));

    // Load + query (lazy CIGAR must be re-attached via reopen_lazy_source).
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "query", idx_path.to_str().unwrap(), "B:0-100"])
        .run();
    assert!(stderr.contains("Loading index"));
    assert!(
        stdout.contains("A\t0\t0\t100\t+\tB"),
        "A not found after reload"
    );
    assert!(
        stdout.contains("cg:Z:100M"),
        "CIGAR not resolved after reload"
    );
}

#[test]
fn command_paf_bgzf_lazy_vs_plain_text_same_result() {
    // Lazy (BGZF) and in-memory (plain text) must produce identical query output.
    use std::fs;
    let temp = tempfile::TempDir::new().unwrap();
    let bgzf_path = temp.path().join("bgzf_vs_plain.paf.gz");
    let plain_path = temp.path().join("bgzf_vs_plain.paf");
    write_bgzf(&bgzf_path, GZ_PAF);
    fs::write(&plain_path, GZ_PAF).unwrap();

    let (lazy_out, _) = PgrCmd::new()
        .args(&["paf", "query", bgzf_path.to_str().unwrap(), "B:0-100"])
        .run();
    let (plain_out, _) = PgrCmd::new()
        .args(&["paf", "query", plain_path.to_str().unwrap(), "B:0-100"])
        .run();
    assert_eq!(
        lazy_out, plain_out,
        "BGZF lazy query output differs from plain text"
    );
}

#[test]
fn command_paf_bgzf_lazy_transitive_bfs() {
    // Transitive BFS over a BGZF index must resolve CIGARs across multiple hops.
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t100\t+\tA\t100\t0\t100\t90\t100\t255\tcg:Z:100M
";
    let temp = tempfile::TempDir::new().unwrap();
    let path = temp.path().join("bgzf_lazy_bfs.paf.gz");
    write_bgzf(&path, paf);
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            path.to_str().unwrap(),
            "B:0-100",
            "--transitive",
        ])
        .run();
    assert!(stderr.contains("mode: lazy"));
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "1-hop A not found");
    assert!(stdout.contains("C\t0\t0\t100\t+\tA"), "2-hop C not found");
}
