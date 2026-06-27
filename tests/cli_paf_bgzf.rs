#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── Gzipped PAF input ───────────────────────────────────────────

/// Write `content` to `path` as a plain gzip file using flate2.
fn write_gz(path: &str, content: &str) {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::fs::File;
    let f = File::create(path).unwrap();
    let mut enc = GzEncoder::new(f, Compression::default());
    std::io::Write::write_all(&mut enc, content.as_bytes()).unwrap();
    enc.finish().unwrap();
}

/// Write `content` to `path` as a BGZF file using noodles-bgzf.
fn write_bgzf(path: &str, content: &str) {
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
    let path = "/tmp/pgr_cli_test_gz_index.paf.gz";
    write_gz(path, GZ_PAF);
    let (_, stderr) = PgrCmd::new().args(&["paf", "index", path]).run();
    assert!(stderr.contains("sequences: 3"));
    assert!(stderr.contains("targets:   1"));
    let _ = std::fs::remove_file(path);
}

#[test]
fn command_paf_index_bgzf() {
    let path = "/tmp/pgr_cli_test_bgzf_index.paf.gz";
    write_bgzf(path, GZ_PAF);
    let (_, stderr) = PgrCmd::new().args(&["paf", "index", path]).run();
    assert!(stderr.contains("sequences: 3"));
    assert!(stderr.contains("targets:   1"));
    let _ = std::fs::remove_file(path);
}

#[test]
fn command_paf_query_gzipped() {
    let path = "/tmp/pgr_cli_test_gz_query.paf.gz";
    write_gz(path, GZ_PAF);
    let (stdout, stderr) = PgrCmd::new().args(&["paf", "query", path, "B:0-100"]).run();
    assert!(stderr.contains("Building index"));
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
    assert!(stdout.contains("C\t0\t0\t50\t+\tB"), "C not found");
    let _ = std::fs::remove_file(path);
}

#[test]
fn command_paf_query_bgzf() {
    let path = "/tmp/pgr_cli_test_bgzf_query.paf.gz";
    write_bgzf(path, GZ_PAF);
    let (stdout, _) = PgrCmd::new().args(&["paf", "query", path, "B:0-100"]).run();
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
    let _ = std::fs::remove_file(path);
}

#[test]
fn command_paf_graph_gzipped() {
    let paf_path = "/tmp/pgr_cli_test_gz_graph.paf.gz";
    let fa_path = "/tmp/pgr_cli_test_gz_graph.fa";
    write_gz(
        paf_path,
        "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n",
    );
    std::fs::write(fa_path, ">A\nACGTACGTAC\n>B\nACGTACGTAC\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "graph", paf_path, "-f", fa_path])
        .run();
    assert!(stdout.contains("S\t1\t"), "missing S line");
    assert!(stdout.contains("P\tA\t"), "missing P line for A");
    let _ = std::fs::remove_file(paf_path);
    let _ = std::fs::remove_file(fa_path);
}

#[test]
fn command_paf_index_gz_save_and_query() {
    let paf_path = "/tmp/pgr_cli_test_gz_persist.paf.gz";
    let idx_path = "/tmp/pgr_cli_test_gz_persist.paf.idx";
    write_gz(paf_path, GZ_PAF);
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", paf_path, "-o", idx_path])
        .run();
    assert!(stderr.contains("saved to"));
    let (stdout, _) = PgrCmd::new()
        .args(&["paf", "query", idx_path, "B:0-100"])
        .run();
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
    let _ = std::fs::remove_file(paf_path);
    let _ = std::fs::remove_file(idx_path);
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
    let path = "/tmp/pgr_cli_test_bgzf_lazy_mode.paf.gz";
    write_bgzf(path, GZ_PAF);
    let (_, stderr) = PgrCmd::new().args(&["paf", "index", path]).run();
    assert!(
        stderr.contains("mode:      lazy"),
        "expected lazy mode banner, got: {stderr}"
    );
    let _ = std::fs::remove_file(path);
}

#[test]
fn command_paf_index_bgzf_plain_gzip_not_lazy() {
    // Plain gzip must NOT trigger lazy mode (no virtual positions).
    let path = "/tmp/pgr_cli_test_gz_not_lazy.paf.gz";
    write_gz(path, GZ_PAF);
    let (_, stderr) = PgrCmd::new().args(&["paf", "index", path]).run();
    assert!(
        !stderr.contains("mode:      lazy"),
        "plain gzip should not report lazy mode, got: {stderr}"
    );
    let _ = std::fs::remove_file(path);
}

#[test]
fn command_paf_query_bgzf_lazy_cigar_resolved() {
    // Querying a BGZF file directly must fetch CIGAR on-demand and
    // produce the same output as the in-memory path.
    let path = "/tmp/pgr_cli_test_bgzf_lazy_query.paf.gz";
    write_bgzf(path, GZ_PAF);
    let (stdout, stderr) = PgrCmd::new().args(&["paf", "query", path, "B:0-100"]).run();
    assert!(
        stderr.contains("mode: lazy"),
        "expected lazy mode banner on query, got: {stderr}"
    );
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "A not found");
    assert!(stdout.contains("C\t0\t0\t50\t+\tB"), "C not found");
    // CIGAR tag must be present (proves lazy fetch worked).
    assert!(stdout.contains("cg:Z:100M"), "CIGAR not resolved");
    let _ = std::fs::remove_file(path);
}

#[test]
fn command_paf_bgzf_lazy_persist_and_reload() {
    // End-to-end: BGZF → index → save .paf.idx → load → query.
    // The persisted index must reopen the BGZF file and resolve CIGAR.
    use std::fs;
    let paf_path = "/tmp/pgr_cli_test_bgzf_lazy_persist.paf.gz";
    let idx_path = "/tmp/pgr_cli_test_bgzf_lazy_persist.paf.idx";
    write_bgzf(paf_path, GZ_PAF);

    // Build + save.
    let (_, stderr) = PgrCmd::new()
        .args(&["paf", "index", paf_path, "-o", idx_path])
        .run();
    assert!(stderr.contains("mode:      lazy"), "lazy mode at build");
    assert!(stderr.contains("saved to"));

    // Load + query (lazy CIGAR must be re-attached via reopen_lazy_source).
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "query", idx_path, "B:0-100"])
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

    let _ = fs::remove_file(paf_path);
    let _ = fs::remove_file(idx_path);
}

#[test]
fn command_paf_bgzf_lazy_vs_plain_text_same_result() {
    // Lazy (BGZF) and in-memory (plain text) must produce identical query output.
    use std::fs;
    let bgzf_path = "/tmp/pgr_cli_test_bgzf_vs_plain.paf.gz";
    let plain_path = "/tmp/pgr_cli_test_bgzf_vs_plain.paf";
    write_bgzf(bgzf_path, GZ_PAF);
    fs::write(plain_path, GZ_PAF).unwrap();

    let (lazy_out, _) = PgrCmd::new()
        .args(&["paf", "query", bgzf_path, "B:0-100"])
        .run();
    let (plain_out, _) = PgrCmd::new()
        .args(&["paf", "query", plain_path, "B:0-100"])
        .run();
    assert_eq!(
        lazy_out, plain_out,
        "BGZF lazy query output differs from plain text"
    );

    let _ = fs::remove_file(bgzf_path);
    let _ = fs::remove_file(plain_path);
}

#[test]
fn command_paf_bgzf_lazy_transitive_bfs() {
    // Transitive BFS over a BGZF index must resolve CIGARs across multiple hops.
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t100\t+\tA\t100\t0\t100\t90\t100\t255\tcg:Z:100M
";
    let path = "/tmp/pgr_cli_test_bgzf_lazy_bfs.paf.gz";
    write_bgzf(path, paf);
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "query", path, "B:0-100", "--transitive"])
        .run();
    assert!(stderr.contains("mode: lazy"));
    assert!(stdout.contains("A\t0\t0\t100\t+\tB"), "1-hop A not found");
    assert!(stdout.contains("C\t0\t0\t100\t+\tA"), "2-hop C not found");
    let _ = std::fs::remove_file(path);
}
