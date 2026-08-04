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

/// Corrupt the contig id of the first occurrence record (0x7f = 127, beyond
/// any single-contig index table). The header stays valid, so a reader that
/// skips per-record validation would only fail (or panic) later.
fn corrupt_first_record_contig(path: &std::path::Path) {
    let mut bytes = std::fs::read(path).unwrap();
    let n_contigs = u32::from_le_bytes(bytes[20..24].try_into().unwrap()) as usize;
    let kmer_bytes = u32::from_le_bytes(bytes[32..36].try_into().unwrap()) as usize;
    let pos_bytes = u32::from_le_bytes(bytes[36..40].try_into().unwrap()) as usize;
    let mut off = 48usize;
    for _ in 0..n_contigs {
        let nb = u32::from_le_bytes(bytes[off..off + 4].try_into().unwrap()) as usize;
        off += 4 + nb + 8;
    }
    let cont_off = off + kmer_bytes + pos_bytes;
    bytes[cont_off] = 0x7f;
    std::fs::write(path, bytes).unwrap();
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

/// A crafted .pgi with an out-of-range record contig id must produce a
/// friendly error, not a panic (Zero Panic).
#[test]
fn command_align_pgi_crafted_index_errors_not_panics() {
    let temp = tempfile::TempDir::new().unwrap();
    let (_, idx) = build_pgi(temp.path(), "ref");
    corrupt_first_record_contig(std::path::Path::new(&idx));

    let out = temp.path().join("out.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&["align", "pgi", &idx, "-o", out.to_str().unwrap()])
        .run();
    assert!(
        stderr.contains("out of range"),
        "crafted index must error, got: {stderr}"
    );
}

/// The sibling index of a `.fa.gz` input must be `ref.fa.pgi` (distinct from
/// a plain `ref.fa`'s `ref.pgi`). Sharing one index would reuse the wrong
/// k-mers when both files exist with the same contig names/lengths but
/// different sequences (silently producing empty or wrong alignments).
#[test]
fn command_align_pgi_gz_sibling_index_distinct() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "ref", &random_seq(400, 71));
    let fa_gz = temp.path().join("ref.fa.gz");
    let mut gz = flate2::write::GzEncoder::new(
        std::fs::File::create(&fa_gz).unwrap(),
        flate2::Compression::default(),
    );
    use std::io::Write;
    write!(gz, ">ref\n{}\n", random_seq(400, 72)).unwrap();
    gz.finish().unwrap();

    // Build the index for the plain file first.
    let out1 = temp.path().join("out1.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &fa,
            "-o",
            out1.to_str().unwrap(),
            "--keep-index",
        ])
        .run();
    assert!(
        stderr.contains("built reference index") && !stderr.contains("out of range"),
        "plain align failed: {stderr}"
    );
    assert!(
        temp.path().join("ref.pgi").is_file(),
        "ref.fa must map to ref.pgi"
    );

    // The gzipped file must get its own sibling index, not reuse ref.pgi.
    let out2 = temp.path().join("out2.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            fa_gz.to_str().unwrap(),
            "-o",
            out2.to_str().unwrap(),
            "--keep-index",
        ])
        .run();
    assert!(
        stderr.contains("ref.fa.pgi"),
        "ref.fa.gz must build its own ref.fa.pgi, got: {stderr}"
    );
}

/// A FASTA edited in place (same contig names/lengths, different sequence)
/// must not silently reuse the stale sibling index; the index is rebuilt
/// when the input's mtime is newer (same convention as the e-kmer cache).
#[test]
fn command_align_pgi_stale_sibling_index_rebuilt() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = temp.path().join("ref.fa");
    std::fs::write(&fa, format!(">ref\n{}\n", random_seq(400, 91))).unwrap();

    let out1 = temp.path().join("out1.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            fa.to_str().unwrap(),
            "-o",
            out1.to_str().unwrap(),
            "--keep-index",
        ])
        .run();
    assert!(
        stderr.contains("built reference index"),
        "build failed: {stderr}"
    );
    assert!(temp.path().join("ref.pgi").is_file());

    // Overwrite with a different sequence, same contig name/length.
    std::thread::sleep(std::time::Duration::from_millis(1100));
    std::fs::write(&fa, format!(">ref\n{}\n", random_seq(400, 92))).unwrap();

    let out2 = temp.path().join("out2.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            fa.to_str().unwrap(),
            "-o",
            out2.to_str().unwrap(),
            "--keep-index",
        ])
        .run();
    assert!(
        stderr.contains("built reference index"),
        "stale index must be rebuilt, got: {stderr}"
    );
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
fn command_align_pgi_output_not_overwrite_sibling_index() {
    // Regression: `-o ref.pgi` must not silently overwrite the sibling index
    // that `ref.fa` maps to. Doing so corrupted the index and broke the next
    // run with a confusing "reading header / failed to fill whole buffer".
    let temp = tempfile::TempDir::new().unwrap();
    let ref_fa = write_fa(temp.path(), "ref", &random_seq(400, 42));

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_fa,
            "-o",
            temp.path().join("ref.pgi").to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("also an input file"),
        "expected -o ref.pgi to be rejected: {stderr}"
    );
    assert!(
        !temp.path().join("ref.pgi").exists(),
        "the sibling index must not be created/corrupted"
    );

    // The `.fa.gz` sibling (`ref.fa.pgi`) is protected too.
    use std::io::Write;
    let fa_gz = temp.path().join("ref.fa.gz");
    let mut gz = flate2::write::GzEncoder::new(
        std::fs::File::create(&fa_gz).unwrap(),
        flate2::Compression::default(),
    );
    write!(gz, ">ref\n{}\n", random_seq(400, 43)).unwrap();
    gz.finish().unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            fa_gz.to_str().unwrap(),
            "-o",
            temp.path().join("ref.fa.pgi").to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("also an input file"),
        "expected -o ref.fa.pgi to be rejected: {stderr}"
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

#[test]
fn command_align_pgi_self_flag_with_query() {
    // --self with the same input passed twice is equivalent to single-input
    // self-alignment: no exact self-identity blocks.
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
            fa.to_str().unwrap(),
            "--self",
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

#[test]
fn command_align_pgi_self_flag_conflicting_query() {
    // --self with a different query input must be rejected.
    let temp = tempfile::TempDir::new().unwrap();
    let ref_fa = write_fa(temp.path(), "ref", &random_seq(400, 42));
    let query_fa = write_fa(temp.path(), "query", &random_seq(400, 7));

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &ref_fa,
            &query_fa,
            "--self",
            "-o",
            temp.path().join("out.psl").to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("--self expects the query"),
        "expected conflicting-query error: {stderr}"
    );
}

#[test]
fn command_align_pgi_default_kmer_conflicts_with_cached_index() {
    // Regression: the sibling-index parameter check only fired when `-k` was
    // given explicitly, so a default run (k=40) silently reused a `k=20`
    // sibling index and reported k=40 semantics with k=20 seeds.
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "genome", &random_seq(400, 42));

    let out = temp.path().join("out.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &fa,
            "-k",
            "20",
            "--keep-index",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "k20 build failed: {stderr}");

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            &fa,
            "--keep-index",
            "-o",
            temp.path().join("default.psl").to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("--kmer 40 conflicts with the cached index"),
        "default k=40 must not reuse the k=20 index: {stderr}"
    );
}

#[test]
fn command_align_pgi_single_ref_seq_on_self_pgi() {
    // Regression: a single-input .pgi self-alignment with only --ref-seq
    // errored ("extension sequences are needed for both sides") because the
    // query side received no sequences. Self mode must reuse the reference
    // extension sequences, matching the direct FASTA output.
    let temp = tempfile::TempDir::new().unwrap();
    let seq = random_seq(400, 42);
    let genome = format!(">genome\n{seq}\n{seq}\n"); // two copies -> repeat blocks
    let fa = temp.path().join("genome.fa");
    fs::write(&fa, &genome).unwrap();
    let pgi = temp.path().join("genome.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            fa.to_str().unwrap(),
            "-o",
            pgi.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "build failed: {stderr}");

    let direct = temp.path().join("direct.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            fa.to_str().unwrap(),
            "-o",
            direct.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "self-align failed: {stderr}");

    let out = temp.path().join("out.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            pgi.to_str().unwrap(),
            "--ref-seq",
            fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "failed: {stderr}");
    assert_eq!(
        fs::read_to_string(&direct).unwrap(),
        fs::read_to_string(&out).unwrap(),
        ".pgi + --ref-seq output differs from direct FASTA self-alignment"
    );

    // The symmetric case (only --query-seq) must behave identically.
    let out = temp.path().join("out_q.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            pgi.to_str().unwrap(),
            "--query-seq",
            fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "failed: {stderr}");
    assert_eq!(
        fs::read_to_string(&direct).unwrap(),
        fs::read_to_string(&out).unwrap(),
        ".pgi + --query-seq output differs from direct FASTA self-alignment"
    );
}

/// A soft-masked (lowercase) region must be skipped by automatic indexing
/// identically whether the genome is a FASTA or a 2bit with stored mask
/// blocks. Regression: `align pgi` read 2bit unmasked, so a masked repeat
/// was indexed and aligned (1 block) while the equivalent FASTA produced 0.
#[test]
fn command_align_pgi_2bit_mask_matches_fasta() {
    let temp = tempfile::TempDir::new().unwrap();
    let pre = random_seq(300, 43);
    let rep = random_seq(200, 44);
    let post = random_seq(300, 45);
    let genome = format!("{pre}{}{post}", rep.to_ascii_lowercase());
    let ref_fa = temp.path().join("ref.fa");
    fs::write(&ref_fa, format!(">ref\n{genome}\n")).unwrap();
    let query_fa = write_fa(temp.path(), "q", &rep);

    // FASTA: the lowercase region is masked out -> no blocks.
    let out_fa = temp.path().join("out_fa.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            ref_fa.to_str().unwrap(),
            &query_fa,
            "-o",
            out_fa.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "failed: {stderr}");
    assert!(
        fs::read_to_string(&out_fa).unwrap().trim().is_empty(),
        "masked FASTA region must not align"
    );

    // Convert to 2bit (keeps the lowercase as mask blocks) and align again.
    let ref_2bit = temp.path().join("ref.2bit");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            ref_fa.to_str().unwrap(),
            "-o",
            ref_2bit.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.is_empty(), "to-2bit failed: {stderr}");
    assert!(ref_2bit.is_file(), "to-2bit must produce the 2bit file");

    let out_2bit = temp.path().join("out_2bit.psl");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "pgi",
            ref_2bit.to_str().unwrap(),
            &query_fa,
            "-o",
            out_2bit.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"), "failed: {stderr}");
    assert!(
        fs::read_to_string(&out_2bit).unwrap().trim().is_empty(),
        "2bit masked region must align exactly like the FASTA (none)"
    );
}

#[test]
fn command_align_pgi_lowercase_copy_has_no_all_zero_blocks() {
    // Regression: the automatic index used to encode lowercase (soft-masked)
    // bases case-insensitively, so a lowercase copy shared seeds with its
    // uppercase twin but the case-sensitive extension DP failed and the chain
    // fell back to an unscored (all-zero) PSL block. The index must now apply
    // FastGA `-M` semantics (skip lowercase), yielding no blocks at all.
    let temp = tempfile::TempDir::new().unwrap();
    let seq = random_seq(300, 42);
    let genome = format!(
        ">genome\n{}{}{}{}{}\n",
        random_seq(300, 43),
        seq,
        random_seq(300, 44),
        seq.to_ascii_lowercase(),
        random_seq(300, 45)
    );
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
    assert!(stderr.contains("wrote"), "align failed: {stderr}");
    let text = fs::read_to_string(&out).unwrap();
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        let scored = f[0].parse::<u32>().unwrap()
            + f[1].parse::<u32>().unwrap()
            + f[2].parse::<u32>().unwrap();
        assert!(
            scored > 0,
            "no all-zero (unscored) block may be emitted: {line}"
        );
    }
}
