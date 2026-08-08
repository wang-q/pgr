#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

/// Return the absolute path to a fixture in `tests/pbit/input`.
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/pbit/input")
        .join(name)
}

/// Read a FASTA file and return the concatenated sequence (all non-header
/// lines joined, uppercased).
fn read_fasta_seq(path: &std::path::Path) -> String {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect::<String>()
        .to_ascii_uppercase()
}

/// Run `pgr pbit create` with the given ref, then `pgr pbit to-fa` and return
/// the extracted sequence for `sample_name`.
fn create_and_extract(
    dir: &std::path::Path,
    ref_fa: &std::path::Path,
    out_pbit: &std::path::Path,
    create_args: &[&str],
    sample_name: &str,
) -> String {
    let mut args = vec!["pbit", "create", "-r", ref_fa.to_str().unwrap()];
    args.extend_from_slice(create_args);
    args.push("-o");
    args.push(out_pbit.to_str().unwrap());
    PgrCmd::new().args(&args).run();

    let out_dir = dir.join("outdir");
    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    let out_fa = out_dir.join(format!("{}.fa", sample_name));
    assert!(
        out_fa.exists(),
        "expected output file: {}",
        out_fa.display()
    );
    read_fasta_seq(&out_fa)
}

// ── Test 1: + strand roundtrip with =/X/I/D CIGAR ──────────────────────

#[test]
fn test_pbit_paf_plus_strand_roundtrip() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_2000.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_2000_plus_strand.fa").to_str().unwrap(),
            "-p",
            fixture("sample_2000_plus_strand.paf").to_str().unwrap(),
        ],
        "sample_2000_plus_strand",
    );
    let expected = read_fasta_seq(&fixture("sample_2000_plus_strand.fa"));
    assert_eq!(got, expected);
}

// ── Test 2: - strand roundtrip (RC semantics) ──────────────────────────

#[test]
fn test_pbit_paf_minus_strand_roundtrip() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_2000.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_2000_minus_strand.fa").to_str().unwrap(),
            "-p",
            fixture("sample_2000_minus_strand.paf").to_str().unwrap(),
        ],
        "sample_2000_minus_strand",
    );
    let expected = read_fasta_seq(&fixture("sample_2000_minus_strand.fa"));
    assert_eq!(got, expected);
}

// ── Test 2b: Identity encoding (segment == reference interval) ─────────

#[test]
fn test_pbit_paf_identity_encoding() {
    // A sample identical to the reference (full-length '=' CIGAR) must be
    // stored as zero-payload Identity segments (v1010) instead of packed
    // CIGAR, and `to-paf` must still rebuild the original chain row.
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let paf = temp.path().join("identical.paf");
    fs::write(
        &paf,
        "chr1\t12289\t0\t12289\t+\tchr1\t12289\t0\t12289\t12289\t12289\t60\tcg:Z:12289=\n",
    )
    .unwrap();

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_12289.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("ref_12289.fa").to_str().unwrap(),
            "-p",
            paf.to_str().unwrap(),
        ],
        "ref_12289",
    );
    let expected = read_fasta_seq(&fixture("ref_12289.fa"));
    assert_eq!(got, expected);

    // All four 4096 bp segments are Identity-encoded, none CIGAR.
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap()])
        .run();
    assert!(
        stdout.contains("Identity 4"),
        "expected 4 Identity segments, got: {}",
        stdout
    );
    assert!(stdout.contains("Cigar 0"), "got: {}", stdout);

    // to-paf rebuilds the single big-chain row from the Identity segments.
    let out_paf = temp.path().join("out.paf");
    PgrCmd::new()
        .args(&[
            "pbit",
            "to-paf",
            out_pbit.to_str().unwrap(),
            "-o",
            out_paf.to_str().unwrap(),
        ])
        .run();
    let content = fs::read_to_string(&out_paf).unwrap();
    let fields: Vec<&str> = content.trim().split('\t').collect();
    assert_eq!(fields.len(), 17, "12 + gi/bi/cg/cs + ms, got: {}", content);
    assert_eq!(fields[4], "+");
    assert_eq!(fields[9], "12289");
    assert_eq!(fields[10], "12289");
    assert!(fields.contains(&"cg:Z:12289="), "got: {}", content);
    assert!(fields.contains(&"cs:Z::12289"), "got: {}", content);
    assert!(fields.contains(&"gi:f:1.000000"), "got: {}", content);
}

// ── Test 2c: mixed Identity + CIGAR segments in one chain ──────────────

#[test]
fn test_pbit_paf_identity_mixed_chain() {
    // A 12289 bp chain with a single SNP at position 5000 (inside segment 1):
    // segments 0/2/3 slice to pure '=' and must be Identity-encoded, segment
    // 1 stays CIGAR. `to-paf` must merge them back into the original chain
    // row (cg and cs rebuilt from the concatenated ops).
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let sample_fa = temp.path().join("sample.fa");
    let paf = temp.path().join("mixed.paf");

    let ref_seq = read_fasta_seq(&fixture("ref_12289.fa"));
    let mut sample_seq = ref_seq.clone();
    let r = ref_seq.as_bytes()[5000] as char;
    let q = if r == 'A' { 'C' } else { 'A' };
    sample_seq.replace_range(5000..5001, &q.to_string());
    fs::write(&sample_fa, format!(">chr1\n{}\n", sample_seq)).unwrap();
    fs::write(
        &paf,
        "chr1\t12289\t0\t12289\t+\tchr1\t12289\t0\t12289\t12288\t12289\t60\tcg:Z:5000=1X7288=\n",
    )
    .unwrap();

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_12289.fa"),
        &out_pbit,
        &[
            "-i",
            sample_fa.to_str().unwrap(),
            "-p",
            paf.to_str().unwrap(),
        ],
        "sample",
    );
    assert_eq!(got, sample_seq);

    // Three pure-match segments are Identity, one carries the X.
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap()])
        .run();
    assert!(stdout.contains("Identity 3"), "got: {}", stdout);
    assert!(stdout.contains("Cigar 1"), "got: {}", stdout);

    // Chain rebuild keeps the single X at the right position.
    let out_paf = temp.path().join("out.paf");
    PgrCmd::new()
        .args(&[
            "pbit",
            "to-paf",
            out_pbit.to_str().unwrap(),
            "-o",
            out_paf.to_str().unwrap(),
        ])
        .run();
    let content = fs::read_to_string(&out_paf).unwrap();
    let fields: Vec<&str> = content.trim().split('\t').collect();
    assert_eq!(fields[9], "12288");
    assert_eq!(fields[10], "12289");
    assert!(fields.contains(&"cg:Z:5000=1X7288="), "got: {}", content);
    let cs_expected = format!("cs:Z::5000*{}{}:7288", r.to_ascii_uppercase(), q);
    assert!(fields.contains(&cs_expected.as_str()), "got: {}", content);
}

// ── Test 2d: to-paf full regression (big chain rebuilt, small verbatim) ─

#[test]
fn test_pbit_paf_to_paf_roundtrip() {
    // v1009/v1010 regression: `to-paf` must rebuild big chains (CIGAR +
    // Identity segments merged, cg/cs/gi/bi/ms recomputed) and pass small
    // chains through verbatim, reproducing the input PAF field-for-field.
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let sample_fa = temp.path().join("sample.fa");
    let in_paf = temp.path().join("in.paf");
    let out_paf = temp.path().join("out.paf");
    let out_dir = temp.path().join("outdir");

    let ref_seq = read_fasta_seq(&fixture("ref_12289.fa"));
    let mut sample_seq = ref_seq.clone();
    let q = if ref_seq.as_bytes()[6000] == b'C' {
        'A'
    } else {
        'C'
    };
    sample_seq.replace_range(6000..6001, &q.to_string());
    fs::write(&sample_fa, format!(">chr1\n{}\n", sample_seq)).unwrap();

    // Row 1: big chain (>= 10 kb, spans all 4 segments) with one SNP.
    // Row 2: small '-' chain, stored verbatim with its own tags.
    let big = "chr1\t12289\t0\t12289\t+\tchr1\t12289\t0\t12289\t12288\t12289\t255\tcg:Z:6000=1X6288=\tms:i:142\n";
    let small = "chr1\t12289\t500\t1800\t-\tchr1\t12289\t0\t1300\t1300\t1300\t60\tcg:Z:1300=\tcs:Z:1300\tms:i:77\n";
    fs::write(&in_paf, format!("{}{}", big, small)).unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_12289.fa").to_str().unwrap(),
            "-i",
            sample_fa.to_str().unwrap(),
            "-p",
            in_paf.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // to-fa stays lossless.
    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();
    assert_eq!(read_fasta_seq(&out_dir.join("sample.fa")), sample_seq);

    // Segment mix: 3 pure-match segments Identity, the SNP segment CIGAR.
    let (stdout, _) = PgrCmd::new()
        .args(&["pbit", "stat", out_pbit.to_str().unwrap()])
        .run();
    assert!(stdout.contains("Identity 3"), "got: {}", stdout);
    assert!(stdout.contains("Cigar 1"), "got: {}", stdout);

    // to-paf: big-chain row rebuilt, small-chain row verbatim.
    PgrCmd::new()
        .args(&[
            "pbit",
            "to-paf",
            out_pbit.to_str().unwrap(),
            "-o",
            out_paf.to_str().unwrap(),
        ])
        .run();
    let content = fs::read_to_string(&out_paf).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 2, "got: {}", content);

    let big_fields: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(
        &big_fields[0..12],
        &[
            "chr1", "12289", "0", "12289", "+", "chr1", "12289", "0", "12289", "12288", "12289",
            "255"
        ],
        "got: {}",
        content
    );
    assert!(big_fields.contains(&"gi:f:0.999919"), "got: {}", content);
    assert!(big_fields.contains(&"bi:f:0.999919"), "got: {}", content);
    assert!(
        big_fields.contains(&"cg:Z:6000=1X6288="),
        "got: {}",
        content
    );
    let cs_expected = format!("cs:Z::6000*{}{}:6288", ref_seq.as_bytes()[6000] as char, q);
    assert!(
        big_fields.contains(&cs_expected.as_str()),
        "got: {}",
        content
    );
    assert!(big_fields.contains(&"ms:i:142"), "got: {}", content);
    assert_eq!(lines[1], small.trim(), "got: {}", content);
}

// ── Test 3: M op split (minimap2 without --eqx) ────────────────────────

#[test]
fn test_pbit_paf_m_op_split() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_2000.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_2000_snp100.fa").to_str().unwrap(),
            "-p",
            fixture("sample_2000_snp100_M.paf").to_str().unwrap(),
        ],
        "sample_2000_snp100",
    );
    let expected = read_fasta_seq(&fixture("sample_2000_snp100.fa"));
    assert_eq!(got, expected);
}

// ── Test 4: mixed mode (CIGAR via create + LZ-diff via append) ─────────

#[test]
fn test_pbit_paf_mixed_mode() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    // Step 1: create with s1 + PAF (CIGAR mode).
    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_snp100.fa").to_str().unwrap(),
            "-p",
            fixture("sample_2000_snp100.paf").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    // Step 2: append s2 without PAF (LZ-diff mode).
    PgrCmd::new()
        .args(&[
            "pbit",
            "append",
            out_pbit.to_str().unwrap(),
            "-i",
            fixture("sample_2000_snp200.fa").to_str().unwrap(),
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    // Verify s1 (CIGAR-encoded).
    let got_s1 = read_fasta_seq(&out_dir.join("sample_2000_snp100.fa"));
    let expected_s1 = read_fasta_seq(&fixture("sample_2000_snp100.fa"));
    assert_eq!(got_s1, expected_s1);

    // Verify s2 (LZ-diff-encoded).
    let got_s2 = read_fasta_seq(&out_dir.join("sample_2000_snp200.fa"));
    let expected_s2 = read_fasta_seq(&fixture("sample_2000_snp200.fa"));
    assert_eq!(got_s2, expected_s2);
}

// ── Test 5: uncovered segment falls back to LZ-diff ────────────────────

#[test]
fn test_pbit_paf_uncovered_fallback() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_5000.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_5000_snp4500.fa").to_str().unwrap(),
            "-p",
            fixture("sample_5000_partial.paf").to_str().unwrap(),
        ],
        "sample_5000_snp4500",
    );
    let expected = read_fasta_seq(&fixture("sample_5000_snp4500.fa"));
    assert_eq!(got, expected);
}

// ── Test 6: empty PAF → all segments fall back to LZ-diff ──────────────

#[test]
fn test_pbit_paf_empty_paf_all_fallback() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_2000.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_2000_snp100.fa").to_str().unwrap(),
            "-p",
            fixture("empty.paf").to_str().unwrap(),
        ],
        "sample_2000_snp100",
    );
    let expected = read_fasta_seq(&fixture("sample_2000_snp100.fa"));
    assert_eq!(got, expected);
}

// ── Test 7: --name TSV with 3 columns (CIGAR + LZ-diff mix) ────────────

#[test]
fn test_pbit_paf_name_tsv_three_columns() {
    let temp = TempDir::new().unwrap();
    let tsv_path = temp.path().join("samples.tsv");
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    // TSV: name<TAB>fasta<TAB>paf (s1 has paf, s2 doesn't).
    let tsv_content = format!(
        "s1\t{}\t{}\ns2\t{}\n",
        fixture("sample_2000_snp100.fa").display(),
        fixture("sample_2000_snp100.paf").display(),
        fixture("sample_2000_snp200.fa").display(),
    );
    fs::write(&tsv_path, tsv_content).unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "--name",
            tsv_path.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    // Verify s1 (CIGAR-encoded).
    let got_s1 = read_fasta_seq(&out_dir.join("s1.fa"));
    let expected_s1 = read_fasta_seq(&fixture("sample_2000_snp100.fa"));
    assert_eq!(got_s1, expected_s1);

    // Verify s2 (LZ-diff-encoded).
    let got_s2 = read_fasta_seq(&out_dir.join("s2.fa"));
    let expected_s2 = read_fasta_seq(&fixture("sample_2000_snp200.fa"));
    assert_eq!(got_s2, expected_s2);
}

// ── Test 8: -i count ≠ --paf count → error ────────────────────────────

#[test]
fn test_pbit_paf_count_mismatch_error() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let (_stdout, stderr) = PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_identical.fa").to_str().unwrap(),
            "-p",
            fixture("sample_2000_identical.paf").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run_fail();

    assert!(
        stderr.contains("count"),
        "expected stderr to mention 'count', got: {}",
        stderr
    );
}

// ── Test 9: --name + --paf mutually exclusive → error ──────────────────

#[test]
fn test_pbit_paf_name_paf_mutex_error() {
    let temp = TempDir::new().unwrap();
    let tsv_path = temp.path().join("samples.tsv");
    let out_pbit = temp.path().join("out.pbit");

    fs::write(
        &tsv_path,
        format!("s1\t{}\n", fixture("sample_2000_identical.fa").display()),
    )
    .unwrap();

    let (_stdout, stderr) = PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "--name",
            tsv_path.to_str().unwrap(),
            "-p",
            fixture("sample_2000_identical.paf").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run_fail();

    assert!(
        stderr.contains("mutually exclusive"),
        "expected stderr to mention 'mutually exclusive', got: {}",
        stderr
    );
}

// ── Test 10: CIGAR target crosses ref segment boundary → fallback ──────

#[test]
fn test_pbit_paf_target_crosses_seg_boundary_fallback() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_8192.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_4096_snp100.fa").to_str().unwrap(),
            "-p",
            fixture("sample_4096_cross_seg.paf").to_str().unwrap(),
        ],
        "sample_4096_snp100",
    );
    let expected = read_fasta_seq(&fixture("sample_4096_snp100.fa"));
    assert_eq!(got, expected);
}

// ── Test 11: rearrangement (sample contig aligned to different ref contig) ─

#[test]
fn test_pbit_paf_rearrangement() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_2000_2contig.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_2000_chr2only.fa").to_str().unwrap(),
            "-p",
            fixture("sample_2000_rearrange.paf").to_str().unwrap(),
        ],
        "sample_2000_chr2only",
    );
    let expected = read_fasta_seq(&fixture("sample_2000_chr2only.fa"));
    assert_eq!(got, expected);
}

// ── Test 12: PAF record without CIGAR tag → skipped (decision 7) ───────

#[test]
fn test_pbit_paf_record_without_cigar_skipped() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_2000.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_2000_snp100.fa").to_str().unwrap(),
            "-p",
            fixture("sample_2000_no_cigar.paf").to_str().unwrap(),
        ],
        "sample_2000_snp100",
    );
    let expected = read_fasta_seq(&fixture("sample_2000_snp100.fa"));
    assert_eq!(got, expected);
}

// ── Test 13: malformed PAF line skipped with warn (decision 8) ─────────

#[test]
fn test_pbit_paf_malformed_line_skipped() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    // Run create and capture stderr to verify the warn was emitted.
    let (_stdout, stderr) = PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_snp100.fa").to_str().unwrap(),
            "-p",
            fixture("sample_2000_malformed.paf").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("skipping invalid PAF line"),
        "expected warn in stderr, got: {}",
        stderr
    );

    // Extract and verify roundtrip.
    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    let got = read_fasta_seq(&out_dir.join("sample_2000_snp100.fa"));
    let expected = read_fasta_seq(&fixture("sample_2000_snp100.fa"));
    assert_eq!(got, expected);
}

// ── Test 14: delta_cache key includes ref_start/ref_end (Bug 1) ────────
//
// Two sample contigs align to different sub-intervals of the SAME ref
// segment with identical CIGAR (`500=`) → packed_data dedup gives them
// the same delta_id, but their ref_start/ref_end differ. Without
// ref_start/ref_end in the cache key, the second decode hits the cache
// and returns the first segment's bytes.

#[test]
fn test_pbit_paf_cache_key_ref_start_ref_end() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");
    let out_dir = temp.path().join("outdir");

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_1000_seed7.fa").to_str().unwrap(),
            "-i",
            fixture("sample_500_A.fa").to_str().unwrap(),
            "-p",
            fixture("sample_500_A.paf").to_str().unwrap(),
            "-i",
            fixture("sample_500_B.fa").to_str().unwrap(),
            "-p",
            fixture("sample_500_B.paf").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    PgrCmd::new()
        .args(&[
            "pbit",
            "to-fa",
            out_pbit.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .run();

    let got_sa = read_fasta_seq(&out_dir.join("sample_500_A.fa"));
    let expected_sa = read_fasta_seq(&fixture("sample_500_A.fa"));
    assert_eq!(got_sa, expected_sa, "sA mismatch");

    let got_sb = read_fasta_seq(&out_dir.join("sample_500_B.fa"));
    let expected_sb = read_fasta_seq(&fixture("sample_500_B.fa"));
    assert_eq!(got_sb, expected_sb, "sB mismatch (cache key bug)");
}

// ── Test 15: - strand multi-segment roundtrip (Bug 2) ──────────────────
//
// ref = 8192 bp → 2 ref segments (segment_size 4096). sample = RC(ref).
// PAF strand='-', CIGAR='8192='. The CIGAR describes RC(query) vs ref;
// forward query [0,4096) maps to RC(query) [4096,8192). Without
// converting forward→RC coords, slice_cigar_by_query assigns the wrong
// CIGAR slice and target range to each segment, corrupting the `=` ops.

#[test]
fn test_pbit_paf_minus_strand_multi_segment() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_8192.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_8192_minus_strand.fa").to_str().unwrap(),
            "-p",
            fixture("sample_8192_minus_strand.paf").to_str().unwrap(),
        ],
        "sample_8192_minus_strand",
    );
    let expected = read_fasta_seq(&fixture("sample_8192_minus_strand.fa"));
    assert_eq!(got, expected, "minus-strand multi-segment mismatch");
}

// ── Test 16: - strand multi-segment with X/I bases (Bug 2 + X/I) ───────
//
// ref = 8192 bp → 2 ref segments (segment_size 4096). RC(sample) = ref with
// a SNP at position 100 and a 2-bp insertion after position 200. PAF
// strand='-', CIGAR='100=1X99=2I7992=' (describes RC(sample) vs ref). The
// X/I bases must be extracted from RC(sample), not forward sample. Segment 1
// (forward [4096,8194)) maps to RC [0,4098) which contains the X and I ops;
// segment 0 (forward [0,4096)) maps to RC [4098,8194) which is pure match.

#[test]
fn test_pbit_paf_minus_strand_with_xi() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_8192.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_8192_minus_xi.fa").to_str().unwrap(),
            "-p",
            fixture("sample_8192_minus_xi.paf").to_str().unwrap(),
        ],
        "sample_8192_minus_xi",
    );
    let expected = read_fasta_seq(&fixture("sample_8192_minus_xi.fa"));
    assert_eq!(got, expected, "minus-strand X/I roundtrip mismatch");
}

// ── Test 17: empty reference contig does not panic (Issue 1) ──────────
//
// ref has chr_empty (0 bp) + chr1 (2000 bp). sample has chr_empty (500 bp)
// + chr1 (ref with SNP at 100). Without the empty-ref guard,
// ref_group_ids[0] panics for chr_empty. With the guard, chr_empty is
// warned+skipped; since v1006 the unmatched contig is stored losslessly via
// content-based LZ-diff / Raw fallback instead, so both contigs must
// roundtrip exactly (no panic, no data loss).

#[test]
fn test_pbit_paf_empty_ref_contig() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_empty_2000.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_empty_2000.fa").to_str().unwrap(),
            "-p",
            fixture("sample_empty_2000.paf").to_str().unwrap(),
        ],
        "sample_empty_2000",
    );
    // Both contigs roundtrip losslessly (chr_empty via Raw/LZ fallback,
    // chr1 via CIGAR encoding).
    let expected = read_fasta_seq(&fixture("sample_empty_2000.fa"));
    assert_eq!(got, expected, "roundtrip mismatch with empty ref contig");
}

#[test]
fn test_pbit_paf_unknown_query_target_name_skipped() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    // The unknown-target record is skipped (falls back to LZ-diff) and the
    // valid record is used for chr1. The command must succeed without panic.
    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            fixture("ref_2000.fa").to_str().unwrap(),
            "-i",
            fixture("sample_2000_snp100.fa").to_str().unwrap(),
            "-p",
            fixture("sample_2000_unknown_target.paf").to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_2000.fa"),
        &out_pbit,
        &[],
        "sample_2000_snp100",
    );
    let expected = read_fasta_seq(&fixture("sample_2000_snp100.fa"));
    assert_eq!(got, expected);
}

#[test]
fn test_pbit_paf_only_comments_fallback_lzdiff() {
    let temp = TempDir::new().unwrap();
    let out_pbit = temp.path().join("out.pbit");

    let got = create_and_extract(
        temp.path(),
        &fixture("ref_2000.fa"),
        &out_pbit,
        &[
            "-i",
            fixture("sample_2000_snp100.fa").to_str().unwrap(),
            "-p",
            fixture("comments_only.paf").to_str().unwrap(),
        ],
        "sample_2000_snp100",
    );
    let expected = read_fasta_seq(&fixture("sample_2000_snp100.fa"));
    assert_eq!(got, expected);
}
