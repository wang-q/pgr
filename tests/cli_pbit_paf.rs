#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use std::io::Write;
use tempfile::TempDir;

/// Write a FASTA file with the given (name, seq) records.
fn write_fasta(path: &std::path::Path, records: &[(&str, &str)]) {
    let mut f = fs::File::create(path).unwrap();
    for (name, seq) in records {
        writeln!(f, ">{}", name).unwrap();
        writeln!(f, "{}", seq).unwrap();
    }
}

/// Generate a deterministic random DNA sequence of the given length.
fn random_dna(len: usize, seed: u64) -> Vec<u8> {
    use rand::rngs::StdRng;
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = StdRng::seed_from_u64(seed);
    (0..len)
        .map(|_| match rng.random_range(0u8..4) {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect()
}

/// Return a different DNA base (for introducing SNPs).
fn snp_of(b: u8) -> u8 {
    match b {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

/// Build a single PAF line with 12 mandatory fields + `cg:Z:` CIGAR tag.
fn make_paf_line(
    qname: &str,
    qlen: u32,
    qs: u32,
    qe: u32,
    strand: &str,
    tname: &str,
    tlen: u32,
    ts: u32,
    te: u32,
    cigar: &str,
) -> String {
    // matches = sum of = ops (approximation; pbit does not use this field).
    // block_length = target span (te - ts). mapq = 60.
    let matches = (qe - qs).min(te - ts);
    let block_len = te - ts;
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t60\tcg:Z:{}",
        qname, qlen, qs, qe, strand, tname, tlen, ts, te, matches, block_len, cigar
    )
}

/// Write PAF lines to a file (one per line).
fn write_paf(path: &std::path::Path, lines: &[String]) {
    let mut f = fs::File::create(path).unwrap();
    for line in lines {
        writeln!(f, "{}", line).unwrap();
    }
}

/// Read a FASTA file produced by `pgr pbit to-fa` and return the concatenated
/// sequence (all non-header lines joined, uppercased).
fn read_extracted_fa(path: &std::path::Path) -> String {
    let content = fs::read_to_string(path).unwrap();
    content
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect::<String>()
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
    read_extracted_fa(&out_fa)
}

// ── Test 1: + strand roundtrip with =/X/I/D CIGAR ──────────────────────

#[test]
fn test_pbit_paf_plus_strand_roundtrip() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let paf_path = temp.path().join("sample.paf");
    let out_pbit = temp.path().join("out.pbit");

    // ref = 2000 bp. Sample (2001 bp) = ref with SNP at 100, 2-bp insertion
    // after position 200, and 1-bp deletion at ref position 299.
    // CIGAR: 100= 1X 99= 2I 99= 1D 1700=
    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);

    let mut sample: Vec<u8> = Vec::with_capacity(2001);
    sample.extend_from_slice(&ref_seq[0..100]); // 100=
    sample.push(snp_of(ref_seq[100])); // 1X
    sample.extend_from_slice(&ref_seq[101..200]); // 99=
    sample.extend_from_slice(b"GT"); // 2I
    sample.extend_from_slice(&ref_seq[200..299]); // 99=
                                                  // ref[299] deleted (1D)
    sample.extend_from_slice(&ref_seq[300..2000]); // 1700=
    write_fasta(
        &sample_fa,
        &[("chr1", std::str::from_utf8(&sample).unwrap())],
    );

    let cigar = "100=1X99=2I99=1D1700=";
    write_paf(
        &paf_path,
        &[make_paf_line(
            "chr1", 2001, 0, 2001, "+", "chr1", 2000, 0, 2000, cigar,
        )],
    );

    let got = create_and_extract(
        temp.path(),
        &ref_fa,
        &out_pbit,
        &[
            "-i",
            sample_fa.to_str().unwrap(),
            "-p",
            paf_path.to_str().unwrap(),
        ],
        "sample",
    );
    let expected: String = sample.iter().map(|&b| b as char).collect();
    assert_eq!(got, expected);
}

// ── Test 2: - strand roundtrip (RC semantics) ──────────────────────────

#[test]
fn test_pbit_paf_minus_strand_roundtrip() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let paf_path = temp.path().join("sample.paf");
    let out_pbit = temp.path().join("out.pbit");

    // ref = 2000 bp. sample = RC(ref). PAF strand='-', CIGAR='2000='
    // (describes RC(sample) vs ref = ref vs ref).
    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);

    let sample: Vec<u8> = pgr::libs::nt::rev_comp(&ref_seq).collect();
    write_fasta(
        &sample_fa,
        &[("chr1", std::str::from_utf8(&sample).unwrap())],
    );

    write_paf(
        &paf_path,
        &[make_paf_line(
            "chr1", 2000, 0, 2000, "-", "chr1", 2000, 0, 2000, "2000=",
        )],
    );

    let got = create_and_extract(
        temp.path(),
        &ref_fa,
        &out_pbit,
        &[
            "-i",
            sample_fa.to_str().unwrap(),
            "-p",
            paf_path.to_str().unwrap(),
        ],
        "sample",
    );
    let expected: String = sample.iter().map(|&b| b as char).collect();
    assert_eq!(got, expected);
}

// ── Test 3: M op split (minimap2 without --eqx) ────────────────────────

#[test]
fn test_pbit_paf_m_op_split() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let paf_path = temp.path().join("sample.paf");
    let out_pbit = temp.path().join("out.pbit");

    // ref = 2000 bp. sample = ref with SNP at position 100.
    // PAF CIGAR = '2000M' (no =/X distinction) → compressor splits M into =/X.
    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);

    let mut sample = ref_seq.clone();
    sample[100] = snp_of(ref_seq[100]);
    write_fasta(
        &sample_fa,
        &[("chr1", std::str::from_utf8(&sample).unwrap())],
    );

    write_paf(
        &paf_path,
        &[make_paf_line(
            "chr1", 2000, 0, 2000, "+", "chr1", 2000, 0, 2000, "2000M",
        )],
    );

    let got = create_and_extract(
        temp.path(),
        &ref_fa,
        &out_pbit,
        &[
            "-i",
            sample_fa.to_str().unwrap(),
            "-p",
            paf_path.to_str().unwrap(),
        ],
        "sample",
    );
    let expected: String = sample.iter().map(|&b| b as char).collect();
    assert_eq!(got, expected);
}

// ── Test 4: mixed mode (CIGAR via create + LZ-diff via append) ─────────

#[test]
fn test_pbit_paf_mixed_mode() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let s1_fa = temp.path().join("s1.fa");
    let s1_paf = temp.path().join("s1.paf");
    let s2_fa = temp.path().join("s2.fa");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);

    // s1: SNP at 100, with PAF (CIGAR mode via `pbit create`).
    let mut s1 = ref_seq.clone();
    s1[100] = snp_of(ref_seq[100]);
    write_fasta(&s1_fa, &[("chr1", std::str::from_utf8(&s1).unwrap())]);
    write_paf(
        &s1_paf,
        &[make_paf_line(
            "chr1",
            2000,
            0,
            2000,
            "+",
            "chr1",
            2000,
            0,
            2000,
            "100=1X1899=",
        )],
    );

    // s2: SNP at 200, no PAF (LZ-diff mode via `pbit append`).
    let mut s2 = ref_seq.clone();
    s2[200] = snp_of(ref_seq[200]);
    write_fasta(&s2_fa, &[("chr1", std::str::from_utf8(&s2).unwrap())]);

    // Step 1: create with s1 + PAF (CIGAR mode).
    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
            "-i",
            s1_fa.to_str().unwrap(),
            "-p",
            s1_paf.to_str().unwrap(),
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
            s2_fa.to_str().unwrap(),
        ])
        .run();

    let out_dir = temp.path().join("outdir");
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
    let got_s1 = read_extracted_fa(&out_dir.join("s1.fa"));
    let expected_s1: String = s1.iter().map(|&b| b as char).collect();
    assert_eq!(got_s1, expected_s1);

    // Verify s2 (LZ-diff-encoded).
    let got_s2 = read_extracted_fa(&out_dir.join("s2.fa"));
    let expected_s2: String = s2.iter().map(|&b| b as char).collect();
    assert_eq!(got_s2, expected_s2);
}

// ── Test 5: uncovered segment falls back to LZ-diff ────────────────────

#[test]
fn test_pbit_paf_uncovered_fallback() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let paf_path = temp.path().join("sample.paf");
    let out_pbit = temp.path().join("out.pbit");

    // ref = 5000 bp → 2 segments (4096 + 904). PAF covers only [0, 4096).
    // Second segment falls back to LZ-diff. SNP at position 4500 (in 2nd seg).
    let ref_seq = random_dna(5000, 42);
    write_fasta(&ref_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);

    let mut sample = ref_seq.clone();
    sample[4500] = snp_of(ref_seq[4500]);
    write_fasta(
        &sample_fa,
        &[("chr1", std::str::from_utf8(&sample).unwrap())],
    );

    // PAF covers [0, 4096) with CIGAR 4096= (first segment matches exactly).
    write_paf(
        &paf_path,
        &[make_paf_line(
            "chr1", 5000, 0, 4096, "+", "chr1", 5000, 0, 4096, "4096=",
        )],
    );

    let got = create_and_extract(
        temp.path(),
        &ref_fa,
        &out_pbit,
        &[
            "-i",
            sample_fa.to_str().unwrap(),
            "-p",
            paf_path.to_str().unwrap(),
        ],
        "sample",
    );
    let expected: String = sample.iter().map(|&b| b as char).collect();
    assert_eq!(got, expected);
}

// ── Test 6: empty PAF → all segments fall back to LZ-diff ──────────────

#[test]
fn test_pbit_paf_empty_paf_all_fallback() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let sample_fa = temp.path().join("sample.fa");
    let paf_path = temp.path().join("empty.paf");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);

    let mut sample = ref_seq.clone();
    sample[100] = snp_of(ref_seq[100]);
    write_fasta(
        &sample_fa,
        &[("chr1", std::str::from_utf8(&sample).unwrap())],
    );

    // Empty PAF file.
    write_paf(&paf_path, &[]);

    let got = create_and_extract(
        temp.path(),
        &ref_fa,
        &out_pbit,
        &[
            "-i",
            sample_fa.to_str().unwrap(),
            "-p",
            paf_path.to_str().unwrap(),
        ],
        "sample",
    );
    let expected: String = sample.iter().map(|&b| b as char).collect();
    assert_eq!(got, expected);
}

// ── Test 7: --name TSV with 3 columns (CIGAR + LZ-diff mix) ────────────

#[test]
fn test_pbit_paf_name_tsv_three_columns() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let s1_fa = temp.path().join("s1.fa");
    let s1_paf = temp.path().join("s1.paf");
    let s2_fa = temp.path().join("s2.fa");
    let tsv_path = temp.path().join("samples.tsv");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);

    // s1: SNP at 100, with PAF (CIGAR mode).
    let mut s1 = ref_seq.clone();
    s1[100] = snp_of(ref_seq[100]);
    write_fasta(&s1_fa, &[("chr1", std::str::from_utf8(&s1).unwrap())]);
    write_paf(
        &s1_paf,
        &[make_paf_line(
            "chr1",
            2000,
            0,
            2000,
            "+",
            "chr1",
            2000,
            0,
            2000,
            "100=1X1899=",
        )],
    );

    // s2: SNP at 200, no PAF (LZ-diff mode).
    let mut s2 = ref_seq.clone();
    s2[200] = snp_of(ref_seq[200]);
    write_fasta(&s2_fa, &[("chr1", std::str::from_utf8(&s2).unwrap())]);

    // TSV: name<TAB>fasta<TAB>paf (s1 has paf, s2 doesn't).
    let tsv_content = format!(
        "s1\t{}\t{}\ns2\t{}\n",
        s1_fa.display(),
        s1_paf.display(),
        s2_fa.display()
    );
    fs::write(&tsv_path, tsv_content).unwrap();

    PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
            "--name",
            tsv_path.to_str().unwrap(),
            "-o",
            out_pbit.to_str().unwrap(),
        ])
        .run();

    let out_dir = temp.path().join("outdir");
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
    let got_s1 = read_extracted_fa(&out_dir.join("s1.fa"));
    let expected_s1: String = s1.iter().map(|&b| b as char).collect();
    assert_eq!(got_s1, expected_s1);

    // Verify s2 (LZ-diff-encoded).
    let got_s2 = read_extracted_fa(&out_dir.join("s2.fa"));
    let expected_s2: String = s2.iter().map(|&b| b as char).collect();
    assert_eq!(got_s2, expected_s2);
}

// ── Test 8: -i count ≠ --paf count → error ────────────────────────────

#[test]
fn test_pbit_paf_count_mismatch_error() {
    let temp = TempDir::new().unwrap();
    let ref_fa = temp.path().join("ref.fa");
    let s1_fa = temp.path().join("s1.fa");
    let s2_fa = temp.path().join("s2.fa");
    let paf_path = temp.path().join("s1.paf");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);
    write_fasta(&s1_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);
    write_fasta(&s2_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);
    write_paf(
        &paf_path,
        &[make_paf_line(
            "chr1", 2000, 0, 2000, "+", "chr1", 2000, 0, 2000, "2000=",
        )],
    );

    let (_stdout, stderr) = PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
            "-i",
            s1_fa.to_str().unwrap(),
            "-i",
            s2_fa.to_str().unwrap(),
            "-p",
            paf_path.to_str().unwrap(),
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
    let ref_fa = temp.path().join("ref.fa");
    let s1_fa = temp.path().join("s1.fa");
    let paf_path = temp.path().join("s1.paf");
    let tsv_path = temp.path().join("samples.tsv");
    let out_pbit = temp.path().join("out.pbit");

    let ref_seq = random_dna(2000, 42);
    write_fasta(&ref_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);
    write_fasta(&s1_fa, &[("chr1", std::str::from_utf8(&ref_seq).unwrap())]);
    write_paf(
        &paf_path,
        &[make_paf_line(
            "chr1", 2000, 0, 2000, "+", "chr1", 2000, 0, 2000, "2000=",
        )],
    );
    fs::write(&tsv_path, format!("s1\t{}\n", s1_fa.display())).unwrap();

    let (_stdout, stderr) = PgrCmd::new()
        .args(&[
            "pbit",
            "create",
            "-r",
            ref_fa.to_str().unwrap(),
            "--name",
            tsv_path.to_str().unwrap(),
            "-p",
            paf_path.to_str().unwrap(),
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
