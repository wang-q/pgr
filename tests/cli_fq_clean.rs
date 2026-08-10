#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::{Read, Write};

fn read_gz(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    let mut dec = flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap());
    dec.read_to_end(&mut out).unwrap();
    out
}

fn write_temp(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

fn write_ref(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

#[test]
fn command_fq_trim_adapter_matches_bbtools_trim_golden() {
    // Byte-level comparison against BBTools 39.38
    // `bbduk.sh ... ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=r trimq=15
    // minlen=60 maxns=0 ftm=5 tossbrokenreads=t ordered=t` on the Lambda
    // golden data (see tests/bbtools/Lambda/README.md).
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("trim.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            "tests/bbtools/Lambda/golden/clumpify.fq.gz",
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/trim.fq.gz")
    );
}

#[test]
fn command_fq_trim_adapter_stats_match_bbtools() {
    // --stats writes the bbduk `stats=` 3-column format; values verified
    // byte-for-byte against BBTools 39.38 on the same input (the #File line
    // carries the input path, so it is reconstructed here).
    let out_dir = tempfile::tempdir().unwrap();
    let stats = out_dir.path().join("trim.stats.txt");

    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            "tests/bbtools/Lambda/golden/clumpify.fq.gz",
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--stats",
            stats.to_str().unwrap(),
            "-o",
            out_dir.path().join("out.fq").to_str().unwrap(),
        ])
        .assert()
        .success();

    let expected = concat!(
        "#File\ttests/bbtools/Lambda/golden/clumpify.fq.gz\n",
        "#Total\t40000\n",
        "#Matched\t767\t1.91750%\n",
        "#Name\tReads\tReadsPct\n",
        "Reverse_adapter\t382\t0.95500%\n",
        "TruSeq_Universal_Adapter\t362\t0.90500%\n",
        "Nextera_LMP_Read2_External_Adapter\t5\t0.01250%\n",
        "PCR_Primers\t3\t0.00750%\n",
        "PhiX_read2_adapter\t3\t0.00750%\n",
        "pcr_dimer\t3\t0.00750%\n",
        "Bisulfite_R1\t2\t0.00500%\n",
        "I5_Adapter_Nextera\t2\t0.00500%\n",
        "I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]501\t2\t0.00500%\n",
        "I7_Nextera_Transposase_2\t1\t0.00250%\n",
        "RNA_PCR_Primer_(RP1)_part_#_15013198\t1\t0.00250%\n",
        "TruSeq_Adapter_Index_1_6\t1\t0.00250%\n",
    );
    assert_eq!(std::fs::read_to_string(&stats).unwrap(), expected);
}

#[test]
fn command_fq_trim_adapter_removes_adapter_and_keeps_clean_read() {
    // A read whose 3' end is a known adapter is trimmed; a clean read is
    // untouched except the ftm multiple-of-5 normalization.
    let input = format!(
        "@r1\n{}AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n+\n{}IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        "A".repeat(70),
        "I".repeat(70)
    );
    let file = write_temp(&input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--no-trim-by-overlap",
            "--no-trim-pair-evenly",
            "--max-ns=-1",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(&out).unwrap();
    let seq = out.lines().nth(1).unwrap();
    assert!(
        !seq.contains("AATGATACGGCGACCACC"),
        "adapter must be trimmed"
    );
    // The hdist=1 table can cut one extra base at the adapter boundary, but
    // the read must keep well above the 60 bp minimum.
    assert!(seq.len() >= 60, "prefix must survive: {seq}");
    assert_eq!(seq.len(), 69, "bbduk-compatible cut position");
}

#[test]
fn command_fq_trim_adapter_parallel_output_matches_single_thread() {
    // The worker pool preserves input order; any thread count must give the
    // same byte output as the golden.
    let out_dir = tempfile::tempdir().unwrap();
    let t1 = out_dir.path().join("t1.fq");
    let t8 = out_dir.path().join("t8.fq");

    for (out, threads) in [(&t1, "1"), (&t8, "8")] {
        PgrCmd::new()
            .args(&[
                "fq",
                "clean",
                "tests/bbtools/Lambda/golden/clumpify.fq.gz",
                "--ref",
                "tests/bbtools/Lambda/illumina_adapters.fa",
                "--parallel",
                threads,
                "-o",
                out.to_str().unwrap(),
            ])
            .assert()
            .success();
    }

    let golden = read_gz("tests/bbtools/Lambda/golden/trim.fq.gz");
    assert_eq!(std::fs::read(&t1).unwrap(), golden);
    assert_eq!(std::fs::read(&t8).unwrap(), golden);
}

#[test]
fn command_fq_trim_adapter_changequality_and_qtrim_empty_edge() {
    // bbduk's default `changequality=t` raises ACGT bases to quality 2 and
    // forces N bases to quality 0; qtrim then clamps a fully-trimmed read to
    // len-1 bases (1 bp for long reads, 0 bp for a 1 bp read). Verified
    // byte-for-byte against BBTools 39.38 with ktrim=f tbo=f tpe=f qtrim=r
    // trimq=15 minlen=0 maxns=-1 ftm=0.
    let input = concat!(
        "@r1\n",
        "AAAAAAAAAAAAAAAAAAAA\n",
        "+\n",
        "!!!!!!!!!!!!!!!!!!!!\n",
        "@r2\n",
        "A\n",
        "+\n",
        "!\n",
    );
    let file = write_temp(input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--no-trim-by-overlap",
            "--no-trim-pair-evenly",
            "--minlen",
            "0",
            "--max-ns=-1",
            "--force-trim-mod",
            "0",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        concat!(
            "@r1\nA\n+\n#\n", // ACGT quality 0 raised to 2, trimmed to 1 bp
            "@r2\n\n+\n\n",   // 1 bp read fully trimmed to empty
        )
    );
}

#[test]
fn command_fq_trim_adapter_no_ref_quality_trim_only() {
    // bbduk `qtrim=r trimq=... minlen=...` without a reference: no k-mer
    // operations, only quality trim + length filter. Verified byte-for-byte
    // against BBTools 39.38 on the same interleaved input.
    let input = concat!(
        "@p1/1\n",
        "ACGTACGTACGTACGTACGT\n",
        "+\n",
        "IIIIIIIIIIIIIIIIIIII\n",
        "@p1/2\n",
        "ACGTACGTACGTACGTACGT\n",
        "+\n",
        "IIIIIIIIIIIIIIIIIIII\n",
        "@p2/1\n",
        "ACGTACGTACGTACGTACGT\n",
        "+\n",
        "IIIIIIIIIIIIIIIIIIII\n",
        "@p2/2\n",
        "ACGTACGTACGT\n",
        "+\n",
        "!!!!!!!!!!!!\n",
    );
    let file = write_temp(input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--no-trim-by-overlap",
            "--no-trim-pair-evenly",
            "--max-ns=-1",
            "--force-trim-mod",
            "0",
            "--trim-quality",
            "15",
            "--minlen",
            "10",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    // p2/2 quality-trims to 1 bp (< minlen), so the whole pair is dropped.
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        concat!(
            "@p1/1\n",
            "ACGTACGTACGTACGTACGT\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIII\n",
            "@p1/2\n",
            "ACGTACGTACGTACGTACGT\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIII\n",
        )
    );
}

#[test]
fn command_fq_trim_adapter_qtrim_rl_trims_both_ends() {
    // Low-quality flanks, high-quality core: qtrim=rl trims both ends
    // (bbduk testOptimal keeps the highest-quality run).
    let input = "\
@r1
ACGTACGTACGTACGT
+
!!!!IIIIIIII!!!!
";
    let file = write_temp(input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--qtrim",
            "rl",
            "--trim-quality",
            "15",
            "--no-trim-by-overlap",
            "--no-trim-pair-evenly",
            "--max-ns=-1",
            "--force-trim-mod",
            "0",
            "--minlen",
            "0",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@r1\nACGTACGT\n+\nIIIIIIII\n"
    );
}

#[test]
fn command_fq_trim_adapter_qtrim_window_trims_low_quality_tail() {
    // qtrim=w uses a sliding window: a low-quality tail triggers trimming
    // once a full window falls below the threshold.
    let input = "\
@r1
ACGTACGTACGTACGT
+
IIIIIIII!!!!!!!!
";
    let file = write_temp(input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--qtrim",
            "w",
            "--qtrim-window",
            "4",
            "--trim-quality",
            "15",
            "--no-trim-by-overlap",
            "--no-trim-pair-evenly",
            "--max-ns=-1",
            "--force-trim-mod",
            "0",
            "--minlen",
            "0",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@r1\nACGTACG\n+\nIIIIIII\n"
    );
}

#[test]
fn command_fq_trim_adapter_polymer_trim_and_gc_filter() {
    // poly-A tail trimmed (trimpolya=4), then a read outside the GC band
    // discarded (mingc/maxgc).
    let input = "\
@r1
ACGTACGTACGTACGTAAAA
+
IIIIIIIIIIIIIIIIIIII
@r2
AAAAAAAAAAAAAAAAAAAA
+
IIIIIIIIIIIIIIIIIIII
@r3
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
";
    let file = write_temp(input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--qtrim",
            "f",
            "--no-trim-by-overlap",
            "--no-trim-pair-evenly",
            "--max-ns=-1",
            "--force-trim-mod",
            "0",
            "--minlen",
            "0",
            "--trim-poly-a",
            "4",
            "--min-gc",
            "0.4",
            "--max-gc",
            "0.6",
            "--no-toss-broken-reads",
            "--no-pair-gc",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("@r1\nACGTACGTACGTACGT\n"), "{text}");
    assert!(!text.contains("@r2"), "{text}");
    assert!(text.contains("@r3"), "{text}");
}

#[test]
fn command_fq_trim_adapter_parallel_out_of_range_is_clap_error() {
    // Regression: an out-of-range --parallel must be rejected with a friendly
    // error before a thread pool is created, not spawn 1025+ worker threads.
    let file = write_temp("@r1\nACGTACGT\n+\nIIIIIIII\n");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--parallel",
            "1000000",
            "-o",
            "stdout",
        ])
        .run_fail();
    assert!(
        stderr.contains("--parallel") || stderr.contains("1..=1024"),
        "stderr: {stderr}"
    );
}

#[test]
fn command_fq_trim_adapter_maq_mbq_and_mcb_filters() {
    // maq=20 keeps r1 (avg Q30), mbq=5 discards r2 (has a Q0 base), mcb=10
    // discards r3 (N run breaks consecutive ACGT).
    let input = "\
@r1
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@r2
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@r3
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII!!!!
@r4
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@r5
ACGTACGTNNNNNNNNNNNN
+
IIIIIIIIIIIIIIIIIIII
@r6
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
";
    let file = write_temp(input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--qtrim",
            "f",
            "--no-trim-by-overlap",
            "--no-trim-pair-evenly",
            "--max-ns=-1",
            "--force-trim-mod",
            "0",
            "--minlen",
            "0",
            "--min-avg-quality",
            "20",
            "--min-base-quality",
            "5",
            "--min-consecutive-bases",
            "10",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.starts_with("@r1\n"), "r1 must survive: {text}");
    assert!(text.contains("@r2\n"), "pair1 must survive: {text}");
    assert!(!text.contains("@r3"), "pair2 must fail mbq: {text}");
    assert!(!text.contains("@r5"), "pair3 must fail mcb: {text}");
}

#[test]
fn command_fq_clean_mask_fully_covered_keeps_hit_regions() {
    // bbduk maskfullycovered: hit k-mer windows stay set, only windows with
    // no match are cleared. A read with a long contiguous adapter run must
    // keep the adapter region masked.
    let input = "\
@r1
GGGGGGGGGGGGGGGGGGGGAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
";
    let file = write_temp(input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "-k",
            "23",
            "--min-k",
            "0",
            "--hamming-distance",
            "0",
            "--mask-kmers",
            "N",
            "--mask-fully-covered",
            "--no-trim-by-overlap",
            "--no-trim-pair-evenly",
            "--no-qtrim",
            "--max-ns=-1",
            "--force-trim-mod",
            "0",
            "--minlen",
            "0",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    let seq = text.lines().nth(1).unwrap();
    assert!(
        seq.ends_with("NNNNNNNNNNNNNNNNNNNNNNNNNNNN"),
        "adapter run must stay masked: {seq}"
    );
}

#[test]
fn command_fq_clean_qtrim_low_trimq_trims_trailing_ns() {
    // bbduk phredToProbError(0)=0.75 makes trailing-N windows (nprob=0.825)
    // negative, trimming them; a missing 0/1 special case kept them.
    let input = "\
@r1
ACGTACGTNNNNNNNN
+
IIIIIIIIIIIIIIII
";
    let file = write_temp(input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--qtrim",
            "r",
            "--trim-quality",
            "0",
            "--no-trim-by-overlap",
            "--no-trim-pair-evenly",
            "--max-ns=-1",
            "--force-trim-mod",
            "0",
            "--minlen",
            "0",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@r1\nACGTACGT\n+\nIIIIIIII\n"
    );
}

#[test]
fn command_fq_clean_force_trim_mod_short_read_no_panic() {
    // ftm=5 on a 4 bp read: bbduk's b0=len-1-len%ftm is negative and the
    // trim clamps to 1 bp; a usize underflow must not panic.
    let input = "\
@r1
ACGT
+
IIII
";
    let file = write_temp(input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--qtrim",
            "f",
            "--no-trim-by-overlap",
            "--no-trim-pair-evenly",
            "--max-ns=-1",
            "--force-trim-mod",
            "5",
            "--minlen",
            "0",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    assert_eq!(std::fs::read_to_string(&out).unwrap(), "@r1\nA\n+\nI\n");
}

#[test]
fn command_fq_clean_kmask_masks_instead_of_trims() {
    // Regression: --mask-kmers must switch the main k-mer operation from
    // right-trimming to masking (bbduk kmask). Adapter k-mers become N and the
    // full read is kept, whereas the default ktrim=right would cut them off.
    let ref_file = write_ref(">adapter\nGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n");
    let input = format!(
        "@r1\nACGTACGTACGTACGTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n+\n{}\n",
        "I".repeat(49)
    );
    let file = write_temp(&input);
    let out_dir = tempfile::tempdir().unwrap();

    let masked = out_dir.path().join("masked.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--ref",
            ref_file.path().to_str().unwrap(),
            "--mask-kmers",
            "N",
            "--qtrim",
            "f",
            "--max-ns=-1",
            "--force-trim-mod",
            "0",
            "--minlen",
            "0",
            "-o",
            masked.to_str().unwrap(),
        ])
        .assert()
        .success();

    let text = std::fs::read_to_string(&masked).unwrap();
    let seq = text.lines().nth(1).unwrap();
    assert_eq!(
        seq.len(),
        49,
        "read must not be trimmed in kmask mode: {seq}"
    );
    assert_eq!(
        seq,
        "ACGTACGTACGTACGT".to_string() + &"N".repeat(33),
        "adapter k-mers must be masked with N"
    );

    // Same input without --mask-kmers keeps the default ktrim=right: the
    // adapter is cut off, not masked.
    let trimmed = out_dir.path().join("trimmed.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--ref",
            ref_file.path().to_str().unwrap(),
            "--qtrim",
            "f",
            "--max-ns=-1",
            "--force-trim-mod",
            "0",
            "--minlen",
            "0",
            "-o",
            trimmed.to_str().unwrap(),
        ])
        .assert()
        .success();

    let seq = std::fs::read_to_string(&trimmed)
        .unwrap()
        .lines()
        .nth(1)
        .unwrap()
        .to_string();
    assert!(seq.len() < 49, "ktrim=right must trim the adapter: {seq}");
    assert!(!seq.contains('N'), "no masking without --mask-kmers: {seq}");
}

#[test]
fn command_fq_clean_kmask_mask_only_options_require_mask_kmers() {
    // --mask-fully-covered / --trim-pad only apply to masking; using them
    // without --mask-kmers is a misconfiguration and must be a friendly error.
    let ref_file = write_ref(">adapter\nGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n");
    let file = write_temp("@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n");

    for extra in [&["--mask-fully-covered"][..], &["--trim-pad", "3"][..]] {
        let (_, stderr) = PgrCmd::new()
            .args(&[
                "fq",
                "clean",
                file.path().to_str().unwrap(),
                "--ref",
                ref_file.path().to_str().unwrap(),
                "-o",
                "stdout",
            ])
            .args(extra)
            .run_fail();
        assert!(stderr.contains("--mask-kmers"), "stderr: {stderr}");
    }
}

#[test]
fn command_fq_clean_kmask_requires_ref() {
    // --mask-kmers only does something with a reference; without one it is a
    // silent no-op, so it must be rejected as documented ("requires --ref").
    let file = write_temp("@r1\nACGTACGT\n+\nIIIIIIII\n");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fq",
            "clean",
            file.path().to_str().unwrap(),
            "--mask-kmers",
            "N",
            "-o",
            "stdout",
        ])
        .run_fail();
    assert!(
        stderr.contains("--mask-kmers") && stderr.contains("--ref"),
        "stderr: {stderr}"
    );
}
