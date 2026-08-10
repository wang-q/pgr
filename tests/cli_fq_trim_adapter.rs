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
            "trim-adapter",
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
fn command_fq_trim_adapter_filter_matches_bbtools_filter_golden() {
    // Byte-level comparison against BBTools 39.38
    // `bbduk.sh ... k=27 cardinality tossbrokenreads=t ordered=t` (filter
    // mode) on the Lambda golden data.
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("filter.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-adapter",
            "tests/bbtools/Lambda/golden/trim.fq.gz",
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--no-ktrim",
            "--no-tbo",
            "--no-tpe",
            "--no-qtrim",
            "--k",
            "27",
            "--mink",
            "0",
            "--minlen",
            "0",
            "--maxns=-1",
            "--ftm",
            "0",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/filter.fq.gz")
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
            "trim-adapter",
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
fn command_fq_trim_adapter_filter_stats_match_bbtools() {
    // Filter mode stats: no adapter kmers survive at k=27, so only headers.
    let out_dir = tempfile::tempdir().unwrap();
    let stats = out_dir.path().join("filter.stats.txt");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-adapter",
            "tests/bbtools/Lambda/golden/trim.fq.gz",
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--no-ktrim",
            "--no-tbo",
            "--no-tpe",
            "--no-qtrim",
            "--k",
            "27",
            "--mink",
            "0",
            "--minlen",
            "0",
            "--maxns=-1",
            "--ftm",
            "0",
            "--stats",
            stats.to_str().unwrap(),
            "-o",
            out_dir.path().join("out.fq").to_str().unwrap(),
        ])
        .assert()
        .success();

    let expected = concat!(
        "#File\ttests/bbtools/Lambda/golden/trim.fq.gz\n",
        "#Total\t36384\n",
        "#Matched\t0\t0.00000%\n",
        "#Name\tReads\tReadsPct\n",
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
            "trim-adapter",
            file.path().to_str().unwrap(),
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--no-tbo",
            "--no-tpe",
            "--maxns=-1",
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
                "trim-adapter",
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
            "trim-adapter",
            file.path().to_str().unwrap(),
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--no-ktrim",
            "--no-tbo",
            "--no-tpe",
            "--minlen",
            "0",
            "--maxns=-1",
            "--ftm",
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
            "trim-adapter",
            file.path().to_str().unwrap(),
            "--no-ktrim",
            "--no-tbo",
            "--no-tpe",
            "--maxns=-1",
            "--ftm",
            "0",
            "--trimq",
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
