#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::Write;

fn write_temp(suffix: &str, content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

#[test]
fn command_fq_trim_q_single_end() {
    // 40 bp read with a 10 bp low-quality tail is trimmed to 30 bp; the short
    // second record is discarded; the header comment is preserved.
    let good = "?".repeat(30);
    let input = format!(
        "@r1 some comment\n{}\n+\n{}{}\n@r2\nACGT\n+\n!!!!\n",
        "A".repeat(40),
        good,
        "!".repeat(10)
    );
    let file = write_temp(".fq", &input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let expected = format!("@r1 some comment\n{}\n+\n{good}\n", "A".repeat(30));
    assert_eq!(std::fs::read_to_string(&out).unwrap(), expected);
}

#[test]
fn command_fq_trim_q_length_threshold_drops_trimmed_reads() {
    // Trimmed to 30 bp but -l 35 discards it.
    let good = "?".repeat(30);
    let input = format!("@r1\n{}\n+\n{}{}\n", "A".repeat(40), good, "!".repeat(10));
    let file = write_temp(".fq", &input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "-l",
            "35",
        ])
        .assert()
        .success();

    assert_eq!(std::fs::read_to_string(&out).unwrap(), "");
}

#[test]
fn command_fq_trim_q_no_fiveprime() {
    // The first window average reaches the threshold despite one bad base, so
    // default sliding cuts the 5' bad base; --no-fiveprime keeps the read.
    let input = format!("@r1\n{}\n+\n{}{}\n", "A".repeat(40), "!", "?".repeat(39));
    let file = write_temp(".fq", &input);
    let out_dir = tempfile::tempdir().unwrap();
    let out_default = out_dir.path().join("out_default.fq");
    let out_no5 = out_dir.path().join("out_no5.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out_default.to_str().unwrap(),
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&out_default).unwrap(),
        format!("@r1\n{}\n+\n{}\n", "A".repeat(39), "?".repeat(39))
    );

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out_no5.to_str().unwrap(),
            "--no-fiveprime",
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&out_no5).unwrap(),
        format!("@r1\n{}\n+\n!{}\n", "A".repeat(40), "?".repeat(39))
    );
}

#[test]
fn command_fq_trim_q_paired_separated_with_singles() {
    // Pair1: both pass (R1 trimmed to 30). Pair2: R1 too short (fails),
    // R2 passes -> R2 goes to singles.
    let dir = tempfile::tempdir().unwrap();
    let r1_path = dir.path().join("R1.fq");
    let r2_path = dir.path().join("R2.fq");
    let good = "?".repeat(30);
    std::fs::write(
        &r1_path,
        format!(
            "@p1/1\n{}\n+\n{}{}\n@p2/1\nACGT\n+\n!!!!\n",
            "A".repeat(40),
            good,
            "!".repeat(10)
        ),
    )
    .unwrap();
    std::fs::write(
        &r2_path,
        format!(
            "@p1/2\n{}\n+\n{}\n@p2/2\n{}\n+\n{}\n",
            "C".repeat(40),
            "?".repeat(40),
            "G".repeat(40),
            "?".repeat(40)
        ),
    )
    .unwrap();

    let out1 = dir.path().join("out1.fq");
    let out2 = dir.path().join("out2.fq");
    let singles = dir.path().join("s.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            r1_path.to_str().unwrap(),
            r2_path.to_str().unwrap(),
            "-o",
            out1.to_str().unwrap(),
            "--outfile-2",
            out2.to_str().unwrap(),
            "--outfile-single",
            singles.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out1).unwrap(),
        format!("@p1/1\n{}\n+\n{good}\n", "A".repeat(30))
    );
    assert_eq!(
        std::fs::read_to_string(&out2).unwrap(),
        format!("@p1/2\n{}\n+\n{}\n", "C".repeat(40), "?".repeat(40))
    );
    assert_eq!(
        std::fs::read_to_string(&singles).unwrap(),
        format!("@p2/2\n{}\n+\n{}\n", "G".repeat(40), "?".repeat(40))
    );
}

#[test]
fn command_fq_trim_q_paired_interleaved() {
    // Without --outfile-2, passing pairs are interleaved into -o.
    let dir = tempfile::tempdir().unwrap();
    let r1_path = dir.path().join("R1.fq");
    let r2_path = dir.path().join("R2.fq");
    let good = "?".repeat(30);
    std::fs::write(
        &r1_path,
        format!(
            "@p1/1\n{}\n+\n{}{}\n@p2/1\nACGT\n+\n!!!!\n",
            "A".repeat(40),
            good,
            "!".repeat(10)
        ),
    )
    .unwrap();
    std::fs::write(
        &r2_path,
        format!(
            "@p1/2\n{}\n+\n{}\n@p2/2\n{}\n+\n{}\n",
            "C".repeat(40),
            "?".repeat(40),
            "G".repeat(40),
            "?".repeat(40)
        ),
    )
    .unwrap();

    let out = dir.path().join("interleaved.fq");
    let singles = dir.path().join("s.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            r1_path.to_str().unwrap(),
            r2_path.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--outfile-single",
            singles.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        format!(
            "@p1/1\n{}\n+\n{good}\n@p1/2\n{}\n+\n{}\n",
            "A".repeat(30),
            "C".repeat(40),
            "?".repeat(40)
        )
    );
    assert_eq!(
        std::fs::read_to_string(&singles).unwrap(),
        format!("@p2/2\n{}\n+\n{}\n", "G".repeat(40), "?".repeat(40))
    );
}

#[test]
fn command_fq_trim_q_mott_method() {
    let good = "?".repeat(30);
    let input = format!("@r1\n{}\n+\n{}{}\n", "A".repeat(40), good, "+".repeat(10));
    let file = write_temp(".fq", &input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--method",
            "mott",
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        format!("@r1\n{}\n+\n{good}\n", "A".repeat(30))
    );
}

#[test]
fn command_fq_trim_q_auto_detects_phred64() {
    // 'Y' (89) under +33 implies Q56 (>54) -> flips to +64, where Q25 < 20
    // trims the trailing '@' (Q0) half; with explicit +33 nothing is trimmed.
    let input = format!(
        "@r1\n{}\n+\n{}{}\n",
        "A".repeat(40),
        "Y".repeat(20),
        "@".repeat(20)
    );
    let file = write_temp(".fq", &input);
    let out_dir = tempfile::tempdir().unwrap();
    let out_auto = out_dir.path().join("out_auto.fq");
    let out_33 = out_dir.path().join("out_33.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out_auto.to_str().unwrap(),
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&out_auto).unwrap(),
        format!("@r1\n{}\n+\n{}\n", "A".repeat(20), "Y".repeat(20))
    );

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out_33.to_str().unwrap(),
            "--quality-base",
            "33",
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&out_33).unwrap(),
        format!(
            "@r1\n{}\n+\n{}{}\n",
            "A".repeat(40),
            "Y".repeat(20),
            "@".repeat(20)
        )
    );
}

#[test]
fn command_fq_trim_q_polyg_right() {
    // 10 trailing Gs are trimmed with --polyg-right 5, then length check passes.
    let input = format!(
        "@r1\n{}{}\n+\n{}\n",
        "A".repeat(30),
        "G".repeat(10),
        "?".repeat(40)
    );
    let file = write_temp(".fq", &input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--polyg-right",
            "5",
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        format!("@r1\n{}\n+\n{}\n", "A".repeat(30), "?".repeat(30))
    );
}

#[test]
fn command_fq_trim_q_gzipped_input() {
    let good = "?".repeat(30);
    let input = format!("@r1\n{}\n+\n{}{}\n", "A".repeat(40), good, "!".repeat(10));
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(input.as_bytes()).unwrap();
    let gz = encoder.finish().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("in.fq.gz");
    std::fs::write(&in_path, &gz).unwrap();
    let out = dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            in_path.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        format!("@r1\n{}\n+\n{good}\n", "A".repeat(30))
    );
}

#[test]
fn command_fq_trim_q_rejects_fasta_input() {
    let file = write_temp(".fa", ">r1\nACGTACGTACGTACGTACGT\n");
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("is not FASTQ"));
}

#[test]
fn command_fq_trim_q_rejects_invalid_quality() {
    // ASCII 32 (space) decodes below 0 for base 33.
    let input = format!("@r1\n{}\n+\n! {}\n", "A".repeat(40), "!".repeat(38));
    let file = write_temp(".fq", &input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("invalid quality character"));
}

#[test]
fn command_fq_trim_q_rejects_outfile2_with_single_input() {
    let file = write_temp(
        ".fq",
        "@r1\nACGTACGTACGTACGTACGTACGT\n+\n!!!!!!!!!!!!!!!!!!!!!!!!\n",
    );
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            file.path().to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--outfile-2",
            out_dir.path().join("out2.fq").to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains(
            "--outfile-2 requires two input files",
        ));
}

#[test]
fn command_fq_trim_q_output_same_as_input_rejected() {
    let temp = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    let path = temp.path().to_str().unwrap().to_string();
    std::fs::write(
        &path,
        "@SEQ\nACGTACGTACGTACGTACGTACGT\n+\n!!!!!!!!!!!!!!!!!!!!!!!!\n",
    )
    .unwrap();

    PgrCmd::new()
        .args(&["fq", "trim-q", &path, "-o", &path])
        .assert()
        .failure()
        .stderr(predicates::str::contains("is also an input file"));
}

#[test]
fn command_fq_trim_q_outputs_distinct() {
    let dir = tempfile::tempdir().unwrap();
    let r1_path = dir.path().join("R1.fq");
    let r2_path = dir.path().join("R2.fq");
    std::fs::write(
        &r1_path,
        format!("@p1/1\n{}\n+\n{}\n", "A".repeat(40), "?".repeat(40)),
    )
    .unwrap();
    std::fs::write(
        &r2_path,
        format!("@p1/2\n{}\n+\n{}\n", "C".repeat(40), "?".repeat(40)),
    )
    .unwrap();
    let out = dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-q",
            r1_path.to_str().unwrap(),
            r2_path.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--outfile-2",
            out.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("output files must be distinct"));
}
