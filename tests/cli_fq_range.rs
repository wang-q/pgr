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
fn command_fq_range_plain_by_name() {
    let input = "@r1 comment\nACGTACGT\n+\n!!!!!!!!\n@r2\nTGCA\n+\n!!!!\n";
    let file = write_temp(".fq", input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "range",
            file.path().to_str().unwrap(),
            "r2",
            "r1",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    // Output preserves query order.
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@r2\nTGCA\n+\n!!!!\n@r1 comment\nACGTACGT\n+\n!!!!!!!!\n"
    );
    // The .loc sidecar was created.
    let loc_path = format!("{}.loc", file.path().display());
    assert!(std::path::Path::new(&loc_path).is_file());
}

#[test]
fn command_fq_range_subsequence() {
    // read1:3-6 cuts both sequence and quality (1-based inclusive).
    let input = "@read1\nACGTACGTACGT\n+\n!!!!IIIIIIII\n";
    let file = write_temp(".fq", input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "range",
            file.path().to_str().unwrap(),
            "read1:3-6",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@read1\nGTAC\n+\n!!II\n"
    );
}

#[test]
fn command_fq_range_paired_end_by_name() {
    let r1 = write_temp(
        ".fq",
        "@read1/1\nACGTACGT\n+\n!!!!!!!!\n@read2/1\nTGCA\n+\n!!!!\n",
    );
    let r2 = write_temp(
        ".fq",
        "@read1/2\nTGCATGCA\n+\n!!!!!!!!\n@read2/2\nACGT\n+\n!!!!\n",
    );
    let out_dir = tempfile::tempdir().unwrap();
    let out1 = out_dir.path().join("r1.out.fq");
    let out2 = out_dir.path().join("r2.out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "range",
            r1.path().to_str().unwrap(),
            "--mate",
            r2.path().to_str().unwrap(),
            "read1",
            "-o",
            out1.to_str().unwrap(),
            "--outfile-2",
            out2.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out1).unwrap(),
        "@read1/1\nACGTACGT\n+\n!!!!!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(&out2).unwrap(),
        "@read1/2\nTGCATGCA\n+\n!!!!!!!!\n"
    );
    // Both .loc sidecars were created.
    assert!(std::path::Path::new(&format!("{}.loc", r1.path().display())).is_file());
    assert!(std::path::Path::new(&format!("{}.loc", r2.path().display())).is_file());
}

#[test]
fn command_fq_range_paired_end_subsequence() {
    let r1 = write_temp(".fq", "@read1/1\nACGTACGTACGT\n+\n!!!!IIIIIIII\n");
    let r2 = write_temp(".fq", "@read1/2\nTGCATGCATGCA\n+\nJJJJJJJJJJJJ\n");
    let out_dir = tempfile::tempdir().unwrap();
    let out1 = out_dir.path().join("r1.out.fq");
    let out2 = out_dir.path().join("r2.out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "range",
            r1.path().to_str().unwrap(),
            "--mate",
            r2.path().to_str().unwrap(),
            "read1:3-6",
            "-o",
            out1.to_str().unwrap(),
            "--outfile-2",
            out2.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out1).unwrap(),
        "@read1/1\nGTAC\n+\n!!II\n"
    );
    assert_eq!(
        std::fs::read_to_string(&out2).unwrap(),
        "@read1/2\nCATG\n+\nJJJJ\n"
    );
}

#[test]
fn command_fq_range_rejects_outfile_2_without_pair() {
    let file = write_temp(".fq", "@r1\nACGT\n+\n!!!!\n");
    let out_dir = tempfile::tempdir().unwrap();
    let out2 = out_dir.path().join("r2.out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "range",
            file.path().to_str().unwrap(),
            "r1",
            "--outfile-2",
            out2.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains(
            "--outfile-2 requires two input files",
        ));
}

#[test]
fn command_fq_range_rejects_plain_gzip() {
    // Plain gzip has no block index, so random access cannot seek by plain
    // offset; only plain text and BGZF are supported.
    let input = "@r1\nACGTACGT\n+\n!!!!!!!!\n@r2\nTGCA\n+\n!!!!\n";
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
            "range",
            in_path.to_str().unwrap(),
            "r1",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains(
            "only plain text and BGZF (.gz) files support range extraction",
        ));
}

#[test]
fn command_fq_range_bgzf_input() {
    let input = "@r1\nACGTACGT\n+\n!!!!!!!!\n@r2\nTGCA\n+\n!!!!\n";
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("in.fq.gz");
    {
        let mut w =
            pgr::libs::bgzf::BgzfWriter::new(std::fs::File::create(&in_path).unwrap()).unwrap();
        w.write_all(input.as_bytes()).unwrap();
        w.finish().unwrap();
    }
    pgr::libs::bgzf::build_gzi_index(in_path.to_str().unwrap()).unwrap();
    let out = dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "range",
            in_path.to_str().unwrap(),
            "r2",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@r2\nTGCA\n+\n!!!!\n"
    );
}

#[test]
fn command_fq_range_interleaved_pair_returned_in_order() {
    // Two reads share the name `read1` (interleaved pair): both are returned.
    let input = "@read1\nACGT\n+\n!!!!\n@read1\nTGCA\n+\n!!!!\n";
    let file = write_temp(".fq", input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "range",
            file.path().to_str().unwrap(),
            "read1",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@read1\nACGT\n+\n!!!!\n@read1\nTGCA\n+\n!!!!\n"
    );
}

#[test]
fn command_fq_range_pair_suffix_matches() {
    // CASAVA-style names read1/1, read1/2: querying `read1` returns both.
    let input = "@read1/1\nACGT\n+\n!!!!\n@read1/2\nTGCA\n+\n!!!!\n";
    let file = write_temp(".fq", input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "range",
            file.path().to_str().unwrap(),
            "read1",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@read1/1\nACGT\n+\n!!!!\n@read1/2\nTGCA\n+\n!!!!\n"
    );
}

#[test]
fn command_fq_range_missing_name_warns() {
    let file = write_temp(".fq", "@r1\nACGT\n+\n!!!!\n");
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "range",
            file.path().to_str().unwrap(),
            "nope",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stderr(predicates::str::contains(
            "not found in the .loc index file",
        ));

    assert_eq!(std::fs::read_to_string(&out).unwrap(), "");
}

#[test]
fn command_fq_range_output_same_as_input_rejected() {
    let temp = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    let path = temp.path().to_str().unwrap().to_string();
    std::fs::write(&path, "@r1\nACGT\n+\n!!!!\n").unwrap();

    PgrCmd::new()
        .args(&["fq", "range", &path, "r1", "-o", &path])
        .assert()
        .failure()
        .stderr(predicates::str::contains("is also an input file"));
}

#[test]
fn command_fq_range_rejects_fasta_input() {
    let file = write_temp(".fa", ">r1\nACGTACGTACGTACGTACGT\n");
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "range",
            file.path().to_str().unwrap(),
            "r1",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("is not FASTQ"));
}
