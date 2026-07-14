#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use pgr::libs::fmt::twobit::TwoBitFile;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

/// Return the absolute path to a fixture in `tests/2bit/input`.
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/2bit/input")
        .join(name)
}

/// Read a FASTA file produced by a CLI command.
fn read_fasta(path: &std::path::Path) -> String {
    fs::read_to_string(path).unwrap()
}

#[test]
fn test_fa_to_2bit() {
    let temp = TempDir::new().unwrap();

    // Basic round-trip.
    let basic_fa = temp.path().join("basic.fa");
    let basic_2bit = temp.path().join("basic.2bit");
    fs::write(&basic_fa, ">seq1\nACGT\n>seq2\nNNNN\n").unwrap();
    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            basic_fa.to_str().unwrap(),
            "-o",
            basic_2bit.to_str().unwrap(),
        ])
        .run();
    let mut tb = TwoBitFile::open(&basic_2bit).unwrap();
    let names = tb.get_sequence_names();
    assert_eq!(names.len(), 2);
    assert!(names.contains(&"seq1".to_string()));
    assert!(names.contains(&"seq2".to_string()));
    assert_eq!(tb.read_sequence("seq1", None, None, false).unwrap(), "ACGT");
    assert_eq!(tb.read_sequence("seq2", None, None, false).unwrap(), "NNNN");

    // Strip version.
    let ver_fa = temp.path().join("ver.fa");
    let ver_2bit = temp.path().join("ver.2bit");
    fs::write(&ver_fa, ">NM_001.1\nACGT\n").unwrap();
    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            ver_fa.to_str().unwrap(),
            "-o",
            ver_2bit.to_str().unwrap(),
            "--strip-version",
        ])
        .run();
    let tb = TwoBitFile::open(&ver_2bit).unwrap();
    assert_eq!(tb.get_sequence_names()[0], "NM_001");

    // Mask preservation.
    let mask_fa = temp.path().join("mask.fa");
    let mask_2bit = temp.path().join("mask.2bit");
    fs::write(&mask_fa, ">seq1\nacgtACGT\n").unwrap();
    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            mask_fa.to_str().unwrap(),
            "-o",
            mask_2bit.to_str().unwrap(),
        ])
        .run();
    let mut tb = TwoBitFile::open(&mask_2bit).unwrap();
    assert_eq!(
        tb.read_sequence("seq1", None, None, false).unwrap(),
        "acgtACGT"
    );
    assert_eq!(
        tb.read_sequence("seq1", None, None, true).unwrap(),
        "ACGTACGT"
    );
}

#[test]
fn test_2bit_to_fa_mask() {
    let temp = TempDir::new().unwrap();
    let masked = temp.path().join("masked.fa");
    let unmasked = temp.path().join("unmasked.fa");

    PgrCmd::new()
        .args(&[
            "2bit",
            "to-fa",
            fixture("mask.2bit").to_str().unwrap(),
            "-o",
            masked.to_str().unwrap(),
        ])
        .run();
    assert!(read_fasta(&masked).contains("acgtACGT"));

    PgrCmd::new()
        .args(&[
            "2bit",
            "to-fa",
            fixture("mask.2bit").to_str().unwrap(),
            "--no-mask",
            "-o",
            unmasked.to_str().unwrap(),
        ])
        .run();
    assert!(read_fasta(&unmasked).contains("ACGTACGT"));
}

#[test]
fn test_2bit_range_basic() {
    let temp = TempDir::new().unwrap();
    let out = temp.path().join("out.fa");
    let out_neg = temp.path().join("out_neg.fa");

    // seq1: ACGTACGT, range 2-5 (1-based) -> CGTA.
    PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            fixture("range.2bit").to_str().unwrap(),
            "seq1:2-5",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(read_fasta(&out).contains(">seq1:2-5\nCGTA"));

    // Negative strand: revcomp(CGTA) = TACG.
    PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            fixture("range.2bit").to_str().unwrap(),
            "seq1(-):2-5",
            "-o",
            out_neg.to_str().unwrap(),
        ])
        .run();
    assert!(read_fasta(&out_neg).contains(">seq1(-):2-5\nTACG"));
}

#[test]
fn test_2bit_range_rgfile() {
    let temp = TempDir::new().unwrap();
    let out = temp.path().join("out.fa");

    PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            fixture("range.2bit").to_str().unwrap(),
            "--rgfile",
            fixture("ranges.txt").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    let content = read_fasta(&out);
    assert!(content.contains(">seq2\nTGCA"));
    assert!(content.contains(">seq1:1-2\nAC"));
}

#[test]
fn test_2bit_size() {
    // Default size.
    let (stdout, _) = PgrCmd::new()
        .args(&["2bit", "size", fixture("basic.2bit").to_str().unwrap()])
        .run();
    assert!(stdout.contains("seq1\t4"));
    assert!(stdout.contains("seq2\t4"));

    // --no-ns flag.
    // seq1: ACGTNNNNACGT (12 bp, 4 Ns) -> 8 non-Ns.
    // seq2: acgt -> 4 non-Ns.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "2bit",
            "size",
            fixture("flags.2bit").to_str().unwrap(),
            "--no-ns",
        ])
        .run();
    assert!(stdout.contains("seq1\t8"));
    assert!(stdout.contains("seq2\t4"));

    // Multiple inputs with disjoint sequence names.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "2bit",
            "size",
            fixture("multi1.2bit").to_str().unwrap(),
            fixture("multi2.2bit").to_str().unwrap(),
        ])
        .run();
    assert!(stdout.contains("seq1\t4"));
    assert!(stdout.contains("seq2\t4"));
}

#[test]
fn test_2bit_some() {
    let temp = TempDir::new().unwrap();
    let out = temp.path().join("out_some.fa");
    let out_inv = temp.path().join("out_some_inv.fa");

    PgrCmd::new()
        .args(&[
            "2bit",
            "some",
            fixture("some.2bit").to_str().unwrap(),
            fixture("list.txt").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    let content = read_fasta(&out);
    assert!(content.contains(">seq1"));
    assert!(content.contains("ACGT"));
    assert!(content.contains(">seq3"));
    assert!(content.contains("NNNN"));
    assert!(!content.contains(">seq2"));

    PgrCmd::new()
        .args(&[
            "2bit",
            "some",
            fixture("some.2bit").to_str().unwrap(),
            fixture("list.txt").to_str().unwrap(),
            "-i",
            "-o",
            out_inv.to_str().unwrap(),
        ])
        .run();
    let content = read_fasta(&out_inv);
    assert!(content.contains(">seq2"));
    assert!(content.contains("TGCA"));
    assert!(!content.contains(">seq1"));
    assert!(!content.contains(">seq3"));
}

#[test]
fn test_2bit_range_seqlist1_file() {
    let temp = TempDir::new().unwrap();
    let output = temp.path().join("out.fa");

    PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            fixture("testMask.2bit").to_str().unwrap(),
            "--rgfile",
            fixture("seqlist1").to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let content = read_fasta(&output);
    assert!(content.contains(">noLower"));
    assert!(content.contains(">startLower"));
    assert!(content.contains(">endLower"));
    assert!(!content.contains(">manyLower"));
}

#[test]
fn test_2bit_masked() {
    let temp = TempDir::new().unwrap();
    let out_mask = temp.path().join("out_mask.txt");
    let out_n = temp.path().join("out_n.txt");

    PgrCmd::new()
        .args(&[
            "2bit",
            "masked",
            fixture("testMask.2bit").to_str().unwrap(),
            "-o",
            out_mask.to_str().unwrap(),
        ])
        .run();
    let content_mask = read_fasta(&out_mask);
    assert!(content_mask.contains("allLower:1-12"));
    assert!(!content_mask.contains("noLower"));

    PgrCmd::new()
        .args(&[
            "2bit",
            "masked",
            fixture("testN.2bit").to_str().unwrap(),
            "--gap",
            "-o",
            out_n.to_str().unwrap(),
        ])
        .run();
    let content_n = read_fasta(&out_n);
    // startN: NANNAANNNAAA, Ns at 1, 3-4, 7-9.
    assert!(content_n.contains("startN:1"));
    assert!(content_n.contains("startN:3-4"));
    assert!(content_n.contains("startN:7-9"));
}

#[test]
fn test_2bit_range_legacy_cases() {
    let temp = TempDir::new().unwrap();

    let test_range = |start: usize, end: usize, expected: &str| {
        let out_path = temp.path().join(format!("out_{}_{}.fa", start, end));
        let range_str = format!("manyLower:{}-{}", start, end);

        PgrCmd::new()
            .args(&[
                "2bit",
                "range",
                fixture("testMask.2bit").to_str().unwrap(),
                &range_str,
                "-o",
                out_path.to_str().unwrap(),
            ])
            .run();

        let content = read_fasta(&out_path);
        assert!(
            content.contains(expected),
            "Failed for {}: expected {}, got {}",
            range_str,
            expected,
            content
        );
    };

    test_range(1, 11, "aCCggTTaCg");
    test_range(2, 10, "CCggTTaC");
    test_range(3, 9, "CggTTa");
    test_range(4, 8, "ggTT");
    test_range(5, 6, "g");
    test_range(5, 7, "gT");
    test_range(6, 7, "T");
    test_range(7, 8, "T");
    test_range(8, 9, "a");
}

#[test]
fn test_2bit_compat_mask_file() {
    let temp = TempDir::new().unwrap();
    let output = temp.path().join("out.fa");

    PgrCmd::new()
        .args(&[
            "2bit",
            "to-fa",
            fixture("testMask.2bit").to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let content = read_fasta(&output);
    assert!(content.contains(">allLower"));
    assert!(content.contains(">endLower"));
    assert!(content.contains(">manyLower"));
    assert!(content.contains(">noLower"));
    assert!(content.contains(">startLower"));
    assert!(content.chars().any(|c| c.is_lowercase()));
}

#[test]
fn test_2bit_compat_n_file() {
    let temp = TempDir::new().unwrap();
    let output = temp.path().join("out.fa");

    PgrCmd::new()
        .args(&[
            "2bit",
            "to-fa",
            fixture("testN.2bit").to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let content = read_fasta(&output);
    assert!(content.contains(">startN"));
    assert!(content.contains("NANNAANNNAAA"));
    assert!(content.contains(">startNonN"));
    assert!(content.contains("ANAANNAAANNN"));
}

#[test]
fn test_2bit_range_complex() {
    let input = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/index/final.contigs.2bit");
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            input.to_str().unwrap(),
            "k81_130",
            "k81_130:11-20",
            "k81_170:304-323",
            "k81_170(-):1-20",
            "k81_158:70001-70020",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains(">k81_130\nAGTTTCAACT"));
    assert!(stdout.contains(">k81_130:11-20\nGGTGAATCAA\n"));
    assert!(stdout.contains(">k81_170:304-323\nAGTTAAAAACCTGATTTATT\n"));
    assert!(stdout.contains(">k81_170(-):1-20\nATTAACCTGTTGTAGGTGTT\n"));
    assert!(stdout.contains(">k81_158:70001-70020\nTGGCTATAACCTAATTTTGT\n"));
}

#[test]
fn test_2bit_range_r_complex() {
    let input = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/index/final.contigs.2bit");
    let rg_file = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/index/sample.rg");

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            input.to_str().unwrap(),
            "-r",
            rg_file.to_str().unwrap(),
        ])
        .run();

    assert_eq!(stdout.lines().count(), 12);
    assert!(stdout.contains(">k81_130:11-20\nGGTGAATCAA\n"));
}

#[test]
fn test_2bit_range_invalid_inverted() {
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            fixture("range.2bit").to_str().unwrap(),
            "seq1:5-2",
        ])
        .run_fail();
    assert!(stderr.contains("range start must not be greater than end"));
}

#[test]
fn test_2bit_range_invalid_zero_coordinate() {
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            fixture("range.2bit").to_str().unwrap(),
            "seq1:0-3",
        ])
        .run_fail();
    assert!(
        stderr.contains("invalid range") || stderr.contains("range coordinates must be positive"),
        "expected invalid range error, got: {}",
        stderr
    );
}

#[test]
fn test_2bit_size_doc_consistency() {
    let tests_pgr = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/pgr");
    let temp = TempDir::new().unwrap();

    for name in ["pseudocat", "pseudopig"] {
        let fa_path = tests_pgr.join(format!("{}.fa", name));
        let twobit_path = tests_pgr.join(format!("{}.2bit", name));
        assert!(fa_path.exists(), "Test file not found: {:?}", fa_path);
        assert!(
            twobit_path.exists(),
            "Test file not found: {:?}",
            twobit_path
        );

        let fa_sizes = temp.path().join(format!("{}.fa.sizes", name));
        let twobit_sizes = temp.path().join(format!("{}.2bit.sizes", name));

        PgrCmd::new()
            .args(&[
                "fa",
                "size",
                fa_path.to_str().unwrap(),
                "-o",
                fa_sizes.to_str().unwrap(),
            ])
            .run();
        PgrCmd::new()
            .args(&[
                "2bit",
                "size",
                twobit_path.to_str().unwrap(),
                "-o",
                twobit_sizes.to_str().unwrap(),
            ])
            .run();

        assert_eq!(
            read_fasta(&fa_sizes),
            read_fasta(&twobit_sizes),
            "pgr fa size and pgr 2bit size output should be identical for {}",
            name
        );
    }
}
