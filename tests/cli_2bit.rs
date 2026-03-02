#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use pgr::libs::twobit::TwoBitFile;
use std::fs;
use tempfile::TempDir;

#[test]
fn test_fa_to_2bit() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let output = temp.path().join("out.2bit");

    fs::write(&input, ">seq1\nACGT\n>seq2\nNNNN\n")?;

    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    assert!(output.exists());

    let mut tb = TwoBitFile::open(&output)?;
    let names = tb.get_sequence_names();
    assert_eq!(names.len(), 2);
    assert!(names.contains(&"seq1".to_string()));
    assert!(names.contains(&"seq2".to_string()));

    let seq1 = tb.read_sequence("seq1", None, None, false)?;
    assert_eq!(seq1, "ACGT");

    let seq2 = tb.read_sequence("seq2", None, None, false)?;
    assert_eq!(seq2, "NNNN");

    Ok(())
}

#[test]
fn test_fa_to_2bit_strip_version() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test_ver.fa");
    let output = temp.path().join("out_ver.2bit");

    fs::write(&input, ">NM_001.1\nACGT\n")?;

    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--strip-version",
        ])
        .run();

    let tb = TwoBitFile::open(&output)?;
    let names = tb.get_sequence_names();
    assert_eq!(names[0], "NM_001");

    Ok(())
}

#[test]
fn test_fa_to_2bit_mask() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test_mask.fa");
    let output = temp.path().join("out_mask.2bit");

    fs::write(&input, ">seq1\nacgtACGT\n")?;

    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let mut tb = TwoBitFile::open(&output)?;
    let seq_masked = tb.read_sequence("seq1", None, None, false)?;
    assert_eq!(seq_masked, "acgtACGT");

    let seq_unmasked = tb.read_sequence("seq1", None, None, true)?;
    assert_eq!(seq_unmasked, "ACGTACGT");

    Ok(())
}

#[test]
fn test_2bit_to_fa_basic() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    let output = temp.path().join("out.fa");

    fs::write(&input, ">seq1\nACGT\n>seq2\nNNNN\n")?;

    // Create 2bit first
    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            input.to_str().unwrap(),
            "-o",
            twobit.to_str().unwrap(),
        ])
        .run();

    // Convert back to FASTA
    PgrCmd::new()
        .args(&[
            "2bit",
            "to-fa",
            twobit.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let content = fs::read_to_string(&output)?;
    // Order might differ, but content should match.
    // >seq1
    // ACGT
    // >seq2
    // NNNN
    assert!(content.contains(">seq1\nACGT"));
    assert!(content.contains(">seq2\nNNNN"));

    Ok(())
}

#[test]
fn test_2bit_range_basic() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    let output = temp.path().join("out.fa");

    fs::write(&input, ">seq1\nACGTACGT\n")?;

    // Create 2bit
    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            input.to_str().unwrap(),
            "-o",
            twobit.to_str().unwrap(),
        ])
        .run();

    // Extract range 2-5 (CGTA) - 1-based
    // 01234567
    // ACGTACGT
    //  CGTA
    PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            twobit.to_str().unwrap(),
            "seq1:2-5",
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let content = fs::read_to_string(&output)?;
    assert!(content.contains(">seq1:2-5\nCGTA"));

    // Extract negative strand
    // seq1:2-5 is CGTA. RevComp: TACG.
    let output_neg = temp.path().join("out_neg.fa");
    PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            twobit.to_str().unwrap(),
            "seq1(-):2-5",
            "-o",
            output_neg.to_str().unwrap(),
        ])
        .run();

    let content_neg = fs::read_to_string(&output_neg)?;
    assert!(content_neg.contains(">seq1(-):2-5\nTACG"));

    Ok(())
}

#[test]
fn test_2bit_range_rgfile() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    let list = temp.path().join("ranges.txt");
    let output = temp.path().join("out.fa");

    fs::write(&input, ">seq1\nACGT\n>seq2\nTGCA\n")?;
    // Request seq2 (entire sequence) and seq1:1-2
    fs::write(&list, "seq2\nseq1:1-2\n")?;

    // Create 2bit
    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            input.to_str().unwrap(),
            "-o",
            twobit.to_str().unwrap(),
        ])
        .run();

    // Extract ranges
    PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            twobit.to_str().unwrap(),
            "-r",
            list.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let content = fs::read_to_string(&output)?;
    assert!(content.contains(">seq2\nTGCA"));
    assert!(content.contains(">seq1:1-2\nAC"));

    Ok(())
}

#[test]
fn test_2bit_to_fa_mask() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    let output_masked = temp.path().join("masked.fa");
    let output_unmasked = temp.path().join("unmasked.fa");

    fs::write(&input, ">seq1\nacgtACGT\n")?;

    // Create 2bit
    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            input.to_str().unwrap(),
            "-o",
            twobit.to_str().unwrap(),
        ])
        .run();

    // Default (masked)
    PgrCmd::new()
        .args(&[
            "2bit",
            "to-fa",
            twobit.to_str().unwrap(),
            "-o",
            output_masked.to_str().unwrap(),
        ])
        .run();

    let content_masked = fs::read_to_string(&output_masked)?;
    assert!(content_masked.contains("acgtACGT"));

    // No mask
    PgrCmd::new()
        .args(&[
            "2bit",
            "to-fa",
            twobit.to_str().unwrap(),
            "--no-mask",
            "-o",
            output_unmasked.to_str().unwrap(),
        ])
        .run();

    let content_unmasked = fs::read_to_string(&output_unmasked)?;
    assert!(content_unmasked.contains("ACGTACGT"));

    Ok(())
}

#[test]
fn test_2bit_size() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");

    fs::write(&input, ">seq1\nACGT\n>seq2\nNNNN\n")?;

    // Create 2bit file first
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    cmd.arg("fa")
        .arg("to-2bit")
        .arg(&input)
        .arg("-o")
        .arg(&twobit);
    cmd.assert().success();

    // Run 2bit size
    let (stdout, _) = PgrCmd::new()
        .args(&["2bit", "size", twobit.to_str().unwrap()])
        .run();
    assert!(stdout.contains("seq1\t4"));
    assert!(stdout.contains("seq2\t4"));

    Ok(())
}

#[test]
fn test_2bit_size_flags() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test_flags.fa");
    let twobit = temp.path().join("test_flags.2bit");

    // seq1: 12 bases, Ns at 4-8 (4 Ns). ACGT NNNN ACGT. Size 12. No-Ns: 8.
    // seq2: 4 bases. acgt. Size 4. No-Ns: 4. Mask: 0-4.
    fs::write(&input, ">seq1\nACGTNNNNACGT\n>seq2\nacgt\n")?;

    // Create 2bit
    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            input.to_str().unwrap(),
            "-o",
            twobit.to_str().unwrap(),
        ])
        .run();

    // Test --no-ns
    let (stdout, _) = PgrCmd::new()
        .args(&["2bit", "size", twobit.to_str().unwrap(), "--no-ns"])
        .run();
    assert!(stdout.contains("seq1\t8"));
    assert!(stdout.contains("seq2\t4"));

    Ok(())
}

#[test]
fn test_2bit_size_multiple() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input1 = temp.path().join("test1.fa");
    let input2 = temp.path().join("test2.fa");
    let twobit1 = temp.path().join("test1.2bit");
    let twobit2 = temp.path().join("test2.2bit");

    fs::write(&input1, ">seq1\nACGT\n")?;
    fs::write(&input2, ">seq2\nTGCA\n")?;

    // Create 2bit files
    for (inp, out) in [(&input1, &twobit1), (&input2, &twobit2)] {
        PgrCmd::new()
            .args(&[
                "fa",
                "to-2bit",
                inp.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ])
            .run();
    }

    // Run 2bit size with multiple inputs
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "2bit",
            "size",
            twobit1.to_str().unwrap(),
            twobit2.to_str().unwrap(),
        ])
        .run();
    assert!(stdout.contains("seq1\t4"));
    assert!(stdout.contains("seq2\t4"));

    Ok(())
}

#[test]
fn test_2bit_range_seqlist1_file() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");
    let list = std::path::Path::new(&manifest_dir).join("tests/2bit/input/seqlist1");

    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");

    PgrCmd::new()
        .args(&[
            "2bit",
            "range",
            input.to_str().unwrap(),
            "--rgfile",
            list.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output)?;

    assert!(output_content.contains(">noLower"));
    assert!(output_content.contains(">startLower"));
    assert!(output_content.contains(">endLower"));
    assert!(!output_content.contains(">manyLower"));

    Ok(())
}

#[test]
fn test_2bit_masked() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input_mask = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");
    let input_n = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testN.2bit");

    let temp = TempDir::new()?;
    let out_mask = temp.path().join("out_mask.txt");
    let out_n = temp.path().join("out_n.txt");

    // 1. testMask.2bit
    PgrCmd::new()
        .args(&[
            "2bit",
            "masked",
            input_mask.to_str().unwrap(),
            "-o",
            out_mask.to_str().unwrap(),
        ])
        .run();

    let content_mask = fs::read_to_string(&out_mask)?;

    // allLower is masked. It has 12 bases.
    assert!(content_mask.contains("allLower:1-12"));
    // noLower should not be in output
    assert!(!content_mask.contains("noLower"));

    // 2. testN.2bit with --gap
    PgrCmd::new()
        .args(&[
            "2bit",
            "masked",
            input_n.to_str().unwrap(),
            "--gap",
            "-o",
            out_n.to_str().unwrap(),
        ])
        .run();

    let content_n = fs::read_to_string(&out_n)?;

    // startN: NANNAANNNAAA
    // Ns at: 1, 3-4, 7-9
    assert!(content_n.contains("startN:1"));
    assert!(content_n.contains("startN:3-4"));
    assert!(content_n.contains("startN:7-9"));

    Ok(())
}

#[test]
fn test_2bit_range_legacy_cases() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");

    let temp = TempDir::new()?;

    // Helper to test range and expected sequence
    let test_range = |start: usize, end: usize, expected: &str| -> anyhow::Result<()> {
        let out_name = format!("out_{}_{}.fa", start, end);
        let out_path = temp.path().join(&out_name);

        let range_str = format!("manyLower:{}-{}", start, end);

        PgrCmd::new()
            .args(&[
                "2bit",
                "range",
                input.to_str().unwrap(),
                &range_str,
                "-o",
                out_path.to_str().unwrap(),
            ])
            .run();

        let content = fs::read_to_string(&out_path)?;
        if !content.contains(expected) {
            anyhow::bail!(
                "Failed for {}: expected {}, got {}",
                range_str,
                expected,
                content
            );
        }
        Ok(())
    };

    // Test cases from original test
    test_range(1, 11, "aCCggTTaCg")?;
    test_range(2, 10, "CCggTTaC")?;
    test_range(3, 9, "CggTTa")?;
    test_range(4, 8, "ggTT")?;
    test_range(5, 6, "g")?;
    test_range(5, 7, "gT")?;
    test_range(6, 7, "T")?;
    test_range(7, 8, "T")?;
    test_range(8, 9, "a")?;

    Ok(())
}

#[test]
fn test_2bit_compat_mask_file() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");

    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");

    PgrCmd::new()
        .args(&[
            "2bit",
            "to-fa",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output)?;

    // Check for sequence names
    assert!(output_content.contains(">allLower"));
    assert!(output_content.contains(">endLower"));
    assert!(output_content.contains(">manyLower"));
    assert!(output_content.contains(">noLower"));
    assert!(output_content.contains(">startLower"));

    // Check masking (lowercase)
    // allLower should be all lowercase
    // We can't easily check full content without reading exact expectation,
    // but we can check if it contains lowercase letters.
    assert!(output_content.chars().any(|c| c.is_lowercase()));

    Ok(())
}

#[test]
fn test_2bit_compat_n_file() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testN.2bit");

    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");

    PgrCmd::new()
        .args(&[
            "2bit",
            "to-fa",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output)?;

    assert!(output_content.contains(">startN"));
    assert!(output_content.contains("NANNAANNNAAA"));

    assert!(output_content.contains(">startNonN"));
    assert!(output_content.contains("ANAANNAAANNN"));

    Ok(())
}

#[test]
fn test_2bit_range_complex() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/index/final.contigs.2bit");

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

    Ok(())
}

#[test]
fn test_2bit_range_r_complex() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/index/final.contigs.2bit");
    let rg_file = std::path::Path::new(&manifest_dir).join("tests/index/sample.rg");

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

    Ok(())
}

#[test]
fn test_2bit_some() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test_some.fa");
    let twobit = temp.path().join("test_some.2bit");
    let list = temp.path().join("list.txt");
    let output = temp.path().join("out_some.fa");
    let output_inv = temp.path().join("out_some_inv.fa");

    // seq1: ACGT (4)
    // seq2: TGCA (4)
    // seq3: NNNN (4)
    fs::write(&input, ">seq1\nACGT\n>seq2\nTGCA\n>seq3\nNNNN\n")?;
    fs::write(&list, "seq1\nseq3\n")?;

    // Create 2bit
    PgrCmd::new()
        .args(&[
            "fa",
            "to-2bit",
            input.to_str().unwrap(),
            "-o",
            twobit.to_str().unwrap(),
        ])
        .run();

    // Test some
    PgrCmd::new()
        .args(&[
            "2bit",
            "some",
            twobit.to_str().unwrap(),
            list.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output)?;
    assert!(output_content.contains(">seq1"));
    assert!(output_content.contains("ACGT"));
    assert!(output_content.contains(">seq3"));
    assert!(output_content.contains("NNNN"));
    assert!(!output_content.contains(">seq2"));

    // Test some invert
    PgrCmd::new()
        .args(&[
            "2bit",
            "some",
            twobit.to_str().unwrap(),
            list.to_str().unwrap(),
            "-i",
            "-o",
            output_inv.to_str().unwrap(),
        ])
        .run();

    let output_inv_content = fs::read_to_string(&output_inv)?;
    assert!(output_inv_content.contains(">seq2"));
    assert!(output_inv_content.contains("TGCA"));
    assert!(!output_inv_content.contains(">seq1"));
    assert!(!output_inv_content.contains(">seq3"));

    Ok(())
}

#[test]
fn test_2bit_size_doc_consistency() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let tests_pgr = std::path::Path::new(&manifest_dir).join("tests/pgr");
    let temp = TempDir::new()?;

    for name in ["pseudocat", "pseudopig"] {
        let fa_path = tests_pgr.join(format!("{}.fa", name));
        let twobit_path = tests_pgr.join(format!("{}.2bit", name));

        // Ensure inputs exist
        assert!(fa_path.exists(), "Test file not found: {:?}", fa_path);
        assert!(
            twobit_path.exists(),
            "Test file not found: {:?}",
            twobit_path
        );

        let fa_sizes_path = temp.path().join(format!("{}.fa.sizes", name));
        let twobit_sizes_path = temp.path().join(format!("{}.2bit.sizes", name));

        // 1. Run pgr fa size
        PgrCmd::new()
            .args(&[
                "fa",
                "size",
                fa_path.to_str().unwrap(),
                "-o",
                fa_sizes_path.to_str().unwrap(),
            ])
            .run();

        // 2. Run pgr 2bit size
        PgrCmd::new()
            .args(&[
                "2bit",
                "size",
                twobit_path.to_str().unwrap(),
                "-o",
                twobit_sizes_path.to_str().unwrap(),
            ])
            .run();

        // 3. Compare
        let content_fa = fs::read_to_string(&fa_sizes_path)?;
        let content_2bit = fs::read_to_string(&twobit_sizes_path)?;

        assert_eq!(
            content_fa, content_2bit,
            "pgr fa size and pgr 2bit size output should be identical for {}",
            name
        );
    }

    Ok(())
}
