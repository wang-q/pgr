#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::Read;

fn read_gz(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    let mut dec = flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap());
    dec.read_to_end(&mut out).unwrap();
    out
}

fn lambda(args: &[&str], out: &str, outu: Option<&str>, ihist: Option<&str>) {
    let mut full: Vec<&str> = vec![
        "fq",
        "merge",
        "tests/bbtools/Lambda/R1.2k.fq.gz",
        "tests/bbtools/Lambda/R2.2k.fq.gz",
        "-o",
        out,
    ];
    if let Some(u) = outu {
        full.push("--outu");
        full.push(u);
    }
    if let Some(h) = ihist {
        full.push("--ihist");
        full.push(h);
    }
    full.extend_from_slice(args);
    PgrCmd::new().args(&full).assert().success();
}

#[test]
fn command_fq_merge_join_matches_bbtools_golden() {
    // `bbmerge.sh ... strict` (net filter on): merged + unmerged + ihist.
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("merged.fq");
    let outu = out_dir.path().join("unmerged.fq");
    let ihist = out_dir.path().join("ihist2.txt");
    lambda(
        &[
            "--strict",
            "--net",
            "tests/bbtools/Lambda/golden/bbmerge.bbnet",
        ],
        out.to_str().unwrap(),
        Some(outu.to_str().unwrap()),
        Some(ihist.to_str().unwrap()),
    );
    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/merge.merged.fq.gz")
    );
    assert_eq!(
        std::fs::read(&outu).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/merge.unmerged.fq.gz")
    );
    assert_eq!(
        std::fs::read(&ihist).unwrap(),
        std::fs::read("tests/bbtools/Lambda/golden/merge.ihist2.txt").unwrap()
    );
}

#[test]
fn command_fq_merge_novector_matches_bbtools_golden() {
    // `bbmerge.sh ... strict makevector=f` (classic filters, no net).
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("merged.fq");
    let outu = out_dir.path().join("unmerged.fq");
    lambda(
        &["--strict", "--no-make-vector"],
        out.to_str().unwrap(),
        Some(outu.to_str().unwrap()),
        None,
    );
    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/merge.novector.merged.fq.gz")
    );
    assert_eq!(
        std::fs::read(&outu).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/merge.novector.unmerged.fq.gz")
    );
}

#[test]
fn command_fq_merge_extend2_rem_matches_bbtools_golden() {
    // BBTools 40.01 `bbmerge-auto.sh ... strict k=81 extend2=80 rem` over
    // the extended reads (anchr merge phase 4): tadpole extension of unmerged
    // pairs + requireExtensionMatch. Merged/unmerged/ihist all byte-identical.
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("merged.fq");
    let outu = out_dir.path().join("unmerged.fq");
    let ihist = out_dir.path().join("ihist.txt");
    PgrCmd::new()
        .args(&[
            "fq",
            "merge",
            "tests/bbtools/Lambda/golden/ext_sub.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--outu",
            outu.to_str().unwrap(),
            "--ihist",
            ihist.to_str().unwrap(),
            "--strict",
            "--no-make-vector",
            "--extend2",
            "80",
            "--rem",
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/merge4.merged.fq.gz")
    );
    assert_eq!(
        std::fs::read(&outu).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/merge4.unmerged.fq.gz")
    );
    assert_eq!(
        std::fs::read(&ihist).unwrap(),
        std::fs::read("tests/bbtools/Lambda/golden/merge4.ihist.txt").unwrap()
    );
}

#[test]
fn command_fq_merge_requires_net_in_make_vector_mode() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "merge",
            "tests/bbtools/Lambda/R1.fq.gz",
            "tests/bbtools/Lambda/R2.fq.gz",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .failure();
}

#[test]
fn command_fq_merge_efilter_keeps_noisy_overlap_that_pfilter_would_reject() {
    // Two high-quality reads with a 40 bp overlap carrying two mismatches:
    // the expected-error filter trips before the probability filter, so the
    // merge is kept (classic mode, strict preset). Without the efilter step
    // the pair would be rejected by pfilter. Verified byte-identical to
    // BBTools 40.01 `BBMerge ... strict makevector=f`.
    let out_dir = tempfile::tempdir().unwrap();
    let in1 = out_dir.path().join("r1.fq");
    let in2 = out_dir.path().join("r2.fq");
    let out = out_dir.path().join("merged.fq");
    let outu = out_dir.path().join("unmerged.fq");
    let r1 = "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTAT";
    let r2 = "ATAGGTTTAGGCAACGGCGTAAGCGTCACTGTCGGAAGGG";
    let q1 = "I".repeat(r1.len());
    let q2 = "I".repeat(r2.len());
    std::fs::write(&in1, format!("@ef1/1\n{r1}\n+\n{q1}\n")).unwrap();
    std::fs::write(&in2, format!("@ef1/2\n{r2}\n+\n{q2}\n")).unwrap();
    PgrCmd::new()
        .args(&[
            "fq",
            "merge",
            in1.to_str().unwrap(),
            in2.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--outu",
            outu.to_str().unwrap(),
            "--strict",
            "--no-make-vector",
        ])
        .assert()
        .success();
    let merged = std::fs::read_to_string(&out).unwrap();
    assert_eq!(merged.lines().count(), 4);
    assert!(merged.starts_with(
        "@ef1/1\nAAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTNCGACAGTGACGCTTNCGCCGTTGCCTAAACCTAT\n"
    ));
    assert!(std::fs::read(&outu).unwrap().is_empty());
}

#[test]
fn command_fq_merge_joins_perfect_overlap() {
    let out_dir = tempfile::tempdir().unwrap();
    let in1 = out_dir.path().join("r1.fq");
    let in2 = out_dir.path().join("r2.fq");
    let out = out_dir.path().join("merged.fq");
    let outu = out_dir.path().join("unmerged.fq");
    // r2 = reverse complement of the last 64 bp of the 96 bp fragment, so
    // RC(r2) overlaps r1's tail by 32 bp (insert size 96).
    std::fs::write(
        &in1,
        concat!(
            "@r1/1\n",
            "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCT\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        ),
    )
    .unwrap();
    std::fs::write(
        &in2,
        concat!(
            "@r1/2\n",
            "GTTTAGGCAACGGCGAAAGCGTCACTGTCGCAAGGGTCGCCGCACATGTCGTTGCCTATATCCC\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        ),
    )
    .unwrap();
    PgrCmd::new()
        .args(&[
            "fq",
            "merge",
            in1.to_str().unwrap(),
            in2.to_str().unwrap(),
            "--no-make-vector",
            "-o",
            out.to_str().unwrap(),
            "--outu",
            outu.to_str().unwrap(),
        ])
        .assert()
        .success();
    let merged = std::fs::read_to_string(&out).unwrap();
    assert!(merged.starts_with("@r1/1\n"));
    assert_eq!(merged.lines().count(), 4);
    assert!(merged.lines().nth(1).unwrap().len() >= 64);
    // The unmerged file is empty: both reads merged.
    assert!(std::fs::read(&outu).unwrap().is_empty());
}

#[test]
fn command_fq_merge_keeps_non_overlapping_pairs_unmerged() {
    let out_dir = tempfile::tempdir().unwrap();
    let in1 = out_dir.path().join("r1.fq");
    let in2 = out_dir.path().join("r2.fq");
    let out = out_dir.path().join("merged.fq");
    let outu = out_dir.path().join("unmerged.fq");
    std::fs::write(
        &in1,
        concat!(
            "@r1/1\n",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        ),
    )
    .unwrap();
    std::fs::write(
        &in2,
        concat!(
            "@r1/2\n",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        ),
    )
    .unwrap();
    PgrCmd::new()
        .args(&[
            "fq",
            "merge",
            in1.to_str().unwrap(),
            in2.to_str().unwrap(),
            "--no-make-vector",
            "-o",
            out.to_str().unwrap(),
            "--outu",
            outu.to_str().unwrap(),
        ])
        .assert()
        .success();
    assert!(std::fs::read(&out).unwrap().is_empty());
    let unmerged = std::fs::read_to_string(&outu).unwrap();
    assert_eq!(unmerged.lines().count(), 8);
}
