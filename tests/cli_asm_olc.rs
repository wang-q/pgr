#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;

/// Parses PAF lines into (query, qlen, qs, qe, strand, target, tlen, ts, te,
/// matches, block, mapq, tags) tuples.
fn parse_paf(data: &str) -> Vec<Vec<String>> {
    data.lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').map(|f| f.to_string()).collect())
        .collect()
}

/// `pgr asm ovlp` reports the exact 10 bp suffix/prefix overlap.
#[test]
fn command_asm_ovlp_finds_dovetail() {
    let out_dir = tempfile::tempdir().unwrap();
    let a = out_dir.path().join("a.fa");
    let b = out_dir.path().join("b.fa");
    fs::write(
        &a,
        ">unitig_1\nTTTTTTTTTTACGTACGTAC\n>unitig_2\nACGTACGTACGGGGGGGGGG\n",
    )
    .unwrap();
    fs::write(&b, ">unitig_1\nGTACGTACGTCCCCCCCCCC\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "asm",
            "ovlp",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--overlap-k",
            "5",
            "--min-overlap",
            "8",
        ])
        .run();
    let paf = parse_paf(&stdout);
    assert!(!paf.is_empty(), "expected at least one overlap");
    // a:unitig_1 suffix 10..20 == b:unitig_1 prefix 0..10 (reverse strand).
    let rev = paf
        .iter()
        .find(|f| {
            f[0] == "a:unitig_1"
                && f[4] == "-"
                && f[5] == "b:unitig_1"
                && f[2] == "10"
                && f[3] == "20"
                && f[7] == "0"
                && f[8] == "10"
        })
        .expect("a:unitig_1 -> b:unitig_1 reverse dovetail");
    assert_eq!(rev[9], "10", "matches = overlap length");
    assert_eq!(rev[10], "10", "block length = overlap length");
    assert_eq!(rev[12], "ov:A:D", "dovetail tag");
}

/// Cross-file names are disambiguated with the file stem prefix.
#[test]
fn command_asm_ovlp_disambiguates_names() {
    let out_dir = tempfile::tempdir().unwrap();
    let a = out_dir.path().join("k21.fa");
    let b = out_dir.path().join("k51.fa");
    fs::write(
        &a,
        ">unitig_1\nTTTTTTTTTTACGTACGTAC\n>unitig_2\nACGTACGTACGGGGGGGGGG\n",
    )
    .unwrap();
    fs::write(&b, ">unitig_1\nGTACGTACGTCCCCCCCCCC\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "asm",
            "ovlp",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--overlap-k",
            "5",
            "--min-overlap",
            "8",
        ])
        .run();
    let paf = parse_paf(&stdout);
    let names: Vec<&str> = paf
        .iter()
        .flat_map(|f| [f[0].as_str(), f[6].as_str()])
        .collect();
    assert!(
        names.contains(&"k21:unitig_1"),
        "expected k21:unitig_1 in {names:?}"
    );
    assert!(
        names.contains(&"k51:unitig_1"),
        "expected k51:unitig_1 in {names:?}"
    );
    assert!(!names.contains(&"unitig_1"), "bare name must be prefixed");
}

/// Determinism: two runs produce identical output.
#[test]
fn command_asm_ovlp_deterministic() {
    let out_dir = tempfile::tempdir().unwrap();
    let a = out_dir.path().join("a.fa");
    fs::write(&a, ">u1\nTTTTTTTTTTACGTACGTAC\n>u2\nACGTACGTACGGGGGGGGGG\n").unwrap();
    let o1 = out_dir.path().join("o1.paf");
    let o2 = out_dir.path().join("o2.paf");
    for o in [&o1, &o2] {
        PgrCmd::new()
            .args(&[
                "asm",
                "ovlp",
                a.to_str().unwrap(),
                "-o",
                o.to_str().unwrap(),
                "--overlap-k",
                "5",
                "--min-overlap",
                "8",
            ])
            .assert()
            .success();
    }
    assert_eq!(fs::read(&o1).unwrap(), fs::read(&o2).unwrap());
}

/// `pgr asm layout` chains ovlp output into one layout with coordinates.
#[test]
fn command_asm_layout_chains_unitigs() {
    let out_dir = tempfile::tempdir().unwrap();
    let ut = out_dir.path().join("ut.fa");
    fs::write(
        &ut,
        ">unitig_1\nAAAAAAAAACGTACGT\n>unitig_2\nACGTACGTCCCCCCCC\n\
         >unitig_3\nCCCCCCCCGGGGGGGG\n",
    )
    .unwrap();
    let ovlp = out_dir.path().join("ovlp.paf");
    let layout = out_dir.path().join("layout.tsv");
    PgrCmd::new()
        .args(&[
            "asm",
            "ovlp",
            ut.to_str().unwrap(),
            "--overlap-k",
            "5",
            "--min-overlap",
            "8",
            "-o",
            ovlp.to_str().unwrap(),
        ])
        .assert()
        .success();
    PgrCmd::new()
        .args(&[
            "asm",
            "layout",
            ovlp.to_str().unwrap(),
            ut.to_str().unwrap(),
            "-o",
            layout.to_str().unwrap(),
        ])
        .assert()
        .success();
    let lines: Vec<String> = fs::read_to_string(&layout)
        .unwrap()
        .lines()
        .map(|l| l.to_string())
        .collect();
    assert_eq!(
        lines,
        vec![
            "contig_1\t0\tut:unitig_1\t+\t0\t16\t0",
            "contig_1\t1\tut:unitig_2\t+\t8\t24\t8",
            "contig_1\t2\tut:unitig_3\t+\t16\t32\t8",
        ]
    );
}

/// Deterministic xorshift PRNG for synthetic test genomes/reads.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

fn rev_comp(s: &[u8]) -> Vec<u8> {
    pgr::libs::nt::rev_comp(s).collect()
}

fn parse_fa_records(data: &[u8]) -> Vec<Vec<u8>> {
    let mut recs = Vec::new();
    let mut cur: Option<Vec<u8>> = None;
    for line in std::str::from_utf8(data).unwrap().lines() {
        if line.starts_with('>') {
            if let Some(c) = cur.take() {
                recs.push(c);
            }
            cur = Some(Vec::new());
        } else if let Some(c) = cur.as_mut() {
            c.extend_from_slice(line.as_bytes());
        }
    }
    if let Some(c) = cur {
        recs.push(c);
    }
    recs
}

/// The full OLC pipeline reconstructs an error-free synthetic genome.
#[test]
fn command_asm_olc_reconstructs_synthetic_genome() {
    let out_dir = tempfile::tempdir().unwrap();
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let bases = *b"ACGT";
    let genome: Vec<u8> = (0..2400).map(|_| bases[rng.below(4)]).collect();
    let read_len = 250usize;
    let mut fa = String::new();
    for i in 0..480 {
        let start = rng.below(genome.len() - read_len + 1);
        let mut r = genome[start..start + read_len].to_vec();
        if rng.below(2) == 1 {
            r = rev_comp(&r);
        }
        fa.push_str(&format!(">r{i}\n{}\n", String::from_utf8(r).unwrap()));
    }
    let reads = out_dir.path().join("reads.fa");
    fs::write(&reads, fa).unwrap();
    let contigs_path = out_dir.path().join("contigs.fa");
    PgrCmd::new()
        .args(&[
            "asm",
            "olc",
            reads.to_str().unwrap(),
            "-o",
            contigs_path.to_str().unwrap(),
            "--kmer",
            "21,51,81",
            "--min-contig-len",
            "100",
        ])
        .assert()
        .success();
    let contigs = parse_fa_records(&fs::read(&contigs_path).unwrap());
    assert!(!contigs.is_empty(), "expected at least one contig");
    let longest = contigs.iter().max_by_key(|c| c.len()).unwrap();
    let rcg = rev_comp(&genome);
    let in_genome = genome
        .windows(longest.len())
        .any(|w| w == longest.as_slice());
    let in_rc = rcg.windows(longest.len()).any(|w| w == longest.as_slice());
    assert!(
        in_genome || in_rc,
        "longest contig is not an exact genome substring"
    );
    assert!(
        longest.len() >= genome.len() * 85 / 100,
        "longest contig covers only {} / {} bp",
        longest.len(),
        genome.len()
    );
}

/// `pgr asm cns` stitches a layout TSV into a FASTA contig.
#[test]
fn command_asm_cns_stitches_layout() {
    let out_dir = tempfile::tempdir().unwrap();
    let ut = out_dir.path().join("ut.fa");
    let layout = out_dir.path().join("layout.tsv");
    let out = out_dir.path().join("contigs.fa");
    fs::write(
        &ut,
        ">unitig_1\nAAAAAAAAACGTACGT\n>unitig_2\nACGTACGTCCCCCCCC\n>unitig_3\nCCCCCCCCGGGGGGGG\n",
    )
    .unwrap();
    fs::write(
        &layout,
        "contig_1\t0\tut:unitig_1\t+\t0\t16\t0\n\
         contig_1\t1\tut:unitig_2\t+\t8\t24\t8\n\
         contig_1\t2\tut:unitig_3\t+\t16\t32\t8\n",
    )
    .unwrap();
    PgrCmd::new()
        .args(&[
            "asm",
            "cns",
            layout.to_str().unwrap(),
            ut.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    let data = fs::read_to_string(&out).unwrap();
    assert!(data.starts_with(">contig_1,len=32,cov=1.5\n"));
    assert!(data.contains("AAAAAAAAACGTACGTCCCCCCCCGGGGGGGG"));
}

/// A layout referencing `contig_0` (id 0) must fail cleanly instead of
/// underflowing `ci - 1` (zero-panic policy).
#[test]
fn command_asm_cns_rejects_contig_zero() {
    let out_dir = tempfile::tempdir().unwrap();
    let ut = out_dir.path().join("ut.fa");
    let layout = out_dir.path().join("layout.tsv");
    let out = out_dir.path().join("contigs.fa");
    fs::write(&ut, ">unitig_1\nAAAAAAAAACGTACGT\n").unwrap();
    fs::write(&layout, "contig_0\t0\tut:unitig_1\t+\t0\t16\t0\n").unwrap();
    PgrCmd::new()
        .args(&[
            "asm",
            "cns",
            layout.to_str().unwrap(),
            ut.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .failure();
}

/// `--keep-dir` writes the intermediate stage files.
#[test]
fn command_asm_olc_keep_dir_writes_stages() {
    let out_dir = tempfile::tempdir().unwrap();
    let keep = out_dir.path().join("stage");
    let mut rng = Rng(0x1234_5678_9ABC_DEF0);
    let bases = *b"ACGT";
    let genome: Vec<u8> = (0..400).map(|_| bases[rng.below(4)]).collect();
    let read_len = 200usize;
    let mut fa = String::new();
    for i in 0..80 {
        let start = rng.below(genome.len() - read_len + 1);
        let mut r = genome[start..start + read_len].to_vec();
        if rng.below(2) == 1 {
            r = rev_comp(&r);
        }
        fa.push_str(&format!(">r{i}\n{}\n", String::from_utf8(r).unwrap()));
    }
    let reads = out_dir.path().join("reads.fa");
    fs::write(&reads, fa).unwrap();
    PgrCmd::new()
        .args(&[
            "asm",
            "olc",
            reads.to_str().unwrap(),
            "-o",
            out_dir.path().join("c.fa").to_str().unwrap(),
            "--kmer",
            "21,31",
            "--min-contig-len",
            "100",
            "--keep-dir",
            keep.to_str().unwrap(),
        ])
        .assert()
        .success();
    for f in ["unitigs.fa", "ovlp.paf", "layout.tsv"] {
        assert!(keep.join(f).is_file(), "missing stage file {f}");
    }
}

/// Reads too shallow for any solid k-mer yield a friendly error.
#[test]
fn command_asm_olc_empty_unitigs_friendly_error() {
    let out_dir = tempfile::tempdir().unwrap();
    let reads = out_dir.path().join("reads.fa");
    fs::write(
        &reads,
        ">r1\nACGTACGTACGTACGTACGT\n>r2\nACGTACGTACGTACGTACGT\n",
    )
    .unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&["asm", "olc", reads.to_str().unwrap(), "--kmer", "21,51"])
        .run_fail();
    assert!(
        stderr.contains("cannot overlap"),
        "expected friendly error, got: {stderr}"
    );
}

/// The driver (with unitig-level contain pre-filtering) and the explicit
/// stage pipeline (unitig -> ovlp -> layout -> cns, unfiltered) produce the
/// same contig sequences on unambiguous graphs (regression guard; on real
/// data the filtered driver additionally merges redundant paths).
#[test]
fn command_asm_olc_matches_stage_pipeline() {
    let out_dir = tempfile::tempdir().unwrap();
    let mut rng = Rng(0xABCD_EF01_2345_6789);
    let bases = *b"ACGT";
    let genome: Vec<u8> = (0..2400).map(|_| bases[rng.below(4)]).collect();
    let read_len = 250usize;
    let mut fa = String::new();
    for i in 0..480 {
        let start = rng.below(genome.len() - read_len + 1);
        let mut r = genome[start..start + read_len].to_vec();
        if rng.below(2) == 1 {
            r = rev_comp(&r);
        }
        fa.push_str(&format!(">r{i}\n{}\n", String::from_utf8(r).unwrap()));
    }
    let reads = out_dir.path().join("reads.fa");
    fs::write(&reads, fa).unwrap();

    // Driver.
    let c1 = out_dir.path().join("c1.fa");
    PgrCmd::new()
        .args(&[
            "asm",
            "olc",
            reads.to_str().unwrap(),
            "-o",
            c1.to_str().unwrap(),
            "--kmer",
            "21,51,81",
            "--min-contig-len",
            "100",
        ])
        .assert()
        .success();
    let mut seqs1 = parse_fa_records(&fs::read(&c1).unwrap());

    // Stage pipeline.
    let mut unitigs = Vec::new();
    for k in ["21", "51", "81"] {
        let u = out_dir.path().join(format!("u{k}.fa"));
        PgrCmd::new()
            .args(&[
                "asm",
                "unitig",
                reads.to_str().unwrap(),
                "-o",
                u.to_str().unwrap(),
                "--kmer",
                k,
                "--min-contig-len",
                "100",
            ])
            .assert()
            .success();
        unitigs.push(u);
    }
    let ovlp = out_dir.path().join("ovlp.paf");
    let layout = out_dir.path().join("layout.tsv");
    let mut args = vec!["asm", "ovlp"];
    for u in &unitigs {
        args.push(u.to_str().unwrap());
    }
    args.extend(["-o", ovlp.to_str().unwrap()]);
    PgrCmd::new().args(&args).assert().success();
    let mut args = vec!["asm", "layout", ovlp.to_str().unwrap()];
    for u in &unitigs {
        args.push(u.to_str().unwrap());
    }
    args.extend(["-o", layout.to_str().unwrap()]);
    PgrCmd::new().args(&args).assert().success();
    let c2 = out_dir.path().join("c2.fa");
    let mut args = vec!["asm", "cns", layout.to_str().unwrap()];
    for u in &unitigs {
        args.push(u.to_str().unwrap());
    }
    args.extend(["-o", c2.to_str().unwrap(), "--min-contig-len", "100"]);
    PgrCmd::new().args(&args).assert().success();
    let mut seqs2 = parse_fa_records(&fs::read(&c2).unwrap());

    seqs1.sort();
    seqs2.sort();
    assert_eq!(seqs1, seqs2, "driver and stage pipeline contigs differ");
}
