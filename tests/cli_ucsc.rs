#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

fn get_path(filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/pgr");
    path.push(filename);
    path
}

fn normalize_chain_output(content: &str) -> String {
    content
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter(|line| !line.trim().is_empty())
        .collect::<Vec<&str>>()
        .join("\n")
}

// 1. Alignment - lavToPsl
#[test]
fn test_lav_to_psl_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("lastz.lav");
    let expected_output = get_path("lastz.psl");
    let output = temp.path().join("out.psl");

    PgrCmd::new()
        .args(&[
            "lav",
            "to-psl",
            input.to_str().unwrap(),
            "--outfile",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();

    let normalize = |s: &str| -> String {
        s.lines()
            .filter(|line| !line.starts_with('#')) // specific to psl header
            .collect::<Vec<&str>>()
            .join("\n")
    };

    assert_eq!(normalize(&output_content), normalize(&expected_content));
}

// 1. Alignment (Prep) - fa size / faToTwoBit
#[test]
fn test_2bit_size_pseudocat() {
    let temp = TempDir::new().unwrap();
    let input = get_path("pseudocat.2bit");
    let expected_output = get_path("pseudocat.sizes");
    let output = temp.path().join("out.sizes");

    PgrCmd::new()
        .args(&[
            "2bit",
            "size",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();

    assert_eq!(output_content, expected_content);
}

// 2. Chain - axtChain
#[test]
fn test_chaining_psl_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("lastz.psl");
    let t_2bit = get_path("pseudocat.2bit");
    let q_2bit = get_path("pseudopig.2bit");
    let expected_output = get_path("pslChain/lastz.raw.chain");
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "psl",
            "chain",
            t_2bit.to_str().unwrap(),
            q_2bit.to_str().unwrap(),
            input.to_str().unwrap(),
            "--outfile",
            output.to_str().unwrap(),
            "--score-scheme",
            "hoxd55",
            "--gap-model",
            "loose",
            "--min-score",
            "1000",
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(
        output_norm.lines().count(),
        expected_norm.lines().count(),
        "Line count mismatch"
    );

    assert_eq!(output_norm, expected_norm);
}

// 2. Chain - chainAntiRepeat
#[test]
fn test_chain_anti_repeat_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("pslChain/lastz.raw.chain");
    let t_2bit = get_path("pseudocat.2bit");
    let q_2bit = get_path("pseudopig.2bit");
    let expected_output = get_path("pslChain/lastz.chain");
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "chain",
            "anti-repeat",
            "--target-2bit",
            t_2bit.to_str().unwrap(),
            "--query-2bit",
            q_2bit.to_str().unwrap(),
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);
}

// 2. Chain - chainMergeSort
#[test]
fn test_chain_sort_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("pslChain/lastz.chain");
    let expected_output = get_path("all.chain");
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "chain",
            "sort",
            input.to_str().unwrap(),
            "--outfile",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);
}

// 2. Chain - chainPreNet
#[test]
fn test_chain_pre_net_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("all.chain");
    let t_sizes = get_path("pseudocat.sizes");
    let q_sizes = get_path("pseudopig.sizes");
    let expected_output = get_path("all.pre.chain");
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "chain",
            "pre-net",
            input.to_str().unwrap(),
            t_sizes.to_str().unwrap(),
            q_sizes.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);
}

// 2. Chain - chainNet
#[test]
fn test_chain_net_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("all.pre.chain");
    let t_sizes = get_path("pseudocat.sizes");
    let q_sizes = get_path("pseudopig.sizes");
    let expected_t_net = get_path("pseudocat.chainnet");
    let expected_q_net = get_path("pseudopig.chainnet");

    let output_t_net = temp.path().join("out.t.net");
    let output_q_net = temp.path().join("out.q.net");

    PgrCmd::new()
        .args(&[
            "chain",
            "net",
            input.to_str().unwrap(),
            t_sizes.to_str().unwrap(),
            q_sizes.to_str().unwrap(),
            output_t_net.to_str().unwrap(),
            output_q_net.to_str().unwrap(),
            "--min-space",
            "1",
            "--min-score",
            "2000",
        ])
        .run();

    let t_net_content = fs::read_to_string(&output_t_net).unwrap();
    let expected_t_content = fs::read_to_string(&expected_t_net).unwrap();
    assert_eq!(t_net_content, expected_t_content);

    let q_net_content = fs::read_to_string(&output_q_net).unwrap();
    let expected_q_content = fs::read_to_string(&expected_q_net).unwrap();
    assert_eq!(q_net_content, expected_q_content);
}

// 2. Chain - netSyntenic
#[test]
fn test_net_syntenic_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("pseudocat.chainnet");
    let expected_output = get_path("noClass.net");
    let output = temp.path().join("out.net");

    PgrCmd::new()
        .args(&[
            "net",
            "syntenic",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();
    assert_eq!(output_content, expected_content);
}

// 2. Chain - netChainSubset
#[test]
fn test_net_subset_lastz() {
    let temp = TempDir::new().unwrap();
    let net_input = get_path("noClass.net");
    let chain_input = get_path("all.chain");
    let expected_output = get_path("subset.chain");
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "net",
            "subset",
            net_input.to_str().unwrap(),
            chain_input.to_str().unwrap(),
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);
}

// 2. Chain - chainStitchId
#[test]
fn test_chain_stitch_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("subset.chain");
    let expected_output = get_path("over.chain");
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "chain",
            "stitch",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);
}

fn normalize_net_output(content: &str) -> String {
    content
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect::<Vec<&str>>()
        .join("\n")
}

#[test]
fn test_net_split_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("noClass.net");
    let output_dir = temp.path().join("net");

    PgrCmd::new()
        .args(&[
            "net",
            "split",
            input.to_str().unwrap(),
            "--outdir",
            output_dir.to_str().unwrap(),
        ])
        .run();

    let output_file = output_dir.join("cat.net");
    assert!(output_file.exists());

    let output_content = fs::read_to_string(&output_file).unwrap();
    let expected_output = get_path("net/cat.net");
    let expected_content = fs::read_to_string(&expected_output).unwrap();

    assert_eq!(
        normalize_net_output(&output_content),
        normalize_net_output(&expected_content)
    );
}

// 3. Axt - netToAxt
#[test]
fn test_net_to_axt_lastz() {
    let temp = TempDir::new().unwrap();
    let net_input = get_path("net/cat.net");
    let chain_input = get_path("all.pre.chain");
    let t_2bit = get_path("pseudocat.2bit");
    let q_2bit = get_path("pseudopig.2bit");
    let output = temp.path().join("cat.axt");

    PgrCmd::new()
        .args(&[
            "net",
            "to-axt",
            net_input.to_str().unwrap(),
            chain_input.to_str().unwrap(),
            t_2bit.to_str().unwrap(),
            q_2bit.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    assert!(output.exists());
    assert!(fs::metadata(&output).unwrap().len() > 0);

    // UCSC netToAxt emits AXT records in net-tree pre-order (parent segments
    // before child fills). The sorted variant is cat.axt (post-axtSort); the
    // raw netToAxt output is cat.tmp.axt.
    let expected_output = get_path("axtNet/cat.tmp.axt");
    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();

    // Normalize newlines
    let output_norm = output_content.replace("\r\n", "\n");
    let expected_norm = expected_content.replace("\r\n", "\n");

    assert_eq!(output_norm, expected_norm);
}

// 3. Axt - axtSort
#[test]
fn test_axt_sort_lastz() {
    let temp = TempDir::new().unwrap();

    // Use the expected output (already sorted) as input to verify idempotency/format handling
    // This avoids dependency on net-to-axt command in this test
    let input_path = get_path("axtNet/cat.axt");
    let output = temp.path().join("cat.axt");

    PgrCmd::new()
        .args(&[
            "axt",
            "sort",
            input_path.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&input_path).unwrap();

    // Normalize newlines for cross-platform comparison
    let output_norm = output_content.replace("\r\n", "\n");
    let expected_norm = expected_content.replace("\r\n", "\n");

    assert_eq!(output_norm, expected_norm);
}

// 3. Axt - axtToMaf
// UCSC `axtToMaf` crashes on macOS arm64 (Trace/BPT trap), so this is a
// pgr self-consistency test: verify the MAF output structure (record count,
// species count, scores) matches the AXT input rather than a UCSC reference.
#[test]
fn test_axt_to_maf_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("axtNet/cat.axt");
    let t_sizes = get_path("pseudocat.sizes");
    let q_sizes = get_path("pseudopig.sizes");
    let output = temp.path().join("cat.maf");

    PgrCmd::new()
        .args(&[
            "axt",
            "to-maf",
            input.to_str().unwrap(),
            "-t",
            t_sizes.to_str().unwrap(),
            "-q",
            q_sizes.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(get_path("axtNet/cat.maf")).unwrap();

    // Normalize newlines, then byte-compare. The expected fixture was
    // generated by pgr itself (UCSC axtToMaf crashes on this platform).
    let output_norm = output_content.replace("\r\n", "\n");
    let expected_norm = expected_content.replace("\r\n", "\n");
    assert_eq!(output_norm, expected_norm);
}

// 4. Chain - chainSplit (by target)
#[test]
fn test_chain_split_target_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("subset.chain");
    let output_dir = temp.path().join("split_target");

    PgrCmd::new()
        .args(&[
            "chain",
            "split",
            input.to_str().unwrap(),
            "--outdir",
            output_dir.to_str().unwrap(),
        ])
        .run();

    let output_file = output_dir.join("cat.chain");
    assert!(output_file.exists());

    let output_content = fs::read_to_string(&output_file).unwrap();
    let expected_content = fs::read_to_string(get_path("chainSplit_target/cat.chain")).unwrap();

    // UCSC chainSplit sorts chains by ID within each file; pgr preserves input
    // order. Compare chain records as a sorted set to ignore ordering.
    assert_eq!(
        normalize_chain_sorted(&output_content),
        normalize_chain_sorted(&expected_content)
    );
}

// 4. Chain - chainSplit (by query)
#[test]
fn test_chain_split_query_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("subset.chain");
    let output_dir = temp.path().join("split_query");

    PgrCmd::new()
        .args(&[
            "chain",
            "split",
            input.to_str().unwrap(),
            "--outdir",
            output_dir.to_str().unwrap(),
            "--by-query",
        ])
        .run();

    for name in &["pig1.chain", "pig2.chain"] {
        let output_file = output_dir.join(name);
        assert!(output_file.exists(), "missing {name}");
        let output_content = fs::read_to_string(&output_file).unwrap();
        let expected_content =
            fs::read_to_string(get_path(&format!("chainSplit_query/{name}"))).unwrap();
        assert_eq!(
            normalize_chain_sorted(&output_content),
            normalize_chain_sorted(&expected_content)
        );
    }
}

// 5. Net - netFilter --syn
#[test]
fn test_net_filter_syn_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("noClass.net");
    let output = temp.path().join("syn.net");

    PgrCmd::new()
        .args(&[
            "net",
            "filter",
            input.to_str().unwrap(),
            "--syn",
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    // Both UCSC netFilter -syn and pgr net filter --syn produce an empty file
    // for this dataset (top score 124204 < default minTopScore 300000).
    let output_content = fs::read_to_string(&output).unwrap();
    assert!(
        output_content.trim().is_empty(),
        "expected empty syn output, got: {output_content}"
    );
}

// 5. Net - netFilter -nonsyn
#[test]
fn test_net_filter_nonsyn_lastz() {
    let temp = TempDir::new().unwrap();
    let input = get_path("noClass.net");
    let output = temp.path().join("nonSyn.net");

    PgrCmd::new()
        .args(&[
            "net",
            "filter",
            input.to_str().unwrap(),
            "--nonsyn",
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(get_path("nonSyn.net")).unwrap();

    // pgr preserves ##matrix/##gapPenalties comments; UCSC netFilter strips
    // them. Strip comment lines before comparing.
    assert_eq!(
        normalize_net_output(&output_content),
        normalize_net_output(&expected_content)
    );
}

/// Parse chain file content into individual chain records, sort them, and
/// join. Used to compare chain files where chain ordering may differ (e.g.
/// UCSC chainSplit sorts by ID, pgr preserves input order).
fn normalize_chain_sorted(content: &str) -> String {
    let mut records: Vec<String> = Vec::new();
    let mut current: Vec<String> = Vec::new();
    for line in content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if line.starts_with("chain") {
            if !current.is_empty() {
                records.push(current.join("\n"));
            }
            current = vec![line.to_string()];
        } else {
            current.push(line.to_string());
        }
    }
    if !current.is_empty() {
        records.push(current.join("\n"));
    }
    records.sort();
    records.join("\n")
}
