#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_nwk_support() {
    // 1. Create target tree
    let mut target_file = NamedTempFile::new().unwrap();
    writeln!(target_file, "((A,B),(C,D));").unwrap();

    // 2. Create replicate trees
    let mut replicates_file = NamedTempFile::new().unwrap();
    writeln!(replicates_file, "((A,B),(C,D));").unwrap();
    writeln!(replicates_file, "((A,B),(C,D));").unwrap();
    writeln!(replicates_file, "((A,C),(B,D));").unwrap(); // different topology

    // 3. Run command (absolute counts)
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "support",
            target_file.path().to_str().unwrap(),
            replicates_file.path().to_str().unwrap(),
        ])
        .run();

    assert!(stdout.contains("((A,B)2,(C,D)2)3;"));

    // 4. Run command (percent)
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "support",
            target_file.path().to_str().unwrap(),
            replicates_file.path().to_str().unwrap(),
            "--percent",
        ])
        .run();

    // 2/3 * 100 = 66
    assert!(stdout.contains("((A,B)66,(C,D)66)100;"));
}
