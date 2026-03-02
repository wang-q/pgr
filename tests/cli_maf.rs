#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_maf_to_fas() {
    let (stdout, _) = PgrCmd::new()
        .args(&["maf", "to-fas", "tests/maf/example.maf"])
        .run();

    assert!(stdout.contains(">S288c.VIII(+):13377-13410"));
    assert!(stdout.contains("TTACTCGTCTTGCGGCCAAAACTCGAAGAAAAAC"));
    assert!(stdout.contains(">Spar.gi_29362578(-):72853-72885"));
    assert!(stdout.contains("TTACCCGTCTTGCGTCCAAAACTCGAA-AAAAAC"));
    assert_eq!(stdout.matches(">").count(), 8); // 2 blocks * 4 sequences
    assert_eq!(stdout.lines().count(), 18);
    assert!(stdout.contains("S288c.VIII"), "name list");
    assert!(stdout.contains(":42072-42168"), "coordinate transformed");
}
