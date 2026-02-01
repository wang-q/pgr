use assert_cmd::Command;
use predicates::prelude::*;

#[test]
fn test_lav_to_psl() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let input = r#"#:lav
s {
    "/path/target.fa" 1 1000
    "/path/query.fa" 1 500
}
h {
    ">target.fa"
    ">query.fa"
}
a {
    s 100
    l 1 1 10 10 95
}
"#;
    
    cmd.arg("lav")
        .arg("to-psl")
        .write_stdin(input)
        .assert()
        .success()
        .stdout(predicate::str::contains("10\t0\t0\t0\t0\t0\t0\t0\t+\tquery.fa\t500\t0\t10\ttarget.fa\t1000\t0\t10\t1\t10,\t0,\t0,"));
        
    Ok(())
}
