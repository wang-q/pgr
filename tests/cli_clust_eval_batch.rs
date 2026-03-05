use assert_cmd::Command;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_workflow_scan_eval_silhouette() -> anyhow::Result<()> {
    // 1. Prepare Tree
    // ((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4);
    // Heights: Leaves=0, (AB)=0.1, (CD)=0.1, Root=0.5
    let mut tree_file = NamedTempFile::new()?;
    write!(tree_file, "((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4);")?;

    // 2. Prepare Matrix
    // A-B: 0.2
    // C-D: 0.2
    // Others: 1.0
    let mut matrix_file = NamedTempFile::new()?;
    write!(
        matrix_file,
        "4
A 0.0 0.2 1.0 1.0
B 0.2 0.0 1.0 1.0
C 1.0 1.0 0.0 0.2
D 1.0 1.0 0.2 0.0
"
    )?;

    // 3. Run nwk cut --scan
    // Scan range: 0.05 (below 0.1), 0.2 (between 0.1 and 0.5), 0.6 (above 0.5)
    // Actually scan logic: start, end, step.
    // 0.0, 0.6, 0.2 -> 0.0, 0.2, 0.4, 0.6
    // 0.0: Cut at 0.0. (AB) height 0.1 > 0. Cut. -> {A},{B},{C},{D}.
    // 0.2: Cut at 0.2. (AB) height 0.1 <= 0.2. Keep. Root height 0.5 > 0.2. Cut. -> {A,B}, {C,D}.
    // 0.4: Cut at 0.4. Same as 0.2.
    // 0.6: Cut at 0.6. Root height 0.5 <= 0.6. Keep. -> {A,B,C,D}.

    let mut cmd_cut = Command::cargo_bin("pgr")?;
    let output_cut = cmd_cut
        .arg("nwk")
        .arg("cut")
        .arg(tree_file.path())
        .arg("--height")
        .arg("0.5") // dummy arg for group requirement, ignored by scan? No, method is required.
        // Actually scan overrides the value?
        // In my code:
        // let method = if matches.contains_id("height") { cut::Method::Height(val) } ...
        // So we need to provide --height <IGNORED_VAL> to satisfy ArgGroup, and --scan controls `val`.
        .arg("--scan")
        .arg("0.0,0.6,0.2")
        .output()?;

    assert!(output_cut.status.success());
    let stdout_cut = String::from_utf8(output_cut.stdout)?;

    // Verify Long Format
    assert!(stdout_cut.starts_with("Group\tClusterID\tSampleID"));
    // Check some content
    // Threshold 0: 4 clusters.
    // Threshold 0.2: 2 clusters.

    // 4. Run clust eval (Batch)
    // Pass stdout_cut to stdin? Or write to file.
    // Let's write to file to be safe.
    let mut partitions_file = NamedTempFile::new()?;
    write!(partitions_file, "{}", stdout_cut)?;

    let mut cmd_eval = Command::cargo_bin("pgr")?;
    let output_eval = cmd_eval
        .arg("clust")
        .arg("eval")
        .arg(partitions_file.path())
        .arg("--format")
        .arg("long")
        .arg("--matrix")
        .arg(matrix_file.path())
        .output()?;

    if !output_eval.status.success() {
        eprintln!("STDERR:\n{}", String::from_utf8_lossy(&output_eval.stderr));
    }
    assert!(output_eval.status.success());
    let stdout_eval = String::from_utf8(output_eval.stdout)?;

    // Output:
    // Group\tsilhouette
    // height=0\t0.000000
    // height=0.2\t0.800000
    // height=0.4\t0.800000
    // height=0.6\t0.000000 (or NaN/0 for single cluster? Sklearn says error if <2 clusters. My impl returns 0.0)

    let lines: Vec<&str> = stdout_eval.lines().collect();
    // Parse output
    // Header
    assert_eq!(lines[0], "Group\tsilhouette");

    // Rows
    for line in &lines[1..] {
        let parts: Vec<&str> = line.split('\t').collect();
        let group = parts[0];
        let sil: f64 = parts[1].parse()?;

        // Parse group "height=val"
        let val_str = group.strip_prefix("height=").unwrap();
        let val: f64 = val_str.parse()?;

        if (val - 0.0).abs() < 1e-6 {
            assert!((sil - 0.0).abs() < 1e-6);
        } else if (val - 0.2).abs() < 1e-6 {
            assert!((sil - 0.8).abs() < 1e-6);
        } else if (val - 0.4).abs() < 1e-6 {
            assert!((sil - 0.8).abs() < 1e-6);
        } else if (val - 0.6).abs() < 1e-6 {
            assert!((sil - 0.0).abs() < 1e-6);
        }
    }

    Ok(())
}
