use assert_cmd::Command;

#[test]
fn command_rename_basic() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("rename")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-n")
        .arg("Homo")
        .arg("-r")
        .arg("Human")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Human"));
    assert!(!stdout.contains("Homo"));
    assert!(stdout.contains("Pan")); // Others preserved

    Ok(())
}

#[test]
fn command_rename_lca() -> anyhow::Result<()> {
    // In catarrhini.nwk, Homo and Pan are children of Hominini.
    // Rename Hominini (LCA of Homo,Pan) to CladeX
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("rename")
        .arg("tests/newick/catarrhini.nwk")
        .arg("--lca")
        .arg("Homo,Pan")
        .arg("-r")
        .arg("CladeX")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("CladeX"));
    assert!(!stdout.contains("Hominini"));

    Ok(())
}

#[test]
fn command_rename_mixed() -> anyhow::Result<()> {
    // ((A,B),C);
    // Rename A -> A1.
    // Rename LCA(A,B) -> AB.

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("rename")
        .arg("stdin")
        .arg("-n")
        .arg("A")
        .arg("-r")
        .arg("A1")
        .arg("-l")
        .arg("A,B")
        .arg("-r")
        .arg("AB")
        .write_stdin("((A,B),C);")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((A1,B)AB,C);"));

    Ok(())
}

#[test]
fn command_replace() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc.replace.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("((Homo,Pan),Gorilla);"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc.replace.tsv")
        .arg("--mode")
        .arg("species")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((A[&&NHX:S=Homo],B[&&NHX:S=Pan]),C[&&NHX:S=Gorilla]);"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc3.replace.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(
        stdout.contains("((Homo[&&NHX:color=red],Pan[&&NHX:color=red]),Gorilla[&&NHX:color=red]);")
    );

    Ok(())
}

#[test]
fn command_replace_comments() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/mixed_comments.replace.tsv")
        .arg("--mode")
        .arg("species")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // A -> Homo [UnstructuredComment]
    // B -> Pan [&&NHX:key=value]
    // C -> Gorilla [tag]
    assert!(stdout.contains("A[&&NHX:S=Homo:UnstructuredComment]"));
    assert!(stdout.contains("B[&&NHX:S=Pan:key=value]"));
    assert!(stdout.contains("C[&&NHX:S=Gorilla:tag]"));

    Ok(())
}

#[test]
fn command_replace_remove() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/replace_remove.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // A -> "" (removed)
    // B -> Pan
    // C -> C (untouched)
    // Original: ((A,B),C);
    // Expected: ((,Pan),C);
    // Note: Newick parser/writer should handle empty names as empty strings,
    // resulting in format like `(:0.1,Pan:0.1)...` or just `(,Pan)...` depending on branch lengths.
    // abc.nwk has no branch lengths.
    assert!(stdout.contains("((,Pan),C);"));

    Ok(())
}

#[test]
fn command_replace_filter() -> anyhow::Result<()> {
    // abc.nwk: ((A,B),C);
    // All are leaves.

    // 1. Skip leaves (should change nothing if all matches are leaves)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc.replace.tsv")
        .arg("--Leaf") // Skip leaves
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // A, B, C are leaves, so they should be skipped. Output remains original (or similar).
    assert!(stdout.contains("((A,B),C);"));

    // 2. Skip internal (should change leaves)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc.replace.tsv")
        .arg("--Internal") // Skip internal
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // A, B, C are leaves, so they should be replaced.
    assert!(stdout.contains("((Homo,Pan),Gorilla);"));

    Ok(())
}

#[test]
fn command_replace_multi() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/forest.nwk")
        .arg("tests/newick/forest.replace.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Check line count (should be 5 trees)
    assert_eq!(stdout.lines().count(), 5);

    // Check replacements in the last tree (line 5)
    // Original: (Homo,(Pan,...
    // Expected: (Human,(Chimp,...
    assert!(stdout.contains("(Human,(Chimp,"));

    Ok(())
}

#[test]
fn command_topo_basic() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/catarrhini.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // catarrhini.nwk has lengths and comments.
    // Default topo removes lengths and comments (properties), keeps labels.
    // Original: ((((Gorilla:16,(Pan:10,Homo:10)Hominini:10)Homininae:15,Pongo:30)Hominidae:15,Hylobates:20):10,(((Macaca:10,Papio:10):20,Cercopithecus:10)Cercopithecinae:25,(Simias:10,Colobus:7)Colobinae:5)Cercopithecidae:10);
    // Expected: ((((Gorilla,(Pan,Homo)Hominini)Homininae,Pongo)Hominidae,Hylobates),(((Macaca,Papio),Cercopithecus)Cercopithecinae,(Simias,Colobus)Colobinae)Cercopithecidae);
    // Note: The root edge length is also removed.

    assert!(stdout.contains("((((Gorilla,(Pan,Homo)Hominini)Homininae,Pongo)Hominidae,Hylobates)"));
    assert!(!stdout.contains(":")); // No lengths

    Ok(())
}

#[test]
fn command_topo_remove_labels() -> anyhow::Result<()> {
    // Test with -I (remove internal labels) and -L (remove leaf labels)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-I")
        .arg("-L")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Should have no labels
    assert!(stdout.contains("((((,(,))"));
    assert!(!stdout.contains("Homo"));
    assert!(!stdout.contains("Hominini"));

    Ok(())
}

#[test]
fn command_topo_keep_bl() -> anyhow::Result<()> {
    // Test --bl (keep branch lengths)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-I")
        .arg("-L")
        .arg("--bl")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains(":16")); // Check for specific length
    assert!(!stdout.contains("Gorilla"));

    Ok(())
}

#[test]
fn command_topo_compat_simple() -> anyhow::Result<()> {
    // simple:newtree.nw
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/newtree.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    let expected = "(FMDV-C,((((((((HRV16,HRV1B)52,(HRV24,HRV85)70)22,(HRV11,(HRV9,(HRV64,HRV94)32)54)1)17,(HRV39,HRV2)92)97,HRV89)62,(HRV78,HRV12)52)100,((((HRV37,HRV3)65,HRV14)89,(HRV52,HRV17)100)75,(HRV93,HRV27)99)83)48,((((POLIO3,((POLIO2,(POLIO1A,COXA18)22)38,COXA17)72)97,COXA1)76,(((ECHO1,COXB2)83,ECHO6)99,(HEV70,HEV68)99)70)64,(COXA14,(COXA6,COXA2))59)100)68);";

    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_topo_compat_multiple() -> anyhow::Result<()> {
    // multiple:forest.nw
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/forest.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    let expected = r#"(Pandion,((Buteo,Aquila,Haliaeetus),(Milvus,Elanus)),Sagittarius,((Micrastur,Falco),(Polyborus,Milvagus)));
((Diomedea,Daption),(Fregata,Phalacrocorax,Sula),(Larus,(Fratercula,Uria)));
(((Ticodendraceae,Betulaceae),Casuarinaceae),(Rhoipteleaceae,Juglandaceae),Myricaceae);
((((Gorilla,(Pan,Homo)Hominini)Homininae,Pongo)Hominidae,Hylobates),(((Macaca,Papio),Cercopithecus)Cercopithecinae,(Simias,Colobus)Colobinae)Cercopithecidae);
(Homo,(Pan,(Gorilla,(Pongo,(Hylobates,(((Cercopithecus,(Macaca,Papio)),Simias),Cebus))))));"#;

    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_topo_compat_rootedge() -> anyhow::Result<()> {
    // rootedge: edged_root.nw
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/edged_root.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    let expected = "(((Btar,Cfam),Hsap),(Mmus,Rnov));";

    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_topo_compat_i() -> anyhow::Result<()> {
    // I:-I newtree.nw
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/newtree.nwk")
        .arg("-I")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    let expected = "(FMDV-C,((((((((HRV16,HRV1B),(HRV24,HRV85)),(HRV11,(HRV9,(HRV64,HRV94)))),(HRV39,HRV2)),HRV89),(HRV78,HRV12)),((((HRV37,HRV3),HRV14),(HRV52,HRV17)),(HRV93,HRV27))),((((POLIO3,((POLIO2,(POLIO1A,COXA18)),COXA17)),COXA1),(((ECHO1,COXB2),ECHO6),(HEV70,HEV68))),(COXA14,(COXA6,COXA2)))));";

    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_topo_compat_l() -> anyhow::Result<()> {
    // L:-L newtree.nw
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/newtree.nwk")
        .arg("-L")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Note: pgr output might differ in empty labels if not handled exactly like newick_utils
    // newick_utils outputs: (,((((((((,)52,(,)70)22,(,(,(,)32)54)1)17,(,)92)97,)62,(,)52)100...
    // It seems it replaces leaf labels with empty string but keeps the node?
    let expected = "(,((((((((,)52,(,)70)22,(,(,(,)32)54)1)17,(,)92)97,)62,(,)52)100,((((,)65,)89,(,)100)75,(,)99)83)48,((((,((,(,)22)38,)72)97,)76,(((,)83,)99,(,)99)70)64,(,(,))59)100)68);";

    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_topo_compat_li() -> anyhow::Result<()> {
    // LI:-LI newtree.nw
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/newtree.nwk")
        .arg("-L")
        .arg("-I")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    let expected = "(,((((((((,),(,)),(,(,(,)))),(,)),),(,)),((((,),),(,)),(,))),((((,((,(,)),)),),(((,),),(,))),(,(,)))));";

    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_topo_compat_bi() -> anyhow::Result<()> {
    // bI:-bI newtree.nw
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/newtree.nwk")
        .arg("-b")
        .arg("-I")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // NOTE: Floating point precision might differ.
    // newick_utils: 2.0799315
    // pgr: ?
    let expected = "(FMDV-C:2.0799315,((((((((HRV16:0.071498,HRV1B:0.082284):0.04546,(HRV24:0.040859,HRV85:0.040089):0.034432):0.023874,(HRV11:0.040805,(HRV9:0.045986,(HRV64:0.048368,HRV94:0.084787):0.018131):0.092702):0.004912):0.018847,(HRV39:0.070769,HRV2:0.039029):0.056213):0.152625,HRV89:0.141183):0.072809,(HRV78:0.230063,HRV12:0.187536):0.069229):0.522696,((((HRV37:0.056416,HRV3:0.111802):0.026307,HRV14:0.031521):0.066208,(HRV52:0.013318,HRV17:0.017873):0.106471):0.052682,(HRV93:0.038271,HRV27:0.0026):0.150076):0.082254):0.091013,((((POLIO3:0,((POLIO2:0,(POLIO1A:0,COXA18:0):0):0,COXA17:0.005726):0.005697):0.051384,COXA1:0.104463):0.058199,(((ECHO1:0,COXB2:0.011614):0.012107,ECHO6:0.005466):0.130995,(HEV70:0.031767,HEV68:0.086627):0.10259):0.062266):0.050449,(COXA14:0.036101,(COXA6:0.011953,COXA2:0.005806):0.016157):0.323718):0.060172):2.0799315);";

    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_topo_compat_bl() -> anyhow::Result<()> {
    // bL:-bL newtree.nw
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/newtree.nwk")
        .arg("-b")
        .arg("-L")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    let expected = "(:2.0799315,((((((((:0.071498,:0.082284)52:0.04546,(:0.040859,:0.040089)70:0.034432)22:0.023874,(:0.040805,(:0.045986,(:0.048368,:0.084787)32:0.018131)54:0.092702)1:0.004912)17:0.018847,(:0.070769,:0.039029)92:0.056213)97:0.152625,:0.141183)62:0.072809,(:0.230063,:0.187536)52:0.069229)100:0.522696,((((:0.056416,:0.111802)65:0.026307,:0.031521)89:0.066208,(:0.013318,:0.017873)100:0.106471)75:0.052682,(:0.038271,:0.0026)99:0.150076)83:0.082254)48:0.091013,((((:0,((:0,(:0,:0)22:0)38:0,:0.005726)72:0.005697)97:0.051384,:0.104463)76:0.058199,(((:0,:0.011614)83:0.012107,:0.005466)99:0.130995,(:0.031767,:0.086627)99:0.10259)70:0.062266)64:0.050449,(:0.036101,(:0.011953,:0.005806):0.016157)59:0.323718)100:0.060172)68:2.0799315);";

    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_topo_compat_bil() -> anyhow::Result<()> {
    // bIL:-bIL newtree.nw
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/newtree.nwk")
        .arg("-b")
        .arg("-I")
        .arg("-L")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    let expected = "(:2.0799315,((((((((:0.071498,:0.082284):0.04546,(:0.040859,:0.040089):0.034432):0.023874,(:0.040805,(:0.045986,(:0.048368,:0.084787):0.018131):0.092702):0.004912):0.018847,(:0.070769,:0.039029):0.056213):0.152625,:0.141183):0.072809,(:0.230063,:0.187536):0.069229):0.522696,((((:0.056416,:0.111802):0.026307,:0.031521):0.066208,(:0.013318,:0.017873):0.106471):0.052682,(:0.038271,:0.0026):0.150076):0.082254):0.091013,((((:0,((:0,(:0,:0):0):0,:0.005726):0.005697):0.051384,:0.104463):0.058199,(((:0,:0.011614):0.012107,:0.005466):0.130995,(:0.031767,:0.086627):0.10259):0.062266):0.050449,(:0.036101,(:0.011953,:0.005806):0.016157):0.323718):0.060172):2.0799315);";

    assert_eq!(stdout.trim(), expected);
    Ok(())
}


#[test]
fn command_order() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abc.nwk")
        .arg("--nd")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("(C,(A,B));"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abc.nwk")
        .arg("--ndr")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((A,B),C);"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abc.nwk")
        .arg("--an")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((A,B),C);"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abc.nwk")
        .arg("--anr")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("(C,(B,A));"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abc.nwk")
        .arg("--anr")
        .arg("--ndr")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((B,A),C);"));

    Ok(())
}

#[test]
fn command_order_list() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abcde.nwk")
        .arg("--list")
        .arg("tests/newick/abcde.list")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("(C:1,(B:1,A:1)D:1)E;"));

    Ok(())
}

#[test]
fn command_order_species() -> anyhow::Result<()> {
    // Create a temporary directory for testing
    let tempdir = tempfile::tempdir()?;
    let temp_path = tempdir.path();

    std::fs::copy("tests/newick/species.nwk", temp_path.join("species.nwk"))?;

    // Generate a list of labels from the tree
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("nwk")
        .arg("label")
        .arg("species.nwk")
        .arg("-o")
        .arg("species.list")
        .current_dir(temp_path)
        .output()?;

    // Order the tree using the generated list
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("species.nwk")
        .arg("--list")
        .arg("species.list")
        .current_dir(temp_path)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Compare the ordered tree with the original one
    // They should be identical as the list was generated from the original order
    let original = std::fs::read_to_string("tests/newick/species.nwk")?;
    assert_eq!(stdout.trim(), original.trim());

    // gene tree
    std::fs::copy("tests/newick/pmxc.nwk", temp_path.join("pmxc.nwk"))?;

    // Order pmxc.nwk using the generated list
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("pmxc.nwk")
        .arg("--list")
        .arg("species.list")
        .current_dir(temp_path)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Read the original pmxc.nwk file
    let original = std::fs::read_to_string("tests/newick/pmxc.nwk")?;

    // The ordered tree should be different from the original one
    assert_ne!(stdout.trim(), original.trim());

    Ok(())
}
