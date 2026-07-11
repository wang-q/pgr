use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

/// Build the clap subcommand for cmp.
pub fn make_subcommand() -> Command {
    Command::new("cmp")
        .about("Compares trees (RF, WRF, KF distances)")
        .after_help(
            r###"
Compare trees using Robinson-Foulds (RF) distance and its variants.

Notes:
* Metrics:
    * RF: Robinson-Foulds distance (Topological difference).
    * WRF: Weighted Robinson-Foulds distance (Branch length difference).
    * KF: Kuhner-Felsenstein (Branch Score) distance.

* Input:
    * One file: Compares all trees in the file against each other (Pairwise).
    * Two files: Compares each tree in file1 against each tree in file2.

* Output:
    * TSV format: Tree1 \t Tree2 \t RF_Dist \t WRF_Dist \t KF_Dist

* IDs are 1-based indices of the trees in the input files.

Examples:
1. Compare all trees in a file:
   pgr nwk cmp trees.nwk

2. Compare trees between two files:
   pgr nwk cmp set1.nwk set2.nwk
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "First input filename (or stdin)",
        ))
        .arg(
            Arg::new("compare_file")
                .num_args(1)
                .index(2)
                .help("Second input filename (optional)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the cmp command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    // 1. Load first file
    let infile = args.get_one::<String>("infile").unwrap();
    let trees1 = Tree::from_file(infile)?;

    // 2. Load second file (if provided) or self-compare against trees1
    let compare_file = args.get_one::<String>("compare_file");
    let trees2_owned: Vec<Tree> = if let Some(f2) = compare_file {
        Tree::from_file(f2)?
    } else {
        Vec::new()
    };
    let trees2: &[Tree] = if compare_file.is_some() {
        &trees2_owned
    } else {
        &trees1
    };

    // 3. Output writer
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    // 4. Compare
    // Header
    writeln!(writer, "Tree1\tTree2\tRF_Dist\tWRF_Dist\tKF_Dist")?;

    // Single-file mode: skip self-comparisons (j == i) and duplicate pairs
    // (j < i) since RF is symmetric. Two-file mode: full cross comparison.
    for (i, t1) in trees1.iter().enumerate() {
        let start_j = if compare_file.is_some() { 0 } else { i + 1 };
        for (j, t2) in trees2.iter().enumerate().skip(start_j) {
            let (rf, wrf, kf) = pgr::libs::phylo::cmp::compute_tree_metrics(t1, t2)?;
            writeln!(writer, "{}\t{}\t{}\t{}\t{}", i + 1, j + 1, rf, wrf, kf)?;
        }
    }

    writer.flush()?;
    Ok(())
}
