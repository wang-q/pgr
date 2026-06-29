use clap::*;
use pgr::libs::phylo::tree::Tree;
use pgr::libs::phylo::TreeComparison;
use std::io::{self, Write};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("cmp")
        .about("Compare trees (RF, WRF, KF distances)")
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
    * TSV format: Tree1_ID \t Tree2_ID \t RF_Dist \t WRF_Dist \t KF_Dist

* IDs are 1-based indices of the trees in the input files.

Examples:
1. Compare all trees in a file:
   pgr nwk cmp trees.nwk

2. Compare trees between two files:
   pgr nwk cmp set1.nwk set2.nwk
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("First input filename (or stdin)"),
        )
        .arg(
            Arg::new("compare_file")
                .num_args(1)
                .index(2)
                .help("Second input filename (optional)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    // 1. Load first file
    let infile = args.get_one::<String>("infile").unwrap();
    let trees1 = Tree::from_file(infile)?;

    // 2. Load second file (if provided) or point to first
    let trees2 = if let Some(f2) = args.get_one::<String>("compare_file") {
        Tree::from_file(f2)?
    } else {
        Vec::new()
    };

    let self_compare = args.get_one::<String>("compare_file").is_none();

    // 3. Output writer
    let outfile = args.get_one::<String>("outfile").unwrap();
    let mut writer: Box<dyn Write> = if outfile == "stdout" {
        Box::new(io::stdout())
    } else {
        Box::new(std::fs::File::create(outfile)?)
    };

    // 4. Compare
    // Header
    writeln!(writer, "Tree1\tTree2\tRF_Dist\tWRF_Dist\tKF_Dist")?;

    if self_compare {
        for (i, t1) in trees1.iter().enumerate() {
            for (j, t2) in trees1.iter().enumerate() {
                let (rf, wrf, kf) = compute_metrics(t1, t2);
                writeln!(writer, "{}\t{}\t{}\t{}\t{}", i + 1, j + 1, rf, wrf, kf)?;
            }
        }
    } else {
        for (i, t1) in trees1.iter().enumerate() {
            for (j, t2) in trees2.iter().enumerate() {
                let (rf, wrf, kf) = compute_metrics(t1, t2);
                writeln!(writer, "{}\t{}\t{}\t{}\t{}", i + 1, j + 1, rf, wrf, kf)?;
            }
        }
    }

    Ok(())
}

fn compute_metrics(t1: &Tree, t2: &Tree) -> (String, String, String) {
    let rf = t1.robinson_foulds(t2);
    let wrf = t1.weighted_robinson_foulds(t2);
    let kf = t1.kuhner_felsenstein(t2);

    let format_float = |v: f64| -> String {
        let s = format!("{:.6}", v);
        let trimmed = s.trim_end_matches('0').trim_end_matches('.');
        if trimmed.is_empty() {
            "0".to_string()
        } else {
            trimmed.to_string()
        }
    };

    match (rf, wrf, kf) {
        (Ok(rf), Ok(wrf), Ok(kf)) => (rf.to_string(), format_float(wrf), format_float(kf)),
        _ => (
            "Error".to_string(),
            "Error".to_string(),
            "Error".to_string(),
        ),
    }
}
