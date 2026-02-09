use clap::*;
use pgr::libs::phylo::tree::Tree;
use pgr::libs::phylo::TreeComparison;
use std::collections::{BTreeMap, HashSet};
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
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename"),
        )
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
    // Check leaves
    let leaves1: HashSet<_> = t1.get_leaf_names().into_iter().flatten().collect();
    let leaves2: HashSet<_> = t2.get_leaf_names().into_iter().flatten().collect();

    if leaves1 != leaves2 {
        return (
            "Error".to_string(),
            "Error".to_string(),
            "Error".to_string(),
        );
    }

    // Build map
    let mut all_leaves: Vec<_> = leaves1.into_iter().collect();
    all_leaves.sort(); // Deterministic order
    let mut leaf_map = BTreeMap::new();
    for (i, name) in all_leaves.iter().enumerate() {
        leaf_map.insert(name.clone(), i);
    }

    // Get splits
    let s1 = t1.get_splits_with_values(&leaf_map);
    let s2 = t2.get_splits_with_values(&leaf_map);

    // RF: symmetric difference count
    let keys1: HashSet<_> = s1.keys().collect();
    let keys2: HashSet<_> = s2.keys().collect();
    let rf = keys1.symmetric_difference(&keys2).count();

    // WRF & KF
    let all_keys: HashSet<_> = s1.keys().chain(s2.keys()).collect();
    let mut wrf = 0.0;
    let mut kf_sq = 0.0;

    for key in all_keys {
        let v1 = s1.get(key).copied().unwrap_or(0.0);
        let v2 = s2.get(key).copied().unwrap_or(0.0);
        let diff = v1 - v2;
        wrf += diff.abs();
        kf_sq += diff.powi(2);
    }

    let format_float = |v: f64| -> String {
        let s = format!("{:.6}", v);
        let trimmed = s.trim_end_matches('0').trim_end_matches('.');
        if trimmed.is_empty() {
            "0".to_string()
        } else {
            trimmed.to_string()
        }
    };

    (
        rf.to_string(),
        format_float(wrf),
        format_float(kf_sq.sqrt()),
    )
}
