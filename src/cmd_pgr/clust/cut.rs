use anyhow::Context;
use clap::{value_parser, Arg, ArgGroup, ArgMatches, Command};
use pgr::libs::clust::tree_cut::{self as cut, RepMode, METHOD_NAMES};
use pgr::libs::phylo::tree::Tree;
use std::io::Write;
/// Build the clap subcommand for cut.
pub fn make_subcommand() -> Command {
    Command::new("cut")
        .about("Cuts a tree into clusters")
        .after_help(
            r###"
Cuts the tree into clusters based on various criteria.

Criteria:
* `--k <K>`: Cut into K clusters (top-down split by height).
* `--height <H>`: Cut at specific height (max distance to leaves).
* `--root-dist <D>`: Cut at specific distance from root.
* `--max-clade <T>`: TreeCluster style (max pairwise distance in clade <= T).
* `--avg-clade <T>`: TreeCluster style (avg pairwise distance in clade <= T).
* `--med-clade <T>`: TreeCluster style (median pairwise distance in clade <= T).
* `--sum-branch <T>`: TreeCluster style (sum of branch lengths in clade <= T).
* `--leaf-dist-max <T>`: TreeCluster style (max distance from cluster root to any leaf <= T).
* `--leaf-dist-min <T>`: TreeCluster style (min distance from cluster root to any leaf <= T).
* `--leaf-dist-avg <T>`: TreeCluster style (avg distance from cluster root to leaves <= T).
* `--max-edge <T>` / `--single-linkage <T>`: Cut branches longer than threshold.
* `--inconsistent <T>`: SciPy style (inconsistent coefficient <= T).
* `--dynamic-tree <N>`: Dynamic Tree Cut (top-down adaptive, N=min cluster size).
* `--dynamic-hybrid <N>`: Hybrid Cut (Dynamic Tree + PAM, N=min cluster size).

Output formats:
    * cluster: Each line contains points of one cluster. The first point is the representative.
    * pair: Each line contains a (representative point, cluster member) pair.

Note:
The representative point is determined by `--rep` (applies to both 'cluster' and 'pair' formats):
* root (default): Member closest to root (alphabetical tie-break).
* medoid: Member with min sum of distances to others (alphabetical tie-break).
* first: Alphabetically first member.

Examples:
1. Cut into 5 clusters:
   pgr clust cut tree.nwk --k 5

2. Cut at height 0.5:
   pgr clust cut tree.nwk --height 0.5

3. Dynamic Tree Cut with min cluster size 20:
   pgr clust cut tree.nwk --dynamic-tree 20

4. Hybrid Cut with PAM (needs distance matrix):
   pgr clust cut tree.nwk --dynamic-hybrid 20 --matrix dist.phy
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input Newick file",
        ))
        .arg(crate::cmd_pgr::args::format_arg())
        .arg(crate::cmd_pgr::args::k_arg())
        .arg(
            Arg::new("height")
                .long("height")
                .value_parser(value_parser!(f64))
                .help("Cut at specific height (max distance to leaves)"),
        )
        .arg(
            Arg::new("root_dist")
                .long("root-dist")
                .value_parser(value_parser!(f64))
                .help("Cut at specific distance from root"),
        )
        .arg(
            Arg::new("max_clade")
                .long("max-clade")
                .value_parser(value_parser!(f64))
                .help("Max pairwise distance in cluster threshold"),
        )
        .arg(
            Arg::new("avg_clade")
                .long("avg-clade")
                .value_parser(value_parser!(f64))
                .help("Average pairwise distance in cluster threshold"),
        )
        .arg(
            Arg::new("med_clade")
                .long("med-clade")
                .value_parser(value_parser!(f64))
                .help("Median pairwise distance in cluster threshold"),
        )
        .arg(
            Arg::new("sum_branch")
                .long("sum-branch")
                .value_parser(value_parser!(f64))
                .help("Sum of branch lengths in cluster threshold"),
        )
        .arg(
            Arg::new("leaf_dist_max")
                .long("leaf-dist-max")
                .value_parser(value_parser!(f64))
                .help("Max distance from cluster root to any leaf"),
        )
        .arg(
            Arg::new("leaf_dist_min")
                .long("leaf-dist-min")
                .value_parser(value_parser!(f64))
                .help("Min distance from cluster root to any leaf"),
        )
        .arg(
            Arg::new("leaf_dist_avg")
                .long("leaf-dist-avg")
                .value_parser(value_parser!(f64))
                .help("Average distance from cluster root to leaves"),
        )
        .arg(
            Arg::new("max_edge")
                .long("max-edge")
                .alias("single-linkage")
                .value_parser(value_parser!(f64))
                .help("Cut branches longer than threshold (Single Linkage)"),
        )
        .arg(crate::cmd_pgr::args::rep_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("inconsistent")
                .long("inconsistent")
                .value_parser(value_parser!(f64))
                .help("Cut by inconsistent coefficient threshold"),
        )
        .arg(crate::cmd_pgr::args::deep_arg())
        .arg(crate::cmd_pgr::args::scan_arg())
        .arg(crate::cmd_pgr::args::stats_out_arg())
        .arg(crate::cmd_pgr::args::support_arg())
        .arg(crate::cmd_pgr::args::dynamic_tree_arg())
        .arg(crate::cmd_pgr::args::dynamic_hybrid_arg())
        .arg(crate::cmd_pgr::args::matrix_arg())
        .arg(crate::cmd_pgr::args::max_pam_dist_arg())
        .arg(crate::cmd_pgr::args::no_pam_dendro_arg())
        .arg(crate::cmd_pgr::args::deep_split_arg())
        .arg(crate::cmd_pgr::args::max_tree_height_arg())
        .group(
            ArgGroup::new("method")
                .args([
                    "k",
                    "height",
                    "root_dist",
                    "max_clade",
                    "avg_clade",
                    "med_clade",
                    "sum_branch",
                    "leaf_dist_max",
                    "leaf_dist_min",
                    "leaf_dist_avg",
                    "max_edge",
                    "inconsistent",
                    "dynamic_tree",
                    "dynamic_hybrid",
                ])
                .required(true),
        )
}
/// Detect which standard cut method was requested.
fn detect_method_name(args: &ArgMatches) -> anyhow::Result<&'static str> {
    METHOD_NAMES
        .iter()
        .find(|&&n| args.contains_id(n))
        .copied()
        .ok_or_else(|| anyhow::anyhow!("no cut method specified"))
}

/// Execute the cut command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .ok_or_else(|| anyhow::anyhow!("missing required argument: infile"))?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    // Remaining arguments have clap default values, so unwrap is safe.
    let format = args.get_one::<String>("clust_format").unwrap();
    let rep_method = args.get_one::<String>("rep").unwrap().as_str();
    let deep = *args.get_one::<usize>("deep").unwrap();

    let mut trees = Tree::from_file(infile)?;
    if trees.len() > 1 {
        anyhow::bail!("Input file contains multiple trees. Only single tree input is supported.");
    }
    if trees.is_empty() {
        anyhow::bail!("Input file contains no tree");
    }

    if let Some(&support_threshold) = args.get_one::<f64>("support") {
        for tree in &mut trees {
            pgr::libs::clust::tree_cut::apply_support_filter(tree, support_threshold);
        }
    }

    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    // Options common to dynamic methods
    let deep_split = args.get_flag("deep_split");
    let max_tree_height = args.get_one::<f64>("max_tree_height").copied();
    let max_pam_dist = args.get_one::<f64>("max_pam_dist").copied();
    let no_pam_dendro = args.get_flag("no_pam_dendro");

    let tree = &trees[0];

    if args.contains_id("scan") {
        return run_scan(
            args,
            tree,
            &mut writer,
            deep,
            max_tree_height,
            deep_split,
            no_pam_dendro,
            max_pam_dist,
        );
    }

    let rep_mode = RepMode::parse(rep_method)?;

    let dynamic_tree = args.get_one::<usize>("dynamic_tree").copied();
    let dynamic_hybrid = args.get_one::<usize>("dynamic_hybrid").copied();

    let matrix = if dynamic_hybrid.is_some() {
        let matrix_file = args
            .get_one::<String>("matrix")
            .ok_or_else(|| anyhow::anyhow!("--matrix is required for dynamic-hybrid"))?;
        Some(pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(
            matrix_file,
        )?)
    } else {
        None
    };

    let (method_name, val) = if dynamic_tree.is_none() && dynamic_hybrid.is_none() {
        let name = detect_method_name(args)?;
        let val = if name == "k" {
            *args
                .get_one::<usize>("k")
                .ok_or_else(|| anyhow::anyhow!("missing --k value"))? as f64
        } else {
            *args
                .get_one::<f64>(name)
                .ok_or_else(|| anyhow::anyhow!("missing --{} value", name))?
        };
        (Some(name), val)
    } else {
        (None, 0.0)
    };

    let dispatch = cut::build_dispatch(
        tree,
        method_name,
        val,
        deep,
        dynamic_tree,
        dynamic_hybrid,
        max_tree_height,
        deep_split,
        no_pam_dendro,
        max_pam_dist,
        matrix,
    )?;

    let (partition, _) = cut::dispatch_cut(tree, dispatch)?;

    let clusters = cut::partition_to_clusters(&partition, tree, rep_mode);
    let output = cut::format_clusters(&clusters, format)?;
    writer.write_all(output.as_bytes())?;

    Ok(())
}

/// Run the `--scan` parameter sweep over a single tree.
///
/// Writes a long-format table (Group, ClusterID, SampleID) to `writer` and,
/// if requested, summary statistics to `--stats-out`.
#[allow(clippy::too_many_arguments)]
fn run_scan(
    args: &ArgMatches,
    tree: &Tree,
    writer: &mut dyn Write,
    deep: usize,
    max_tree_height: Option<f64>,
    deep_split: bool,
    no_pam_dendro: bool,
    max_pam_dist: Option<f64>,
) -> anyhow::Result<()> {
    if args.contains_id("dynamic_hybrid") {
        anyhow::bail!("--scan is not supported with --dynamic-hybrid");
    }

    let scan_str = args.get_one::<String>("scan").unwrap();
    let (start, end, step) = parse_scan_range(scan_str)?;

    let mut stats_writer = init_stats_writer(args)?;
    writer.write_all(b"Group\tClusterID\tSampleID\n")?;

    let dynamic_tree = args.get_one::<usize>("dynamic_tree").copied();
    let method_name = if dynamic_tree.is_none() {
        Some(detect_method_name(args)?)
    } else {
        None
    };

    let n_steps = compute_n_steps(start, end, step)?;
    for i in 0..=n_steps {
        let val = start + (i as f64) * step;
        if val > end + 1e-9 {
            break;
        }

        let dispatch = if let Some(min_size) = dynamic_tree {
            build_dynamic_tree_dispatch(
                tree,
                val,
                deep,
                min_size,
                max_tree_height,
                deep_split,
                no_pam_dendro,
                max_pam_dist,
            )?
        } else {
            // detect_method_name already verified that a method is present.
            let name = method_name.unwrap();
            cut::build_dispatch(
                tree,
                Some(name),
                val,
                deep,
                None,
                None,
                max_tree_height,
                deep_split,
                no_pam_dendro,
                max_pam_dist,
                None,
            )?
        };

        let (partition, method_name) = cut::dispatch_cut(tree, dispatch)?;
        let group_label = format!("{}={}", method_name, val);

        if let Some(w) = &mut stats_writer {
            let (n_clusters, n_single, n_non_single, max_size) = partition.get_stats();
            w.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\t{}\n",
                group_label, n_clusters, n_single, n_non_single, max_size
            ))?;
        }

        let rows = cut::format_scan_rows(&partition, tree, &group_label);
        writer.write_all(rows.as_bytes())?;
    }

    writer.flush()?;
    Ok(())
}

/// Parse the `--scan` argument of the form `start,end,step`.
fn parse_scan_range(scan_str: &str) -> anyhow::Result<(f64, f64, f64)> {
    let parts: Vec<&str> = scan_str.split(',').collect();
    if parts.len() != 3 {
        anyhow::bail!("--scan format must be start,end,step");
    }
    let start: f64 = parts[0].parse()?;
    let end: f64 = parts[1].parse()?;
    let step: f64 = parts[2].parse()?;

    if step <= 0.0 {
        anyhow::bail!("Scan step must be positive");
    }
    Ok((start, end, step))
}

/// Compute the number of scan steps using integer arithmetic to avoid
/// floating-point drift.
fn compute_n_steps(start: f64, end: f64, step: f64) -> anyhow::Result<i64> {
    let n_steps_f = ((end - start) / step).round();
    if !n_steps_f.is_finite() || n_steps_f < 0.0 || n_steps_f > i64::MAX as f64 {
        anyhow::bail!(
            "scan range too large: start={}, end={}, step={}",
            start,
            end,
            step
        );
    }
    Ok(n_steps_f as i64)
}

/// Open the `--stats-out` writer and write its header.
fn init_stats_writer(args: &ArgMatches) -> anyhow::Result<Option<Box<dyn Write>>> {
    if let Some(stats_file) = args.get_one::<String>("stats_out") {
        let mut w = Box::new(
            pgr::writer(stats_file)
                .with_context(|| format!("Failed to open writer for {}", stats_file))?,
        );
        w.write_all(b"Group\tClusters\tSingletons\tNon-Singletons\tMaxSize\n")?;
        Ok(Some(w))
    } else {
        Ok(None)
    }
}

/// Build a dispatch for dynamic-tree scan values, validating that the value is
/// a non-negative integer.
#[allow(clippy::too_many_arguments)]
fn build_dynamic_tree_dispatch(
    tree: &Tree,
    val: f64,
    deep: usize,
    min_size: usize,
    max_tree_height: Option<f64>,
    deep_split: bool,
    no_pam_dendro: bool,
    max_pam_dist: Option<f64>,
) -> anyhow::Result<cut::CutDispatch> {
    if !val.is_finite() || val < 0.0 || val > usize::MAX as f64 {
        anyhow::bail!("scan value out of range: {}", val);
    }
    if val.fract() != 0.0 {
        anyhow::bail!("scan value must be integer for dynamic-tree: {}", val);
    }
    cut::build_dispatch(
        tree,
        None,
        val,
        deep,
        Some(min_size),
        None,
        max_tree_height,
        deep_split,
        no_pam_dendro,
        max_pam_dist,
        None,
    )
}
