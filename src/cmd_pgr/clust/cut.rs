use anyhow::Context;
use clap::{builder::PossibleValue, value_parser, Arg, ArgAction, ArgGroup, ArgMatches, Command};
use pgr::libs::clust::tree_cut::dynamic::DynamicTreeOptions;
use pgr::libs::clust::tree_cut::hybrid::HybridOptions;
use pgr::libs::clust::tree_cut::{self as cut, CutDispatch, RepMode, METHOD_NAMES};
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
        .arg(
            Arg::new("rep")
                .long("rep")
                .value_parser([
                    PossibleValue::new("root"),
                    PossibleValue::new("first"),
                    PossibleValue::new("medoid"),
                ])
                .default_value("root")
                .help("Representative selection method"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("inconsistent")
                .long("inconsistent")
                .value_parser(value_parser!(f64))
                .help("Cut by inconsistent coefficient threshold"),
        )
        .arg(
            Arg::new("deep")
                .long("deep")
                .value_parser(value_parser!(usize))
                .default_value("2")
                .help("Depth for inconsistent coefficient calculation (default: 2)"),
        )
        .arg(
            Arg::new("scan")
                .long("scan")
                .help("Scan thresholds (format: start,end,step)"),
        )
        .arg(
            Arg::new("stats_out")
                .long("stats-out")
                .help("Output statistics to a separate file (useful when format is 'long')"),
        )
        .arg(
            Arg::new("support")
                .long("support")
                .value_parser(value_parser!(f64))
                .help("Branch support threshold (edges with support < S will be treated as infinite length)"),
        )
        .arg(
            Arg::new("dynamic_tree")
                .long("dynamic-tree")
                .value_parser(value_parser!(usize))
                .help("Use dynamic tree cut method (value: min cluster size)"),
        )
        .arg(
            Arg::new("dynamic_hybrid")
                .long("dynamic-hybrid")
                .value_parser(value_parser!(usize))
                .help("Use dynamic hybrid cut method (value: min cluster size)"),
        )
        .arg(crate::cmd_pgr::args::matrix_arg())
        .arg(
            Arg::new("max_pam_dist")
                .long("max-pam-dist")
                .value_parser(value_parser!(f64))
                .help("Maximum distance to medoid for PAM reassignment"),
        )
        .arg(
            Arg::new("no_pam_dendro")
                .long("no-pam-dendro")
                .action(ArgAction::SetTrue)
                .help("Disable dendrogram respect in PAM stage (allow assigning to clusters across high branches)"),
        )
        .arg(
            Arg::new("deep_split")
                .long("deep-split")
                .action(ArgAction::SetTrue)
                .help("Enable deep split for dynamic tree cut (default: false)"),
        )
        .arg(
            Arg::new("max_tree_height")
                .long("max-tree-height")
                .value_parser(value_parser!(f64))
                .help("Maximum joining height for dynamic tree cut (default: 99% of tree height)"),
        )
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
/// Execute the cut command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
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

    if let Some(scan_str) = args.get_one::<String>("scan") {
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

        let mut stats_writer: Option<Box<dyn Write>> =
            if let Some(stats_file) = args.get_one::<String>("stats_out") {
                let mut w = Box::new(
                    pgr::writer(stats_file)
                        .with_context(|| format!("Failed to open writer for {}", stats_file))?,
                );
                w.write_all(b"Group\tClusters\tSingletons\tNon-Singletons\tMaxSize\n")?;
                Some(w)
            } else {
                None
            };

        writer.write_all(b"Group\tClusterID\tSampleID\n")?;

        let tree = &trees[0];

        // Pre-calculate leaf depths for scanning if needed
        let leaf_depths_scan = if args.contains_id("leaf_dist_max")
            || args.contains_id("leaf_dist_min")
            || args.contains_id("leaf_dist_avg")
        {
            Some(pgr::libs::phylo::tree::stat::get_leaf_depth_stats(tree))
        } else {
            None
        };

        // Use integer step count + multiplication to avoid floating-point
        // accumulation drift (e.g. 0.1 added 1000 times != 100.0).
        let n_steps_f = ((end - start) / step).round();
        if !n_steps_f.is_finite() || n_steps_f < 0.0 || n_steps_f > i64::MAX as f64 {
            anyhow::bail!(
                "scan range too large: start={}, end={}, step={}",
                start,
                end,
                step
            );
        }
        let n_steps = n_steps_f as i64;
        for i in 0..=n_steps {
            let val = start + (i as f64) * step;
            if val > end + 1e-9 {
                break;
            }

            let dispatch = if args.contains_id("dynamic_tree") {
                if !val.is_finite() || val < 0.0 || val > usize::MAX as f64 {
                    anyhow::bail!("scan value out of range: {}", val);
                }
                if val.fract() != 0.0 {
                    anyhow::bail!("scan value must be integer for dynamic-tree: {}", val);
                }
                CutDispatch::DynamicTree(DynamicTreeOptions {
                    min_module_size: val as usize,
                    deep_split: args.get_flag("deep_split"),
                    max_tree_height: args.get_one::<f64>("max_tree_height").copied(),
                })
            } else {
                let method_name = METHOD_NAMES
                    .iter()
                    .find(|&&n| args.contains_id(n))
                    .copied()
                    .ok_or_else(|| anyhow::anyhow!("no cut method specified"))?;
                CutDispatch::Standard {
                    name: method_name,
                    val,
                    deep,
                    leaf_depths: leaf_depths_scan,
                }
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
        return Ok(());
    }

    let rep_mode = RepMode::parse(rep_method)
        .map_err(|e| anyhow::anyhow!("invalid --rep '{}': {}", rep_method, e))?;

    for tree in trees.iter() {
        let dispatch = if let Some(&min_cluster_size) = args.get_one::<usize>("dynamic_tree") {
            CutDispatch::DynamicTree(DynamicTreeOptions {
                min_module_size: min_cluster_size,
                deep_split: args.get_flag("deep_split"),
                max_tree_height: args.get_one::<f64>("max_tree_height").copied(),
            })
        } else if let Some(&min_cluster_size) = args.get_one::<usize>("dynamic_hybrid") {
            let matrix_file = args
                .get_one::<String>("matrix")
                .ok_or_else(|| anyhow::anyhow!("--matrix is required for dynamic-hybrid"))?;
            let dist_matrix = pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(matrix_file)?;

            CutDispatch::DynamicHybrid(HybridOptions {
                min_cluster_size,
                dist_matrix,
                cut_height: args.get_one::<f64>("max_tree_height").copied(),
                deep_split: if args.get_flag("deep_split") { 1 } else { 0 },
                max_core_scatter: None,
                min_gap: None,
                pam_stage: true, // Default to true
                pam_respects_dendro: !args.get_flag("no_pam_dendro"),
                max_pam_dist: args.get_one::<f64>("max_pam_dist").copied(),
                respect_small_clusters: true, // Default to true to match R
            })
        } else {
            let method_name = METHOD_NAMES
                .iter()
                .find(|&&n| args.contains_id(n))
                .copied()
                .ok_or_else(|| anyhow::anyhow!("no cut method specified"))?;
            let val = if method_name == "k" {
                *args
                    .get_one::<usize>("k")
                    .ok_or_else(|| anyhow::anyhow!("missing --k value"))? as f64
            } else {
                *args
                    .get_one::<f64>(method_name)
                    .ok_or_else(|| anyhow::anyhow!("missing --{} value", method_name))?
            };
            let leaf_depths = if method_name.starts_with("leaf_dist_") {
                Some(pgr::libs::phylo::tree::stat::get_leaf_depth_stats(tree))
            } else {
                None
            };
            CutDispatch::Standard {
                name: method_name,
                val,
                deep,
                leaf_depths,
            }
        };

        let (partition, _) = cut::dispatch_cut(tree, dispatch)?;

        let clusters = cut::partition_to_clusters(&partition, tree, rep_mode);
        let output = cut::format_clusters(&clusters, format).map_err(|e| anyhow::anyhow!(e))?;
        writer.write_all(output.as_bytes())?;
    }

    Ok(())
}
