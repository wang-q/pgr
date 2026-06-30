use clap::{Arg, ArgMatches, Command};
use pgr::libs::clust::eval::{
    load_batch_partitions, load_partition, remove_singletons, run_batch, run_single, Coordinates,
    DistanceMatrix, EvalTarget, PartitionFormat, TreeDistance,
};
use pgr::libs::pairmat::NamedMatrix;
use pgr::libs::phylo::tree::Tree;

pub fn make_subcommand() -> Command {
    Command::new("eval")
        .about("Evaluates clustering quality")
        .after_help(
            r###"Calculates clustering evaluation metrics.

Modes:
1. External Evaluation (Partition vs Partition):
   Compares two clustering partitions (e.g., ground truth vs result).
   Metrics: ARI, AMI, V-Measure.

2. Internal Evaluation (Partition + Matrix):
   Evaluates a single partition using a distance matrix.
   Metrics: Silhouette Coefficient.

3. Batch Evaluation (Long Format):
   Evaluates multiple partitions (e.g. from parameter scan) against a ground truth or using internal metrics.
   Input file must be in 'long' format (Group, Cluster, Sample).

Examples:
1. Compare result with ground truth:
   $ pgr clust eval result.tsv --other other.tsv -o eval.tsv

2. Evaluate result using distance matrix:
   $ pgr clust eval result.tsv --matrix dist.phy

3. Batch evaluation of scan results:
   $ pgr clust eval scan.tsv --format long --matrix dist.phy
"###,
        )
        .arg(
            Arg::new("p1")
                .required(true)
                .index(1)
                .help("Partition file"),
        )
        .arg(
            Arg::new("other")
                .long("other")
                .alias("truth")
                .num_args(1)
                .help("Other partition file (for external evaluation)"),
        )
        .arg(
            Arg::new("matrix")
                .long("matrix")
                .num_args(1)
                .help("Distance matrix file (for internal evaluation: Silhouette)"),
        )
        .arg(
            Arg::new("tree")
                .long("tree")
                .num_args(1)
                .help("Tree file (for internal evaluation: Silhouette, using patristic distance)"),
        )
        .arg(
            Arg::new("coords")
                .long("coords")
                .num_args(1)
                .help("Coordinate matrix file (for internal evaluation: Davies-Bouldin)"),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .value_parser(["cluster", "pair", "long"])
                .default_value("pair")
                .help("Input format for partition files"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("no-singletons")
                .long("no-singletons")
                .action(clap::ArgAction::SetTrue)
                .help("Exclude true singletons (from Reference/Ground Truth) from evaluation"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let p1_path = matches.get_one::<String>("p1").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(matches);

    let format_str = matches.get_one::<String>("format").unwrap();
    let format: PartitionFormat = format_str.parse().expect("Invalid format");

    let mut writer = pgr::writer(outfile)?;

    let remove_singletons_flag = matches.get_flag("no-singletons");

    if format == PartitionFormat::Long {
        // Batch Mode
        let batches = load_batch_partitions(p1_path)?;

        // Prepare resources (I/O stays in cmd layer)
        let p2 = if let Some(p2_path) = matches.get_one::<String>("other") {
            let mut truth = load_partition(p2_path, PartitionFormat::Pair)?;
            if remove_singletons_flag {
                remove_singletons(&mut truth);
            }
            Some(truth)
        } else {
            None
        };

        let dist_provider: Option<Box<dyn DistanceMatrix>> =
            if let Some(matrix_path) = matches.get_one::<String>("matrix") {
                Some(Box::new(NamedMatrix::from_relaxed_phylip(matrix_path)?))
            } else if let Some(tree_path) = matches.get_one::<String>("tree") {
                let trees = Tree::from_file(tree_path)?;
                if trees.len() != 1 {
                    anyhow::bail!("Tree file must contain exactly one tree.");
                }
                Some(Box::new(TreeDistance::new(
                    trees.into_iter().next().unwrap(),
                )))
            } else {
                None
            };

        let coords = if let Some(coords_path) = matches.get_one::<String>("coords") {
            Some(Coordinates::from_path(coords_path)?)
        } else {
            None
        };

        if p2.is_none() && dist_provider.is_none() && coords.is_none() {
            anyhow::bail!(
                "Batch mode requires at least one evaluation target: --other/--truth, --matrix, --tree, or --coords."
            );
        }

        let mut targets: Vec<EvalTarget<'_>> = vec![];
        if let Some(ref truth) = p2 {
            targets.push(EvalTarget::External(truth));
        }
        if let Some(ref d) = dist_provider {
            targets.push(EvalTarget::Matrix(&**d));
        }
        if let Some(ref c) = coords {
            targets.push(EvalTarget::Coords(c));
        }

        run_batch(batches, &targets, &mut writer)?;
        return Ok(());
    }

    // Single Mode
    let p1 = load_partition(p1_path, format)?;

    if let Some(p2_path) = matches.get_one::<String>("other") {
        let mut p2 = load_partition(p2_path, format)?;
        if remove_singletons_flag {
            remove_singletons(&mut p2);
        }
        run_single(&p1, EvalTarget::External(&p2), &mut writer)?;
    } else if let Some(matrix_path) = matches.get_one::<String>("matrix") {
        let matrix = NamedMatrix::from_relaxed_phylip(matrix_path)?;
        run_single(&p1, EvalTarget::Matrix(&matrix), &mut writer)?;
    } else if let Some(tree_path) = matches.get_one::<String>("tree") {
        let trees = Tree::from_file(tree_path)?;
        if trees.len() != 1 {
            anyhow::bail!("Tree file must contain exactly one tree.");
        }
        let dist = TreeDistance::new(trees.into_iter().next().unwrap());
        run_single(&p1, EvalTarget::Matrix(&dist), &mut writer)?;
    } else if let Some(coords_path) = matches.get_one::<String>("coords") {
        let coords = Coordinates::from_path(coords_path)?;
        run_single(&p1, EvalTarget::Coords(&coords), &mut writer)?;
    } else {
        anyhow::bail!(
            "Either --other/--truth (for external eval), --matrix, --tree, or --coords (for internal eval) must be provided."
        );
    }

    Ok(())
}
