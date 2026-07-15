use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::clust::eval::{
    load_batch_partitions, load_partition, remove_singletons, run_batch, run_single, Coordinates,
    DistanceMatrix, EvalTarget, PartitionFormat, TreeDistance,
};
use pgr::libs::pairmat::NamedMatrix;
use pgr::libs::phylo::tree::Tree;
use std::io::Write;
/// Build the clap subcommand for eval.
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
        .arg(crate::cmd_pgr::args::other_partition_arg())
        .arg(crate::cmd_pgr::args::matrix_arg())
        .arg(crate::cmd_pgr::args::tree_arg())
        .arg(crate::cmd_pgr::args::coords_arg())
        .arg(crate::cmd_pgr::args::clust_input_format_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(crate::cmd_pgr::args::no_singletons_arg())
}
/// Execute the eval command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let p1_path = args
        .get_one::<String>("p1")
        .ok_or_else(|| anyhow::anyhow!("missing required argument: p1"))?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    // clust_input_format has a clap default value, so unwrap is safe.
    let format_str = args.get_one::<String>("clust_input_format").unwrap();
    let format: PartitionFormat = match format_str.parse() {
        Ok(f) => f,
        Err(e) => anyhow::bail!("Invalid format: {}", e),
    };

    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    let remove_singletons_flag = args.get_flag("no_singletons");

    if format == PartitionFormat::Long {
        // Batch Mode
        let batches = load_batch_partitions(p1_path)?;

        // Prepare resources (I/O stays in cmd layer)
        let p2 = if let Some(p2_path) = args.get_one::<String>("other") {
            let mut truth = load_partition(p2_path, PartitionFormat::Pair)?;
            if remove_singletons_flag {
                remove_singletons(&mut truth);
            }
            Some(truth)
        } else {
            None
        };

        let dist_provider: Option<Box<dyn DistanceMatrix>> =
            if let Some(matrix_path) = args.get_one::<String>("matrix") {
                Some(Box::new(NamedMatrix::from_relaxed_phylip(matrix_path)?))
            } else if let Some(tree_path) = args.get_one::<String>("tree") {
                let trees = Tree::from_file(tree_path)?;
                if trees.len() != 1 {
                    anyhow::bail!("Tree file must contain exactly one tree.");
                }
                let tree = trees
                    .into_iter()
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("expected exactly one tree"))?;
                Some(Box::new(TreeDistance::new(tree)))
            } else {
                None
            };

        let coords = if let Some(coords_path) = args.get_one::<String>("coords") {
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
        writer.flush()?;
        return Ok(());
    }

    // Single Mode
    let p1 = load_partition(p1_path, format)?;

    if let Some(p2_path) = args.get_one::<String>("other") {
        let mut p2 = load_partition(p2_path, format)?;
        if remove_singletons_flag {
            remove_singletons(&mut p2);
        }
        run_single(&p1, EvalTarget::External(&p2), &mut writer)?;
    } else if let Some(matrix_path) = args.get_one::<String>("matrix") {
        let matrix = NamedMatrix::from_relaxed_phylip(matrix_path)?;
        run_single(&p1, EvalTarget::Matrix(&matrix), &mut writer)?;
    } else if let Some(tree_path) = args.get_one::<String>("tree") {
        let trees = Tree::from_file(tree_path)?;
        if trees.len() != 1 {
            anyhow::bail!("Tree file must contain exactly one tree.");
        }
        let tree = trees
            .into_iter()
            .next()
            .ok_or_else(|| anyhow::anyhow!("expected exactly one tree"))?;
        let dist = TreeDistance::new(tree);
        run_single(&p1, EvalTarget::Matrix(&dist), &mut writer)?;
    } else if let Some(coords_path) = args.get_one::<String>("coords") {
        let coords = Coordinates::from_path(coords_path)?;
        run_single(&p1, EvalTarget::Coords(&coords), &mut writer)?;
    } else {
        anyhow::bail!(
            "Either --other/--truth (for external eval), --matrix, --tree, or --coords (for internal eval) must be provided."
        );
    }

    writer.flush()?;
    Ok(())
}
