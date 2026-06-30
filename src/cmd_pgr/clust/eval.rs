use clap::{Arg, ArgMatches, Command};
use pgr::libs::clust::eval::format::{
    coord_metric_values, distance_metric_values, external_metric_values, format_metrics_row,
    COORD_METRIC_NAMES, DISTANCE_METRIC_NAMES, EXTERNAL_METRIC_NAMES,
};
use pgr::libs::clust::eval::{
    evaluate, load_batch_partitions, load_partition, remove_singletons, Coordinates,
    DistanceMatrix, PartitionFormat, TreeDistance,
};
use pgr::libs::pairmat::NamedMatrix;
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

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

        // Prepare resources
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

        // Write Header
        let mut header = vec!["Group"];
        if p2.is_some() {
            header.extend_from_slice(EXTERNAL_METRIC_NAMES);
        }
        if dist_provider.is_some() {
            header.extend_from_slice(DISTANCE_METRIC_NAMES);
        }
        if coords.is_some() {
            header.extend_from_slice(COORD_METRIC_NAMES);
        }
        writeln!(writer, "{}", header.join("\t"))?;

        // Process batches
        for (group, p1) in batches {
            let mut row = vec![group];

            if let Some(ref truth) = p2 {
                let metrics = evaluate(&p1, truth);
                row.push(format_metrics_row(&external_metric_values(&metrics)));
            }

            if let Some(ref d) = dist_provider {
                let values = distance_metric_values(&p1, d.as_ref());
                row.push(format_metrics_row(&values));
            }

            if let Some(ref c) = coords {
                let values = coord_metric_values(&p1, c);
                row.push(format_metrics_row(&values));
            }

            writeln!(writer, "{}", row.join("\t"))?;
        }

        return Ok(());
    }

    // Single Mode
    let p1 = load_partition(p1_path, format)?;

    if let Some(p2_path) = matches.get_one::<String>("other") {
        let mut p2 = load_partition(p2_path, format)?;
        if remove_singletons_flag {
            remove_singletons(&mut p2);
        }
        let metrics = evaluate(&p1, &p2);

        writeln!(writer, "{}", EXTERNAL_METRIC_NAMES.join("\t"))?;
        writeln!(
            writer,
            "{}",
            format_metrics_row(&external_metric_values(&metrics))
        )?;
    } else if let Some(matrix_path) = matches.get_one::<String>("matrix") {
        let matrix = NamedMatrix::from_relaxed_phylip(matrix_path)?;
        let values = distance_metric_values(&p1, &matrix);

        writeln!(writer, "{}", DISTANCE_METRIC_NAMES.join("\t"))?;
        writeln!(writer, "{}", format_metrics_row(&values))?;
    } else if let Some(tree_path) = matches.get_one::<String>("tree") {
        let trees = Tree::from_file(tree_path)?;
        if trees.len() != 1 {
            anyhow::bail!("Tree file must contain exactly one tree.");
        }
        let dist = TreeDistance::new(trees.into_iter().next().unwrap());
        let values = distance_metric_values(&p1, &dist);

        writeln!(writer, "{}", DISTANCE_METRIC_NAMES.join("\t"))?;
        writeln!(writer, "{}", format_metrics_row(&values))?;
    } else if let Some(coords_path) = matches.get_one::<String>("coords") {
        let coords = Coordinates::from_path(coords_path)?;
        let values = coord_metric_values(&p1, &coords);

        writeln!(writer, "{}", COORD_METRIC_NAMES.join("\t"))?;
        writeln!(writer, "{}", format_metrics_row(&values))?;
    } else {
        anyhow::bail!(
            "Either --other/--truth (for external eval), --matrix, --tree, or --coords (for internal eval) must be provided."
        );
    }

    Ok(())
}
