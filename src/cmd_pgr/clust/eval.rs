use clap::{Arg, ArgMatches, Command};
use pgr::libs::clust::eval::{
    davies_bouldin_score, evaluate, load_batch_partitions, load_partition, silhouette_score,
    Coordinates, PartitionFormat,
};
use pgr::libs::pairmat::NamedMatrix;
use std::fs::File;
use std::io::{self, Write};

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
   $ pgr clust eval result.tsv truth.tsv -o eval.tsv

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
            Arg::new("p2")
                .required(false)
                .index(2)
                .help("Second partition file (Ground Truth). If p1 is batch/long, p2 is assumed to be in 'pair' format."),
        )
        .arg(
            Arg::new("matrix")
                .long("matrix")
                .num_args(1)
                .help("Distance matrix file (for internal evaluation: Silhouette)"),
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
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let p1_path = matches.get_one::<String>("p1").unwrap();
    let outfile = matches.get_one::<String>("outfile").unwrap();

    let format_str = matches.get_one::<String>("format").unwrap();
    let format: PartitionFormat = format_str.parse().expect("Invalid format");

    let mut writer: Box<dyn Write> = if outfile == "stdout" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(outfile)?)
    };

    if format == PartitionFormat::Long {
        // Batch Mode
        let batches = load_batch_partitions(p1_path)?;

        // Prepare resources
        let p2 = if let Some(p2_path) = matches.get_one::<String>("p2") {
            Some(load_partition(p2_path, PartitionFormat::Pair)?)
        } else {
            None
        };

        let matrix = if let Some(matrix_path) = matches.get_one::<String>("matrix") {
            Some(NamedMatrix::from_relaxed_phylip(matrix_path))
        } else {
            None
        };

        let coords = if let Some(coords_path) = matches.get_one::<String>("coords") {
            Some(Coordinates::from_path(coords_path)?)
        } else {
            None
        };

        if p2.is_none() && matrix.is_none() && coords.is_none() {
            anyhow::bail!(
                "Batch mode requires at least one evaluation target: <p2>, --matrix, or --coords."
            );
        }

        // Write Header
        let mut header = vec!["Group"];
        if p2.is_some() {
            header.extend_from_slice(&[
                "ari",
                "ami",
                "homogeneity",
                "completeness",
                "v_measure",
                "fmi",
                "nmi",
                "mi",
            ]);
        }
        if matrix.is_some() {
            header.push("silhouette");
        }
        if coords.is_some() {
            header.push("davies_bouldin");
        }
        writeln!(writer, "{}", header.join("\t"))?;

        // Process batches
        for (group, p1) in batches {
            let mut row = vec![group];

            if let Some(ref truth) = p2 {
                let metrics = evaluate(&p1, truth);
                row.push(format!("{:.6}", metrics.ari));
                row.push(format!("{:.6}", metrics.ami));
                row.push(format!("{:.6}", metrics.homogeneity));
                row.push(format!("{:.6}", metrics.completeness));
                row.push(format!("{:.6}", metrics.v_measure));
                row.push(format!("{:.6}", metrics.fmi));
                row.push(format!("{:.6}", metrics.nmi));
                row.push(format!("{:.6}", metrics.mi));
            }

            if let Some(ref m) = matrix {
                let score = silhouette_score(&p1, m);
                row.push(format!("{:.6}", score));
            }

            if let Some(ref c) = coords {
                let score = davies_bouldin_score(&p1, c);
                row.push(format!("{:.6}", score));
            }

            writeln!(writer, "{}", row.join("\t"))?;
        }

        return Ok(());
    }

    // Single Mode
    let p1 = load_partition(p1_path, format)?;

    if let Some(p2_path) = matches.get_one::<String>("p2") {
        let p2 = load_partition(p2_path, format)?;
        let metrics = evaluate(&p1, &p2);

        writeln!(
            writer,
            "ari\tami\thomogeneity\tcompleteness\tv_measure\tfmi\tnmi\tmi"
        )?;
        writeln!(
            writer,
            "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            metrics.ari,
            metrics.ami,
            metrics.homogeneity,
            metrics.completeness,
            metrics.v_measure,
            metrics.fmi,
            metrics.nmi,
            metrics.mi
        )?;
    } else if let Some(matrix_path) = matches.get_one::<String>("matrix") {
        let matrix = NamedMatrix::from_relaxed_phylip(matrix_path);
        let score = silhouette_score(&p1, &matrix);

        writeln!(writer, "silhouette")?;
        writeln!(writer, "{:.6}", score)?;
    } else if let Some(coords_path) = matches.get_one::<String>("coords") {
        let coords = Coordinates::from_path(coords_path)?;
        let score = davies_bouldin_score(&p1, &coords);

        writeln!(writer, "davies_bouldin")?;
        writeln!(writer, "{:.6}", score)?;
    } else {
        anyhow::bail!(
            "Either <p2> (for external eval), --matrix, or --coords (for internal eval) must be provided."
        );
    }

    Ok(())
}
