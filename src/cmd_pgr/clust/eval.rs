use clap::{Arg, ArgMatches, Command};
use pgr::libs::clust::eval::{
    davies_bouldin_score, evaluate, load_partition, silhouette_score, Coordinates, PartitionFormat,
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

Examples:
1. Compare result with ground truth:
   $ pgr clust eval result.tsv truth.tsv -o eval.tsv

2. Evaluate result using distance matrix:
   $ pgr clust eval result.tsv --matrix dist.phy
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
                .help("Second partition file (for external evaluation)"),
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
                .value_parser(["cluster", "pair"])
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

    let p1 = load_partition(p1_path, format)?;

    let mut writer: Box<dyn Write> = if outfile == "stdout" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(outfile)?)
    };

    if let Some(p2_path) = matches.get_one::<String>("p2") {
        let p2 = load_partition(p2_path, format)?;
        let metrics = evaluate(&p1, &p2);

        writeln!(writer, "ari\tami\thomogeneity\tcompleteness\tv_measure")?;
        writeln!(
            writer,
            "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            metrics.ari, metrics.ami, metrics.homogeneity, metrics.completeness, metrics.v_measure
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
