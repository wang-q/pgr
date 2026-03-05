use clap::{Arg, ArgMatches, Command};
use pgr::libs::clust::eval::{evaluate, load_partition, PartitionFormat};
use std::fs::File;
use std::io::{self, Write};

pub fn make_subcommand() -> Command {
    Command::new("eval")
        .about("Evaluates clustering quality by comparing two partitions")
        .after_help(
            r###"Calculates clustering evaluation metrics (ARI, AMI, V-Measure).

This command compares two clustering partitions (e.g., ground truth vs result, or result A vs result B).
It supports two input formats:
* Cluster-based: Each line is a cluster, items separated by tabs.
* Pair-based: Each line is "Representative/Label <tab> Member".

Metrics:
* ARI (Adjusted Rand Index): [-1, 1], 1=perfect, 0=random.
* AMI (Adjusted Mutual Information): [0, 1], 1=perfect, 0=random.
* V-Measure: Harmonic mean of Homogeneity and Completeness.

Examples:
1. Compare clustering result with ground truth (default pair format):
   $ pgr clust eval result.tsv truth.tsv -o eval.tsv

2. Compare two clustering results in cluster format:
   $ pgr clust eval method1.tsv method2.tsv --format cluster
"###,
        )
        .arg(
            Arg::new("p1")
                .required(true)
                .index(1)
                .help("First partition file"),
        )
        .arg(
            Arg::new("p2")
                .required(true)
                .index(2)
                .help("Second partition file"),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .value_parser(["cluster", "pair"])
                .default_value("pair")
                .help("Input format for both partition files"),
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
    let p2_path = matches.get_one::<String>("p2").unwrap();
    let outfile = matches.get_one::<String>("outfile").unwrap();

    let format_str = matches.get_one::<String>("format").unwrap();
    let format: PartitionFormat = format_str.parse().expect("Invalid format");

    let p1 = load_partition(p1_path, format)?;
    let p2 = load_partition(p2_path, format)?;

    let metrics = evaluate(&p1, &p2);

    let mut writer: Box<dyn Write> = if outfile == "stdout" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(outfile)?)
    };

    writeln!(writer, "ari\tami\thomogeneity\tcompleteness\tv_measure")?;
    writeln!(
        writer,
        "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
        metrics.ari, metrics.ami, metrics.homogeneity, metrics.completeness, metrics.v_measure
    )?;

    Ok(())
}
