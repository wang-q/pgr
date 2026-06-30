use clap::*;
use std::io::Write;

pub fn make_subcommand() -> Command {
    Command::new("compare")
        .about("Compare two distance matrices")
        .after_help(
            r###"
Compare two PHYLIP distance matrices and calculate similarity metrics.

Methods:
    * all:       Calculate all metrics below
    * pearson:   Pearson correlation coefficient (-1 to 1)
    * spearman:  Spearman rank correlation (-1 to 1)
    * mae:       Mean absolute error
    * cosine:    Cosine similarity (-1 to 1)
    * jaccard:   Weighted Jaccard similarity (0 to 1)
    * euclid:    Euclidean distance

Examples:
    # Compare using Pearson correlation
    pgr mat compare matrix1.phy matrix2.phy --method pearson

    # Compare using multiple methods
    pgr mat compare matrix1.phy matrix2.phy --method pearson,cosine,jaccard
"###,
        )
        .arg(
            Arg::new("matrix1")
                .required(true)
                .index(1)
                .help("First PHYLIP matrix file"),
        )
        .arg(
            Arg::new("matrix2")
                .required(true)
                .index(2)
                .help("Second PHYLIP matrix file"),
        )
        .arg(
            Arg::new("method")
                .long("method")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("all"),
                    builder::PossibleValue::new("pearson"),
                    builder::PossibleValue::new("spearman"),
                    builder::PossibleValue::new("mae"),
                    builder::PossibleValue::new("cosine"),
                    builder::PossibleValue::new("jaccard"),
                    builder::PossibleValue::new("euclid"),
                ])
                .default_value("pearson")
                .help("Comparison method(s), comma-separated"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let matrix1_file = args.get_one::<String>("matrix1").unwrap();
    let matrix2_file = args.get_one::<String>("matrix2").unwrap();
    let methods = if args.get_one::<String>("method").unwrap() == "all" {
        "pearson,spearman,mae,cosine,jaccard,euclid"
    } else {
        args.get_one::<String>("method").unwrap()
    };
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap())?;

    // Load matrices
    let matrix1 = pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(matrix1_file)?;
    let matrix2 = pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(matrix2_file)?;

    // Get common sequence names
    let names1 = matrix1.get_names();
    let names2 = matrix2.get_names();
    let common_names: Vec<_> = names1.iter().filter(|name| names2.contains(name)).collect();

    // Report sequence counts
    log::info!(
        "Sequences in matrices: {} and {}",
        names1.len(),
        names2.len()
    );
    log::info!("Common sequences: {}", common_names.len());

    if common_names.is_empty() {
        return Err(anyhow::anyhow!(
            "No common sequence names found between matrices"
        ));
    }

    // Extract values for comparison
    let mut values1 = Vec::with_capacity(common_names.len() * (common_names.len() - 1) / 2);
    let mut values2 = Vec::with_capacity(common_names.len() * (common_names.len() - 1) / 2);

    for i in 0..common_names.len() {
        for j in 0..i {
            if let (Some(v1), Some(v2)) = (
                matrix1.get_by_name(common_names[i], common_names[j]),
                matrix2.get_by_name(common_names[i], common_names[j]),
            ) {
                values1.push(v1);
                values2.push(v2);
            }
        }
    }

    // Write header
    writer.write_all(b"Method\tScore\n")?;

    // Calculate and output metrics
    for method in methods.split(',') {
        let result = match method {
            "pearson" => pgr::libs::linalg::pearson_correlation(&values1, &values2),
            "spearman" => pgr::libs::linalg::spearman_correlation(&values1, &values2),
            "mae" => pgr::libs::linalg::mean_absolute_error(&values1, &values2),
            "cosine" => pgr::libs::linalg::cosine_similarity(&values1, &values2),
            "jaccard" => pgr::libs::linalg::weighted_jaccard_similarity(&values1, &values2),
            "euclid" => pgr::libs::linalg::euclidean_distance(&values1, &values2),
            _ => anyhow::bail!("unknown method: {}", method),
        };
        writer.write_fmt(format_args!("{}\t{:.6}\n", method, result))?;
    }

    Ok(())
}
