use clap::{Arg, ArgMatches, Command};
use std::io::Write;

pub fn make_subcommand() -> Command {
    Command::new("compare")
        .about("Compares two distance matrices")
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
1. Compare using Pearson correlation:
   pgr mat compare matrix1.phy matrix2.phy --method pearson

2. Compare using multiple methods:
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
        .arg(crate::cmd_pgr::args::mat_method_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let matrix1_file = args.get_one::<String>("matrix1").unwrap();
    let matrix2_file = args.get_one::<String>("matrix2").unwrap();
    let methods = if args.get_one::<String>("mat_method").unwrap() == "all" {
        "pearson,spearman,mae,cosine,jaccard,euclid"
    } else {
        args.get_one::<String>("mat_method").unwrap()
    };
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    // Load matrices
    let matrix1 = pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(matrix1_file)?;
    let matrix2 = pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(matrix2_file)?;

    // Report sequence counts
    log::info!(
        "Sequences in matrices: {} and {}",
        matrix1.size(),
        matrix2.size()
    );

    // Extract paired values from common upper triangle
    let (common_names, values1, values2) =
        pgr::libs::pairmat::extract_common_upper_triangle(&matrix1, &matrix2)?;

    log::info!("Common sequences: {}", common_names.len());

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
