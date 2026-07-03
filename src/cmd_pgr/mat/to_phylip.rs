use clap::{ArgMatches, Command};
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("to-phylip")
        .about("Converts pairwise distances to a phylip distance matrix")
        .after_help(
            r###"
Input format:
    * Tab-separated values (TSV)
    * Three columns: name1, name2, distance

Examples:
    1. Convert pairwise distances to PHYLIP matrix:
       pgr mat to-phylip input.tsv -o output.phy
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input file containing pairwise distances",
        ))
        .arg(crate::cmd_pgr::args::same_arg("0.0"))
        .arg(crate::cmd_pgr::args::missing_arg("1.0"))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let opt_same = *args.get_one::<f32>("same").unwrap();
    let opt_missing = *args.get_one::<f32>("missing").unwrap();
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    //----------------------------
    // Ops
    //----------------------------
    // Load matrix from pairwise distances
    let matrix = pgr::libs::pairmat::NamedMatrix::from_pair_scores(infile, opt_same, opt_missing)?;
    let names = matrix.get_names();
    let size = matrix.size();

    // Write sequence count
    writer.write_fmt(format_args!("{:>4}\n", size))?;

    // Output full matrix
    for (i, name) in names.iter().enumerate().take(size) {
        writer.write_fmt(format_args!("{}", name))?;
        for j in 0..size {
            writer.write_fmt(format_args!("\t{}", matrix.get(i, j)))?;
        }
        writer.write_fmt(format_args!("\n"))?;
    }

    Ok(())
}
