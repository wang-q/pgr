use clap::{ArgMatches, Command};
/// Build the clap subcommand for subset.
pub fn make_subcommand() -> Command {
    Command::new("subset")
        .about("Extracts a submatrix from a PHYLIP matrix using a list of names")
        .after_help(
            r###"
Input formats:
    * Matrix: PHYLIP distance matrix (full or lower-triangular)
    * List: One name per line
    * Empty lines and lines starting with '#' in the list file are ignored

Examples:
    1. Create a full submatrix:
       pgr mat subset input.phy list.txt -o output.phy

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PHYLIP matrix file",
        ))
        .arg(crate::cmd_pgr::args::fa_name_list_arg(true))
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the subset command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let list_file = args.get_one::<String>("name_list").unwrap();
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let wanted_names = pgr::libs::io::read_names::<Vec<String>>(list_file)?;

    // Load and process matrix
    let matrix = pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(infile)?;

    let missing = pgr::libs::pairmat::write_subset(&matrix, &wanted_names, &mut writer)?;
    for name in &missing {
        log::warn!("Name not found in matrix: {}", name);
    }

    Ok(())
}
