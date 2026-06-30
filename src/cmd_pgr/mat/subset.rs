use clap::*;

pub fn make_subcommand() -> Command {
    Command::new("subset")
        .about("Extract a submatrix from a PHYLIP matrix using a list of names")
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
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input PHYLIP matrix file"),
        )
        .arg(
            Arg::new("list")
                .required(true)
                .index(2)
                .help("File containing sequence names to extract"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let list_file = args.get_one::<String>("list").unwrap();
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let wanted_names = pgr::libs::io::read_names::<Vec<String>>(list_file)?;

    //----------------------------
    // Load and process matrix
    //----------------------------
    let matrix = pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(infile)?;

    let missing = pgr::libs::pairmat::write_subset(&matrix, &wanted_names, &mut writer)?;
    for name in &missing {
        log::warn!("Name not found in matrix: {}", name);
    }

    Ok(())
}
