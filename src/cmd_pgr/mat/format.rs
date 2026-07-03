use clap::*;

pub fn make_subcommand() -> Command {
    Command::new("format")
        .about("Convert between different PHYLIP matrix formats")
        .after_help(
            r###"
Convert a PHYLIP matrix between different formats while preserving all distance values.

Input format:
    * PHYLIP distance matrix (full or lower-triangular)
    * Optional first line: number of sequences
    * Each line: sequence name followed by distances

Output formats:
    * full
        - Full square matrix
        - Tab-separated values
        - Original sequence names preserved
    * lower
        - Lower triangular matrix
        - Tab-separated values
        - Original sequence names preserved
    * strict
        - Standard PHYLIP format
        - Names truncated to 10 characters
        - Names left-aligned with space padding
        - Distances in fixed-width format (6 decimal places)
        - Space-separated values

Examples:
    1. Create a full matrix:
       pgr mat format input.phy -o output.phy

    2. Create a lower-triangular matrix:
       pgr mat format input.phy --format lower -o output.phy

    3. Create a strict PHYLIP matrix:
       pgr mat format input.phy --format strict -o output.phy
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PHYLIP matrix file",
        ))
        .arg(
            Arg::new("format")
                .long("format")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("full"),
                    builder::PossibleValue::new("lower"),
                    builder::PossibleValue::new("strict"),
                ])
                .default_value("full")
                .help("Output format"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let opt_mode = args.get_one::<String>("format").unwrap();
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    //----------------------------
    // Ops
    //----------------------------
    let matrix = pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(infile)?;
    let fmt = pgr::libs::pairmat::MatrixFormat::from_mode(opt_mode)?;

    pgr::libs::pairmat::write_phylip_matrix(&matrix, fmt, &mut writer)?;

    Ok(())
}
