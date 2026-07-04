use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::fmt::axt::{write_axt, AxtReader};

/// Build the clap subcommand for sort.
pub fn make_subcommand() -> Command {
    Command::new("sort")
        .about("Sorts axt files")
        .after_help(
            r###"
Sorts axt files by target, query, or score.

Examples:
  # Sort by target (default)
  pgr axt sort in.axt -o out.axt

  # Sort by query
  pgr axt sort in.axt --by-query -o out.axt

  # Sort by score (descending)
  pgr axt sort in.axt --by-score -o out.axt
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input AXT file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(crate::cmd_pgr::args::by_query_arg(
            "Sort by query position, not target",
        ))
        .arg(
            Arg::new("by_score")
                .long("by-score")
                .action(ArgAction::SetTrue)
                .help("Sort by score"),
        )
        .arg(
            Arg::new("renumber")
                .long("renumber")
                .action(ArgAction::SetTrue)
                .help("Renumber AXT records"),
        )
}
/// Execute the sort command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let by_query = args.get_flag("by_query");
    let by_score = args.get_flag("by_score");
    let renumber = args.get_flag("renumber");

    let reader = pgr::reader(input)?;
    let mut writer = pgr::writer(output)?;

    let mut axt_reader = AxtReader::new(reader);

    let mut axts = Vec::new();
    for result in axt_reader.by_ref() {
        axts.push(result?);
    }

    for header in &axt_reader.headers {
        writeln!(writer, "{}", header)?;
    }

    let by = if by_score {
        pgr::libs::fmt::axt::AxtSortBy::Score
    } else if by_query {
        pgr::libs::fmt::axt::AxtSortBy::Query
    } else {
        pgr::libs::fmt::axt::AxtSortBy::Target
    };

    pgr::libs::fmt::axt::sort_axts(&mut axts, by, renumber);

    for axt in &axts {
        write_axt(&mut writer, axt)?;
    }

    Ok(())
}
