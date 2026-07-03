use clap::*;
use pgr::libs::fmt::axt::{write_axt, AxtReader};

// Create subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("sort")
        .about("Sort axt files")
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
        .arg(
            Arg::new("by_query")
                .long("by-query")
                .action(ArgAction::SetTrue)
                .help("Sort by query position, not target"),
        )
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
        let axt = result?;
        axts.push(axt);
    }

    for header in &axt_reader.headers {
        writeln!(writer, "{}", header)?;
    }

    if by_score {
        // Sort by score. Assuming higher score is better (descending).
        axts.sort_by_key(|b| std::cmp::Reverse(b.score.unwrap_or(0)));
    } else if by_query {
        axts.sort_by(|a, b| a.q_name.cmp(&b.q_name).then(a.q_start.cmp(&b.q_start)));
    } else {
        // Sort by target (default)
        axts.sort_by(|a, b| a.t_name.cmp(&b.t_name).then(a.t_start.cmp(&b.t_start)));
    }

    for (i, axt) in axts.iter_mut().enumerate() {
        if renumber {
            axt.id = i as u64;
        }
        write_axt(&mut writer, axt)?;
    }

    Ok(())
}
