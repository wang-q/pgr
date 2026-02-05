use clap::*;
use pgr::libs::axt::{write_axt, AxtReader};

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
  pgr axt sort in.axt --query -o out.axt

  # Sort by score (descending)
  pgr axt sort in.axt --by-score -o out.axt
"###,
        )
        .arg(
            Arg::new("input")
                .index(1)
                .default_value("stdin")
                .help("Input AXT file"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Output AXT file")
                .num_args(1)
                .value_name("FILE")
                .default_value("stdout"),
        )
        .arg(
            Arg::new("query")
                .long("query")
                .action(ArgAction::SetTrue)
                .help("Sort by query position, not target"),
        )
        .arg(
            Arg::new("by-score")
                .long("by-score")
                .action(ArgAction::SetTrue)
                .help("Sort by score"),
        )
        .arg(
            Arg::new("renumber")
                .long("renumber")
                .short('r')
                .action(ArgAction::SetTrue)
                .help("Renumber AXT records"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = args.get_one::<String>("input").unwrap();
    let output = args.get_one::<String>("output").unwrap();
    let by_query = args.get_flag("query");
    let by_score = args.get_flag("by-score");
    let renumber = args.get_flag("renumber");

    let reader = intspan::reader(input);
    let mut writer = intspan::writer(output);

    let axt_reader = AxtReader::new(reader);

    let mut axts = Vec::new();
    for result in axt_reader {
        let axt = result?;
        axts.push(axt);
    }

    if by_score {
        // Sort by score. Assuming higher score is better (descending).
        axts.sort_by(|a, b| b.score.unwrap_or(0).cmp(&a.score.unwrap_or(0)));
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
