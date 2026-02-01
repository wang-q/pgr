use clap::*;
use pgr::libs::axt::{AxtReader, write_axt};

// Create subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("sort")
        .about("Sort axt files")
        .arg(
            Arg::new("input")
                .index(1)
                .default_value("stdin")
                .help("Input axt file (or stdin if not specified)")
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Output axt file (or stdout if not specified)")
                .num_args(1)
                .default_value("stdout")
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
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = args.get_one::<String>("input").unwrap();
    let output = args.get_one::<String>("output").unwrap();
    let by_query = args.get_flag("query");
    let by_score = args.get_flag("by-score");

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
        axts.sort_by(|a, b| {
            a.q_name.cmp(&b.q_name)
                .then(a.q_start.cmp(&b.q_start))
        });
    } else {
        // Sort by target (default)
        axts.sort_by(|a, b| {
            a.t_name.cmp(&b.t_name)
                .then(a.t_start.cmp(&b.t_start))
        });
    }

    for axt in axts {
        write_axt(&mut writer, &axt)?;
    }

    Ok(())
}
