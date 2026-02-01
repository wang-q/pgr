use clap::*;
use pgr::libs::axt::{AxtReader, write_axt};
use std::io::BufWriter;

// Create subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("sort")
        .about("Sort axt files")
        .arg(
            Arg::new("input")
                .help("Input axt file")
                .index(1)
                .required(true),
        )
        .arg(
            Arg::new("output")
                .help("Output axt file")
                .index(2)
                .required(true),
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

    let reader = std::fs::File::open(input)?;
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

    let writer = std::fs::File::create(output)?;
    let mut buf_writer = BufWriter::new(writer);
    
    for axt in axts {
        write_axt(&mut buf_writer, &axt)?;
    }

    Ok(())
}
