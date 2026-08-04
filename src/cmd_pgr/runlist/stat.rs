//! `pgr runlist stat` — coverage stats of a runlist against chromosome sizes.

use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for stat.
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Coverage on chromosomes for runlists")
        .after_help(
            r###"
Prints per-chromosome coverage as TSV (`key\tchr\tchrLength\tsize\tcoverage`
plus an `all` row). `--all` keeps only the whole-genome stats.

Examples:
1. Per-chromosome stats:
   pgr runlist stat chr.sizes in.json -o stat.tsv
"###,
        )
        .arg(
            Arg::new("chr.sizes")
                .required(true)
                .index(1)
                .help("Sets the input file to use"),
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(2)
                .help("Sets the input file to use"),
        )
        .arg(
            Arg::new("all")
                .long("all")
                .action(ArgAction::SetTrue)
                .help("Only write whole genome stats"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the stat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let sizes = pgr::read_sizes::<i32>(args.get_one::<String>("chr.sizes").unwrap())?;
    let json = pgr::libs::runlist::read_json(args.get_one::<String>("infile").unwrap())?;
    let set_of = pgr::libs::runlist::json_to_sets(&json)?;
    let is_multi = set_of.len() > 1 || !set_of.contains_key("__single__");
    let is_all = args.get_flag("all");

    let mut lines: Vec<String> = Vec::new();
    let mut header = "key\tchr\tchrLength\tsize\tcoverage".to_string();
    if is_multi {
        if is_all {
            header = header.replacen("chr\t", "", 1);
        }
        lines.push(header);
        for (name, set) in &set_of {
            lines.push(pgr::libs::runlist::stat_lines(
                set,
                &sizes,
                is_all,
                Some(name),
            )?);
        }
    } else {
        header = header.replacen("key\t", "", 1);
        if is_all {
            header = header.replacen("chr\t", "", 1);
        }
        lines.push(header);
        lines.push(pgr::libs::runlist::stat_lines(
            &set_of["__single__"],
            &sizes,
            is_all,
            None,
        )?);
    }
    let mut w = pgr::writer(outfile)?;
    for line in &lines {
        writeln!(w, "{}", line)?;
    }
    w.flush()?;
    Ok(())
}
