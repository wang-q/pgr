use clap::*;
use pgr::libs::phylo::reader;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("indent")
        .about("Indent the Newick file")
        .after_help(
            r#"
* Set `--text` to something other than whitespaces will result in an invalid Newick file
    * Use `--text ".   "` can produce visual guide lines

* Set `--text` to empty ("") will remove indentation

"#,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input filename. [stdin] for standard input"),
        )
        .arg(
            Arg::new("text")
                .long("text")
                .short('t')
                .num_args(1)
                .default_value("  ")
                .help("Use this text instead of the default two spaces"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
        .arg(
            Arg::new("compact")
                .long("compact")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Compact output (remove indentation)"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    let compact = args.get_flag("compact");
    let text = if compact {
        ""
    } else {
        args.get_one::<String>("text").unwrap()
    };

    let infile = args.get_one::<String>("infile").unwrap();
    let trees = reader::from_file(infile);

    for tree in trees {
        let out_string = tree.to_newick_with_format(text);
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}

