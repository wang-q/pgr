use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

/// Build the clap subcommand for indent.
pub fn make_subcommand() -> Command {
    Command::new("indent")
        .about("Formats Newick trees with indentation")
        .after_help(
            r###"
Re-formats the Newick tree, making structure more clear.

Notes:
* By default, prints the input tree indented with two spaces ("  ").
* The default output is valid Newick.
* Use `--compact` to remove all indentation (output single line).
* Using non-whitespace characters for `--text` may result in invalid Newick.

Examples:
1. Default indentation:
   pgr nwk indent tests/newick/catarrhini.nwk

2. Compact output (remove indentation):
   pgr nwk indent tests/newick/catarrhini.nwk --compact

3. Indent with visual guides (NOT valid Newick):
   pgr nwk indent tests/newick/catarrhini.nwk --text ".   "
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(
            Arg::new("text")
                .long("text")
                .num_args(1)
                .default_value("  ")
                .help("Use this text instead of the default two spaces"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("compact")
                .long("compact")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Compact output (remove indentation)"),
        )
}

/// Execute the indent command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    let compact = args.get_flag("compact");
    let text = if compact {
        ""
    } else {
        args.get_one::<String>("text").unwrap()
    };

    let infile = args.get_one::<String>("infile").unwrap();
    let trees = Tree::from_file(infile)?;

    for tree in trees {
        let out_string = tree.to_newick_with_format(text);
        writer.write_all((out_string + "\n").as_ref())?;
    }

    writer.flush()?;
    Ok(())
}
