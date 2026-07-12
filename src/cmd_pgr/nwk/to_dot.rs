use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

/// Build the clap subcommand for to-dot.
pub fn make_subcommand() -> Command {
    Command::new("to-dot")
        .about("Converts Newick trees to Graphviz DOT format")
        .after_help(
            r###"
Convert Newick trees to Graphviz DOT format for visualization.

Examples:
1. Convert to DOT:
   pgr nwk to-dot tests/newick/catarrhini.nwk

2. Save to file:
   pgr nwk to-dot tests/newick/catarrhini.nwk -o tree.dot

3. Create an image (requires Graphviz installed):
   pgr nwk to-dot tests/newick/catarrhini.nwk | dot -Tpng -o tree.png
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the to-dot command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    let infile = args
        .get_one::<String>("infile")
        .ok_or_else(|| anyhow::anyhow!("missing required argument: infile"))?;
    let trees = Tree::from_file(infile)?;

    for tree in trees {
        let out_string = tree.to_dot();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    writer.flush()?;
    Ok(())
}
