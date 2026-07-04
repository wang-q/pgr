use clap::{ArgMatches, Command};
use pgr::libs::phylo::tree::io::to_forest;
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

/// Build the clap subcommand for to-forest.
pub fn make_subcommand() -> Command {
    Command::new("to-forest")
        .about("Converts Newick trees to raw LaTeX Forest code")
        .after_help(
            r###"
Convert Newick trees to raw LaTeX Forest code.

This command is designed for manually modifying the generated Forest code.

Notes:
* Styles are stored in the comments of each node
* Drawing a cladogram by default
* Set `--bl` to draw a phylogenetic tree

Examples:
1. Convert to Forest code:
   pgr nwk to-forest tests/newick/catarrhini.nwk

2. Convert with branch lengths:
   pgr nwk to-forest tests/newick/catarrhini.nwk --bl

3. Save to file:
   pgr nwk to-forest tests/newick/catarrhini.nwk -o forest.tex
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::bl_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the to-forest command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let is_bl = args.get_flag("bl");

    let infile = args.get_one::<String>("infile").unwrap();

    let tree = Tree::from_file(infile)?
        .into_iter()
        .next()
        .unwrap_or(Tree::new());

    let height = if is_bl {
        tree.get_root()
            .map(|r| tree.get_height(r, true))
            .unwrap_or(0.0)
    } else {
        0.0
    };

    let out_string = to_forest(&tree, height);

    writer.write_all((out_string + "\n").as_ref())?;

    Ok(())
}
