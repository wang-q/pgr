use clap::*;
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("to-dot")
        .about("Convert Newick trees to Graphviz DOT format")
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
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input filename. [stdin] for standard input"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let infile = args.get_one::<String>("infile").unwrap();
    let trees = Tree::from_file(infile)?;

    for tree in trees {
        let out_string = tree.to_dot();
        writer.write_all((out_string + "\n").as_ref())?;
    }

    Ok(())
}
