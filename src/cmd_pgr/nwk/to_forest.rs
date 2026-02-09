use clap::*;
use pgr::libs::phylo::tree::io::to_forest;
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("to-forest")
        .about("Convert Newick trees to raw LaTeX Forest code")
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
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input filename. [stdin] for standard input"),
        )
        .arg(
            Arg::new("bl")
                .long("bl")
                .action(ArgAction::SetTrue)
                .help("With branch lengths"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
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
