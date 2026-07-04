use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

/// Build the clap subcommand for to-svg.
pub fn make_subcommand() -> Command {
    Command::new("to-svg")
        .about("Converts Newick trees to SVG format")
        .after_help(
            r###"
Convert Newick trees to SVG format for visualization.

Notes:
* Automatically draws a phylogram if branch lengths are present, otherwise a cladogram
* Underscore `_` in names will be replaced as spaces " "
* Default styles match the LaTeX Forest template (grey branches, black dots)
* Scale bar is drawn in phylogram mode

Examples:
1. Convert to SVG:
   pgr nwk to-svg tests/newick/catarrhini.nwk -o tree.svg

2. Custom width and spacing:
   pgr nwk to-svg tests/newick/catarrhini.nwk -w 1200 -v 30 -o tree.svg
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(
            Arg::new("width")
                .short('w')
                .long("width")
                .num_args(1)
                .default_value("800")
                .value_parser(clap::value_parser!(f64))
                .help("SVG width in pixels"),
        )
        .arg(
            Arg::new("vskip")
                .short('v')
                .long("vskip")
                .num_args(1)
                .default_value("20")
                .value_parser(clap::value_parser!(f64))
                .help("Vertical spacing between leaf nodes in pixels"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the to-svg command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    let width: f64 = *args.get_one::<f64>("width").unwrap();
    let vskip: f64 = *args.get_one::<f64>("vskip").unwrap();

    let infile = args.get_one::<String>("infile").unwrap();

    let tree = Tree::from_file(infile)?
        .into_iter()
        .next()
        .ok_or_else(|| anyhow::anyhow!("no trees found in {}", infile))?;

    // Auto-detect: if any node has a branch length, draw phylogram
    let has_bl = if let Some(root) = tree.get_root() {
        let ids = tree
            .preorder(&root)
            .map_err(|e| anyhow::anyhow!("preorder traversal failed: {}", e))?;
        ids.iter().any(|&id| {
            tree.get_node(id)
                .map(|n| n.length.is_some())
                .unwrap_or(false)
        })
    } else {
        false
    };
    let height = if has_bl {
        tree.get_root()
            .map(|r| tree.get_height(r, true))
            .unwrap_or(0.0)
    } else {
        0.0
    };

    let out_string = pgr::libs::phylo::tree::io::to_svg(&tree, height, vskip, width);

    writer.write_all(out_string.as_ref())?;

    Ok(())
}
