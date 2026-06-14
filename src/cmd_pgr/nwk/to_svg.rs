use clap::*;
use pgr::libs::phylo::tree::Tree;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("to-svg")
        .about("Convert Newick trees to SVG format")
        .after_help(
            r###"
Convert Newick trees to SVG format for visualization.

Notes:
* Draws a cladogram by default
* Set `--bl` to draw a phylogenetic tree with scaled branch lengths
* Underscore `_` in names will be replaced as spaces " "
* Default styles match the LaTeX Forest template (grey branches, black dots)
* Scale bar is drawn in phylogram mode

Examples:
1. Convert to SVG (cladogram):
   pgr nwk to-svg tests/newick/catarrhini.nwk -o tree.svg

2. Convert with branch lengths (phylogram):
   pgr nwk to-svg tests/newick/catarrhini.nwk --bl -o tree.svg

3. Custom width and spacing:
   pgr nwk to-svg tests/newick/catarrhini.nwk -w 1200 -v 30 -o tree.svg
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
                .help("Draw a phylogram with scaled branch lengths instead of a cladogram"),
        )
        .arg(
            Arg::new("width")
                .short('w')
                .long("width")
                .num_args(1)
                .default_value("800")
                .help("SVG width in pixels"),
        )
        .arg(
            Arg::new("vskip")
                .short('v')
                .long("vskip")
                .num_args(1)
                .default_value("20")
                .help("Vertical spacing between leaf nodes in pixels"),
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
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());
    let is_bl = args.get_flag("bl");

    let width: f64 = args
        .get_one::<String>("width")
        .unwrap()
        .parse()
        .unwrap_or(800.0);
    let vskip: f64 = args
        .get_one::<String>("vskip")
        .unwrap()
        .parse()
        .unwrap_or(20.0);

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

    let out_string = pgr::libs::phylo::tree::io::to_svg(&tree, height, vskip, width);

    writer.write_all(out_string.as_ref())?;

    Ok(())
}
