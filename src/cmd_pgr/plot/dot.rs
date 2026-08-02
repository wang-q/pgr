use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::paf::parser::parse_paf;
use pgr::libs::plot::dot::{render_dot_svg, DotOpts, PlotRange};
use std::io::Write;

/// Build the clap subcommand for dot.
pub fn make_subcommand() -> Command {
    Command::new("dot")
        .about("Plots a dot plot from PAF alignments")
        .after_help(
            r###"
This command draws a static collinear plot (dot plot) of PAF alignments.

* Input is a PAF file, .paf.gz is supported. Reads from stdin if input file is 'stdin'.

* The two axes are the target and query sequences, laid out contig by contig
  in first-appearance order.
* Each alignment is drawn as a line segment colored by identity on a blue
  scale: --min-identity (lightest) to --max-identity (deepest).
* Alignments below --min-identity, shorter than --min-len, or less identical
  than --min-identity are skipped; at most --max-align alignments are drawn
  (longest first).
* With --range, only alignments overlapping the target-side region are drawn;
  the query axis auto-focuses on the aligned regions (local zoom-in).

* Output is an SVG file, rendered with no external dependencies. Convert to
  PDF/PNG with rsvg-convert, inkscape, or cairosvg if needed.

* Examples
    pgr plot dot input.paf -o dot.svg

    pgr plot dot input.paf.gz | rsvg-convert -f pdf -o dot.pdf

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("min_len")
                .long("min-len")
                .num_args(1)
                .value_parser(value_parser!(u32))
                .default_value("100")
                .help("Minimum alignment block length"),
        )
        .arg(
            Arg::new("min_identity")
                .long("min-identity")
                .num_args(1)
                .value_parser(value_parser!(f64))
                .default_value("0.7")
                .help("Minimum identity to plot; lightest color of the scale"),
        )
        .arg(
            Arg::new("identity_max")
                .long("max-identity")
                .num_args(1)
                .value_parser(value_parser!(f64))
                .default_value("1.0")
                .help("Identity at which the color scale saturates"),
        )
        .arg(
            Arg::new("max_align")
                .long("max-align")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .default_value("100000")
                .help("Maximum number of alignments to plot (0 = all)"),
        )
        .arg(
            Arg::new("width")
                .long("width")
                .num_args(1)
                .value_parser(value_parser!(u32))
                .default_value("1200")
                .help("Plot width in pixels; height is scaled automatically"),
        )
        .arg(
            Arg::new("range")
                .long("range")
                .num_args(1)
                .help("Target-side region to zoom into (chr:start-end, 1-based)"),
        )
}

/// Execute the dot command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let opts = DotOpts {
        min_len: *args.get_one::<u32>("min_len").unwrap(),
        min_identity: *args.get_one::<f64>("min_identity").unwrap(),
        identity_max: *args.get_one::<f64>("identity_max").unwrap(),
        max_align: *args.get_one::<usize>("max_align").unwrap(),
        width: *args.get_one::<u32>("width").unwrap(),
        range: args
            .get_one::<String>("range")
            .map(|s| s.parse::<PlotRange>())
            .transpose()
            .context("failed to parse --range")?,
    };

    let reader = pgr::reader(infile)?;
    let records = parse_paf(reader)?;
    let svg = render_dot_svg(&records, &opts)?;

    let mut writer = pgr::writer(outfile)?;
    writer
        .write_all(svg.as_bytes())
        .context("failed to write SVG output")?;

    Ok(())
}
