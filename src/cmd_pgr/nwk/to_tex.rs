use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::phylo::tree::io::{compute_scale_bar, to_forest};
use pgr::libs::phylo::tree::Tree;
use std::io::{Read, Write};

/// Build the clap subcommand for to-tex.
pub fn make_subcommand() -> Command {
    Command::new("to-tex")
        .about("Converts Newick trees to a full LaTeX document")
        .after_help(
            r###"
Convert Newick trees to a full LaTeX document.

Notes:
* Styles are stored in the comments of each node
* Drawing a cladogram by default
* Set `--bl` to draw a phylogenetic tree
* Underscore `_` is a control character in LaTeX
  * All `_`s in names, labels and comments will be replaced as spaces " "
* To compile the .tex files to pdf, you need LaTeX or Tectonic
  * `XeLaTeX` and `latexmk` for compiling unicode .tex
  * `xeCJK` package for East Asian characters
  * `Forest` package for drawing trees

Examples:
1. Generate LaTeX file:
   pgr nwk to-tex tests/newick/catarrhini.nwk -o tree.tex

2. Compile with Tectonic (recommended):
   tectonic tree.tex

3. Compile with Latexmk:
   latexmk -xelatex tree.tex
   latexmk -c tree.tex
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::bl_arg())
        .arg(
            Arg::new("forest")
                .long("forest")
                .action(ArgAction::SetTrue)
                .help("Treat input as a file containing pre-generated Forest code (pass-through mode)"),
        )
        .arg(
            Arg::new("no_default_style")
                .long("no-default-style")
                .action(ArgAction::SetTrue)
                .help("Skip default font settings in the template to allow custom styles"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the to-tex command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let is_bl = args.get_flag("bl");
    let is_style = args.get_flag("no_default_style");

    let infile = args.get_one::<String>("infile").unwrap();

    let out_string = if args.get_flag("forest") {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
        let mut s = String::new();
        reader
            .read_to_string(&mut s)
            .with_context(|| format!("Failed to read from {}", infile))?;

        s
    } else {
        let tree = Tree::from_file(infile)?
            .into_iter()
            .next()
            .ok_or_else(|| anyhow::anyhow!("no trees found in {}", infile))?;

        let height = if is_bl {
            tree.get_root()
                .map(|r| tree.get_height(r, true))
                .unwrap_or(0.0)
        } else {
            0.0
        };
        let mut s = to_forest(&tree, height);

        // a bar of unit length
        if is_bl && height > 0.0 {
            let (scale, bar_mm) = compute_scale_bar(height);

            // Draw scale bar
            s += "\\draw[-, grey, line width=1pt]";
            s += " ($(current bounding box.south east)+(-10mm,-2mm)$)";
            s += &format!(" --++ (-{}mm,0mm)", bar_mm);
            s += &format!(" node[midway, below]{{\\scriptsize{{{}}}}};\n", scale);
        }

        s
    };

    static FILE_TEMPLATE: &str = include_str!("../../assets/template.tex");
    let mut template = FILE_TEMPLATE.to_string();

    {
        // Section forest
        let begin = template
            .find("%FOREST_BEGIN")
            .ok_or_else(|| anyhow::anyhow!("template marker %FOREST_BEGIN missing"))?;
        let end = template
            .find("%FOREST_END")
            .ok_or_else(|| anyhow::anyhow!("template marker %FOREST_END missing"))?;
        anyhow::ensure!(begin < end, "template markers %FOREST out of order");
        template.replace_range(begin..end, &out_string);
    }

    if !is_style {
        let default_font = r#"\setmainfont{NotoSans}[
    Extension      = .ttf,
    UprightFont    = *-Regular,
    BoldFont       = *-Bold,
    ItalicFont     = *-Italic,
    BoldItalicFont = *-BoldItalic
]
"#;

        // Section style
        let begin = template
            .find("%STYLE_BEGIN")
            .ok_or_else(|| anyhow::anyhow!("template marker %STYLE_BEGIN missing"))?;
        let end = template
            .find("%STYLE_END")
            .ok_or_else(|| anyhow::anyhow!("template marker %STYLE_END missing"))?;
        anyhow::ensure!(begin < end, "template markers %STYLE out of order");
        template.replace_range(begin..end, default_font);
    }

    writer.write_all(template.as_ref())?;

    writer.flush()?;
    Ok(())
}
