use clap::*;
use pgr::libs::phylo::tree::io::to_forest;
use pgr::libs::phylo::tree::Tree;
use std::io::Read;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("to-tex")
        .about("Convert Newick trees to a full LaTeX document")
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
            Arg::new("forest")
                .long("forest")
                .action(ArgAction::SetTrue)
                .help("Treat input as a file containing pre-generated Forest code (pass-through mode)"),
        )
        .arg(
            Arg::new("style")
                .long("style")
                .short('s')
                .action(ArgAction::SetTrue)
                .help("Skip default font settings in the template to allow custom styles"),
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
    let is_style = args.get_flag("style");

    let infile = args.get_one::<String>("infile").unwrap();

    let out_string = if args.get_flag("forest") {
        let mut reader = pgr::reader(infile);
        let mut s = String::new();
        reader.read_to_string(&mut s).expect("Read error");

        s
    } else {
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
        // eprintln!("height = {:#?}", tree.height().unwrap());
        // eprintln!("diameter = {:#?}", tree.diameter().unwrap());

        let mut s = to_forest(&tree, height);

        // a bar of unit length
        if is_bl {
            // Determine scale dynamically
            let target_scale = height / 5.0; // Target ~20% of tree height
            let magnitude = target_scale.log10().floor();
            let base = 10.0_f64.powf(magnitude);

            // Find the best step: 1x, 2x, or 5x of the base magnitude
            // e.g., if target is 0.035 (base 0.01), candidates are 0.01, 0.02, 0.05
            // we want the largest one <= target.
            let scale = [1.0, 2.0, 5.0]
                .iter()
                .map(|&x| base * x)
                .filter(|&x| x <= target_scale)
                .last()
                .unwrap_or(base); // Fallback to base (1x) if even 1x is too big? Should not happen if logic is sound.

            // Calculate actual length in millimeters
            let bar_mm = (scale * 100.0 / height).round() as i32;
            
            // If the bar is too small to see (< 5mm), don't draw it or warn?
            // Current logic just draws it.

            // Draw scale bar
            s += "\\draw[-, grey, line width=1pt]";
            s += " ($(current bounding box.south east)+(-10mm,-2mm)$)";
            s += &format!(" --++ (-{}mm,0mm)", bar_mm);
            s += &format!(" node[midway, below]{{\\scriptsize{{{}}}}};\n", scale);
        }

        s
    };

    static FILE_TEMPLATE: &str = include_str!("../../../docs/template.tex");
    let mut template = FILE_TEMPLATE.to_string();

    {
        // Section forest
        let begin = template.find("%FOREST_BEGIN").unwrap();
        let end = template.find("%FOREST_END").unwrap();
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
        let begin = template.find("%STYLE_BEGIN").unwrap();
        let end = template.find("%STYLE_END").unwrap();
        template.replace_range(begin..end, default_font);
    }

    writer.write_all(template.as_ref())?;

    Ok(())
}
