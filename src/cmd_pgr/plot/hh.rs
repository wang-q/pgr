use crate::cmd_pgr::plot::common::{context_get_str, render_and_write, replace_section};
use clap::*;
use indexmap::IndexMap;
use pgr::libs::plot::histogram::{calc_density, calc_hist, create_table};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("hh")
        .about("Histo-heatmap")
        .after_help(
            r###"
* Input file is a tab-separated file
    * First column: X values
    * Second column: Group names (optional)
    * Header line is required

* The output will be a LaTeX file containing a heatmap
    * Showing the distribution of X values across groups
    * Using colors from white to red
    * With a color bar showing the scale
    * For single group data, a dummy group will be added to meet the matrix
      plot requirements

* To convert the .tex files to pdf
    * Install tectonic (https://tectonic-typesetting.github.io)
    * It will automatically handle all required LaTeX packages

* Examples
    pgr plot hh input.tsv -o output.tex

    pgr plot hh input.tsv  |
        tectonic - &&
        mv texput.pdf hh.pdf

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
        .arg(Arg::new("xl").long("xl").num_args(1).help("X label"))
        .arg(Arg::new("yl").long("yl").num_args(1).help("Y label"))
        .arg(
            Arg::new("col")
                .long("col")
                .short('c')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .default_value("1")
                .help("Which column to count"),
        )
        .arg(
            Arg::new("group")
                .long("group")
                .short('g')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("The group column"),
        )
        .arg(
            Arg::new("bins")
                .long("bins")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .default_value("40")
                .help("Number of bins"),
        )
        .arg(Arg::new("xmm").long("xmm").num_args(1).help("X min,max"))
        .arg(
            Arg::new("unit")
                .long("unit")
                .num_args(1)
                .default_value("0.5,1.5")
                .help("Cell width,height"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();

    // Optional arguments with defaults
    let col = args.get_one::<usize>("col").unwrap();
    let group = args.get_one::<usize>("group");
    let bins = args.get_one::<usize>("bins").unwrap();

    // Optional labels
    let xlabel = args.get_one::<String>("xl").map(|s| s.to_string());
    let ylabel = args.get_one::<String>("yl").map(|s| s.to_string());

    // Parse X min,max if provided
    let xmm = args.get_one::<String>("xmm").and_then(|s| {
        let parts: Vec<f64> = s.split(',').filter_map(|x| x.trim().parse().ok()).collect();
        if parts.len() == 2 {
            Some((parts[0], parts[1]))
        } else {
            None
        }
    });

    // Parse unit
    let unit = args
        .get_one::<String>("unit")
        .map(|s| {
            let parts: Vec<f64> = s.split(',').filter_map(|x| x.trim().parse().ok()).collect();
            (parts[0], parts[1])
        })
        .unwrap();

    //----------------------------
    // Table section
    //----------------------------
    let (data, col_name, group_name) = load_data(infile, *col, group)?;

    // Calculate histogram for each group
    let (hist_data, bin_edges) = calc_hist(&data, *bins, xmm)?;
    let density_data = calc_density(&hist_data);
    let table = create_table(&density_data);

    //----------------------------
    // Axis section
    //----------------------------
    // Use column names if labels are not specified
    let xlabel = xlabel.unwrap_or(col_name);
    let ylabel = ylabel.unwrap_or(group_name);

    // Width unit per bin
    let width = (*bins as f64) * unit.0;
    // Calculate height, 1 unit per group, minimum 2
    let height = (density_data.len() as f64).max(2.0) * unit.1;

    // Y
    let ygroups: Vec<_> = density_data.keys().cloned().collect();
    let yticks = (0..=density_data.len().max(2))
        .map(|i| i as f64 - 0.5)
        .collect::<Vec<_>>(); // Generate ticks from -0.5 to n-0.5

    // Find the longest group name for label width
    let label_len = ygroups.iter().map(|s| s.len()).max().unwrap_or(0).max(3);

    // X
    let xticks = (0..=*bins)
        .filter_map(|i| {
            if i % 5 == 0 {
                Some(i as f64 - 0.5)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    let xtick_labels = bin_edges
        .iter()
        .enumerate()
        .filter_map(|(i, &edge)| {
            if i % 5 == 0 {
                Some(format!("{}", edge))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    //----------------------------
    // Context
    //----------------------------
    let mut context = tera::Context::new();

    context.insert("outfile", args.get_one::<String>("outfile").unwrap());
    context.insert("table", &table);
    context.insert("xlabel", &xlabel);
    context.insert("ylabel", &ylabel);
    context.insert("width", &width);
    context.insert("height", &height);
    context.insert("xticks", &xticks);
    context.insert("xtick_labels", &xtick_labels);
    context.insert("ygroups", &ygroups);
    context.insert("yticks", &yticks);
    context.insert("label_len", &label_len);

    gen_hh(&context)?;

    Ok(())
}

fn load_data(
    infile: &str,
    col: usize,
    group: Option<&usize>,
) -> anyhow::Result<(IndexMap<String, Vec<f64>>, String, String)> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(infile)?;

    let headers = rdr.headers()?.clone();
    let mut data: IndexMap<String, Vec<f64>> = IndexMap::new();

    // Get column headers
    let xlabel = headers[col - 1].to_string();
    let ylabel = match group {
        Some(g) => headers[*g - 1].to_string(),
        None => String::new(),
    };

    for result in rdr.records() {
        let record = result?;

        if let Ok(val) = record[col - 1].parse::<f64>() {
            // Get group name, use "default" if group column not specified
            let group_name = match group {
                Some(g) => record[*g - 1].to_string(),
                None => "default".to_string(),
            };

            // Add value to corresponding group
            data.entry(group_name).or_default().push(val);
        }
    }

    Ok((data, xlabel, ylabel))
}

fn gen_hh(context: &tera::Context) -> anyhow::Result<()> {
    let outfile = context_get_str(context, "outfile")?;
    let mut writer = pgr::writer(outfile)?;

    static FILE_TEMPLATE: &str = include_str!("../../../docs/heatmap.tex");
    let mut template = FILE_TEMPLATE.to_string();

    let out_string = r###"%
width={{ width }}cm,
height={{ height }}cm,
xlabel={ {{ xlabel }} },
ylabel={ {{ ylabel }} },
extra x ticks={ {{ xticks | join(sep=", ") }} },
extra x tick labels={ {{ xtick_labels | join(sep=", ") }} },
yticklabels={ {{ ygroups | join(sep=", ") }} },
extra y ticks={ {{ yticks | join(sep=", ") }} },
y tick label style={
    text width={{ label_len }}ex,
},
    "###;
    replace_section(&mut template, "%AXIS_BEGIN", "%AXIS_END", out_string)?;

    let table = context_get_str(context, "table")?;
    replace_section(&mut template, "%TABLE_BEGIN", "%TABLE_END", table)?;

    render_and_write(&template, context, &mut writer)?;
    Ok(())
}
