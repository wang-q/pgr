use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::plot::histogram::{
    calc_density, calc_hist, compute_hh_axis, create_table, load_data, render_hh_tex,
};

/// Build the clap subcommand for hh.
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
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("xlabel")
                .long("xlabel")
                .num_args(1)
                .help("X label"),
        )
        .arg(
            Arg::new("ylabel")
                .long("ylabel")
                .num_args(1)
                .help("Y label"),
        )
        .arg(
            Arg::new("column")
                .long("column")
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
        .arg(
            Arg::new("xmin_max")
                .long("xmin-max")
                .num_args(1)
                .help("X min,max"),
        )
        .arg(
            Arg::new("unit")
                .long("unit")
                .num_args(1)
                .default_value("0.5,1.5")
                .help("Cell width,height"),
        )
}

/// Execute the hh command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();

    // Optional arguments with defaults
    let col = args.get_one::<usize>("column").unwrap();
    let group = args.get_one::<usize>("group");
    let bins = args.get_one::<usize>("bins").unwrap();

    // Optional labels
    let xlabel = args.get_one::<String>("xlabel").map(|s| s.to_string());
    let ylabel = args.get_one::<String>("ylabel").map(|s| s.to_string());

    // Parse X min,max if provided
    let xmm = if let Some(s) = args.get_one::<String>("xmin_max") {
        let parts: Vec<f64> = s
            .split(',')
            .map(|x| x.trim().parse::<f64>())
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| anyhow::anyhow!("invalid --xmin_max value: {}", e))?;
        if parts.len() != 2 {
            anyhow::bail!("--xmin_max must be two comma-separated floats, got: {}", s);
        }
        Some((parts[0], parts[1]))
    } else {
        None
    };

    // Parse unit
    let unit_str = args
        .get_one::<String>("unit")
        .ok_or_else(|| anyhow::anyhow!("--unit required"))?;
    let parts: Vec<f64> = unit_str
        .split(',')
        .map(|x| x.trim().parse::<f64>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| anyhow::anyhow!("invalid --unit value: {}", e))?;
    let unit: (f64, f64) = if parts.len() == 2 {
        (parts[0], parts[1])
    } else {
        anyhow::bail!(
            "--unit must be two comma-separated floats, got: {}",
            unit_str
        );
    };

    // Table section
    let (data, col_name, group_name) = load_data(infile, *col, group)?;

    // Calculate histogram for each group
    let (hist_data, bin_edges) = calc_hist(&data, *bins, xmm)?;
    let density_data = calc_density(&hist_data);
    let table = create_table(&density_data);

    // Axis section
    let xlabel = xlabel.unwrap_or(col_name);
    let ylabel = ylabel.unwrap_or(group_name);

    let axis = compute_hh_axis(&density_data, *bins, &bin_edges, unit);

    // Context
    let mut context = tera::Context::new();

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer = pgr::writer(outfile)?;
    context.insert("table", &table);
    context.insert("xlabel", &xlabel);
    context.insert("ylabel", &ylabel);
    context.insert("width", &axis.width);
    context.insert("height", &axis.height);
    context.insert("xticks", &axis.xticks);
    context.insert("xtick_labels", &axis.xtick_labels);
    context.insert("ygroups", &axis.ygroups);
    context.insert("yticks", &axis.yticks);
    context.insert("label_len", &axis.label_len);

    render_hh_tex(&context, &mut writer)?;

    Ok(())
}
