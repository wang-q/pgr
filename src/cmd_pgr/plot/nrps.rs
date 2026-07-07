use anyhow::{anyhow, Context};
use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::plot::nrps::parse_nrps;
use std::io::Read;

/// Build the clap subcommand for nrps.
pub fn make_subcommand() -> Command {
    Command::new("nrps")
        .about("Plots an NRPS structure diagram")
        .after_help(
            r###"
* Input file is a tab-separated file
    * First column: Domain type (A, C, E, CE, T, Te, R, M)
    * Second column: Text (optional, amino acid/name)
    * Third column: Color (optional)

* Colors
    * black: RGB(26,25,25)
    * grey: RGB(129,130,132)
    * red: RGB(188,36,46)
    * brown: RGB(121,37,0)
    * green: RGB(32,128,108)
    * purple: RGB(160,90,150)
    * blue: RGB(0,103,149)

* Examples
    pgr plot nrps input.tsv -o output.tex

    pgr plot nrps input.tsv |
        tectonic - &&
        mv texput.pdf nrps.pdf

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required())
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("legend")
                .long("legend")
                .action(ArgAction::SetTrue)
                .help("Include legend in the output"),
        )
        .arg(crate::cmd_pgr::args::color_arg(
            Some("grey"),
            "Default color for modules",
        ))
}

/// Execute the nrps command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let default_color = args.get_one::<String>("color").unwrap().clone();
    let is_legend = args.get_flag("legend");

    // Read TSV file
    let mut content = String::new();
    pgr::reader(infile)
        .with_context(|| format!("Failed to open reader for {}", infile))?
        .read_to_string(&mut content)
        .with_context(|| format!("Failed to read {}", infile))?;
    let nrps_data = parse_nrps(&content, &default_color)?;

    // Generate all modules
    let mut all_tex = String::new();
    for (module_name, domains) in &nrps_data.modules {
        let info = nrps_data
            .module_info
            .get(module_name)
            .ok_or_else(|| anyhow!("missing module info: {}", module_name))?;
        let module_tex = pgr::libs::plot::nrps::gen_module(info, domains)?;
        all_tex.push_str(&module_tex);
        all_tex.push('\n');
    }

    // Context
    let mut context = tera::Context::new();

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    context.insert("outfile", outfile);
    context.insert("all_tex", &all_tex);
    context.insert("is_legend", &is_legend);
    context.insert("default_color", &default_color);

    pgr::libs::plot::nrps::gen_nrps(&context)?;

    Ok(())
}
