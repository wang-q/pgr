use anyhow::{anyhow, Context, Result};
use clap::*;
use pgr::libs::plot::common::{context_get_str, render_and_write, replace_section};
use pgr::libs::plot::nrps::parse_nrps;
use std::collections::HashMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("nrps")
        .about("NRPS structure diagram")
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
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input filename. [stdin] for standard input"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("legend")
                .long("legend")
                .action(ArgAction::SetTrue)
                .help("Include legend in the output"),
        )
        .arg(
            Arg::new("color")
                .long("color")
                .short('c')
                .num_args(1)
                .default_value("grey")
                .help("Default color for modules"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args
        .get_one::<String>("infile")
        .ok_or_else(|| anyhow!("missing infile"))?;
    let default_color = args
        .get_one::<String>("color")
        .ok_or_else(|| anyhow!("missing color"))?
        .clone();
    let is_legend = args.get_flag("legend");

    //----------------------------
    // Read TSV file
    //----------------------------
    let content = std::fs::read_to_string(infile)?;
    let nrps_data = parse_nrps(&content, &default_color)?;

    // Generate all modules
    let mut all_tex = String::new();
    for (module_name, domains) in &nrps_data.modules {
        let info = nrps_data
            .module_info
            .get(module_name)
            .ok_or_else(|| anyhow!("missing module info: {}", module_name))?;
        let module_tex = gen_module(info, domains)?;
        all_tex.push_str(&module_tex);
        all_tex.push('\n');
    }

    //----------------------------
    // Context
    //----------------------------
    let mut context = tera::Context::new();

    let outfile = args
        .get_one::<String>("outfile")
        .ok_or_else(|| anyhow!("missing outfile"))?;
    context.insert("outfile", outfile);
    context.insert("all_tex", &all_tex);
    context.insert("is_legend", &is_legend);
    context.insert("default_color", &default_color);

    gen_nrps(&context)?;

    Ok(())
}

fn gen_module(
    info: &HashMap<String, String>,
    domains: &Vec<HashMap<String, String>>,
) -> Result<String> {
    let mut context = tera::Context::new();
    context.insert("info", info);
    context.insert("domains", domains);
    let last_domain = domains
        .last()
        .ok_or_else(|| anyhow!("empty domains in gen_module"))?;
    let first_domain = domains
        .first()
        .ok_or_else(|| anyhow!("empty domains in gen_module"))?;
    context.insert("last_domain", last_domain);
    context.insert("first_domain", first_domain);

    let template = r###"
    \begin{scope}[shift={([shift={({{ first_domain.dx_before }}cm,0)}]{{ info.prev }}.east)}]
{% for domain in domains -%}
        \node[{{ domain.type }}, {{ domain.color }}] ({{ domain.id }}) at ({{ domain.pos }}cm,0) {};
{% if domain.text != "" -%}
        \node[text=white,anchor=center,align=left] at ({{ domain.id }}) { {{ domain.text }}};
{% endif -%}
{% endfor -%}
        \begin{scope}[on background layer]
            \draw[{{ first_domain.color }}, line width=0.5mm, yshift=-1cm]
                let \p1 = ({{ first_domain.id }}), \p2 = ({{ last_domain.id }}) in
                (\x1,0) -- (\x2,0)
                node[midway, below, text={{ first_domain.color }}] { {{ info.id }} };
            \draw[{{ first_domain.color }}, line width=2mm]
                let \p1 = ({{ last_domain.id }}) in
                (-{{ first_domain.dx_before }}cm,0) -- (\x1 + {{ last_domain.dx_after }}cm,0)
                coordinate ({{ info.id }});
        \end{scope}
    \end{scope}"###;

    let mut tera = tera::Tera::default();
    tera.add_raw_templates(vec![("t", template)])
        .context("failed to register nrps module template")?;

    let rendered = tera
        .render("t", &context)
        .context("failed to render nrps module template")?;
    Ok(rendered)
}

fn gen_nrps(context: &tera::Context) -> Result<()> {
    let outfile = context_get_str(context, "outfile")?;
    let all_tex = context_get_str(context, "all_tex")?;
    let mut writer = pgr::writer(outfile)?;

    static FILE_TEMPLATE: &str = include_str!("../../../docs/nrps.tex");
    let mut template = FILE_TEMPLATE.to_string();

    // Section color
    let default_color = context_get_str(context, "default_color")?;
    let color_section = format!(
        r###"%
        draw={},
        fill={},
        text=white,
        "###,
        default_color, default_color
    );
    replace_section(
        &mut template,
        "%COLOR_BEGIN%",
        "%COLOR_END%",
        &color_section,
    )?;

    // Section module
    replace_section(&mut template, "%MODULE_BEGIN%", "%MODULE_END%", all_tex)?;

    // Section legend
    let is_legend = context
        .get("is_legend")
        .ok_or_else(|| anyhow!("missing context key: is_legend"))?
        .as_bool()
        .ok_or_else(|| anyhow!("context key is_legend is not a bool"))?;
    if !is_legend {
        replace_section(&mut template, "%LEGEND_BEGIN%", "%LEGEND_END%", "")?;
    }

    render_and_write(&template, context, &mut writer)?;
    Ok(())
}
