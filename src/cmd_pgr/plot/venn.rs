use crate::cmd_pgr::args::outfile_arg;
use anyhow::{anyhow, Context, Result};
use clap::*;
use indexmap::IndexMap;
use std::path::Path;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("venn")
        .about("Plot Venn diagram for 2-4 sets")
        .after_help(
            r###"
Generates a LaTeX file for a Venn diagram representing the intersections of sets.

Notes:
* Input files should contain lists of items (one per line).
* Supports 2, 3, or 4 sets.
* Output is a standalone LaTeX file using TikZ.

Examples:
1. Two sets:
   pgr plot venn list1.txt list2.txt -o venn.tex

2. Three sets:
   pgr plot venn list1.txt list2.txt list3.txt -o venn.tex
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..=4)
                .index(1)
                .help("Input list files (2-4 files)"),
        )
        .arg(outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .ok_or_else(|| anyhow!("missing infiles"))?
        .map(|s| s.to_string())
        .collect();

    //----------------------------
    // Ops
    //----------------------------
    let mut ints_of: IndexMap<String, _> = indexmap::IndexMap::new();
    let mut all_elems = indexmap::IndexSet::new();

    for (i, file) in infiles.iter().enumerate() {
        // Get filename as label
        let mut basename = Path::new(file)
            .file_name()
            .ok_or_else(|| anyhow!("invalid filename: {}", file))?
            .to_str()
            .ok_or_else(|| anyhow!("invalid UTF-8 in filename: {}", file))?
            .split('.')
            .next()
            .ok_or_else(|| anyhow!("empty filename after splitting: {}", file))?
            .to_string();

        // Handle duplicate names
        if ints_of.contains_key(&basename) {
            basename = format!("cat{}", i + 1);
        }

        // Read file content and convert to IntSpan
        let vec = intspan::read_first_column(file);
        let mut ints = intspan::IntSpan::new();

        for e in &vec {
            all_elems.insert(e.clone());
            let idx = all_elems
                .get_index_of(e)
                .ok_or_else(|| anyhow!("element not found after insert: {}", e))?;
            ints.add_n(idx as i32);
        }
        ints_of.insert(basename, ints);
    }

    let mut excls = Vec::new(); // sizes of exclusive elements
    let mut inter = Vec::new(); // sizes of intersections

    match ints_of.len() {
        2 => {
            let set_a = ints_of
                .get_index(0)
                .ok_or_else(|| anyhow!("missing set 0"))?
                .1;
            let set_b = ints_of
                .get_index(1)
                .ok_or_else(|| anyhow!("missing set 1"))?
                .1;

            // A ∩ B
            let i_ab = set_a.intersect(set_b).size();
            inter.push(i_ab);

            // A - B
            excls.push(set_a.diff(set_b).size());
            // B - A
            excls.push(set_b.diff(set_a).size());
        }
        3 => {
            let set_a = ints_of
                .get_index(0)
                .ok_or_else(|| anyhow!("missing set 0"))?
                .1;
            let set_b = ints_of
                .get_index(1)
                .ok_or_else(|| anyhow!("missing set 1"))?
                .1;
            let set_c = ints_of
                .get_index(2)
                .ok_or_else(|| anyhow!("missing set 2"))?
                .1;

            // A ∩ B ∩ C
            let i_abc = set_a.intersect(set_b).intersect(set_c);

            // Binary intersections minus triple intersection
            let sets_arr = [set_a, set_b, set_c];
            for i in 0..2 {
                for j in (i + 1)..=2 {
                    let intersection = sets_arr[i].intersect(sets_arr[j]).diff(&i_abc).size();
                    inter.push(intersection);
                }
            }

            inter.push(i_abc.size());

            // A - B - C
            excls.push(set_a.diff(set_b).diff(set_c).size());
            // B - A - C
            excls.push(set_b.diff(set_a).diff(set_c).size());
            // C - A - B
            excls.push(set_c.diff(set_a).diff(set_b).size());
        }
        4 => {
            let set_a = ints_of
                .get_index(0)
                .ok_or_else(|| anyhow!("missing set 0"))?
                .1;
            let set_b = ints_of
                .get_index(1)
                .ok_or_else(|| anyhow!("missing set 1"))?
                .1;
            let set_c = ints_of
                .get_index(2)
                .ok_or_else(|| anyhow!("missing set 2"))?
                .1;
            let set_d = ints_of
                .get_index(3)
                .ok_or_else(|| anyhow!("missing set 3"))?
                .1;

            // Quadruple intersection
            let i_abcd = set_a.intersect(set_b).intersect(set_c).intersect(set_d);

            // Binary intersections
            let sets_arr = [set_a, set_b, set_c, set_d];
            for i in 0..3 {
                for j in (i + 1)..=3 {
                    let mut i_temp = sets_arr[i].intersect(sets_arr[j]);
                    // Subtract all higher-order intersections containing these two sets
                    for (k, _) in sets_arr.iter().enumerate() {
                        if k != i && k != j {
                            i_temp.subtract(sets_arr[k]);
                        }
                    }
                    inter.push(i_temp.size());
                }
            }

            // Triple intersections
            for i in 0..2 {
                for j in (i + 1)..3 {
                    for k in (j + 1)..=3 {
                        let i_temp = sets_arr[i]
                            .intersect(sets_arr[j])
                            .intersect(sets_arr[k])
                            .diff(&i_abcd);
                        inter.push(i_temp.size());
                    }
                }
            }

            // Quadruple intersection
            inter.push(i_abcd.size());

            // Exclusive elements
            excls.push(set_a.diff(set_b).diff(set_c).diff(set_d).size());
            excls.push(set_b.diff(set_a).diff(set_c).diff(set_d).size());
            excls.push(set_c.diff(set_a).diff(set_b).diff(set_d).size());
            excls.push(set_d.diff(set_a).diff(set_b).diff(set_c).size());
        }
        _ => {}
    }

    //----------------------------
    // Context
    //----------------------------
    let mut context = tera::Context::new();

    let outfile = args
        .get_one::<String>("outfile")
        .ok_or_else(|| anyhow!("missing outfile"))?;
    context.insert("outfile", outfile);
    context.insert("label", &ints_of.keys().collect::<Vec<&String>>());
    context.insert("excls", &excls);
    context.insert("inter", &inter);

    if ints_of.len() == 2 {
        gen_venn_2(&context)?;
    } else if ints_of.len() == 3 {
        gen_venn_3(&context)?;
    } else if ints_of.len() == 4 {
        gen_venn_4(&context)?;
    }

    Ok(())
}

fn gen_venn_2(context: &tera::Context) -> Result<()> {
    let outfile = context_get_str(context, "outfile")?;
    let mut writer = pgr::writer(outfile)?;

    static FILE_TEMPLATE: &str = include_str!("../../../docs/venn.tex");
    let mut template = FILE_TEMPLATE.to_string();

    let out_string = r###"
% Basic parameters for circles
\def\radius{2cm}
\def\overlap{1.2cm}

% Draw two circles
\begin{scope}[opacity=0.5]
    \fill[indianred1] (-\overlap,0) circle (\radius);
    \fill[deepskyblue] (\overlap,0) circle (\radius);
\end{scope}

% Add circle edges
\draw[grey, thick] (-\overlap,0) circle (\radius);
\draw[grey, thick] (\overlap,0) circle (\radius);

% Add labels
\node[text centered] at (-2.8, -1.8) { {{ label.0 }} };
\node[text centered] at (2.8,  -1.8) { {{ label.1 }} };

% Add numbers
\node[text centered] at (-2,    0) { {{ excls.0 }} };
\node[text centered] at (2,     0) { {{ excls.1 }} };
\node[text centered] at (0,     0) { {{ inter.0 }} };
    "###;

    {
        // Section venn
        let begin = template
            .find("%VENN_BEGIN")
            .ok_or_else(|| anyhow!("venn template anchor %VENN_BEGIN not found"))?;
        let end = template
            .find("%VENN_END")
            .ok_or_else(|| anyhow!("venn template anchor %VENN_END not found"))?;
        template.replace_range(begin..end, out_string);
    }

    let mut tera = tera::Tera::default();
    tera.add_raw_templates(vec![("t", template)])
        .context("failed to register venn template")?;

    let rendered = tera
        .render("t", context)
        .context("failed to render venn template")?;
    writer.write_all(rendered.as_ref())?;

    Ok(())
}

fn gen_venn_3(context: &tera::Context) -> Result<()> {
    let outfile = context_get_str(context, "outfile")?;
    let mut writer = pgr::writer(outfile)?;

    static FILE_TEMPLATE: &str = include_str!("../../../docs/venn.tex");
    let mut template = FILE_TEMPLATE.to_string();

    let out_string = r###"
% Basic parameters for circles
\def\radius{2cm}
\def\xshift{1.2cm}
\def\yshift{2.08cm}

% Draw three circles
\begin{scope}[opacity=0.5]
    \fill[indianred1] (-\xshift,0) circle (\radius);
    \fill[deepskyblue] (0,\yshift) circle (\radius);
    \fill[palegreen] (\xshift,0) circle (\radius);
\end{scope}

% Add circle edges
\draw[grey, thick] (-\xshift,0) circle (\radius);
\draw[grey, thick] (0,\yshift) circle (\radius);
\draw[grey, thick] (\xshift,0) circle (\radius);

% Add labels
\node[text centered] at (-2.8, -1.8) { {{ label.0 }} };
\node[text centered] at (-1.8,  3.8) { {{ label.1 }} };
\node[text centered] at (2.8,  -1.8) { {{ label.2 }} };

% Add numbers for exclusive regions
\node[text centered] at (-2,   -0.2) { {{ excls.0 }} };
\node[text centered] at (0,     2.8) { {{ excls.1 }} };
\node[text centered] at (2,    -0.2) { {{ excls.2 }} };

% Add numbers for binary intersections
\node[text centered] at (-1.2,  1.2) { {{ inter.0 }} }; % AB
\node[text centered] at (0,    -0.7) { {{ inter.1 }} }; % AC
\node[text centered] at (1.2,   1.2) { {{ inter.2 }} }; % BC

% Add number for triple intersection
\node[text centered] at (0,     0.6) { {{ inter.3 }} }; % ABC
    "###;

    {
        // Section venn
        let begin = template
            .find("%VENN_BEGIN")
            .ok_or_else(|| anyhow!("venn template anchor %VENN_BEGIN not found"))?;
        let end = template
            .find("%VENN_END")
            .ok_or_else(|| anyhow!("venn template anchor %VENN_END not found"))?;
        template.replace_range(begin..end, out_string);
    }

    let mut tera = tera::Tera::default();
    tera.add_raw_templates(vec![("t", template)])
        .context("failed to register venn template")?;

    let rendered = tera
        .render("t", context)
        .context("failed to render venn template")?;
    writer.write_all(rendered.as_ref())?;

    Ok(())
}

fn gen_venn_4(context: &tera::Context) -> Result<()> {
    let outfile = context_get_str(context, "outfile")?;
    let mut writer = pgr::writer(outfile)?;

    static FILE_TEMPLATE: &str = include_str!("../../../docs/venn.tex");
    let mut template = FILE_TEMPLATE.to_string();

    let out_string = r###"
% Basic parameters for ellipses
\def\xradius{3.5cm}
\def\yradius{2cm}
\def\yradiusB{1.8cm}
\def\xshift{1.5cm}
\def\xshiftB{0.7cm}
\def\yshift{1cm}

% Draw four ellipses
\begin{scope}[opacity=0.5]
    \fill[indianred1] (-\xshift,0) ellipse [x radius=\xradius, y radius=\yradius, rotate=-45];    % A
    \fill[deepskyblue] (-\xshiftB,\yshift) ellipse [x radius=\xradius, y radius=\yradiusB, rotate=-45];    % B
    \fill[palegreen] (\xshiftB,\yshift) ellipse [x radius=\xradius, y radius=\yradiusB, rotate=45];      % C
    \fill[burlywood] (\xshift,0) ellipse [x radius=\xradius, y radius=\yradius, rotate=45];    % D
\end{scope}

% Add ellipse edges
\draw[grey, thick] (-\xshift,0) ellipse [x radius=\xradius, y radius=\yradius, rotate=-45];
\draw[grey, thick] (-\xshiftB,\yshift) ellipse [x radius=\xradius, y radius=\yradiusB, rotate=-45];
\draw[grey, thick] (\xshiftB,\yshift) ellipse [x radius=\xradius, y radius=\yradiusB, rotate=45];
\draw[grey, thick] (\xshift,0) ellipse [x radius=\xradius, y radius=\yradius, rotate=45];

% Add labels
\node[text centered] at (-2.2, -2.6) { {{ label.0 }} }; % A
\node[text centered] at (-3.8,  3.6) { {{ label.1 }} }; % B
\node[text centered] at (3.8,   3.6) { {{ label.2 }} }; % C
\node[text centered] at (2.2,  -2.6) { {{ label.3 }} }; % D

% Add numbers for exclusive regions
\node[text centered] at (-3.2,    0) { {{ excls.0 }} }; % A
\node[text centered] at (-2,    3.2) { {{ excls.1 }} }; % B
\node[text centered] at (2,     3.2) { {{ excls.2 }} }; % C
\node[text centered] at (3.2,     0) { {{ excls.3 }} }; % D

%%% Add numbers for binary intersections
\node[text centered] at (-2.2,  1.5) { {{ inter.0 }} }; % AB
\node[text centered] at (-1.7, -1.0) { {{ inter.1 }} }; % AC
\node[text centered] at (0,    -2.2) { {{ inter.2 }} }; % AD
\node[text centered] at (0,     2.0) { {{ inter.3 }} }; % BC
\node[text centered] at (1.7,  -1.0) { {{ inter.4 }} }; % BD
\node[text centered] at (2.2,   1.5) { {{ inter.5 }} }; % CD

% Add numbers for triple intersections
\node[text centered] at (-1.1,  0.7) { {{ inter.6 }} }; % ABC
\node[text centered] at (0.9,  -1.5) { {{ inter.7 }} }; % ABD
\node[text centered] at (-0.9, -1.5) { {{ inter.8 }} }; % ACD
\node[text centered] at (1.1,   0.7) { {{ inter.9 }} }; % BCD

% Add number for quadruple intersection
\node[text centered] at (0,    -0.2) { {{ inter.10 }} }; % ABCD
    "###;

    {
        // Section venn
        let begin = template
            .find("%VENN_BEGIN")
            .ok_or_else(|| anyhow!("venn template anchor %VENN_BEGIN not found"))?;
        let end = template
            .find("%VENN_END")
            .ok_or_else(|| anyhow!("venn template anchor %VENN_END not found"))?;
        template.replace_range(begin..end, out_string);
    }

    let mut tera = tera::Tera::default();
    tera.add_raw_templates(vec![("t", template)])
        .context("failed to register venn template")?;

    let rendered = tera
        .render("t", context)
        .context("failed to render venn template")?;
    writer.write_all(rendered.as_ref())?;

    Ok(())
}

// Helper: get a string value from tera::Context, replacing the common
// `context.get(k).unwrap().as_str().unwrap()` pattern with a friendly error.
fn context_get_str<'a>(context: &'a tera::Context, key: &str) -> Result<&'a str> {
    context
        .get(key)
        .ok_or_else(|| anyhow!("missing context key: {}", key))?
        .as_str()
        .ok_or_else(|| anyhow!("context key {} is not a string", key))
}
