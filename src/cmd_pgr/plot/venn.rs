use crate::cmd_pgr::args::{infiles_arg_with_numargs, outfile_arg};
use anyhow::{anyhow, Result};
use clap::{ArgMatches, Command};
use pgr::libs::plot::common::{context_get_str, render_and_write, replace_section};
use pgr::libs::plot::venn::{venn_sets_2, venn_sets_3, venn_sets_4};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("venn")
        .about("Plots Venn diagram for 2-4 sets")
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
        .arg(infiles_arg_with_numargs(
            "Input list files (2-4 files)",
            1..=4,
        ))
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
    let ints_of = pgr::libs::plot::venn::build_venn_sets_from_files(&infiles)?;

    let get_set = |i: usize| -> anyhow::Result<&intspan::IntSpan> {
        Ok(ints_of
            .get_index(i)
            .ok_or_else(|| anyhow!("missing set {}", i))?
            .1)
    };

    let (excls, inter) = match ints_of.len() {
        2 => {
            let r = venn_sets_2(get_set(0)?, get_set(1)?);
            (r.excls, r.inter)
        }
        3 => {
            let r = venn_sets_3(get_set(0)?, get_set(1)?, get_set(2)?);
            (r.excls, r.inter)
        }
        4 => {
            let r = venn_sets_4(get_set(0)?, get_set(1)?, get_set(2)?, get_set(3)?);
            (r.excls, r.inter)
        }
        _ => (Vec::new(), Vec::new()),
    };

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

    let out_string: &str = match ints_of.len() {
        2 => VENN_2,
        3 => VENN_3,
        4 => VENN_4,
        _ => return Ok(()),
    };
    gen_venn(&context, out_string)?;

    Ok(())
}

const VENN_2: &str = r###"
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

const VENN_3: &str = r###"
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

const VENN_4: &str = r###"
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

fn gen_venn(context: &tera::Context, out_string: &str) -> Result<()> {
    let outfile = context_get_str(context, "outfile")?;
    let mut writer = pgr::writer(outfile)?;

    static FILE_TEMPLATE: &str = include_str!("../../assets/venn.tex");
    let mut template = FILE_TEMPLATE.to_string();

    replace_section(&mut template, "%VENN_BEGIN", "%VENN_END", out_string)?;

    render_and_write(&template, context, &mut writer)?;
    Ok(())
}
