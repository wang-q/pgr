use anyhow::anyhow;
use indexmap::IndexMap;
use intspan::IntSpan;
use std::path::Path;

use super::common::{context_get_str, render_and_write, replace_section};

/// Build named IntSpan sets from a list of input files (one item per line).
/// Each set is labeled by the file's basename (extension stripped); duplicate
/// labels are renamed to `cat{i}` (1-based).
pub fn build_venn_sets_from_files(infiles: &[String]) -> anyhow::Result<IndexMap<String, IntSpan>> {
    let mut ints_of: IndexMap<String, IntSpan> = IndexMap::new();
    let mut all_elems = indexmap::IndexSet::new();

    for (i, file) in infiles.iter().enumerate() {
        let mut basename = Path::new(file)
            .file_name()
            .ok_or_else(|| anyhow!("invalid filename: {}", file))?
            .to_str()
            .ok_or_else(|| anyhow!("invalid UTF-8 in filename: {}", file))?
            .split('.')
            .next()
            .ok_or_else(|| anyhow!("empty filename after splitting: {}", file))?
            .to_string();

        if ints_of.contains_key(&basename) {
            basename = format!("cat{}", i + 1);
        }

        let vec = crate::libs::io::read_names::<Vec<String>>(file)?;
        let mut ints = IntSpan::new();

        for e in &vec {
            all_elems.insert(e.clone());
            let idx = all_elems
                .get_index_of(e)
                .ok_or_else(|| anyhow!("element not found after insert: {}", e))?;
            ints.add_n(idx as i32);
        }
        ints_of.insert(basename, ints);
    }

    Ok(ints_of)
}

/// Result of a Venn set-operation computation: exclusive element counts per set
/// and intersection counts ordered from lowest-order to highest-order intersections.
pub struct VennResult {
    /// Sizes of elements exclusive to each set (A only, B only, ...).
    pub excls: Vec<i32>,
    /// Sizes of intersections, ordered binary, then triple, ..., then the n-fold intersection.
    pub inter: Vec<i32>,
}

/// Compute Venn counts for 2 sets.
pub fn venn_sets_2(a: &IntSpan, b: &IntSpan) -> VennResult {
    let mut excls = Vec::new();
    let mut inter = Vec::new();

    // A ∩ B
    let i_ab = a.intersect(b).size();
    inter.push(i_ab);

    // A - B
    excls.push(a.diff(b).size());
    // B - A
    excls.push(b.diff(a).size());

    VennResult { excls, inter }
}

/// Compute Venn counts for 3 sets.
pub fn venn_sets_3(a: &IntSpan, b: &IntSpan, c: &IntSpan) -> VennResult {
    let mut excls = Vec::new();
    let mut inter = Vec::new();

    // A ∩ B ∩ C
    let i_abc = a.intersect(b).intersect(c);

    // Binary intersections minus triple intersection
    let sets_arr = [a, b, c];
    for i in 0..2 {
        for j in (i + 1)..=2 {
            let intersection = sets_arr[i].intersect(sets_arr[j]).diff(&i_abc).size();
            inter.push(intersection);
        }
    }

    inter.push(i_abc.size());

    // A - B - C
    excls.push(a.diff(b).diff(c).size());
    // B - A - C
    excls.push(b.diff(a).diff(c).size());
    // C - A - B
    excls.push(c.diff(a).diff(b).size());

    VennResult { excls, inter }
}

/// Compute Venn counts for 4 sets.
pub fn venn_sets_4(a: &IntSpan, b: &IntSpan, c: &IntSpan, d: &IntSpan) -> VennResult {
    let mut excls = Vec::new();
    let mut inter = Vec::new();

    // Quadruple intersection
    let i_abcd = a.intersect(b).intersect(c).intersect(d);

    // Binary intersections
    let sets_arr = [a, b, c, d];
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
    excls.push(a.diff(b).diff(c).diff(d).size());
    excls.push(b.diff(a).diff(c).diff(d).size());
    excls.push(c.diff(a).diff(b).diff(d).size());
    excls.push(d.diff(a).diff(b).diff(c).size());

    VennResult { excls, inter }
}

/// LaTeX/TikZ snippet for the 2-set Venn diagram body.
pub const VENN_2: &str = r###"
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

/// LaTeX/TikZ snippet for the 3-set Venn diagram body.
pub const VENN_3: &str = r###"
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

/// LaTeX/TikZ snippet for the 4-set Venn diagram body.
pub const VENN_4: &str = r###"
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

/// Render the Venn diagram LaTeX file by substituting the body template
/// selected by `n` (2, 3, or 4) into the file template, then rendering via Tera.
/// For unsupported `n`, returns `Ok(())` without writing.
pub fn gen_venn(context: &tera::Context, n: usize) -> anyhow::Result<()> {
    let out_string = match n {
        2 => VENN_2,
        3 => VENN_3,
        4 => VENN_4,
        _ => return Ok(()),
    };

    let outfile = context_get_str(context, "outfile")?;
    let mut writer = crate::writer(outfile)?;

    static FILE_TEMPLATE: &str = include_str!("../../assets/venn.tex");
    let mut template = FILE_TEMPLATE.to_string();

    replace_section(&mut template, "%VENN_BEGIN", "%VENN_END", out_string)?;

    render_and_write(&template, context, &mut writer)?;
    Ok(())
}
