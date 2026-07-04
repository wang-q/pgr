use crate::cmd_pgr::args::{infiles_arg_with_numargs, outfile_arg};
use anyhow::anyhow;
use clap::{ArgMatches, Command};
use pgr::libs::plot::venn::{venn_sets_2, venn_sets_3, venn_sets_4};

/// Build the clap subcommand for venn.
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
            2..=4,
        ))
        .arg(outfile_arg())
}

/// Execute the venn command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .ok_or_else(|| anyhow!("missing infiles"))?
        .map(|s| s.to_string())
        .collect();

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

    // Context
    let mut context = tera::Context::new();

    let outfile = args
        .get_one::<String>("outfile")
        .ok_or_else(|| anyhow!("missing outfile"))?;
    context.insert("outfile", outfile);
    context.insert("label", &ints_of.keys().collect::<Vec<&String>>());
    context.insert("excls", &excls);
    context.insert("inter", &inter);

    pgr::libs::plot::venn::gen_venn(&context, ints_of.len())?;

    Ok(())
}
