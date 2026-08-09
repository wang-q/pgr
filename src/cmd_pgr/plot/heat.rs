use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for heat.
pub fn make_subcommand() -> Command {
    Command::new("heat")
        .about("Plots a GC-content x coverage heatmap as LaTeX")
        .after_help(
            r###"
Plots the GC-content x k-mer coverage heatmap (KatGC heat plot equivalent)
from a `.kgc` matrix (`pgr kmer gc`) as a standalone LaTeX document.

* Compile with tectonic to get a PDF

Examples:
1. Plot from a GC matrix:
   pgr kmer gc reads.fq.gz -k 21 -o reads.kgc
   pgr plot heat reads.kgc -o heat.tex
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input .kgc GC matrix file",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the heat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let text =
        std::fs::read_to_string(infile).with_context(|| format!("failed to read {infile}"))?;
    let mut rows = Vec::new();
    let mut zmax = 0u64;
    for line in text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 3 {
            continue;
        }
        let gc = cols[0].trim_end_matches(".5").parse::<usize>().ok();
        let kf = cols[1].trim_end_matches(".5").parse::<usize>().ok();
        let count = cols[2].parse::<u64>().ok();
        if let (Some(i), Some(a), Some(c)) = (gc, kf, count) {
            rows.push((i, a, c));
            zmax = zmax.max(c);
        }
    }
    anyhow::ensure!(!rows.is_empty(), "no matrix rows found in {infile}");

    let hm = pgr::libs::plot::heat::heatmap_from_kgc(&rows, zmax);
    let mut w = pgr::writer(outfile)?;
    pgr::libs::plot::heat::render_heat(&mut w, &hm)?;
    w.flush()?;
    log::info!("==> Wrote GC x coverage heatmap to {}", outfile);
    Ok(())
}
