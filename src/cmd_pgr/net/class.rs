use anyhow::Context;
use clap::{ArgMatches, Command};
use std::collections::HashMap;
use std::io::Write;

use pgr::libs::chain::net::{collect_stats_gap, read_nets, Stats};
/// Build the clap subcommand for class.
pub fn make_subcommand() -> Command {
    Command::new("class")
        .about("Shows stats of net")
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input net file (or stdin if 'stdin')",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the class command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("infile").unwrap();

    let reader = pgr::reader(input_path)
        .with_context(|| format!("Failed to open reader for {}", input_path))?;

    let chroms = read_nets(reader)?;

    let mut stats_map: HashMap<String, Stats> = HashMap::new();
    let mut total_bases = 0;

    for chrom in chroms {
        total_bases += chrom.size;
        collect_stats_gap(&chrom.root, &mut stats_map);
    }

    // Print results
    // We want to sort by bases desc? Or just alphabetical?
    // UCSC usually sorts by bases or hierarchy.
    // Let's sort by bases desc.

    let mut results: Vec<(String, u64, u64)> = stats_map
        .into_iter()
        .map(|(k, v)| (k, v.count, v.bases))
        .collect();

    results.sort_by_key(|b| std::cmp::Reverse(b.2)); // Sort by bases descending

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    writeln!(
        writer,
        "{:<20} {:>10} {:>15} {:>10}",
        "Class", "Count", "Bases", "%"
    )?;
    writeln!(
        writer,
        "{:<20} {:>10} {:>15} {:>10}",
        "-----", "-----", "-----", "-"
    )?;

    for (class, count, bases) in results {
        let pct = if total_bases > 0 {
            (bases as f64 / total_bases as f64) * 100.0
        } else {
            0.0
        };
        writeln!(
            writer,
            "{:<20} {:>10} {:>15} {:>10.2}%",
            class, count, bases, pct
        )?;
    }

    writeln!(
        writer,
        "\nTotal bases covered by nets/gaps: {}",
        total_bases
    )?;

    writer.flush()?;
    Ok(())
}
