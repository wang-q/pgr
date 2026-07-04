use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use itertools::Itertools;
use std::io::Write;

/// Build the clap subcommand for link.
pub fn make_subcommand() -> Command {
    Command::new("link")
        .about("Outputs bi/multi-lateral range links from block FA files")
        .after_help(
            r###"
Outputs bi/multi-lateral range links from block FA files.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* By default, the tool outputs multi-lateral links (all ranges in a block)
* Use `--pair` to output bilateral (pairwise) links
* Use `--best` to output best-to-best bilateral links based on sequence similarity

Examples:
1. Output multi-lateral links:
   pgr fas link tests/fas/example.fas

2. Output bilateral (pairwise) links:
   pgr fas link tests/fas/example.fas --pair

3. Output best-to-best bilateral links:
   pgr fas link tests/fas/example.fas --best

4. Output results to a file:
   pgr fas link tests/fas/example.fas -o output.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(
            Arg::new("pair")
                .long("pair")
                .action(ArgAction::SetTrue)
                .help("Output bilateral (pairwise) links"),
        )
        .arg(
            Arg::new("best")
                .long("best")
                .action(ArgAction::SetTrue)
                .help("Output best-to-best bilateral links"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the link command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args)).with_context(|| {
        format!(
            "Failed to open writer for {}",
            crate::cmd_pgr::args::get_outfile(args)
        )
    })?;
    let is_pair = args.get_flag("pair");
    let is_best = args.get_flag("best");

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            let headers: Vec<String> = block
                .entries
                .iter()
                .map(|entry| entry.range().to_string())
                .collect();

            if is_pair {
                // Output bilateral (pairwise) links
                for (i, j) in (0..headers.len()).tuple_combinations() {
                    writer.write_all(format!("{}\t{}\n", headers[i], headers[j]).as_ref())?;
                }
            } else if is_best {
                let best_pair = pgr::libs::fmt::fas::find_best_pairs(&block.entries)?;
                for (i, j) in best_pair {
                    writer.write_all(format!("{}\t{}\n", headers[i], headers[j]).as_ref())?;
                }
            } else {
                // Output multi-lateral links
                writer.write_all(format!("{}\n", headers.join("\t")).as_ref())?;
            }
        }
    }

    Ok(())
}
