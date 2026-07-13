use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for filter.
pub fn make_subcommand() -> Command {
    Command::new("filter")
        .about("Filters blocks and optionally formats sequences")
        .after_help(
            r###"
Filters blocks in block FA files based on species name and sequence length.
It can also format sequences by converting them to uppercase or removing dashes.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* If `--name` is not specified, the first species in each block is used as the default
* Sequences can be filtered based on length using `--min-len` (greater than or equal) and `--max-len` (less than or equal)
* Sequences can be formatted using `--upper` (convert to uppercase) and `--dash` (remove dashes)

Examples:
1. Filter blocks for a specific species:
   pgr fas filter tests/fas/example.fas --name S288c

2. Filter blocks with sequences >= 100 bp:
   pgr fas filter tests/fas/example.fas --min-len 100

3. Filter blocks with sequences <= 200 bp:
   pgr fas filter tests/fas/example.fas --max-len 200

4. Convert sequences to uppercase and remove dashes:
   pgr fas filter tests/fas/example.fas --upper --dash

5. Output results to a file:
   pgr fas filter tests/fas/example.fas -o output.fas

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(
            crate::cmd_pgr::args::fas_name_arg("Filter blocks based on this species"),
        )
        .arg(crate::cmd_pgr::args::min_len_arg())
        .arg(crate::cmd_pgr::args::max_len_arg())
        .arg(crate::cmd_pgr::args::upper_arg())
        .arg(crate::cmd_pgr::args::dash_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the filter command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let opt_name: &str = args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("");
    let opt_min = args.get_one::<usize>("min_len").copied();
    let opt_max = args.get_one::<usize>("max_len").copied();

    let is_upper = args.get_flag("upper");
    let is_dash = args.get_flag("dash");

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            // Determine the index of the species
            if block.entries.is_empty() {
                continue;
            }
            let idx = if !opt_name.is_empty() {
                match block.names.iter().position(|x| x == opt_name) {
                    Some(i) => i,
                    None => continue,
                }
            } else {
                0
            };

            let idx_seq = block.entries[idx].seq();

            if let Some(min) = opt_min {
                if idx_seq.len() < min {
                    continue;
                }
            }

            if let Some(max) = opt_max {
                if idx_seq.len() > max {
                    continue;
                }
            }

            for entry in &block.entries {
                let out_seq = pgr::libs::fasta::filter::format_sequence(
                    entry.seq(),
                    is_dash,
                    false,
                    is_upper,
                );

                let out_entry =
                    pgr::libs::fmt::fas::FasEntry::from(entry.range(), out_seq.as_bytes());
                writer.write_all(out_entry.to_string().as_ref())?;
            }

            // end of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    writer.flush()?;
    Ok(())
}
