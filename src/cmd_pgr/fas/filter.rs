use clap::*;
use std::io::Write;

// Create clap subcommand arguments
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
   pgr fas filter tests/fasr/example.fas --name S288c

2. Filter blocks with sequences >= 100 bp:
   pgr fas filter tests/fasr/example.fas --min-len 100

3. Filter blocks with sequences <= 200 bp:
   pgr fas filter tests/fasr/example.fas --max-len 200

4. Convert sequences to uppercase and remove dashes:
   pgr fas filter tests/fasr/example.fas --upper --dash

5. Output results to a file:
   pgr fas filter tests/fasr/example.fas -o output.fas

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

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let opt_name = &args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("")
        .to_string();
    let opt_ge = args
        .get_one::<usize>("min_len")
        .copied()
        .unwrap_or(usize::MAX);
    let opt_le = args
        .get_one::<usize>("max_len")
        .copied()
        .unwrap_or(usize::MAX);

    let is_upper = args.get_flag("upper");
    let is_dash = args.get_flag("dash");

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        'BLOCK: while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
            // Determine the index of the species
            let idx = if !opt_name.is_empty() {
                if !block.names.contains(opt_name) {
                    continue 'BLOCK;
                }
                block.names.iter().position(|x| x == opt_name).unwrap()
            } else {
                0
            };

            let idx_seq = block.entries[idx].seq();

            // --ge
            if opt_ge != usize::MAX && idx_seq.len() < opt_ge {
                continue 'BLOCK;
            }

            // --le
            if opt_le != usize::MAX && idx_seq.len() > opt_le {
                continue 'BLOCK;
            }

            for entry in &block.entries {
                let out_seq = pgr::libs::fasta::filter::format_sequence(
                    entry.seq(),
                    is_dash,
                    false,
                    is_upper,
                );

                //----------------------------
                // Output
                //----------------------------
                let out_entry =
                    pgr::libs::fmt::fas::FasEntry::from(entry.range(), out_seq.as_bytes());
                writer.write_all(out_entry.to_string().as_ref())?;
            }

            // end of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    Ok(())
}
