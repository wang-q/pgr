use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for stat.
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Calculates basic statistics of block FA files")
        .after_help(
            r###"
Calculates basic statistics of block FA files (length, comparable, difference, etc.).

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Get statistics for block FA files:
   pgr fas stat tests/fas/example.fas

2. Statistics treating the last sequence as an outgroup:
   pgr fas stat tests/fas/example.fas --outgroup

3. Output results to a file:
   pgr fas stat tests/fas/example.fas -o output.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::outgroup_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the stat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let has_outgroup = args.get_flag("outgroup");

    let field_names = [
        "target",
        "length",
        "comparable",
        "difference",
        "gap",
        "ambiguous",
        "D",
        "indel",
    ];

    // Operating
    writer.write_all(format!("{}\n", field_names.join("\t")).as_ref())?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            if block.entries.is_empty() {
                continue;
            }
            let target = block.entries.first().unwrap().range().to_string();

            let mut seqs: Vec<&[u8]> = vec![];
            for entry in &block.entries {
                seqs.push(entry.seq());
            }

            if has_outgroup {
                if seqs.len() < 2 {
                    anyhow::bail!(
                        "block has only {} entries, cannot apply --outgroup",
                        seqs.len()
                    );
                }
                seqs.pop();
            }

            let (length, comparable, difference, gap, ambiguous, mean_d) =
                pgr::libs::alignment::alignment_stat(&seqs)?;

            let mut indel_ints = intspan::IntSpan::new();
            for seq in seqs {
                indel_ints.merge(&pgr::libs::alignment::indel_intspan(seq));
            }

            writer.write_all(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    target,
                    length,
                    comparable,
                    difference,
                    gap,
                    ambiguous,
                    mean_d,
                    indel_ints.span_size(),
                )
                .as_ref(),
            )?;
        }
    }

    writer.flush()?;
    Ok(())
}
