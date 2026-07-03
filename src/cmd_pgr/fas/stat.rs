use clap::{ArgMatches, Command};
use std::io::Write;

// Create clap subcommand arguments
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
   pgr fas stat tests/fas/part1.fas

2. Statistics treating the last sequence as an outgroup:
   pgr fas stat tests/fas/part1.fas --outgroup

3. Output results to a file:
   pgr fas stat tests/fas/part1.fas -o output.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::outgroup_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
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

    //----------------------------
    // Operating
    //----------------------------
    writer.write_all(format!("{}\n", field_names.join("\t")).as_ref())?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
            let target = block
                .entries
                .first()
                .ok_or_else(|| anyhow::anyhow!("empty block"))?
                .range()
                .to_string();

            let mut seqs: Vec<&[u8]> = vec![];
            for entry in &block.entries {
                seqs.push(entry.seq().as_ref());
            }

            if has_outgroup {
                seqs.pop();
            }

            // let (length, comparable, difference, gap, ambiguous, mean_d) = alignment_stat(&seqs);
            let result = pgr::libs::alignment::alignment_stat(&seqs);

            let mut indel_ints = intspan::IntSpan::new();
            for seq in seqs {
                indel_ints.merge(&pgr::libs::alignment::indel_intspan(seq));
            }

            writer.write_all(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    target,
                    result.0,
                    result.1,
                    result.2,
                    result.3,
                    result.4,
                    result.5,
                    indel_ints.span_size(),
                )
                .as_ref(),
            )?;
        }
    }

    Ok(())
}
