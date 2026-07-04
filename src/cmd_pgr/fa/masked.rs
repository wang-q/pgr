use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for masked.
pub fn make_subcommand() -> Command {
    Command::new("masked")
        .about("Identifies masked regions in FASTA file(s)")
        .after_help(
            r###"
This command identifies masked regions in one or more FASTA files. Masked regions can be:
- Lowercase letters
- Regions of N/n

The output is a list of regions in the format:
    seq_name:start-end        # For regions spanning multiple positions
    seq_name:position        # For single positions

Notes:
* Coordinates are 1-based, inclusive
* Supports both plain text and gzipped (.gz) files
* Adjacent masked positions are merged into a single region

Examples:
1. Identify masked regions (lowercase and N/n):
   pgr fa masked input.fa -o masked_regions.txt

2. Identify only N/n gap regions:
   pgr fa masked input.fa --gap -o gap_regions.txt

3. Process multiple input files:
   pgr fa masked input1.fa input2.fa -o masked_regions.txt

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("FASTA"))
        .arg(crate::cmd_pgr::args::gap_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the masked command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let is_gap = args.get_flag("gap");
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut fa_in = pgr::libs::fmt::fa::reader(infile)?;

        for result in fa_in.records() {
            let record = result?;
            let name = String::from_utf8(record.name().into())?;
            let seq = record.sequence();

            for (begin, end) in pgr::libs::fmt::fa::find_masked_regions(&seq[..], is_gap) {
                writer.write_all(out_line(&name, begin, end).as_ref())?;
            }
        }
    }

    Ok(())
}

// Generate the output line for a masked region
fn out_line(name: &str, begin: usize, end: usize) -> String {
    if begin == end {
        format!("{}:{}\n", name, begin + 1)
    } else {
        format!("{}:{}-{}\n", name, begin + 1, end + 1)
    }
}
