use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::paf::parser::parse_paf_line;
use pgr::libs::paf::validate::ValidationReport;
use std::io::BufRead;
/// Build the clap subcommand for validate.
pub fn make_subcommand() -> Command {
    Command::new("validate")
        .about("Checks PAF end coordinates against the cg:Z: CIGAR tag")
        .after_help(
            r###"
For each PAF record, reconstructs the expected query and target end positions
from the cg:Z: CIGAR tag (matches + mismatches + insertion/deletion bases) and
compares them with the declared coordinates. Ends that disagree are reported
as query/target invalid.

Notes:
* Records without a usable cg:Z: tag are counted and skipped (not fatal)
* Malformed cg:Z: tags are counted and skipped (not fatal)
* Supports both plain text and gzipped (.gz) files
* Reads PAF from stdin if input file is 'stdin'

Examples:
1. Validate a PAF file and print the report to screen:
   pgr paf validate alignments.paf

2. Save the report to a file:
   pgr paf validate alignments.paf -o report.txt
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PAF file (or 'stdin' for piped input)",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the validate command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, [infile.as_str()])?;

    let reader =
        pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
    let mut report = ValidationReport::default();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let rec = parse_paf_line(&line).with_context(|| format!("Invalid PAF line: {line}"))?;
        report.validate(&rec)?;
    }

    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    report.write_report(&mut writer)?;

    Ok(())
}
