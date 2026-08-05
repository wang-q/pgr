use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

use super::expand::open_aln;
use pgr::libs::onepack::expand::record_to_psl;

/// Build the clap subcommand for to-psl.
pub fn make_subcommand() -> Command {
    Command::new("to-psl")
        .about("Converts a FastGA .1aln file to PSL format")
        .after_help(
            r###"
Expands each alignment record in a FastGA `.1aln` (ONEcode trace-point) file
back into base-level aligned columns and emits a PSL record per alignment.

The `.1aln` header stores only the source genome file references, so the two
source genomes must be supplied with --ref-seq and --query-seq (FASTA or
gzipped FASTA). The reference genome is the `a` side and becomes the query
(`q`) of the PSL; the other genome is the `b` side and becomes the target
(`t`).

A `+-` strand means the target (`b`) sequence was stored reverse-complemented
in the `.1aln`; the emitted target coordinates are in forward orientation.

Notes:
* Requires --ref-seq and --query-seq (the source genomes).
* Reads a single .1aln file; does not support gzip or stdin (the ONEcode
  container requires random access to the footer offset at EOF).

Examples:
1. Convert:
   pgr 1aln to-psl mg1655-sakai.1aln --ref-seq mg1655.fa.gz \
       --query-seq sakai.fa.gz -o out.psl

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input .1aln file",
        ))
        .arg(
            clap::Arg::new("ref_seq")
                .long("ref-seq")
                .required(true)
                .num_args(1)
                .help("Reference (a side) genome FASTA(.gz)"),
        )
        .arg(
            clap::Arg::new("query_seq")
                .long("query-seq")
                .required(true)
                .num_args(1)
                .help("Query (b side) genome FASTA(.gz)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the to-psl command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let (mut aln, genomes) = open_aln(args)?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {outfile}"))?;

    while let Some(rec) = aln.next_record()? {
        let psl = record_to_psl(&rec, aln.tspace, &genomes, &aln)?;
        psl.write_to(&mut writer)?;
    }
    writer.flush()?;
    Ok(())
}
