use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

use super::expand::{open_aln, record_to_paf, write_paf};

/// Build the clap subcommand for to-paf.
pub fn make_subcommand() -> Command {
    Command::new("to-paf")
        .about("Converts a FastGA .1aln file to PAF format")
        .after_help(
            r###"
Expands each alignment record in a FastGA `.1aln` (ONEcode trace-point) file
back into base-level aligned columns and emits a PAF record per alignment.

The `.1aln` header stores only the source genome file references, so the two
source genomes must be supplied with --ref-seq and --query-seq (FASTA or
gzipped FASTA). The reference genome is the `a` side and becomes the query
column of the PAF; the other genome is the `b` side and becomes the target.

Custom PAF tags:
* `dv:f:` - identity (matches / aligned bases)
* `df:i:` - number of differences (substitutions + indels)
* `cg:Z:` - X-CIGAR (only with --cigar)

A `-` strand means the `b`-side sequence was stored reverse-complemented in
the `.1aln`; its PAF coordinates are given in forward orientation.

Notes:
* Requires --ref-seq and --query-seq (the source genomes).
* Reads a single .1aln file; does not support gzip or stdin (the ONEcode
  container requires random access to the footer offset at EOF).

Examples:
1. Convert with default tags (no CIGAR):
   pgr 1aln to-paf mg1655-sakai.1aln --ref-seq mg1655.fa.gz \
       --query-seq sakai.fa.gz -o out.paf
2. Include the cg:Z CIGAR tag:
   pgr 1aln to-paf mg1655-sakai.1aln --ref-seq mg1655.fa.gz \
       --query-seq sakai.fa.gz --cigar -o out.paf

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
        .arg(
            clap::Arg::new("cigar")
                .long("cigar")
                .action(clap::ArgAction::SetTrue)
                .help("Emit the cg:Z CIGAR tag"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the to-paf command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let (mut aln, genomes) = open_aln(args)?;
    let with_cigar = args.get_flag("cigar");
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {outfile}"))?;

    while let Some(rec) = aln.next_record()? {
        let paf = record_to_paf(&rec, aln.tspace, &genomes, &aln, with_cigar)?;
        write_paf(&mut writer, &paf)?;
    }
    writer.flush()?;
    Ok(())
}
