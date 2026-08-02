//! `pgr sd decompose` — elementary SD decomposition from cluster FASTA.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for decompose.
pub fn make_subcommand() -> Command {
    Command::new("decompose")
        .about("Decomposes a cluster FASTA into elementary SD fragments")
        .after_help(
            r###"
Reads one cluster FASTA (headers in `{species}#{chrom}{strand}#{start}#{end}`
form, as produced by `pgr sd cluster`) and writes elementary SD BED rows:
`species<TAB>chrom<TAB>begin<TAB>end<TAB>set_id<TAB>length<TAB>score<TAB>strand`.

Shared k-mers (present in >= 2 sequences of the cluster) seed fragments that
are merged with a gap tolerance of 50 bp; fragments shorter than 100 bp are
dropped (T2T-CHM13 style, see notes/design/sd.md §4.5).

Examples:
1. Decompose one cluster:
   pgr sd decompose cluster_1.fa -o cluster_1.elem.bed
"###,
        )
        .arg(
            Arg::new("infile")
                .index(1)
                .required(true)
                .help("Cluster FASTA file"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the decompose command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let reader = pgr::reader(infile)?;
    let mut writer = pgr::writer(outfile)?;
    pgr::libs::sd::decompose::decompose_fasta(reader, &mut writer)?;
    Ok(())
}
