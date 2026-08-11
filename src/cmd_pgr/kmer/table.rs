use clap::{ArgMatches, Command};

/// Build the clap subcommand for table.
pub fn make_subcommand() -> Command {
    Command::new("table")
        .about("Builds a k-mer count table (.pkt)")
        .after_help(
            r###"
Builds a canonical k-mer count table from FASTA/FASTQ sequences. Counts
accumulate across all input files (FastK `-t1` semantics: every k-mer is
kept, even singletons).

* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Build a table from reads:
   pgr kmer table reads.fq.gz -k 21 -o reads.pkt
2. Merge several genomes into one table:
   pgr kmer table a.fa b.fa.gz -k 17 -o all.pkt
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("FASTA/FASTQ"))
        .arg(crate::cmd_pgr::args::kmer_arg_with_default("17"))
        .arg(crate::cmd_pgr::args::outfile_arg_required())
}

/// Execute the table command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles = crate::cmd_pgr::args::collect_infiles(args);
    let k = *args.get_one::<usize>("kmer").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, infiles.iter().copied())?;

    let mut seqs = Vec::new();
    for f in infiles {
        seqs.extend(super::read_seqs(f)?);
    }
    let table = pgr::libs::kmer::count::build_table(&seqs, k)?;
    pgr::libs::kmer::count::save(&table, std::path::Path::new(outfile))?;
    log::info!(
        "==> Wrote {} unique {}-mers to {}",
        table.counts.len(),
        k,
        outfile
    );
    Ok(())
}
