use clap::{Arg, ArgAction, ArgMatches, Command};

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
* Default counter is the direct radix path; --supermer switches to the
  FastK-style super-mer two-stage counter (identical output, faster only on
  read data with high coverage and moderate k; see notes/design/kmer.md)

Examples:
1. Build a table from reads:
   pgr kmer table reads.fq.gz -k 21 -o reads.pkt
2. Merge several genomes into one table:
   pgr kmer table a.fa b.fa.gz -k 17 -o all.pkt
3. Use the super-mer counter explicitly:
   pgr kmer table reads.fq.gz -k 31 --supermer -o reads.pkt
4. Override the super-mer minimizer length:
   pgr kmer table reads.fq.gz -k 31 --supermer --minimizer 8 -o reads.pkt
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("FASTA/FASTQ"))
        .arg(crate::cmd_pgr::args::kmer_arg_with_default("17"))
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("supermer")
                .long("supermer")
                .action(ArgAction::SetTrue)
                .help("Uses the FastK-style super-mer two-stage counter"),
        )
        .arg(
            Arg::new("minimizer")
                .long("minimizer")
                .value_parser(clap::value_parser!(usize))
                .requires("supermer")
                .help("Minimizer length for --supermer (default: min(12, max(5, ceil(k/4))))"),
        )
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
    let table = if args.get_flag("supermer") {
        match args.get_one::<usize>("minimizer") {
            Some(&m) => pgr::libs::kmer::supermer::build_table_with_m(&seqs, k, m)?,
            None => pgr::libs::kmer::supermer::build_table(&seqs, k)?,
        }
    } else {
        pgr::libs::kmer::count::build_table(&seqs, k)?
    };
    pgr::libs::kmer::count::save(&table, std::path::Path::new(outfile))?;
    log::info!(
        "==> Wrote {} unique {}-mers to {}",
        table.counts.len(),
        k,
        outfile
    );
    Ok(())
}
