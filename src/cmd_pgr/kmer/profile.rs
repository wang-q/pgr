use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for profile.
pub fn make_subcommand() -> Command {
    Command::new("profile")
        .about("Generates per-sequence k-mer profiles (.pkp)")
        .after_help(
            r###"
Generates one k-mer count profile per sequence (read or chromosome) and
writes them to a `.pkp` file. For every k-mer position of every input
sequence the profile records one count; the counts are looked up from a
k-mer table, either built on the fly from the input or reused via --table.

* Without --table (self): the input sequences are counted first and each
  position reports how many times its k-mer occurs in the input dataset;
  repeated regions therefore show high values (FastK `-p` semantics).
* With --table (relative): each position reports the count of its k-mer in
  the given table; positions whose k-mer is absent from the table report 0
  (FastK `-p:<table>` semantics). This is a lookup against an external
  table, not a comparison between profiles.

Both modes write the same `.pkp` format; only the source of the counts
differs.

* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Self profile of a genome:
   pgr kmer profile genome.fa -k 17 -o genome.pkp
2. Reads profile relative to a repeat table:
   pgr kmer profile reads.fq.gz -t lib.pkt -o reads.pkp
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTA/FASTQ file to process",
        ))
        .arg(table_arg())
        .arg(kmer_arg())
        .arg(crate::cmd_pgr::args::outfile_arg_required())
}

/// Optional `-t/--table` argument for profile and hist.
pub(super) fn table_arg() -> Arg {
    Arg::new("table")
        .long("table")
        .short('t')
        .num_args(1)
        .help("Reuse a k-mer table (.pkt); k is read from the table")
}

/// Optional `-k/--kmer` argument; required unless --table is given.
pub(super) fn kmer_arg() -> Arg {
    Arg::new("kmer")
        .long("kmer")
        .short('k')
        .num_args(1)
        .value_parser(clap::value_parser!(usize))
        .help("K-mer size (required unless --table is given)")
}

/// Execute the profile command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();
    let k = super::resolve_k(
        args.get_one::<usize>("kmer"),
        args.get_one::<String>("table").map(String::as_str),
    )?;
    let table_path = args.get_one::<String>("table").map(String::as_str);
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        [infile.as_str()].into_iter().chain(table_path),
    )?;

    let seqs = super::read_seqs(infile)?;
    let profiles = if let Some(t) = table_path {
        let table = pgr::libs::kmer::count::load(std::path::Path::new(t), k)?;
        pgr::libs::kmer::profile::relative_profiles(&seqs, k, &table)
    } else {
        let table = pgr::libs::kmer::count::build_table(&seqs, k)?;
        pgr::libs::kmer::profile::self_profiles(&seqs, k, &table)
    };
    pgr::libs::kmer::profile::save_profiles(std::path::Path::new(outfile), k, &profiles)?;
    log::info!("==> Wrote {} profiles to {}", profiles.len(), outfile);
    Ok(())
}
