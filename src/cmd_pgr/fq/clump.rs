use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::fq::clump::clump;
use std::io::Write;

/// Build the clap subcommand for clump.
pub fn make_subcommand() -> Command {
    Command::new("clump")
        .about("Sorts reads by k-mer signature (clumpify-compatible)")
        .after_help(
            r###"
This command sorts interleaved paired reads by the pivot k-mer of R1,
reproducing the BBTools `clumpify.sh` default output order byte for byte.
The sorting clusters reads that share k-mers, which speeds up the k-mer
steps that follow in a read-cleaning pipeline.

Notes:
* Paired input must be interleaved; mates are kept together
* Deterministic for a given k-mer size and seed
* Supports both plain text and gzipped (.gz) files

Examples:
1. Sort reads with the BBTools-compatible defaults:
   pgr fq clump in.fq.gz -o clumped.fq.gz

2. Reproduce a BBTools run with a different seed:
   pgr fq clump in.fq.gz -o out.fq --seed 2
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input FASTQ file(s): 1 interleaved or 2 paired (R1, R2)",
            1..=2,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("31")
                .value_parser(value_parser!(usize))
                .help("K-mer size (2..=31)"),
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(u64))
                .help("Comparator seed"),
        )
}

/// Execute the clump command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let k = *args.get_one::<usize>("kmer").unwrap();
    let seed = *args.get_one::<u64>("seed").unwrap();
    if !(2..=31).contains(&k) {
        anyhow::bail!("--kmer must be in 2..=31, got {}", k);
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, infiles.iter().map(String::as_str))?;
    let mut out =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    clump(&infiles, &mut out, k, seed)?;
    out.flush()?;
    Ok(())
}
