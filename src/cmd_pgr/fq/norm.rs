use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::fq::norm::{norm, NormOptions};
use std::io::Write;

/// Build the clap subcommand for norm.
pub fn make_subcommand() -> Command {
    Command::new("norm")
        .about("Filters reads by k-mer depth (bbnorm-style cutoff)")
        .after_help(
            r###"
This command removes reads whose k-mer coverage is below a minimum depth,
following the BBTools 39.38 `bbnorm.sh passes=1 bits=16 min=<n>
target=9999999` read decision logic. The k-mer counts are exact (canonical
table) instead of the approximate `bits=16` hash table used by bbnorm.

Notes:
* Paired reads are kept only when both mates pass the depth threshold
* Input is one interleaved FASTQ or two files (R1, R2)
* Supports both plain text and gzipped (.gz) files

Examples:
1. Keep reads with at least one k-mer at depth 3:
   pgr fq norm reads.fq.gz -k 31 --min 3 -o out.fq
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
                .help("K-mer size"),
        )
        .arg(
            Arg::new("min")
                .long("min")
                .num_args(1)
                .default_value("3")
                .value_parser(value_parser!(usize))
                .help("Minimum k-mer depth cutoff"),
        )
}

/// Execute the norm command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let opts = NormOptions {
        k: *args.get_one::<usize>("kmer").unwrap(),
        min_depth: *args.get_one::<usize>("min").unwrap(),
    };
    if !(2..=31).contains(&opts.k) {
        anyhow::bail!("--kmer must be in 2..=31, got {}", opts.k);
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, infiles.iter().map(String::as_str))?;
    let mut out =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    norm(&infiles, &mut out, &opts)?;
    out.flush()?;
    Ok(())
}
