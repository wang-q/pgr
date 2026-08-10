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
* --mem sets the in-memory count budget (KMG, default 2g); data estimated to
  exceed it is counted via external hash buckets with bounded memory
* Supports both plain text and gzipped (.gz) files

Examples:
1. Keep reads with at least one k-mer at depth 3:
   pgr fq norm reads.fq.gz -k 31 --min 3 -o out.fq

2. Bound memory to 1 GiB (external bucket path for larger data):
   pgr fq norm R1.fq.gz R2.fq.gz -k 31 --min 3 --mem 1g -o out.fq
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
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("auto")
                .help("Worker threads (default: logical CPU count)"),
        )
        .arg(
            Arg::new("mem")
                .long("mem")
                .num_args(1)
                .default_value("2g")
                .help("In-memory count budget (KMG; default 2g)"),
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
        mem: Some(pgr::libs::sys::parse_mem_size(
            args.get_one::<String>("mem").unwrap(),
        )?),
    };
    let parallel =
        crate::cmd_pgr::args::parse_parallel_auto(args.get_one::<String>("parallel").unwrap())?;
    if !(2..=31).contains(&opts.k) {
        anyhow::bail!("--kmer must be in 2..=31, got {}", opts.k);
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, infiles.iter().map(String::as_str))?;
    let mut out =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    norm(&infiles, &mut out, &opts, parallel)?;
    out.flush()?;
    Ok(())
}
