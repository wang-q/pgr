use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::fq::trim_adapter::{trim_adapter, AdapterTrimOptions};
use std::io::Write;

/// Build the clap subcommand for trim-adapter.
pub fn make_subcommand() -> Command {
    Command::new("trim-adapter")
        .about("Trims adapters by k-mer matching (bbduk-compatible)")
        .after_help(
            r###"
This command removes adapter/contaminant sequences by matching read k-mers
against a reference, then quality-trims and length-filters the reads. It
reproduces BBTools 39.38 `bbduk.sh` output byte for byte for the anchr trim
pipeline parameters (`ordered=t`, deterministic).

Notes:
* `ktrim=r` right-trims at the first matching reference k-mer (`--ktrim`)
* `--mink` enables short k-mer matching at read ends (adapters shorter than k)
* `--hdist` stores single-substitution reference variants
* `--tbo` trims implied adapters from mate overlap; `--tpe` equalizes mates
* `--qtrim r` right-trims below `--trimq`; `--minlen` drops short reads
* `--maxns` drops reads with too many N bases; `--ftm` right-trims lengths to
  a multiple; `--toss-broken-reads` drops pairs where one mate fails
* Input is one interleaved FASTQ or two files (R1, R2)
* --threads controls the worker pool (default: logical CPU count); output
  order is preserved regardless of thread count
* Supports both plain text and gzipped (.gz) files

Examples:
1. Adapter trim with the anchr pipeline defaults:
   pgr fq trim-adapter R1.fq.gz R2.fq.gz --ref illumina_adapters.fa \
       -o out.fq

2. K-mer filtering mode (match and discard, bbduk filter step):
   pgr fq trim-adapter in.fq --ref illumina_adapters.fa --k 27 \
       --no-ktrim -o out.fq
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input FASTQ file(s): 1 interleaved or 2 paired (R1, R2)",
            1..=2,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("ref")
                .long("ref")
                .num_args(1)
                .help("Reference FASTA of adapters/contaminants"),
        )
        .arg(
            Arg::new("k")
                .long("k")
                .short('k')
                .num_args(1)
                .default_value("23")
                .value_parser(value_parser!(usize))
                .help("K-mer size"),
        )
        .arg(
            Arg::new("mink")
                .long("mink")
                .num_args(1)
                .default_value("11")
                .value_parser(value_parser!(usize))
                .help("Minimum short k-mer size at read ends"),
        )
        .arg(
            Arg::new("hdist")
                .long("hdist")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Reference hamming distance"),
        )
        .arg(
            Arg::new("no_ktrim")
                .long("no-ktrim")
                .action(clap::ArgAction::SetTrue)
                .help("Disable k-mer trimming (filtering mode)"),
        )
        .arg(
            Arg::new("no_tbo")
                .long("no-tbo")
                .action(clap::ArgAction::SetTrue)
                .help("Disable overlap trimming"),
        )
        .arg(
            Arg::new("no_tpe")
                .long("no-tpe")
                .action(clap::ArgAction::SetTrue)
                .help("Disable even pair trimming"),
        )
        .arg(
            Arg::new("no_qtrim")
                .long("no-qtrim")
                .action(clap::ArgAction::SetTrue)
                .help("Disable quality trimming"),
        )
        .arg(
            Arg::new("trimq")
                .long("trimq")
                .num_args(1)
                .default_value("15")
                .value_parser(value_parser!(u8))
                .help("Quality threshold for qtrim"),
        )
        .arg(
            Arg::new("minlen")
                .long("minlen")
                .num_args(1)
                .default_value("60")
                .value_parser(value_parser!(usize))
                .help("Minimum kept read length"),
        )
        .arg(
            Arg::new("maxns")
                .long("maxns")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(i64))
                .help("Maximum allowed N bases"),
        )
        .arg(
            Arg::new("ftm")
                .long("ftm")
                .num_args(1)
                .default_value("5")
                .value_parser(value_parser!(usize))
                .help("Right-trim lengths to a multiple (0 disables)"),
        )
        .arg(
            Arg::new("no_toss")
                .long("no-toss-broken-reads")
                .action(clap::ArgAction::SetTrue)
                .help("Keep surviving mates of discarded reads"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .short('p')
                .num_args(1)
                .default_value("auto")
                .help("Worker threads (default: logical CPU count)"),
        )
}

/// Execute the trim-adapter command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let ref_file = args.get_one::<String>("ref").unwrap();
    let threads = match args.get_one::<String>("threads").unwrap().as_str() {
        "auto" => pgr::libs::sys::logical_cpus(),
        s => s
            .parse::<usize>()
            .map_err(|_| anyhow::anyhow!("invalid --threads: {}", s))?,
    };
    let opts = AdapterTrimOptions {
        k: *args.get_one::<usize>("k").unwrap(),
        mink: *args.get_one::<usize>("mink").unwrap(),
        hdist: *args.get_one::<usize>("hdist").unwrap(),
        ktrim_right: !args.get_flag("no_ktrim"),
        tbo: !args.get_flag("no_tbo"),
        tpe: !args.get_flag("no_tpe"),
        qtrim_right: !args.get_flag("no_qtrim"),
        trimq: *args.get_one::<u8>("trimq").unwrap(),
        minlen: *args.get_one::<usize>("minlen").unwrap(),
        maxns: *args.get_one::<i64>("maxns").unwrap(),
        ftm: *args.get_one::<usize>("ftm").unwrap(),
        toss_broken_reads: !args.get_flag("no_toss"),
        ref_file: ref_file.clone(),
        quality_base: 33,
        max_bad_kmers: 0,
    };
    if !(2..=31).contains(&opts.k) {
        anyhow::bail!("--k must be in 2..=31, got {}", opts.k);
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, infiles.iter().map(String::as_str))?;
    let mut out =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    trim_adapter(&infiles, &mut out, &opts, threads)?;
    out.flush()?;
    Ok(())
}
