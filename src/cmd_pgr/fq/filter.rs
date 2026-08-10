use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::fq::trim_adapter::{trim_adapter, AdapterTrimOptions};
use std::io::Write;

/// Build the clap subcommand for filter.
pub fn make_subcommand() -> Command {
    Command::new("filter")
        .about("Filters reads matching reference k-mers (bbduk kfilter)")
        .after_help(
            r###"
Discards reads containing k-mers matching a reference (adapters,
contaminants, spike-ins). Reproduces the second BBTools 39.38 `bbduk.sh`
call of the anchr trim pipeline (`k=<matchk> cardinality`) byte for byte
(`ordered=t`, deterministic).

Notes:
* A read is discarded when more than zero k-mers match the reference
  (bbduk `minkmerhits=1`); surviving mates follow `--toss-broken-reads`
* Defaults match the bbduk filter call: k=27, mink=0, hdist=0, minlen=10
* For k-mer *trimming* (the first bbduk call) use `pgr fq clean`; for
  sickle-style pure quality trimming use `pgr fq trim-qual`
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Filter adapter/artifact matches:
   pgr fq filter in.fq --ref adapters.fa -o out.fq

2. With per-reference match statistics:
   pgr fq filter in.fq --ref adapters.fa -k 27 --stats R.filter.stats.txt \
       -o out.fq
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
                .required(true)
                .help("Reference FASTA of contaminants/adapters"),
        )
        .arg(
            Arg::new("k")
                .long("k")
                .short('k')
                .num_args(1)
                .default_value("27")
                .value_parser(value_parser!(usize))
                .help("K-mer size"),
        )
        .arg(
            Arg::new("min_k")
                .long("min-k")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Minimum short k-mer size at read ends"),
        )
        .arg(
            Arg::new("hamming_distance")
                .long("hamming-distance")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Reference hamming distance (bbduk: hdist)"),
        )
        .arg(
            Arg::new("minlen")
                .long("minlen")
                .num_args(1)
                .default_value("10")
                .value_parser(value_parser!(usize))
                .help("Minimum kept read length"),
        )
        .arg(
            Arg::new("max_ns")
                .long("max-ns")
                .num_args(1)
                .default_value("-1")
                .value_parser(value_parser!(i64))
                .help("Maximum allowed N bases; negative disables"),
        )
        .arg(
            Arg::new("no_toss")
                .long("no-toss-broken-reads")
                .action(clap::ArgAction::SetTrue)
                .help("Keep surviving mates of discarded reads (bbduk: removeifeitherbad=f)"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("auto")
                .help("Worker threads (bbduk: threads)"),
        )
        .arg(
            Arg::new("stats")
                .long("stats")
                .num_args(1)
                .help("Write per-reference match statistics"),
        )
}

/// Execute the filter command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let ref_file = args.get_one::<String>("ref").cloned();
    let parallel =
        crate::cmd_pgr::args::parse_parallel_auto(args.get_one::<String>("parallel").unwrap())?;
    let opts = AdapterTrimOptions {
        k: *args.get_one::<usize>("k").unwrap(),
        mink: *args.get_one::<usize>("min_k").unwrap(),
        hdist: *args.get_one::<usize>("hamming_distance").unwrap(),
        ktrim_right: false,
        tbo: false,
        tpe: false,
        qtrim_right: false,
        qtrim_left: false,
        qtrim_window: 0,
        trimq: 0.0,
        minlen: *args.get_one::<usize>("minlen").unwrap(),
        maxns: *args.get_one::<i64>("max_ns").unwrap(),
        ftm: 0,
        force_trim_left: 0,
        force_trim_right: 0,
        force_trim_right2: 0,
        toss_broken_reads: !args.get_flag("no_toss"),
        ref_file: ref_file.clone(),
        quality_base: 33,
        max_bad_kmers: 0,
        trim_poly_a: 0,
        trim_poly_g_left: 0,
        trim_poly_g_right: 0,
        filter_poly_g: 0,
        trim_poly_c_left: 0,
        trim_poly_c_right: 0,
        filter_poly_c: 0,
        max_non_poly: 1,
        min_avg_quality: 0.0,
        min_avg_quality_bases: 0,
        min_base_quality: 0,
        max_n_rate: 1.0,
        min_len_fraction: 0.0,
        max_length: usize::MAX,
        min_consecutive_bases: 0,
        min_gc: 0.0,
        max_gc: 1.0,
        use_pair_gc: true,
        kmask_symbol: None,
        kmask_lowercase: false,
        kmask_fully_covered: false,
        trim_pad: 0,
        stats: args.get_one::<String>("stats").cloned(),
    };
    if !(2..=31).contains(&opts.k) {
        anyhow::bail!("--k must be in 2..=31, got {}", opts.k);
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, infiles.iter().map(String::as_str))?;
    let mut out =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    trim_adapter(&infiles, &mut out, &opts, parallel)?;
    out.flush()?;
    Ok(())
}
