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
* Without --ref, no k-mer operations run and the command only quality-trims
  and filters (bbduk `qtrim=r minlen=...` without a reference)
* `ktrim=r` right-trims at the first matching reference k-mer (`--ktrim`)
* `--mink` enables short k-mer matching at read ends (adapters shorter than k)
* `--hdist` stores single-substitution reference variants
* `--tbo` trims implied adapters from mate overlap; `--tpe` equalizes mates
* `--qtrim r` right-trims below `--trimq`; `--minlen` drops short reads
* `--maxns` drops reads with too many N bases; `--ftm` right-trims lengths to
  a multiple; `--toss-broken-reads` drops pairs where one mate fails
* Input is one interleaved FASTQ or two files (R1, R2)
* --parallel controls the worker pool (default: logical CPU count); output
  order is preserved regardless of thread count
* --stats writes per-reference match statistics in the bbduk `stats=`
  tab-separated format
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
                .help("Reference FASTA of adapters/contaminants (omit for quality-trim-only)"),
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
            Arg::new("qtrim")
                .long("qtrim")
                .num_args(1)
                .default_value("r")
                .value_parser(["r", "l", "rl", "w", "f"])
                .help("Quality trim mode: r, l, rl, w (sliding window), or f"),
        )
        .arg(
            Arg::new("qtrim_window")
                .long("qtrim-window")
                .num_args(1)
                .default_value("4")
                .value_parser(value_parser!(usize))
                .help("Window size for qtrim=w"),
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
            Arg::new("forcetrim_left")
                .long("forcetrim-left")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim bases left of this position (0 disables)"),
        )
        .arg(
            Arg::new("forcetrim_right")
                .long("forcetrim-right")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim bases right of this position (0 disables)"),
        )
        .arg(
            Arg::new("forcetrim_right2")
                .long("forcetrim-right2")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim this many bases on the right end (0 disables)"),
        )
        .arg(
            Arg::new("trim_poly_a")
                .long("trim-poly-a")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim poly-A/T tails of at least this length"),
        )
        .arg(
            Arg::new("trim_poly_g_left")
                .long("trim-poly-g-left")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim poly-G prefixes of at least this length"),
        )
        .arg(
            Arg::new("trim_poly_g_right")
                .long("trim-poly-g-right")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim poly-G tails of at least this length"),
        )
        .arg(
            Arg::new("filter_poly_g")
                .long("filter-poly-g")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Discard reads with a poly-G prefix of at least this length"),
        )
        .arg(
            Arg::new("trim_poly_c_left")
                .long("trim-poly-c-left")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim poly-C prefixes of at least this length"),
        )
        .arg(
            Arg::new("trim_poly_c_right")
                .long("trim-poly-c-right")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim poly-C tails of at least this length"),
        )
        .arg(
            Arg::new("filter_poly_c")
                .long("filter-poly-c")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Discard reads with a poly-C prefix of at least this length"),
        )
        .arg(
            Arg::new("max_non_poly")
                .long("max-non-poly")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Allowed non-polymer bases inside a polymer run"),
        )
        .arg(
            Arg::new("maq")
                .long("maq")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(f64))
                .help("Discard reads with average quality below this"),
        )
        .arg(
            Arg::new("maqb")
                .long("maqb")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Use only this many leading bases for maq"),
        )
        .arg(
            Arg::new("mbq")
                .long("mbq")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(u8))
                .help("Discard reads with any base below this quality"),
        )
        .arg(
            Arg::new("maxnrate")
                .long("maxnrate")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(f64))
                .help("Discard reads with more than this fraction of Ns"),
        )
        .arg(
            Arg::new("minlen_fraction")
                .long("minlen-fraction")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(f64))
                .help("Minimum read length as a fraction of the original (mlf)"),
        )
        .arg(
            Arg::new("maxlength")
                .long("maxlength")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Discard reads longer than this (0 disables)"),
        )
        .arg(
            Arg::new("mcb")
                .long("mcb")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Discard reads without this many consecutive ACGT bases"),
        )
        .arg(
            Arg::new("mingc")
                .long("mingc")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(f64))
                .help("Discard reads with GC content below this"),
        )
        .arg(
            Arg::new("maxgc")
                .long("maxgc")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(f64))
                .help("Discard reads with GC content above this"),
        )
        .arg(
            Arg::new("no_pair_gc")
                .long("no-pair-gc")
                .action(clap::ArgAction::SetTrue)
                .help("Check GC per read instead of the pair average"),
        )
        .arg(
            Arg::new("kmask")
                .long("kmask")
                .num_args(1)
                .help("Mask matching k-mers: a symbol, 'lc' for lowercase, or 't' for N"),
        )
        .arg(
            Arg::new("mask_fully_covered")
                .long("mask-fully-covered")
                .action(clap::ArgAction::SetTrue)
                .help("Only mask bases fully covered by matching k-mers"),
        )
        .arg(
            Arg::new("trim_pad")
                .long("trim-pad")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Extra bases to mask around matching k-mers"),
        )
        .arg(
            Arg::new("no_toss")
                .long("no-toss-broken-reads")
                .action(clap::ArgAction::SetTrue)
                .help("Keep surviving mates of discarded reads"),
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
            Arg::new("stats")
                .long("stats")
                .num_args(1)
                .help("Write per-reference match statistics (bbduk stats=)"),
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
    let ref_file = args.get_one::<String>("ref").cloned();
    let parallel = match args.get_one::<String>("parallel").unwrap().as_str() {
        "auto" => pgr::libs::sys::logical_cpus(),
        s => s
            .parse::<usize>()
            .map_err(|_| anyhow::anyhow!("invalid --threads: {}", s))?,
    };
    let no_qtrim = args.get_flag("no_qtrim");
    let qtrim = args.get_one::<String>("qtrim").unwrap().as_str();
    let (qtrim_left, qtrim_right, qtrim_window) = if no_qtrim || qtrim == "f" {
        (false, false, 0)
    } else {
        match qtrim {
            "l" => (true, false, 0),
            "rl" => (true, true, 0),
            "w" => (false, true, *args.get_one::<usize>("qtrim_window").unwrap()),
            _ => (false, true, 0),
        }
    };
    let kmask = args.get_one::<String>("kmask").map(String::as_str);
    let (kmask_symbol, kmask_lowercase) = match kmask {
        Some("lc") => (None, true),
        Some("t") | Some("true") => (Some(b'N'), false),
        Some(c) if c.len() == 1 => (Some(c.as_bytes()[0]), false),
        Some(other) => anyhow::bail!("invalid --kmask value: {other}"),
        None => (None, false),
    };
    let maxlength = *args.get_one::<usize>("maxlength").unwrap();
    let opts = AdapterTrimOptions {
        k: *args.get_one::<usize>("k").unwrap(),
        mink: *args.get_one::<usize>("mink").unwrap(),
        hdist: *args.get_one::<usize>("hdist").unwrap(),
        ktrim_right: !args.get_flag("no_ktrim"),
        tbo: !args.get_flag("no_tbo"),
        tpe: !args.get_flag("no_tpe"),
        qtrim_right,
        qtrim_left,
        qtrim_window,
        trimq: *args.get_one::<u8>("trimq").unwrap(),
        minlen: *args.get_one::<usize>("minlen").unwrap(),
        maxns: *args.get_one::<i64>("maxns").unwrap(),
        ftm: *args.get_one::<usize>("ftm").unwrap(),
        force_trim_left: *args.get_one::<usize>("forcetrim_left").unwrap(),
        force_trim_right: *args.get_one::<usize>("forcetrim_right").unwrap(),
        force_trim_right2: *args.get_one::<usize>("forcetrim_right2").unwrap(),
        toss_broken_reads: !args.get_flag("no_toss"),
        ref_file: ref_file.clone(),
        quality_base: 33,
        max_bad_kmers: 0,
        trim_poly_a: *args.get_one::<usize>("trim_poly_a").unwrap(),
        trim_poly_g_left: *args.get_one::<usize>("trim_poly_g_left").unwrap(),
        trim_poly_g_right: *args.get_one::<usize>("trim_poly_g_right").unwrap(),
        filter_poly_g: *args.get_one::<usize>("filter_poly_g").unwrap(),
        trim_poly_c_left: *args.get_one::<usize>("trim_poly_c_left").unwrap(),
        trim_poly_c_right: *args.get_one::<usize>("trim_poly_c_right").unwrap(),
        filter_poly_c: *args.get_one::<usize>("filter_poly_c").unwrap(),
        max_non_poly: *args.get_one::<usize>("max_non_poly").unwrap(),
        min_avg_quality: *args.get_one::<f64>("maq").unwrap(),
        min_avg_quality_bases: *args.get_one::<usize>("maqb").unwrap(),
        min_base_quality: *args.get_one::<u8>("mbq").unwrap(),
        max_n_rate: *args.get_one::<f64>("maxnrate").unwrap(),
        min_len_fraction: *args.get_one::<f64>("minlen_fraction").unwrap(),
        max_length: if maxlength == 0 {
            usize::MAX
        } else {
            maxlength
        },
        min_consecutive_bases: *args.get_one::<usize>("mcb").unwrap(),
        min_gc: *args.get_one::<f64>("mingc").unwrap(),
        max_gc: *args.get_one::<f64>("maxgc").unwrap(),
        use_pair_gc: !args.get_flag("no_pair_gc"),
        kmask_symbol,
        kmask_lowercase,
        kmask_fully_covered: args.get_flag("mask_fully_covered"),
        trim_pad: *args.get_one::<usize>("trim_pad").unwrap(),
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
