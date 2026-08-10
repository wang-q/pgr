use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::fq::trim_adapter::{trim_adapter, AdapterTrimOptions};
use std::io::Write;

/// Build the clap subcommand for clean.
pub fn make_subcommand() -> Command {
    Command::new("clean")
        .about("Cleans reads: adapter k-mer trimming, quality and composition filtering")
        .after_help(
            r###"
This command cleans reads in one pass: adapter/contaminant k-mer trimming,
quality trimming, polymer and GC filtering, and masking. It reproduces the
first BBTools 39.38 `bbduk.sh` call of the anchr trim pipeline (the `ktrim`
pass) byte for byte (`ordered=t`, deterministic).

Notes:
* Without --ref, no k-mer operations run and the command only quality-trims
  and filters (bbduk `qtrim=r minlen=...` without a reference)
* k-mer right-trimming is on by default (bbduk `ktrim=r`); `--min-k` enables
  short k-mer matching at read ends (adapters shorter than k)
* `--hamming-distance` stores single-substitution reference variants
* `--trim-by-overlap` trims implied adapters from mate overlap;
  `--trim-pair-evenly` equalizes mates
* `--qtrim r` right-trims below `--trim-quality`; `--minlen` drops short reads
* `--max-ns` drops reads with too many N bases; `--force-trim-mod`
  right-trims lengths to a multiple; `--toss-broken-reads` drops pairs where
  one mate fails
* For k-mer contaminant filtering (bbduk `kfilter`, the second bbduk call of
  the pipeline) use `pgr fq filter` instead; for sickle-style pure quality
  trimming use `pgr fq trim-qual`.
* Options renamed from their bbduk counterparts show the bbduk name in
  parentheses (e.g. `--min-k` = `mink`); identical options are not annotated
* Input is one interleaved FASTQ or two files (R1, R2)
* --parallel controls the worker pool (default: logical CPU count); output
  order is preserved regardless of thread count
* --stats writes per-reference match statistics in the bbduk `stats=`
  tab-separated format
* Supports both plain text and gzipped (.gz) files

Examples:
1. Adapter trim with the anchr pipeline defaults:
   pgr fq clean R1.fq.gz R2.fq.gz --ref illumina_adapters.fa \
       -o out.fq

2. K-mer filtering mode (match and discard, bbduk filter step):
   pgr fq clean in.fq --ref illumina_adapters.fa -o out.fq
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
            Arg::new("min_k")
                .long("min-k")
                .num_args(1)
                .default_value("11")
                .value_parser(value_parser!(usize))
                .help("Minimum short k-mer size at read ends"),
        )
        .arg(
            Arg::new("hamming_distance")
                .long("hamming-distance")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Reference hamming distance (bbduk: hdist)"),
        )
        .arg(
            Arg::new("no_trim_by_overlap")
                .long("no-trim-by-overlap")
                .action(clap::ArgAction::SetTrue)
                .help("Disable overlap trimming (bbduk: tbo=f)"),
        )
        .arg(
            Arg::new("no_trim_pair_evenly")
                .long("no-trim-pair-evenly")
                .action(clap::ArgAction::SetTrue)
                .help("Disable even pair trimming (bbduk: tpe=f)"),
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
                .help("Quality trim mode: r, l, rl, w, or f"),
        )
        .arg(
            Arg::new("qtrim_window")
                .long("qtrim-window")
                .num_args(1)
                .default_value("4")
                .value_parser(value_parser!(usize))
                .help("Window size for qtrim=w (bbduk: qtrim=w,N)"),
        )
        .arg(
            Arg::new("trim_quality")
                .long("trim-quality")
                .num_args(1)
                .default_value("15")
                .value_parser(value_parser!(u8))
                .help("Quality threshold for qtrim (bbduk: trimq)"),
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
            Arg::new("max_ns")
                .long("max-ns")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(i64))
                .help("Maximum allowed N bases; negative disables"),
        )
        .arg(
            Arg::new("force_trim_mod")
                .long("force-trim-mod")
                .num_args(1)
                .default_value("5")
                .value_parser(value_parser!(usize))
                .help("Right-trim lengths to a multiple (bbduk: ftm)"),
        )
        .arg(
            Arg::new("force_trim_left")
                .long("force-trim-left")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim bases left of this position"),
        )
        .arg(
            Arg::new("force_trim_right")
                .long("force-trim-right")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim bases right of this position"),
        )
        .arg(
            Arg::new("force_trim_right2")
                .long("force-trim-right2")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim this many bases on the right end"),
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
                .help("Discard reads with a poly-G prefix"),
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
                .help("Discard reads with a poly-C prefix"),
        )
        .arg(
            Arg::new("max_non_poly")
                .long("max-non-poly")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Allowed non-polymer bases in a polymer run"),
        )
        .arg(
            Arg::new("min_avg_quality")
                .long("min-avg-quality")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(f64))
                .help("Discard reads with average quality below this (bbduk: maq)"),
        )
        .arg(
            Arg::new("min_avg_quality_bases")
                .long("min-avg-quality-bases")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Use only this many leading bases for min-avg-quality (bbduk: maqb)"),
        )
        .arg(
            Arg::new("min_base_quality")
                .long("min-base-quality")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(u8))
                .help("Discard reads with any base below this quality (bbduk: mbq)"),
        )
        .arg(
            Arg::new("max_n_rate")
                .long("max-n-rate")
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
                .help("Minimum read length as a fraction of the original (bbduk: mlf)"),
        )
        .arg(
            Arg::new("maxlength")
                .long("maxlength")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Discard reads longer than this"),
        )
        .arg(
            Arg::new("min_consecutive_bases")
                .long("min-consecutive-bases")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Discard reads without this many consecutive ACGT bases (bbduk: mcb)"),
        )
        .arg(
            Arg::new("min_gc")
                .long("min-gc")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(f64))
                .help("Discard reads with GC content below this"),
        )
        .arg(
            Arg::new("max_gc")
                .long("max-gc")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(f64))
                .help("Discard reads with GC content above this"),
        )
        .arg(
            Arg::new("no_pair_gc")
                .long("no-pair-gc")
                .action(clap::ArgAction::SetTrue)
                .help("Check GC per read instead of the pair average (bbduk: gcpairs=f)"),
        )
        .arg(
            Arg::new("mask_kmers")
                .long("mask-kmers")
                .num_args(1)
                .help("Mask matching k-mers: a symbol, 'lc', or 't' (bbduk: kmask)"),
        )
        .arg(
            Arg::new("mask_fully_covered")
                .long("mask-fully-covered")
                .action(clap::ArgAction::SetTrue)
                .help("Only mask bases fully covered by k-mers"),
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

/// Execute the clean command.
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
    let kmask = args.get_one::<String>("mask_kmers").map(String::as_str);
    let (kmask_symbol, kmask_lowercase) = match kmask {
        Some("lc") => (None, true),
        Some("t") | Some("true") => (Some(b'N'), false),
        Some(c) if c.len() == 1 => (Some(c.as_bytes()[0]), false),
        Some(other) => anyhow::bail!("invalid --mask-kmers value: {other}"),
        None => (None, false),
    };
    let maxlength = *args.get_one::<usize>("maxlength").unwrap();
    let opts = AdapterTrimOptions {
        k: *args.get_one::<usize>("k").unwrap(),
        mink: *args.get_one::<usize>("min_k").unwrap(),
        hdist: *args.get_one::<usize>("hamming_distance").unwrap(),
        ktrim_right: true,
        tbo: !args.get_flag("no_trim_by_overlap"),
        tpe: !args.get_flag("no_trim_pair_evenly"),
        qtrim_right,
        qtrim_left,
        qtrim_window,
        trimq: *args.get_one::<u8>("trim_quality").unwrap(),
        minlen: *args.get_one::<usize>("minlen").unwrap(),
        maxns: *args.get_one::<i64>("max_ns").unwrap(),
        ftm: *args.get_one::<usize>("force_trim_mod").unwrap(),
        force_trim_left: *args.get_one::<usize>("force_trim_left").unwrap(),
        force_trim_right: *args.get_one::<usize>("force_trim_right").unwrap(),
        force_trim_right2: *args.get_one::<usize>("force_trim_right2").unwrap(),
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
        min_avg_quality: *args.get_one::<f64>("min_avg_quality").unwrap(),
        min_avg_quality_bases: *args.get_one::<usize>("min_avg_quality_bases").unwrap(),
        min_base_quality: *args.get_one::<u8>("min_base_quality").unwrap(),
        max_n_rate: *args.get_one::<f64>("max_n_rate").unwrap(),
        min_len_fraction: *args.get_one::<f64>("minlen_fraction").unwrap(),
        max_length: if maxlength == 0 {
            usize::MAX
        } else {
            maxlength
        },
        min_consecutive_bases: *args.get_one::<usize>("min_consecutive_bases").unwrap(),
        min_gc: *args.get_one::<f64>("min_gc").unwrap(),
        max_gc: *args.get_one::<f64>("max_gc").unwrap(),
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
