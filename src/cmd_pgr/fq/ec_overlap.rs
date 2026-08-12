use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::fq::bbnet::CellNet;
use pgr::libs::fq::merge::{merge, write_ihist, MergeOptions, Preset};
use std::io::Write;

/// Build the clap subcommand for ec-overlap.
pub fn make_subcommand() -> Command {
    Command::new("ec-overlap")
        .about("Error-corrects paired reads by overlap without joining (bbmerge ecco)")
        .after_help(
            r###"
This command error-corrects paired-end reads using the evidence of their
overlapping region, without joining the pair, reproducing the BBTools
`bbmerge.sh ... ecco` mode (anchr merge phase 1). It is the overlap-based
counterpart of `pgr fq ec-kmer` (tadpole ecc): this command needs paired
reads with a true overlap, while `ec-kmer` corrects from the k-mer graph.

Notes:
* Input is 1 interleaved FASTQ file or 2 paired files (R1, R2)
* Corrected pairs are written to the output by default (bbmerge `mix`);
  pass `--no-mix` to send only the corrected pairs to the output and the
  untouched pairs to `--outu`
* `--strict`/`--vstrict` apply the bbmerge strict/vstrict parameter sets;
  explicit options override the preset values
* `--ihist` writes the insert-size histogram in the bbmerge `ihist` format
* By default the BBMerge overlap net (bbmerge.bbnet) filters merges, so
  `--net FILE` is required unless `--no-make-vector` is given
* Processing is ordered and deterministic (equivalent to `threads=1`)
* Supports both plain text and gzipped (.gz) files

Examples:
1. Error-correct by overlap, keeping all pairs (anchr merge phase 1):
   pgr fq ec-overlap R1.fq.gz R2.fq.gz -o ecco.fq.gz --vstrict \
       --net bbmerge.bbnet --ihist ihist.merge1.txt

2. Only corrected pairs to the output, the rest to outu:
   pgr fq ec-overlap in.fq.gz -o ecco.fq.gz --outu rest.fq.gz \
       --no-mix --no-make-vector
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input FASTQ file(s): 1 interleaved or 2 paired (R1, R2)",
            1..=2,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("outu")
                .long("outu")
                .num_args(1)
                .help("Output file for untouched pairs (when --no-mix)"),
        )
        .arg(
            Arg::new("ihist")
                .long("ihist")
                .num_args(1)
                .help("Write the insert-size histogram to this file"),
        )
        .arg(
            Arg::new("no_mix")
                .long("no-mix")
                .action(clap::ArgAction::SetTrue)
                .help("Send only corrected pairs to the output (bbmerge: mix=f)"),
        )
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(clap::ArgAction::SetTrue)
                .help("Apply the bbmerge strict parameter set"),
        )
        .arg(
            Arg::new("vstrict")
                .long("vstrict")
                .action(clap::ArgAction::SetTrue)
                .help("Apply the bbmerge vstrict parameter set"),
        )
        .arg(
            Arg::new("min_overlap")
                .long("min-overlap")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimum overlap (bbmerge: minoverlap)"),
        )
        .arg(
            Arg::new("min_overlap0")
                .long("min-overlap0")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimum overlap for pre-screening (bbmerge: minoverlap0)"),
        )
        .arg(
            Arg::new("min_insert")
                .long("min-insert")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimum insert size (bbmerge: mininsert)"),
        )
        .arg(
            Arg::new("min_insert0")
                .long("min-insert0")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimum insert size for pre-screening (bbmerge: mininsert0)"),
        )
        .arg(
            Arg::new("max_ratio")
                .long("max-ratio")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help("Maximum error ratio (bbmerge: maxratio)"),
        )
        .arg(
            Arg::new("ratio_margin")
                .long("ratio-margin")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help("Ratio margin (bbmerge: ratiomargin)"),
        )
        .arg(
            Arg::new("ratio_offset")
                .long("ratio-offset")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help("Ratio offset (bbmerge: ratiooffset)"),
        )
        .arg(
            Arg::new("min_second_ratio")
                .long("min-second-ratio")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help("Minimum ratio for the second-best overlap (bbmerge: minsecondratio)"),
        )
        .arg(
            Arg::new("ratio_reduction")
                .long("ratio-reduction")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Overlap reduction for ratio mode (bbmerge: ratiominoverlapreduction)"),
        )
        .arg(
            Arg::new("min_entropy")
                .long("min-entropy")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimum entropy score (bbmerge: minentropy)"),
        )
        .arg(
            Arg::new("efilter")
                .long("efilter")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help("Expected-error filter ratio; 0 disables it"),
        )
        .arg(
            Arg::new("pfilter")
                .long("pfilter")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help("Probability filter; 0 disables it"),
        )
        .arg(
            Arg::new("no_make_vector")
                .long("no-make-vector")
                .action(clap::ArgAction::SetTrue)
                .help("Disable the BBMerge MAKE_VECTOR behavior (ratio maxratio 0.7)"),
        )
        .arg(
            Arg::new("net")
                .long("net")
                .num_args(1)
                .help("BBMerge overlap-filter net file (bbmerge.bbnet)"),
        )
}

/// Execute the ec-overlap command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let outu = args.get_one::<String>("outu").cloned();
    let ihist = args.get_one::<String>("ihist").cloned();

    if args.get_flag("strict") && args.get_flag("vstrict") {
        anyhow::bail!("--strict and --vstrict are mutually exclusive");
    }
    let preset = if args.get_flag("vstrict") {
        Preset::VStrict
    } else if args.get_flag("strict") {
        Preset::Strict
    } else {
        Preset::Normal
    };
    let mut opts = MergeOptions::from_preset(preset);
    opts.ecco = true;
    opts.mix = !args.get_flag("no_mix");
    if let Some(x) = args.get_one::<usize>("min_overlap") {
        opts.min_overlap = *x;
    }
    if let Some(x) = args.get_one::<usize>("min_overlap0") {
        opts.min_overlap0 = *x;
    }
    if let Some(x) = args.get_one::<usize>("min_insert") {
        opts.min_insert = *x;
    }
    if let Some(x) = args.get_one::<usize>("min_insert0") {
        opts.min_insert0 = Some(*x);
    }
    if let Some(x) = args.get_one::<f32>("max_ratio") {
        opts.max_ratio = *x;
    }
    if let Some(x) = args.get_one::<f32>("ratio_margin") {
        opts.ratio_margin = *x;
    }
    if let Some(x) = args.get_one::<f32>("ratio_offset") {
        opts.ratio_offset = *x;
    }
    if let Some(x) = args.get_one::<f32>("min_second_ratio") {
        opts.min_second_ratio = *x;
    }
    if let Some(x) = args.get_one::<usize>("ratio_reduction") {
        opts.ratio_reduction = *x;
    }
    if let Some(x) = args.get_one::<usize>("min_entropy") {
        opts.min_entropy = *x;
    }
    if let Some(x) = args.get_one::<f32>("efilter") {
        // `0` disables the filter, matching `fq merge` and the help text.
        opts.efilter = (*x > 0.0).then_some(*x);
    }
    if let Some(x) = args.get_one::<f32>("pfilter") {
        opts.pfilter = *x;
    }
    if args.get_flag("no_make_vector") {
        opts.make_vector = false;
    }
    if let Some(net) = args.get_one::<String>("net") {
        opts.net = Some(CellNet::load(net).context("failed to load overlap net")?);
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, infiles.iter().map(String::as_str))?;
    if let Some(p) = &outu {
        crate::cmd_pgr::args::ensure_outfile_distinct(p, infiles.iter().map(String::as_str))?;
    }
    if let Some(p) = &ihist {
        crate::cmd_pgr::args::ensure_outfile_distinct(p, infiles.iter().map(String::as_str))?;
    }

    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    let mut outu_w = outu
        .as_ref()
        .map(|p| pgr::libs::io::writer(p))
        .transpose()?;
    let stats = merge(&infiles, &mut out, outu_w.as_mut(), &opts)?;
    out.flush()?;
    if let Some(path) = &ihist {
        let mut w =
            pgr::libs::io::writer(path).with_context(|| format!("failed to open ihist {path}"))?;
        write_ihist(&mut w, &stats)?;
        w.flush()?;
    }
    eprintln!(
        "Pairs: {}  Joined: {}  Bases corrected: {}",
        stats.pairs, stats.joined, stats.errors_corrected
    );
    Ok(())
}
