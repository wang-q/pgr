use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::fq::bbnet::CellNet;
use pgr::libs::fq::merge::{merge, write_ihist, MergeOptions, Preset};
use std::io::Write;

/// Build the clap subcommand for merge.
pub fn make_subcommand() -> Command {
    Command::new("merge")
        .about("Merges overlapping paired-end reads (bbmerge-compatible)")
        .after_help(
            r###"
This command merges overlapping paired-end reads into single reads and/or
error-corrects them by overlap, reproducing the BBTools `bbmerge.sh` /
`bbmerge-auto.sh` overlap pipeline. `--ecco` corrects pairs without joining
them (anchr merge phase 1); the default mode joins overlapping pairs and
writes the unmerged pairs to `--outu` (anchr merge phase 4). `--extend2`
with `--rem` reproduces the bbmerge-auto `extend2=N rem` mode: unmerged
pairs are extended along a k-mer graph (k=81) and the overlap is retried,
requiring the extended overlap to match the unextended one.

Notes:
* Input is 1 interleaved FASTQ file or 2 paired files (R1, R2)
* `--strict`/`--vstrict` apply the bbmerge strict/vstrict parameter sets;
  explicit options override the preset values
* `--ecco` defaults to `--mix` (all reads are written to the output, like
  `bbmerge.sh ... ecco`); pass `--no-mix` to send only the corrected pairs
  to the main output and the rest to `--outu`
* `--ihist` writes the insert-size histogram in the bbmerge `ihist` format
* By default the BBMerge overlap net (bbmerge.bbnet) filters merges, so
  `--net FILE` is required unless `--no-make-vector` is given
* Processing is ordered and deterministic (equivalent to `threads=1`)
* Supports both plain text and gzipped (.gz) files

Examples:
1. Error-correct by overlap, keeping all pairs (anchr phase 1):
   pgr fq merge R1.fq.gz R2.fq.gz -o ecco.fq.gz --ecco --mix --vstrict \
       --net bbmerge.bbnet --ihist ihist.merge1.txt

2. Merge overlapping pairs, unmerged to outu (anchr phase 4):
   pgr fq merge in.fq.gz -o merged.fq.gz --outu unmerged.fq.gz \
       --strict --no-make-vector --ihist ihist.merge.txt

3. Tune the overlap parameters explicitly:
   pgr fq merge R1.fq R2.fq -o out.fq --min-overlap 11 --max-ratio 0.075

4. Merge with tadpole extension retry (anchr merge phase 4):
   pgr fq merge in.fq.gz -o merged.fq.gz --outu unmerged.fq.gz \
       --strict --no-make-vector --extend2 80 --rem --ihist ihist.merge.txt
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
                .help("Output file for unmerged read pairs"),
        )
        .arg(
            Arg::new("ihist")
                .long("ihist")
                .num_args(1)
                .help("Write the insert-size histogram to this file"),
        )
        .arg(
            Arg::new("ecco")
                .long("ecco")
                .action(clap::ArgAction::SetTrue)
                .help("Error-correct pairs by overlap without joining"),
        )
        .arg(
            Arg::new("mix")
                .long("mix")
                .action(clap::ArgAction::SetTrue)
                .help("Also write unmerged pairs to the main output (bbmerge: mix)"),
        )
        .arg(
            Arg::new("no_mix")
                .long("no-mix")
                .action(clap::ArgAction::SetTrue)
                .help("Do not auto-mix when --ecco is set (bbmerge: mix=f)"),
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
        .arg(
            Arg::new("extend2")
                .long("extend2")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Extend unmerged pairs by up to this many bases and retry (bbmerge-auto: extend2)"),
        )
        .arg(
            Arg::new("rem")
                .long("rem")
                .action(clap::ArgAction::SetTrue)
                .help("Require the extended overlap to match the unextended one (bbmerge-auto: rem)"),
        )
}

/// Execute the merge command.
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
        opts.efilter = (*x > 0.0).then_some(*x);
    }
    if let Some(x) = args.get_one::<f32>("pfilter") {
        opts.pfilter = *x;
    }
    opts.ecco = args.get_flag("ecco");
    opts.mix = args.get_flag("mix") || (opts.ecco && !args.get_flag("no_mix"));
    if args.get_flag("no_make_vector") {
        opts.make_vector = false;
    }
    if let Some(path) = args.get_one::<String>("net") {
        opts.net = Some(CellNet::load(path)?);
    }
    if let Some(x) = args.get_one::<usize>("extend2") {
        opts.extend2 = *x;
    }
    opts.rem = args.get_flag("rem");
    if opts.extend2 > 0 {
        // bbmerge-auto forces MAKE_VECTOR=false whenever a tadpole (extend2 /
        // eccTadpole / kfilter) is active.
        opts.make_vector = false;
    }
    if opts.make_vector && opts.net.is_none() {
        anyhow::bail!(
            "make-vector mode (the default) requires --net with a bbmerge.bbnet \
             file; use --no-make-vector for the classic overlap filters"
        );
    }

    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    let mut outu_w = match &outu {
        Some(path) => Some(
            pgr::libs::io::writer(path).with_context(|| format!("failed to open output {path}"))?,
        ),
        None => None,
    };
    let stats = merge(&infiles, &mut out, outu_w.as_mut(), &opts)?;
    out.flush()?;
    if let Some(path) = &ihist {
        let mut w =
            pgr::libs::io::writer(path).with_context(|| format!("failed to open output {path}"))?;
        write_ihist(&mut w, &stats)?;
        w.flush()?;
    }
    Ok(())
}
