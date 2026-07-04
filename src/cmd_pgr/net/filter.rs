use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::chain::net::{filter_chrom, prune_gap, read_nets, FilterCriteria};
/// Build the clap subcommand for filter.
pub fn make_subcommand() -> Command {
    Command::new("filter")
        .about("Filters out parts of net")
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input net file (or stdin if 'stdin')",
        ))
        .arg(crate::cmd_pgr::args::min_score_arg_optional(
            "Restrict to those scoring at least N",
        ))
        .arg(crate::cmd_pgr::args::max_score_arg_optional(
            "Restrict to those scoring less than N",
        ))
        .arg(
            Arg::new("min_gap")
                .long("min-gap")
                .num_args(1)
                .value_parser(clap::value_parser!(u64))
                .help("Restrict to those with gap size (tSize) >= minSize"),
        )
        .arg(
            Arg::new("min_ali")
                .long("min-ali")
                .num_args(1)
                .value_parser(clap::value_parser!(u64))
                .help("Restrict to those with at least given bases aligning"),
        )
        .arg(
            Arg::new("max_ali")
                .long("max-ali")
                .num_args(1)
                .value_parser(clap::value_parser!(u64))
                .help("Restrict to those with at most given bases aligning"),
        )
        .arg(
            Arg::new("min_size_t")
                .long("min-size-t")
                .num_args(1)
                .value_parser(clap::value_parser!(u64))
                .help("Restrict to those at least this big on target"),
        )
        .arg(
            Arg::new("min_size_q")
                .long("min-size-q")
                .num_args(1)
                .value_parser(clap::value_parser!(u64))
                .help("Restrict to those at least this big on query"),
        )
        .arg(
            Arg::new("target_names")
                .long("target-names")
                .num_args(1)
                .help("Restrict target side sequence to those named (comma separated)"),
        )
        .arg(
            Arg::new("not_target_names")
                .long("not-target-names")
                .num_args(1)
                .help("Restrict target side sequence to those not named (comma separated)"),
        )
        .arg(
            Arg::new("query_names")
                .long("query-names")
                .num_args(1)
                .help("Restrict query side sequence to those named (comma separated)"),
        )
        .arg(
            Arg::new("not_query_names")
                .long("not-query-names")
                .num_args(1)
                .help("Restrict query side sequence to those not named (comma separated)"),
        )
        .arg(crate::cmd_pgr::args::net_type_arg(
            ArgAction::Append,
            "Restrict to given type, maybe repeated",
        ))
        .arg(crate::cmd_pgr::args::syn_arg(
            "Do filtering based on synteny (tuned for human/mouse)",
        ))
        .arg(
            Arg::new("nonsyn")
                .long("nonsyn")
                .action(ArgAction::SetTrue)
                .help("Do inverse filtering based on synteny"),
        )
        .arg(
            Arg::new("fill_only")
                .long("fill-only")
                .action(ArgAction::SetTrue)
                .help("Only pass fills, not gaps"),
        )
        .arg(
            Arg::new("gap_only")
                .long("gap-only")
                .action(ArgAction::SetTrue)
                .help("Only pass gaps, not fills"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the filter command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("infile").unwrap();

    if args.get_flag("fill_only") && args.get_flag("gap_only") {
        anyhow::bail!("--fill-only and --gap-only are mutually exclusive");
    }
    if args.get_flag("syn") && args.get_flag("nonsyn") {
        anyhow::bail!("--syn and --nonsyn are mutually exclusive");
    }

    let mut criteria = FilterCriteria::default();

    if let Some(v) = args.get_one::<f64>("min_score") {
        criteria.min_score = Some(*v);
    }
    if let Some(v) = args.get_one::<f64>("max_score") {
        criteria.max_score = Some(*v);
    }
    if let Some(v) = args.get_one::<u64>("min_gap") {
        criteria.min_gap = Some(*v);
    }
    if let Some(v) = args.get_one::<u64>("min_ali") {
        criteria.min_ali = Some(*v);
    }
    if let Some(v) = args.get_one::<u64>("max_ali") {
        criteria.max_ali = Some(*v);
    }
    if let Some(v) = args.get_one::<u64>("min_size_t") {
        criteria.min_size_t = Some(*v);
    }
    if let Some(v) = args.get_one::<u64>("min_size_q") {
        criteria.min_size_q = Some(*v);
    }

    if let Some(s) = args.get_one::<String>("target_names") {
        criteria.t_names = Some(s.split(',').map(|s| s.to_string()).collect());
    }
    if let Some(s) = args.get_one::<String>("not_target_names") {
        criteria.not_t_names = Some(s.split(',').map(|s| s.to_string()).collect());
    }
    if let Some(s) = args.get_one::<String>("query_names") {
        criteria.q_names = Some(s.split(',').map(|s| s.to_string()).collect());
    }
    if let Some(s) = args.get_one::<String>("not_query_names") {
        criteria.not_q_names = Some(s.split(',').map(|s| s.to_string()).collect());
    }
    if let Some(vals) = args.get_many::<String>("type") {
        criteria.types = Some(vals.map(|s| s.to_string()).collect());
    }

    criteria.do_syn = args.get_flag("syn");
    criteria.do_nonsyn = args.get_flag("nonsyn");
    criteria.fill_only = args.get_flag("fill_only");
    criteria.gap_only = args.get_flag("gap_only");

    let reader = pgr::reader(input_path)
        .with_context(|| format!("Failed to open reader for {}", input_path))?;
    let chroms = read_nets(reader)?;

    let out_path = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(out_path).with_context(|| format!("Failed to open writer for {}", out_path))?;

    for chrom in chroms {
        if !filter_chrom(&chrom, &criteria) {
            continue;
        }

        prune_gap(&chrom.root, &criteria);

        if !chrom.root.borrow().fills.is_empty() {
            chrom.write(&mut writer)?;
        }
    }

    Ok(())
}
