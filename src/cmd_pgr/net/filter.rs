use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::chain::net::{filter_chrom, prune_gap, read_nets, FilterCriteria};

pub fn make_subcommand() -> Command {
    Command::new("filter")
        .about("Filter out parts of net")
        .arg(
            Arg::new("input")
                .index(1)
                .required(true)
                .help("Input net file (or stdin if 'stdin')"),
        )
        .arg(
            Arg::new("min_score")
                .long("min-score")
                .value_parser(clap::value_parser!(f64))
                .help("Restrict to those scoring at least N"),
        )
        .arg(
            Arg::new("max_score")
                .long("max-score")
                .value_parser(clap::value_parser!(f64))
                .help("Restrict to those scoring less than N"),
        )
        .arg(
            Arg::new("min_gap")
                .long("min-gap")
                .value_parser(clap::value_parser!(u64))
                .help("Restrict to those with gap size (tSize) >= minSize"),
        )
        .arg(
            Arg::new("min_ali")
                .long("min-ali")
                .value_parser(clap::value_parser!(u64))
                .help("Restrict to those with at least given bases aligning"),
        )
        .arg(
            Arg::new("max_ali")
                .long("max-ali")
                .value_parser(clap::value_parser!(u64))
                .help("Restrict to those with at most given bases aligning"),
        )
        .arg(
            Arg::new("min_size_t")
                .long("min-size-t")
                .value_parser(clap::value_parser!(u64))
                .help("Restrict to those at least this big on target"),
        )
        .arg(
            Arg::new("min_size_q")
                .long("min-size-q")
                .value_parser(clap::value_parser!(u64))
                .help("Restrict to those at least this big on query"),
        )
        .arg(
            Arg::new("t")
                .long("t")
                .help("Restrict target side sequence to those named (comma separated)"),
        )
        .arg(
            Arg::new("not_t")
                .long("not-t")
                .help("Restrict target side sequence to those not named (comma separated)"),
        )
        .arg(
            Arg::new("q")
                .long("q")
                .help("Restrict query side sequence to those named (comma separated)"),
        )
        .arg(
            Arg::new("not_q")
                .long("not-q")
                .help("Restrict query side sequence to those not named (comma separated)"),
        )
        .arg(
            Arg::new("type")
                .long("type")
                .action(ArgAction::Append)
                .help("Restrict to given type, maybe repeated"),
        )
        .arg(
            Arg::new("syn")
                .long("syn")
                .action(ArgAction::SetTrue)
                .help("Do filtering based on synteny (tuned for human/mouse)"),
        )
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

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("input").unwrap();

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

    if let Some(s) = args.get_one::<String>("t") {
        criteria.t_names = Some(s.split(',').map(|s| s.to_string()).collect());
    }
    if let Some(s) = args.get_one::<String>("not_t") {
        criteria.not_t_names = Some(s.split(',').map(|s| s.to_string()).collect());
    }
    if let Some(s) = args.get_one::<String>("q") {
        criteria.q_names = Some(s.split(',').map(|s| s.to_string()).collect());
    }
    if let Some(s) = args.get_one::<String>("not_q") {
        criteria.not_q_names = Some(s.split(',').map(|s| s.to_string()).collect());
    }
    if let Some(vals) = args.get_many::<String>("type") {
        criteria.types = Some(vals.map(|s| s.to_string()).collect());
    }

    criteria.do_syn = args.get_flag("syn");
    criteria.do_nonsyn = args.get_flag("nonsyn");
    criteria.fill_only = args.get_flag("fill_only");
    criteria.gap_only = args.get_flag("gap_only");

    let reader = pgr::reader(input_path)?;
    let chroms = read_nets(reader)?;

    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

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
