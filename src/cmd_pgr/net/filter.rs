use clap::{Arg, ArgAction, ArgMatches, Command};
use std::cell::{Ref, RefCell};
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use std::rc::Rc;

use pgr::libs::net::{read_nets, Chrom, Fill, Gap};

pub fn make_subcommand() -> Command {
    Command::new("filter")
        .about("Filter out parts of net")
        .arg(
            Arg::new("input")
                .index(1)
                .required(true)
                .help("Input net file (or stdin if '-')"),
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
}

struct FilterCriteria {
    min_score: Option<f64>,
    max_score: Option<f64>,
    min_gap: Option<u64>,
    min_ali: Option<u64>,
    max_ali: Option<u64>,
    min_size_t: Option<u64>,
    min_size_q: Option<u64>,
    t_names: Option<HashSet<String>>,
    not_t_names: Option<HashSet<String>>,
    q_names: Option<HashSet<String>>,
    not_q_names: Option<HashSet<String>>,
    types: Option<HashSet<String>>,

    // Synteny specific
    do_syn: bool,
    do_nonsyn: bool,
    min_top_score: f64,
    min_syn_score: f64,
    min_syn_size: f64, // tSize in original? "Min syntenic block size". Assuming tSize based on logic.
    min_syn_ali: u64,
    max_far: i64,

    fill_only: bool,
    gap_only: bool,
}

impl Default for FilterCriteria {
    fn default() -> Self {
        Self {
            min_score: None,
            max_score: None,
            min_gap: None,
            min_ali: None,
            max_ali: None,
            min_size_t: None,
            min_size_q: None,
            t_names: None,
            not_t_names: None,
            q_names: None,
            not_q_names: None,
            types: None,
            do_syn: false,
            do_nonsyn: false,
            min_top_score: 300000.0,
            min_syn_score: 200000.0,
            min_syn_size: 20000.0,
            min_syn_ali: 10000,
            max_far: 200000,
            fill_only: false,
            gap_only: false,
        }
    }
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

    let reader: Box<dyn io::BufRead> = if input_path == "-" {
        Box::new(BufReader::new(io::stdin()))
    } else {
        Box::new(BufReader::new(File::open(input_path)?))
    };

    let chroms = read_nets(reader)?;

    let mut writer = BufWriter::new(io::stdout());

    for chrom in chroms {
        if !filter_chrom(&chrom, &criteria) {
            continue;
        }

        // Prune the tree
        // We need to mutate the chrom. The read_nets returns owned Chroms.
        // We can modify the tree in place.

        prune_gap(&chrom.root, &criteria);

        // If after pruning, the root gap has no fills, we might still want to write the net header?
        // netFilter.c: "if ((net->fillList = cnPrune(net->fillList)) != NULL) chainNetWrite(net, f);"
        // So if the list is empty, it doesn't write anything.

        if !chrom.root.borrow().fills.is_empty() {
            chrom.write(&mut writer)?;
        }
    }

    Ok(())
}

fn filter_chrom(chrom: &Chrom, c: &FilterCriteria) -> bool {
    if let Some(names) = &c.t_names {
        if !names.contains(&chrom.name) {
            return false;
        }
    }
    if let Some(names) = &c.not_t_names {
        if names.contains(&chrom.name) {
            return false;
        }
    }
    true
}

fn syn_filter(fill: &Fill, c: &FilterCriteria) -> bool {
    if fill.class.is_empty() {
        // errAbort("No type field...");
        // For CLI tool, maybe return false or print warning?
        // UCSC aborts. Let's return false for safety/simplicity or assume not syntenic.
        return false;
    }
    let t_size = fill.end - fill.start;

    if fill.score >= c.min_syn_score
        && (t_size as f64) >= c.min_syn_size
        && fill.ali >= c.min_syn_ali
    {
        return true;
    }
    if fill.class == "top" {
        return fill.score >= c.min_top_score;
    }
    if fill.class == "nonSyn" {
        return false;
    }
    if fill.q_far.unwrap_or(0) > c.max_far {
        return false;
    }
    true
}

fn filter_one(fill: &Fill, c: &FilterCriteria) -> bool {
    if let Some(names) = &c.q_names {
        if !names.contains(&fill.o_chrom) {
            return false;
        }
    }
    if let Some(names) = &c.not_q_names {
        if names.contains(&fill.o_chrom) {
            return false;
        }
    }
    if let Some(types) = &c.types {
        if !types.contains(&fill.class) {
            return false;
        }
    }

    // In UCSC netFilter, if fill->chainId (it's a fill, not a gap wrapper?)
    // UCSC fills always have chainId? Wait.
    // In UCSC code: if (fill->chainId) { ... checks ... } else { ... checks for gap wrapper? ... }
    // In our Rust struct, Fill always has chain_id. Gap is separate.
    // So this is always true for Fill.

    if c.gap_only {
        return false;
    }

    if let Some(min_q) = c.min_size_q {
        let q_size = fill.o_end - fill.o_start;
        if q_size < min_q {
            return false;
        }
    }
    if let Some(min_t) = c.min_size_t {
        let t_size = fill.end - fill.start;
        if t_size < min_t {
            return false;
        }
    }

    if let Some(min_s) = c.min_score {
        if fill.score < min_s {
            return false;
        }
    }
    if let Some(max_s) = c.max_score {
        if fill.score > max_s {
            return false;
        }
    }

    if let Some(min_a) = c.min_ali {
        if fill.ali < min_a {
            return false;
        }
    }
    if let Some(max_a) = c.max_ali {
        if fill.ali > max_a {
            return false;
        }
    }

    // Skip range checks for now unless requested

    if c.do_syn && !syn_filter(fill, c) {
        return false;
    }
    if c.do_nonsyn && syn_filter(fill, c) {
        return false;
    }

    true
}

fn prune_gap(gap: &Rc<RefCell<Gap>>, c: &FilterCriteria) {
    let mut gap_mut = gap.borrow_mut();

    // Filter fills
    // We want to keep fills that pass filterOne
    // And if they pass, we recurse into their children

    let mut new_fills = Vec::new();

    for fill_rc in &gap_mut.fills {
        let keep = {
            let fill: Ref<Fill> = fill_rc.borrow();
            filter_one(&fill, c)
        };

        if keep {
            // Recurse
            prune_fill(fill_rc, c);
            new_fills.push(fill_rc.clone());
        }
    }

    gap_mut.fills = new_fills;
}

fn prune_fill(fill: &Rc<RefCell<Fill>>, c: &FilterCriteria) {
    let mut fill_mut = fill.borrow_mut();

    if c.fill_only {
        // If fillOnly is set, we don't want to pass gaps?
        // UCSC: if (fillOnly) return FALSE; (inside filterOne "else" block, which handles "gap" wrappers?)
        // Wait, UCSC net structure is Fill -> children (Fill list).
        // It seems UCSC's net structure is a bit different or I'm misinterpreting "gap" wrapper.
        // In UCSC, cnFill can be a "gap" if chainId is 0?
        // In our Rust struct, we have explicit Gap and Fill types.
        // So `fill_only` means we don't want children gaps?
        // But `Gap` in our struct is the container of children Fills.
        // If we remove Gaps, we remove all children Fills too.
        // Let's assume `fill_only` implies checking logic on Gaps.
    }

    // Filter gaps in this fill
    // Our Fill has `gaps: Vec<Rc<RefCell<Gap>>>`
    // We should filter these gaps.

    let mut new_gaps = Vec::new();
    for gap_rc in &fill_mut.gaps {
        // Check if gap passes criteria?
        // UCSC filterOne: else { if (fillOnly) return FALSE; if (fill->tSize < minGap) return FALSE; }
        // This "else" block corresponds to when `fill->chainId` is 0.
        // In our case, `Gap` corresponds to this.

        let keep = {
            let gap: Ref<Gap> = gap_rc.borrow();
            if c.fill_only {
                false
            } else if let Some(min_g) = c.min_gap {
                (gap.end - gap.start) >= min_g
            } else {
                true
            }
        };

        if keep {
            prune_gap(gap_rc, c);
            new_gaps.push(gap_rc.clone());
        }
    }
    fill_mut.gaps = new_gaps;
}
