//! `pgr runlist statop` — cross-set coverage stats.

use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;
use std::path::Path;

/// Build the clap subcommand for statop.
pub fn make_subcommand() -> Command {
    Command::new("statop")
        .about("Coverage on chromosomes for one runlist crossed another")
        .after_help(
            r###"
Prints CSV stats comparing `infile1` (possibly multi) against `infile2`
(single): `key,chr,chrLength,size,<base>Length,<base>Size,c1,c2,ratio`,
where `<base>` is the stem of `infile2` (or `--base`). `--all` keeps only
the whole-genome stats.

Examples:
1. Intersection stats:
   pgr runlist statop chr.sizes a.json b.json -o statop.csv
"###,
        )
        .arg(
            Arg::new("chr.sizes")
                .required(true)
                .index(1)
                .help("Sets the input file to use"),
        )
        .arg(
            Arg::new("infile1")
                .required(true)
                .index(2)
                .help("Sets the input file to use"),
        )
        .arg(
            Arg::new("infile2")
                .required(true)
                .index(3)
                .help("Sets the input file to use"),
        )
        .arg(
            Arg::new("all")
                .long("all")
                .action(ArgAction::SetTrue)
                .help("Only write whole genome stats"),
        )
        .arg(
            Arg::new("op")
                .long("op")
                .num_args(1)
                .default_value("intersect")
                .value_parser(["intersect", "union", "diff", "xor"])
                .help("operations: intersect, union, diff or xor"),
        )
        .arg(
            Arg::new("base")
                .long("base")
                .num_args(1)
                .help("basename of infile2"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the statop command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let sizes = pgr::read_sizes::<i32>(args.get_one::<String>("chr.sizes").unwrap())?;
    let json1 = pgr::libs::runlist::read_json(args.get_one::<String>("infile1").unwrap())?;
    let s1_of = pgr::libs::runlist::json_to_sets(&json1)?;
    let json2 = pgr::libs::runlist::read_json(args.get_one::<String>("infile2").unwrap())?;
    let s2 = pgr::libs::runlist::json_to_set(&json2)?;
    let is_multi = s1_of.len() > 1 || !s1_of.contains_key("__single__");
    let is_all = args.get_flag("all");
    let base = match args.get_one::<String>("base") {
        Some(b) => b.clone(),
        None => Path::new(args.get_one::<String>("infile2").unwrap())
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or_default()
            .to_string(),
    };
    let op = match args.get_one::<String>("op").unwrap().as_str() {
        "intersect" => pgr::libs::runlist::CompareOp::Intersect,
        "union" => pgr::libs::runlist::CompareOp::Union,
        "diff" => pgr::libs::runlist::CompareOp::Diff,
        "xor" => pgr::libs::runlist::CompareOp::Xor,
        _ => unreachable!("invalid statop op"),
    };

    // Apply the operation per set, filling missing chromosomes from sizes.
    let chrs: Vec<String> = sizes.keys().cloned().collect();
    let mut res_of: std::collections::BTreeMap<
        String,
        std::collections::BTreeMap<String, pgr::libs::ds::IntSpan>,
    > = std::collections::BTreeMap::new();
    let mut filled_of: std::collections::BTreeMap<
        String,
        std::collections::BTreeMap<String, pgr::libs::ds::IntSpan>,
    > = std::collections::BTreeMap::new();
    for (name, s1) in &s1_of {
        // Fill missing chromosomes with empty sets so per-chr rows and the
        // whole-genome row cover every chromosome (spanr `fill_up_m` parity).
        let mut filled: std::collections::BTreeMap<String, pgr::libs::ds::IntSpan> = s1.clone();
        for chr in &chrs {
            filled.entry(chr.clone()).or_default();
        }
        let mut set: std::collections::BTreeMap<String, pgr::libs::ds::IntSpan> =
            std::collections::BTreeMap::new();
        for chr in &chrs {
            let a = filled[chr].clone();
            let b = s2.get(chr).cloned().unwrap_or_default();
            set.insert(
                chr.clone(),
                match op {
                    pgr::libs::runlist::CompareOp::Intersect => a.intersect(&b),
                    pgr::libs::runlist::CompareOp::Union => a.union(&b),
                    pgr::libs::runlist::CompareOp::Diff => a.diff(&b),
                    pgr::libs::runlist::CompareOp::Xor => a.xor(&b),
                },
            );
        }
        res_of.insert(name.clone(), set);
        filled_of.insert(name.clone(), filled);
    }

    let mut lines: Vec<String> = Vec::new();
    let mut header = format!(
        "key,chr,chrLength,size,{}Length,{}Size,c1,c2,ratio",
        base, base
    );
    if is_multi {
        if is_all {
            header = header.replacen("chr,", "", 1);
        }
        lines.push(header);
        for (name, s1) in &filled_of {
            lines.push(pgr::libs::runlist::statop_lines(
                s1,
                &sizes,
                &s2,
                &res_of[name],
                is_all,
                Some(name),
            )?);
        }
    } else {
        header = header.replacen("key,", "", 1);
        if is_all {
            header = header.replacen("chr,", "", 1);
        }
        lines.push(header);
        lines.push(pgr::libs::runlist::statop_lines(
            &filled_of["__single__"],
            &sizes,
            &s2,
            &res_of["__single__"],
            is_all,
            None,
        )?);
    }
    let mut w = pgr::writer(outfile)?;
    for line in &lines {
        writeln!(w, "{}", line)?;
    }
    w.flush()?;
    Ok(())
}
