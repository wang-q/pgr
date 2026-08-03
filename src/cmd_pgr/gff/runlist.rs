//! `pgr gff runlist` — convert GFF to a per-chromosome runlist JSON.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for runlist.
pub fn make_subcommand() -> Command {
    Command::new("runlist")
        .about("Converts GFF files to per-chromosome runlists")
        .after_help(
            r###"
Reads GFF records (start/end 1-based inclusive) and writes the per-chromosome
union as a runlist JSON ready for `pgr fa mask`. `--tag` restricts to one
feature type (third column). Migrated from the external `spanr gff` command.

Examples:
1. All features:
   pgr gff runlist in.gff -o out.json
2. Only genes:
   pgr gff runlist in.gff --tag gene -o out.json
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("tag")
                .long("tag")
                .num_args(1)
                .help("primary tag (the third field)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the runlist command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let tag = args
        .get_one::<String>("tag")
        .map(|s| s.as_str())
        .unwrap_or("");
    let mut set: std::collections::BTreeMap<String, pgr::libs::ds::IntSpan> =
        std::collections::BTreeMap::new();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile)?;
        for record in pgr::libs::fmt::gff::read_records(reader)? {
            if !tag.is_empty() && record.ty != tag {
                continue;
            }
            set.entry(record.seqid)
                .or_default()
                .add_pair(record.start as i32, record.end as i32);
        }
    }
    let json = pgr::libs::ds::intspan::set2json(&set);
    pgr::libs::ds::intspan::write_json(outfile, &json)?;
    Ok(())
}
