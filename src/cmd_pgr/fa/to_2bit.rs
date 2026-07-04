use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::fmt::twobit::TwoBitWriter;
use std::collections::HashSet;
/// Build the clap subcommand for to-2bit.
pub fn make_subcommand() -> Command {
    Command::new("to-2bit")
        .about("Converts FASTA to 2bit")
        .after_help(
            r###"
This command converts FASTA files to 2bit format.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
  # Convert FASTA to 2bit
  pgr fa to-2bit in.fa -o out.2bit
  pgr fa to-2bit in1.fa in2.fa -o out.2bit
  pgr fa to-2bit in.fa -o out.2bit --no-mask
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_help(
            "Input FASTA files",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(crate::cmd_pgr::args::no_mask_arg())
        .arg(
            Arg::new("strip_version")
                .long("strip-version")
                .action(ArgAction::SetTrue)
                .help("Strip version number from sequence names (e.g. NM_001.1 -> NM_001)"),
        )
        .arg(
            Arg::new("ignore_dups")
                .long("ignore-dups")
                .action(ArgAction::SetTrue)
                .help("Ignore duplicate sequence names (keep first)"),
        )
        .arg(crate::cmd_pgr::args::name_prefix_arg(None))
}
/// Execute the to-2bit command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles = args.get_many::<String>("infiles").unwrap();
    let output = crate::cmd_pgr::args::get_outfile(args);
    let no_mask = args.get_flag("no_mask");
    let strip_version = args.get_flag("strip_version");
    let ignore_dups = args.get_flag("ignore_dups");
    let name_prefix = args
        .get_one::<String>("name_prefix")
        .map(|s| s.as_str())
        .unwrap_or("");

    let mut seen_names = HashSet::new();
    let mut data = Vec::new();

    for infile in infiles {
        let reader = pgr::reader(infile)?;
        let reader = bio::io::fasta::Reader::new(reader);
        for result in reader.records() {
            let record = result?;
            let mut name = record.id().to_string();

            if strip_version {
                if let Some(idx) = name.rfind('.') {
                    if name[idx + 1..].chars().all(|c| c.is_ascii_digit()) {
                        name.truncate(idx);
                    }
                }
            }

            if !name_prefix.is_empty() {
                name = format!("{}{}", name_prefix, name);
            }

            if seen_names.contains(&name) {
                if ignore_dups {
                    continue;
                } else {
                    anyhow::bail!("Duplicate sequence name: {}", name);
                }
            }

            seen_names.insert(name.clone());

            let seq = std::str::from_utf8(record.seq())?.to_string();
            data.push((name, seq));
        }
    }

    let refs: Vec<(&str, &str)> = data.iter().map(|(n, s)| (n.as_str(), s.as_str())).collect();

    let mut writer = pgr::writer(output)?;
    let mut tb_writer = TwoBitWriter::new(&mut writer);

    tb_writer.write(&refs, !no_mask)?;

    Ok(())
}
