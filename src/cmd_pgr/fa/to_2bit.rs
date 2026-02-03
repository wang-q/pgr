use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::twobit::TwoBitWriter;
use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;

pub fn make_subcommand() -> Command {
    Command::new("to-2bit")
        .about("Convert FASTA to 2bit")
        .after_help(
            r###"
Examples:
  # Convert FASTA to 2bit
  pgr fa to-2bit in.fa -o out.2bit
  pgr fa to-2bit in1.fa in2.fa -o out.2bit
  pgr fa to-2bit in.fa -o out.2bit --no-mask
"###,
        )
        .arg(
            Arg::new("infiles")
                .help("Input FASTA files")
                .required(true)
                .num_args(1..)
                .index(1),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output 2bit file")
                .required(true),
        )
        .arg(
            Arg::new("no_mask")
                .long("no-mask")
                .action(ArgAction::SetTrue)
                .help("Do not create mask blocks from lowercase letters"),
        )
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
        .arg(
            Arg::new("name_prefix")
                .long("name-prefix")
                .value_name("STR")
                .help("Add prefix to sequence names"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles = args.get_many::<String>("infiles").unwrap();
    let output = args.get_one::<String>("output").unwrap();
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
        let reader = bio::io::fasta::Reader::from_file(infile)?;
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

    let f = File::create(output)?;
    let mut writer = BufWriter::new(f);
    let mut tb_writer = TwoBitWriter::new(&mut writer);

    tb_writer.write(&refs, !no_mask)?;

    Ok(())
}
