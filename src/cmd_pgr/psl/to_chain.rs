use clap::{Arg, ArgMatches, Command};
use pgr::libs::fmt::psl::Psl;
use std::io::BufRead;
use std::str::FromStr;
/// Build the clap subcommand for to-chain.
pub fn make_subcommand() -> Command {
    Command::new("to-chain")
        .about("Converts PSL to Chain format")
        .after_help(
            r###"
Examples:
  # Convert PSL to Chain
  pgr psl to-chain in.psl -o out.chain
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input PSL file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("fix_strand")
                .long("fix-strand")
                .action(clap::ArgAction::SetTrue)
                .help("Fix '-' target strand by reverse complementing the record"),
        )
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(clap::ArgAction::SetTrue)
                .help("Fail on parse errors instead of skipping malformed lines"),
        )
}
/// Execute the to-chain command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let fix_strand = args.get_flag("fix_strand");
    let strict = args.get_flag("strict");

    let reader = pgr::reader(input)?;
    let mut writer = pgr::writer(output)?;

    let mut chain_id = 1;

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut psl = match Psl::from_str(&line) {
            Ok(p) => p,
            Err(e) => {
                if strict {
                    anyhow::bail!("failed to parse psl line: {}: {}", line, e);
                }
                log::warn!("skipping unparseable psl line: {}: {}", line, e);
                continue;
            }
        };

        // Handle negative target strand
        let t_strand_char = psl.strand.chars().nth(1).unwrap_or('+');
        if t_strand_char == '-' {
            if fix_strand {
                psl.rc();
            } else {
                // In strict mode we might abort, but for now maybe just warn or skip?
                // UCSC pslToChain aborts by default.
                // Let's abort to match behavior, or maybe just skip?
                // "errAbort" in C.
                anyhow::bail!("PSL record has '-' for target strand. Use --fix-strand to fix.");
            }
        }

        psl.write_chain(&mut writer, chain_id)?;

        chain_id += 1;
    }

    Ok(())
}
