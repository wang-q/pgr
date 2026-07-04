use clap::{Arg, ArgMatches, Command};
use pgr::libs::fmt::psl::parse_or_warn;
use std::io::BufRead;
/// Build the clap subcommand for to-chain.
pub fn make_subcommand() -> Command {
    Command::new("to-chain")
        .about("Converts PSL to Chain format")
        .after_help(
            r###"
Examples:
1. Convert PSL to Chain:
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
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        // Skip PSL header lines (psLayout version 3, column names, separator)
        if line.starts_with("psLayout") || line.starts_with("match") || line.starts_with("------") {
            continue;
        }

        let mut psl = match parse_or_warn(&line, strict)? {
            Some(p) => p,
            None => continue,
        };

        // Handle negative target strand
        let strand_bytes = psl.strand.as_bytes();
        if strand_bytes.len() < 2 {
            anyhow::bail!("malformed PSL strand (expected 2 chars): {}", psl.strand);
        }
        let t_strand_char = strand_bytes[1] as char;
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
