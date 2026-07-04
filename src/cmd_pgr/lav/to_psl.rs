use clap::{Arg, ArgMatches, Command};
use pgr::libs::fmt::lav::{LavReader, LavStanza};
/// Build the clap subcommand for to-psl.
pub fn make_subcommand() -> Command {
    Command::new("to-psl")
        .about("Converts from lav to psl format")
        .after_help(
            r###"
Convert blastz lav to psl format.

Examples:
1. Convert lav to psl:
   pgr lav to-psl in.lav -o out.psl
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input LAV file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("target_strand")
                .long("target-strand")
                .help("Set the target strand (default is no strand)"),
        )
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(clap::ArgAction::SetTrue)
                .help("Fail on unknown LAV stanzas instead of warning and skipping"),
        )
}
/// Execute the to-psl command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let target_strand = args.get_one::<String>("target_strand");
    let strict = args.get_flag("strict");

    let reader = pgr::reader(input)?;
    let mut writer = pgr::writer(output)?;

    let mut lav_reader = LavReader::new(reader);

    let mut t_size: Option<u32> = None;
    let mut q_size: Option<u32> = None;
    let mut t_name: Option<String> = None;
    let mut q_name: Option<String> = None;
    let mut strand: Option<String> = None;

    while let Some(stanza) = lav_reader.next_stanza()? {
        match stanza {
            LavStanza::Sizes {
                t_size: t,
                q_size: q,
            } => {
                t_size =
                    Some(u32::try_from(t).map_err(|_| anyhow::anyhow!("invalid t_size: {}", t))?);
                q_size =
                    Some(u32::try_from(q).map_err(|_| anyhow::anyhow!("invalid q_size: {}", q))?);
            }
            LavStanza::Header {
                t_name: t,
                q_name: q,
                is_rc,
            } => {
                t_name = Some(t);
                q_name = Some(q);
                strand = Some(if is_rc {
                    "-".to_string()
                } else {
                    "+".to_string()
                });
            }
            LavStanza::Alignment { blocks } => {
                if blocks.is_empty() {
                    continue;
                }

                let t_size = t_size.ok_or_else(|| {
                    anyhow::anyhow!("Alignment stanza encountered before Sizes stanza")
                })?;
                let q_size = q_size.ok_or_else(|| {
                    anyhow::anyhow!("Alignment stanza encountered before Sizes stanza")
                })?;
                let t_name = t_name.as_deref().ok_or_else(|| {
                    anyhow::anyhow!("Alignment stanza encountered before Header stanza")
                })?;
                let q_name = q_name.as_deref().ok_or_else(|| {
                    anyhow::anyhow!("Alignment stanza encountered before Header stanza")
                })?;
                let strand = strand.as_deref().ok_or_else(|| {
                    anyhow::anyhow!("Alignment stanza encountered before Header stanza")
                })?;

                let mut psl = pgr::libs::fmt::lav::blocks_to_psl(
                    &blocks, t_size, q_size, t_name, q_name, strand,
                );

                if let Some(ts) = target_strand {
                    // Append target strand if provided
                    if psl.strand.len() == 1 {
                        let ts_char = ts
                            .chars()
                            .next()
                            .ok_or_else(|| anyhow::anyhow!("--target-strand cannot be empty"))?;
                        psl.strand.push(ts_char);
                    }
                }

                psl.write_to(&mut writer)?;
            }
            other => {
                if strict {
                    anyhow::bail!("unknown lav stanza: {:?}", other);
                }
                log::warn!("skipping unknown lav stanza: {:?}", other);
            }
        }
    }

    Ok(())
}
