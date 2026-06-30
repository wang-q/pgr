use clap::{Arg, ArgMatches, Command};
use pgr::libs::fmt::lav::{LavReader, LavStanza};

pub fn make_subcommand() -> Command {
    Command::new("to-psl")
        .about("Convert from lav to psl format")
        .after_help(
            r###"
Convert blastz lav to psl format.

Examples:
  # Convert lav to psl
  pgr lav to-psl in.lav -o out.psl
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input LAV file"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("target_strand")
                .long("target-strand")
                .value_name("STRAND")
                .help("Set the target strand (default is no strand)"),
        )
        .arg(
            Arg::new("score_file")
                .long("score-file")
                .value_name("FILE")
                .help("Output lav scores to side file"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let target_strand = args.get_one::<String>("target_strand");
    // let score_file = args.get_one::<String>("score_file"); // TODO: Implement score file output

    let reader = pgr::reader(input)?;
    let mut writer = pgr::writer(output)?;

    //----------------------------
    // Ops
    //----------------------------
    let mut lav_reader = LavReader::new(reader);

    // State
    let mut t_size = 0;
    let mut q_size = 0;
    let mut t_name = String::new();
    let mut q_name = String::new();
    let mut strand = String::from("+");

    while let Some(stanza) = lav_reader.next_stanza()? {
        match stanza {
            LavStanza::Sizes {
                t_size: t,
                q_size: q,
            } => {
                t_size = t as u32;
                q_size = q as u32;
            }
            LavStanza::Header {
                t_name: t,
                q_name: q,
                is_rc,
            } => {
                t_name = t;
                q_name = q;
                strand = if is_rc {
                    "-".to_string()
                } else {
                    "+".to_string()
                };
            }
            LavStanza::Alignment { blocks } => {
                if blocks.is_empty() {
                    continue;
                }

                let mut psl = pgr::libs::fmt::lav::blocks_to_psl(
                    &blocks, t_size, q_size, &t_name, &q_name, &strand,
                );

                if let Some(ts) = target_strand {
                    // Append target strand if provided
                    if psl.strand.len() == 1 {
                        psl.strand.push(ts.chars().next().unwrap_or('+'));
                    }
                }

                psl.write_to(&mut writer)?;
            }
            _ => {}
        }
    }

    Ok(())
}
