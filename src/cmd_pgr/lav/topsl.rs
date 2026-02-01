use pgr::libs::lav::{LavReader, LavStanza, blocks_to_psl};
use clap::{Arg, ArgMatches, Command};

pub fn make_subcommand() -> Command {
    Command::new("topsl")
        .about("Convert LAV to PSL format")
        .after_help(
            r###"
* <input> is the path to a LAV file, .lav.gz is supported
    * input == stdin means reading from STDIN

"###,
        )
        .arg(
            Arg::new("input")
                .index(1)
                .default_value("stdin")
                .help("Input LAV file (or stdin if not specified)")
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Output PSL file (or stdout if not specified)")
                .num_args(1)
                .default_value("stdout")
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let input = args.get_one::<String>("input").unwrap();
    let output = args.get_one::<String>("output").unwrap();

    let reader = intspan::reader(input);
    let mut writer = intspan::writer(output);

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
            LavStanza::Sizes { t_size: t, q_size: q } => {
                t_size = t as u32;
                q_size = q as u32;
            }
            LavStanza::Header { t_name: t, q_name: q, is_rc } => {
                t_name = t;
                q_name = q;
                strand = if is_rc { "-".to_string() } else { "+".to_string() };
            }
            LavStanza::Alignment { blocks } => {
                if blocks.is_empty() { continue; }
                
                let psl = blocks_to_psl(&blocks, t_size, q_size, &t_name, &q_name, &strand);
                psl.write_to(&mut writer)?;
            }
            _ => {}
        }
    }

    Ok(())
}


