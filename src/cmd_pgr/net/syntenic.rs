use clap::{ArgMatches, Command};
use pgr::libs::chain::net::{classify_syntenic, read_nets, write_net};

pub fn make_subcommand() -> Command {
    Command::new("syntenic")
        .about("Adds synteny info to net")
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input net file",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(crate::cmd_pgr::args::min_score_arg("0.0"))
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let in_file = matches.get_one::<String>("infile").unwrap();
    let out_file = crate::cmd_pgr::args::get_outfile(matches);
    let min_score = *matches.get_one::<f64>("min_score").unwrap();

    let reader = pgr::reader(in_file)?;
    let nets = read_nets(reader)?;

    classify_syntenic(&nets);

    let mut writer = pgr::writer(out_file)?;
    for net in &nets {
        write_net(net, &mut writer, false, min_score, 0)?;
    }

    Ok(())
}
