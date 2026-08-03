use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for multiz.
pub fn make_subcommand() -> Command {
    Command::new("multiz")
        .about("Merges block FA files using multiz-like DP on reference")
        .after_help(
            r###"
Merge multiple block FA files in the shared reference coordinate system using a multiz-like banded DP.

Notes:
* Takes two or more .fas inputs that share a reference name.
* Automatically derives windows from reference coverage with radius padding.
* Merges every window covered by at least one input and keeps the union of
  species across inputs.

Examples:
1. Merge with default radius:
   pgr fas multiz -r S288c tests/fas/S288cvsRM11_1a.slice.fas tests/fas/S288cvsSpar.slice.fas

2. Merge with a larger radius and minimum width:
   pgr fas multiz -r S288c --radius 30 --min-width 1000 tests/fas/S288cvsRM11_1a.slice.fas tests/fas/S288cvsYJM789.slice.fas tests/fas/S288cvsSpar.slice.fas

3. Write merged blocks to a file:
   pgr fas multiz -r S288c tests/fas/S288cvsRM11_1a.slice.fas tests/fas/S288cvsSpar.slice.fas -o merged.fas
"###,
        )
        .arg(
            Arg::new("ref_name")
                .short('r')
                .long("ref-name")
                .num_args(1)
                .required(true)
                .help("Reference sequence name present in all inputs"),
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input block FA file(s) to merge",
            2..,
        ))
        .arg(
            Arg::new("radius")
                .long("radius")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("30")
                .help("Banded DP radius around the reference diagonal"),
        )
        .arg(
            Arg::new("min_width")
                .long("min-width")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("1")
                .help("Minimum window width to consider for merging"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the multiz command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let ref_name = args.get_one::<String>("ref_name").unwrap().to_string();
    let radius = *args.get_one::<usize>("radius").unwrap();
    let min_width = *args.get_one::<usize>("min_width").unwrap();

    let cfg = pgr::libs::fas_multiz::FasMultizConfig {
        ref_name: ref_name.clone(),
        radius,
        min_width,
    };

    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();

    let blocks = pgr::libs::fas_multiz::merge_fas_files_auto_windows(&ref_name, &infiles, &cfg)?;

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    for block in blocks {
        for entry in &block.entries {
            writer.write_all(entry.to_string().as_ref())?;
        }
        writer.write_all(b"\n")?;
    }

    writer.flush()?;
    Ok(())
}
