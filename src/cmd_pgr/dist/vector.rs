use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::BufRead;

use pgr::libs::clust::feature::FeatureVector;
use pgr::libs::linalg;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("vector")
        .about("Calculates similarity/distance between vectors")
        .after_help(
            r###"
This command calculates pairwise similarity/distance between vectors in input file(s).

modes:
    * euclidean distance
        * --mode euclid
    * euclidean distance to similarity
        * --mode euclid --sim
    * binary euclidean distance
        * --mode euclid --binary
    * binary euclidean distance to dissimilarity
        * --mode euclid --binary --sim --dis

    * cosine similarity, -1 -- 1
        * --mode cosine
    * cosine distance, 0 -- 2
        * --mode cosine --dis
    * binary cosine similarity
        * --mode cosine --binary
    * binary cosine similarity
        * --mode cosine --binary --dis

    * jaccard index
        * --mode jaccard --binary
    * weighted jaccard similarity
        * --mode jaccard

"###,
        )
        .arg(crate::cmd_pgr::args::pair_infiles_arg())
        .arg(crate::cmd_pgr::args::mode_arg(
            "euclid",
            &["euclid", "cosine", "jaccard"],
            "Mode of calculation",
        ))
        .arg(
            Arg::new("binary")
                .long("binary")
                .action(ArgAction::SetTrue)
                .help("Treat values in list as binary (0 or 1)"),
        )
        .arg(crate::cmd_pgr::args::sim_arg())
        .arg(
            Arg::new("dis")
                .long("dis")
                .action(ArgAction::SetTrue)
                .help("Convert to dissimilarity"),
        )
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let opt_mode = args.get_one::<String>("mode").unwrap();

    let is_bin = args.get_flag("binary");
    let is_sim = args.get_flag("sim");
    let is_dis = args.get_flag("dis");

    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    let infiles = crate::cmd_pgr::args::collect_infiles(args);

    let (sender, writer_thread) = pgr::libs::par::spawn_writer_and_pool(
        crate::cmd_pgr::args::get_outfile(args),
        opt_parallel,
    )?;

    //----------------------------
    // Ops
    //----------------------------
    let (entries1, entries2) =
        pgr::libs::par::load_two_sets(&infiles, false, |paths| load_file(&paths[0], is_bin))?;

    pgr::libs::par::par_run_pairs(
        &entries1,
        &entries2,
        &sender,
        |e1, e2| match linalg::vector_score(e1.list(), e2.list(), opt_mode, is_sim, is_dis) {
            Ok(score) => {
                let line = format!("{}\t{}\t{:.4}\n", e1.name(), e2.name(), score);
                Some(line)
            }
            Err(e) => {
                log::error!("{}", e);
                None
            }
        },
    );

    // Drop the sender to signal the writer thread to exit
    drop(sender);
    // Wait for the writer thread to finish
    writer_thread
        .join()
        .map_err(|_| anyhow::anyhow!("writer thread panicked"))?;

    Ok(())
}

fn load_file(infile: &str, is_bin: bool) -> anyhow::Result<Vec<FeatureVector>> {
    let mut entries = vec![];
    let reader = pgr::reader(infile)?;
    'LINE: for line in reader.lines().map_while(Result::ok) {
        let mut entry = FeatureVector::parse(&line);
        if entry.name().is_empty() {
            continue 'LINE;
        }
        if is_bin {
            let bin_list = entry
                .list()
                .iter()
                .map(|e| if *e > 0.0 { 1.0 } else { 0.0 })
                .collect::<Vec<f32>>();
            entry = FeatureVector::from(entry.name(), &bin_list);
        }
        entries.push(entry);
    }
    Ok(entries)
}
