use clap::Command;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("hv")
        .about("Estimate distances between DNA/protein files using hypervectors")
        .after_help(
            r###"
This command calculates pairwise distances between files in FA file(s) using minimizers and hypervectors.

* The outputs are printed to stdout in the following format:
    <file1> <file2> <total1> <total2> <inter> <union> <mash_distance> <jaccard_index> <containment_index>

* Minimizers and Hash Algorithms are the same as `pgr dist seq`

* Input Modes:
    * For a single sequence file: Merge all sequences within the file into a single hypervector.
      Note that comparing this set to itself (self-comparison) is not meaningful,
      as the distance will always be 0 and the similarity will always be 1.
    * For two sequence files: Merge all sequences within each file into a single hypervector,
      and calculate distances between the two hypervectors.
    * When --list is set:
      - For each file listed in the list file, merge all sequences within that file
        into a single hypervector, and calculate distances between these hypervectors.
      - The merging does not span across multiple files listed in the list file.

Examples:
1. Merge all sequences in a file and compare to another:
   pgr dist hv file1.fa file2.fa

2. Use Mod-Minimizer for DNA sequences (canonical k-mers):
   pgr dist hv file1.fa file2.fa --hasher mod -k 21 -w 5

3. Treat input as a list file and calculate distances:
   pgr dist hv list.txt --list

4. Use 4 threads for parallel processing:
   pgr dist hv input.fa --parallel 4

5. Perform six-frame translation on a FA file and match to another
    pgr fa six-frame input.fa |
        pgr dist hv stdin match.fa

"###,
        )
        .arg(crate::cmd_pgr::args::pair_infiles_arg())
        .arg(crate::cmd_pgr::args::hasher_arg())
        .arg(crate::cmd_pgr::args::kmer_arg())
        .arg(crate::cmd_pgr::args::window_arg())
        .arg(
            clap::Arg::new("dim")
                .long("dim")
                .short('d')
                .num_args(1)
                .default_value("4096")
                .value_parser(clap::value_parser!(usize))
                .help("The dimension size should be a multiple of 32."),
        )
        .arg(crate::cmd_pgr::args::sim_arg())
        .arg(crate::cmd_pgr::args::list_arg())
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &clap::ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let opt_hasher = args.get_one::<String>("hasher").unwrap();
    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_window = *args.get_one::<usize>("window").unwrap();
    let opt_dim = *args.get_one::<usize>("dim").unwrap();

    let is_sim = args.get_flag("sim");
    let is_list = args.get_flag("list_files");
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    let infiles = crate::cmd_pgr::args::collect_infiles(args);

    let (sender, writer_thread) = pgr::libs::par::spawn_writer_and_pool(
        crate::cmd_pgr::args::get_outfile(args),
        opt_parallel,
    )?;

    //----------------------------
    // Ops
    //----------------------------
    let (entries1, entries2) = pgr::libs::par::load_two_sets(&infiles, is_list, |paths| {
        pgr::libs::par::load_entries(paths, |p| {
            let entry =
                pgr::libs::hv::load_hv_from_fasta(p, opt_hasher, opt_kmer, opt_window, opt_dim)?;
            Ok(vec![entry])
        })
    })?;

    pgr::libs::par::par_run_pairs(&entries1, &entries2, &sender, |e1, e2| {
        let d = pgr::libs::hv::calc_distances(&e1.set, &e2.set, opt_kmer);

        let dist = if is_sim { 1.0 - d.mash } else { d.mash };

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
            e1.name, e2.name, d.card1, d.card2, d.inter, d.union, dist, d.jaccard, d.containment
        );
        Some(line)
    });

    // Drop the sender to signal the writer thread to exit
    drop(sender);
    // Wait for the writer thread to finish
    writer_thread.join().unwrap();

    Ok(())
}
