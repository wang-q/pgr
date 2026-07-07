use clap::{ArgMatches, Command};

/// Build the clap subcommand for seq.
pub fn make_subcommand() -> Command {
    Command::new("seq")
        .about("Estimates sequence distances using minimizers")
        .after_help(
            r###"
This command calculates pairwise distances between sequences in FA file(s) using minimizers.

* The outputs are printed to stdout in the following format:
    <sequence1> <sequence2> <mash_distance> <jaccard_index> <containment_index>
* With --merge
    <file1> <file2> <total1> <total2> <inter> <union> <mash_distance> <jaccard_index> <containment_index>

* Minimizers
    Given a $(k + w - 1)$-mer, consider the $w$ contained $k$-mers. The (rightmost) $k$-mer with
    minimal hash (for some given hash function) is the minimizer.

* We use minimizers here to sample kmers
    * For proteins, the length is short, so the window size can be set as: `-k 7 -w 2`
    * DNA: `-k 21 -w 5`
    * Increasing the window size speeds up processing

* Hash Algorithms (--hasher):
    * The `--hasher` parameter selects the hash algorithm used for minimizer calculation.
    * Available options:
        - `rapid`: RapidHash (default)
        - `fx`: FxHash
        - `murmur`: MurmurHash3
    * Note: The `mod` option is not a hash algorithm but a special mode for DNA sequences.

* Mod-Minimizer (--hasher mod):
    * It generates canonical k-mers, meaning that a sequence and its reverse complement
      are generating the same k-mer set.

* To get accurate pairwise sequence identities, use clustalo
  https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity

* Input Modes:
    * By default (--list-files is false):
        * Single file: Treat the file as a sequence file and calculate pairwise distances
          for all sequences within it.
        * Two files: Treat both files as sequence files and calculate pairwise distances
          between sequences from the two files.
    * When --list-files is set:
        * Single file: Treat the file as a list file (each line is a path to a sequence file)
          and calculate pairwise distances for all sequences in the listed files.
        * Two files: Treat both files as list files and calculate pairwise distances
          between sequences from the two list files.

* --merge Behavior:
  - By default (--merge is false):
    * Distances are calculated between individual sequences.
  - When --merge is set:
    * For a single sequence file: Merge all sequences within the file into a single set
      of minimizers. Note that comparing this set to itself (self-comparison) is not
      meaningful, as the distance will always be 0 and the similarity will always be 1.
    * For two sequence files: Merge all sequences within each file into a single set,
      and calculate distances between the two sets.
    * When --list-files is set, --merge operates on each sequence file individually:
      - For each file listed in the list file, merge all sequences within that file
        into a single set, and calculate distances between these sets.
      - The merging does not span across multiple files listed in the list file.

Examples:
1. Calculate distances with default parameters:
   pgr dist seq input.fa

2. Use Mod-Minimizer for DNA sequences (canonical k-mers):
   pgr dist seq input.fa --hasher mod -k 21 -w 5

3. Compare two FA files:
   pgr dist seq file1.fa file2.fa

4. Merge all sequences in a file and compare to another:
   pgr dist seq file1.fa file2.fa --merge

5. Treat input as a list file and calculate distances:
   pgr dist seq list.txt --list-files

6. Use 4 threads for parallel processing:
   pgr dist seq input.fa --parallel 4

"###,
        )
        .arg(crate::cmd_pgr::args::pair_infiles_arg())
        .arg(crate::cmd_pgr::args::hasher_arg())
        .arg(crate::cmd_pgr::args::kmer_arg())
        .arg(crate::cmd_pgr::args::window_arg())
        .arg(crate::cmd_pgr::args::sim_arg())
        .arg(
            clap::Arg::new("zero")
                .long("zero")
                .action(clap::ArgAction::SetTrue)
                .help("Also write results with zero Jaccard index"),
        )
        .arg(
            clap::Arg::new("merge")
                .long("merge")
                .action(clap::ArgAction::SetTrue)
                .help("Merge all sequences within a file into a single set for comparison"),
        )
        .arg(crate::cmd_pgr::args::list_arg())
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the seq command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opt_hasher = args.get_one::<String>("hasher").unwrap();
    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_window = *args.get_one::<usize>("window").unwrap();
    anyhow::ensure!(opt_kmer > 0, "--kmer must be positive: {}", opt_kmer);
    anyhow::ensure!(opt_window > 0, "--window must be positive: {}", opt_window);

    let is_sim = args.get_flag("sim");
    let is_zero = args.get_flag("zero");
    let is_merge = args.get_flag("merge");
    let is_list = args.get_flag("list_files");
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    let infiles = crate::cmd_pgr::args::collect_infiles(args);

    let (sender, writer_thread) = pgr::libs::par::spawn_writer_and_pool(
        crate::cmd_pgr::args::get_outfile(args),
        opt_parallel,
    )?;

    let (entries1, entries2) = pgr::libs::par::load_two_sets(&infiles, is_list, |paths| {
        pgr::libs::par::load_entries(paths, |p| {
            pgr::libs::hash::load_minimizers(p, opt_hasher, opt_kmer, opt_window, is_merge)
        })
    })?;

    pgr::libs::par::par_run_pairs(&entries1, &entries2, &sender, |e1, e2| {
        let d = pgr::libs::hash::set_distances(&e1.set, &e2.set, opt_kmer);

        if !is_zero && d.jaccard == 0. {
            return None;
        }

        let dist = if is_sim {
            pgr::libs::hash::mash_to_sim(d.mash)
        } else {
            d.mash
        };

        let line = if is_merge {
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
                e1.name,
                e2.name,
                d.total1,
                d.total2,
                d.inter,
                d.union,
                dist,
                d.jaccard,
                d.containment
            )
        } else {
            format!(
                "{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
                e1.name, e2.name, dist, d.jaccard, d.containment
            )
        };
        Some(line)
    });

    // Drop the sender to signal the writer thread to exit
    drop(sender);
    // Wait for the writer thread to finish
    writer_thread.join().map_err(|e| {
        let msg = e
            .downcast_ref::<String>()
            .map(|s| s.as_str())
            .or_else(|| e.downcast_ref::<&str>().copied())
            .unwrap_or("<non-string panic payload>");
        anyhow::anyhow!("writer thread panicked: {}", msg)
    })?;

    Ok(())
}
