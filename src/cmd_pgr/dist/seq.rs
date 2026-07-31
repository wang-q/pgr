use clap::{parser::ValueSource, ArgMatches, Command};

/// Build the clap subcommand for seq.
pub fn make_subcommand() -> Command {
    Command::new("seq")
        .about("Estimates sequence distances using minimizers or syncmers")
        .after_help(
            r###"
This command calculates pairwise distances between sequences in FA file(s) using
minimizers or closed syncmers.

* The outputs are printed to stdout in the following format:
    <sequence1> <sequence2> <mash_distance> <jaccard_index> <containment_index>
* With --merge
    <file1> <file2> <total1> <total2> <inter> <union> <mash_distance> <jaccard_index> <containment_index>

* Samplers (--sampler):
    * `minimizer` (default): Given a $(k + w - 1)$-mer, consider the $w$ contained $k$-mers.
      The (rightmost) $k$-mer with minimal hash is the minimizer.
    * `syncmer`: Closed syncmers per Edgar (2021). A window of `w` s-mers is emitted iff
      its minimal s-mer hash falls at the first or last position. This gives a sparse but
      complete cover and localizes indel/rearrangement perturbation better than minimizers.
      `-k` is the s-mer size; `-w` is the number of s-mers per window (span `k + w - 1`).

* We use these samplers to sample kmers
    * For proteins, the length is short, so the window size can be set as: `-k 7 -w 2`
    * DNA (minimizer): `-k 21 -w 5`
    * DNA (syncmer): `-k 8 -w 55` (syng defaults; applied automatically when not set)
    * Protein (syncmer): `-k 7 -w 5` (applied automatically with --protein when not set;
      k=7 keeps random match prob negligible, w=5 gives ~33% density for short sequences)
    * Increasing the window size speeds up processing

* Hash Algorithms (--hasher):
    * The `--hasher` parameter selects the hash algorithm used for minimizer calculation.
    * Available options:
        - `rapid`: RapidHash (default)
        - `fx`: FxHash
        - `murmur`: MurmurHash3
    * Note: The `mod` option is not a hash algorithm but a special mode for DNA sequences.
    * For `--sampler syncmer`, `--hasher` is ignored: DNA uses a 2-bit canonical rolling
      hash and protein uses RapidHash on s-mer bytes.

* Mod-Minimizer (--hasher mod):
    * It generates canonical k-mers, meaning that a sequence and its reverse complement
      are generating the same k-mer set.

* --protein:
    * Declares that the input is protein sequence; affects all samplers.
    * With `--sampler syncmer`: switches to the protein s-mer byte-hash path
      (no reverse complement); DNA uses the 2-bit canonical rolling hash.
    * With `--sampler minimizer`: `--hasher mod` is rejected (mod-minimizer
      requires DNA reverse complement); the byte hashers (rapid/fx/murmur) hash
      raw bytes and work for both DNA and protein.

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

3. Use closed syncmers for DNA (syng defaults applied automatically):
   pgr dist seq input.fa --sampler syncmer

4. Use closed syncmers for proteins (defaults smer=7, window=5 applied automatically):
   pgr dist seq proteins.fa --sampler syncmer --protein

5. Compare two FA files:
   pgr dist seq file1.fa file2.fa

6. Merge all sequences in a file and compare to another:
   pgr dist seq file1.fa file2.fa --merge

7. Treat input as a list file and calculate distances:
   pgr dist seq list.txt --list-files

8. Use 4 threads for parallel processing:
   pgr dist seq input.fa --parallel 4

"###,
        )
        .arg(crate::cmd_pgr::args::pair_infiles_arg())
        .arg(crate::cmd_pgr::args::sampler_arg())
        .arg(crate::cmd_pgr::args::hasher_arg())
        .arg(crate::cmd_pgr::args::kmer_arg())
        .arg(crate::cmd_pgr::args::window_arg())
        .arg(crate::cmd_pgr::args::protein_arg())
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

/// Resolve `--kmer` and `--window`, applying syncmer defaults when `--sampler
/// syncmer` is used without explicit `-k`/`-w`:
/// * DNA (`--protein` off): syng defaults smer=8, window=55.
/// * Protein (`--protein` on): smer=7, window=5 — k=7 keeps random match prob
///   negligible (20^7 ≈ 1.3e9) and matches the protein k=7 convention; w=5
///   gives ~33% density so short proteins still yield enough syncmers.
fn resolve_kmer_window(args: &ArgMatches, opt_sampler: &str, is_protein: bool) -> (usize, usize) {
    let kmer_cli = matches!(args.value_source("kmer"), Some(ValueSource::CommandLine));
    let window_cli = matches!(args.value_source("window"), Some(ValueSource::CommandLine));
    let (def_k, def_w) = if is_protein { (7, 5) } else { (8, 55) };
    let default_k = if opt_sampler == "syncmer" && !kmer_cli {
        def_k
    } else {
        *args.get_one::<usize>("kmer").unwrap()
    };
    let default_w = if opt_sampler == "syncmer" && !window_cli {
        def_w
    } else {
        *args.get_one::<usize>("window").unwrap()
    };
    (default_k, default_w)
}

/// Execute the seq command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opt_sampler = args.get_one::<String>("sampler").unwrap();
    let opt_hasher = args.get_one::<String>("hasher").unwrap();
    let is_protein = args.get_flag("protein");
    let (opt_kmer, opt_window) = resolve_kmer_window(args, opt_sampler, is_protein);
    anyhow::ensure!(opt_kmer > 0, "--kmer must be positive: {}", opt_kmer);
    anyhow::ensure!(opt_window > 0, "--window must be positive: {}", opt_window);

    let is_sim = args.get_flag("sim");
    let is_zero = args.get_flag("zero");
    let is_merge = args.get_flag("merge");
    let is_list = args.get_flag("list_files");
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    // mod-minimizer relies on DNA reverse complement and is meaningless on protein.
    if is_protein && opt_sampler == "minimizer" && opt_hasher == "mod" {
        anyhow::bail!(
            "--hasher mod is DNA-only (canonical reverse complement) and cannot be used with --protein"
        );
    }

    let infiles = crate::cmd_pgr::args::collect_infiles(args);

    let (sender, writer_thread) = pgr::libs::par::spawn_writer_and_pool(
        crate::cmd_pgr::args::get_outfile(args),
        opt_parallel,
    )?;

    let (entries1, entries2) = match opt_sampler.as_str() {
        "syncmer" => pgr::libs::par::load_two_sets(&infiles, is_list, |paths| {
            pgr::libs::par::load_entries(paths, |p| {
                pgr::libs::hash::load_syncmers(p, opt_kmer, opt_window, is_protein, is_merge)
            })
        })?,
        _ => pgr::libs::par::load_two_sets(&infiles, is_list, |paths| {
            pgr::libs::par::load_entries(paths, |p| {
                pgr::libs::hash::load_minimizers(p, opt_hasher, opt_kmer, opt_window, is_merge)
            })
        })?,
    };

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
