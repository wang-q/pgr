use clap::{value_parser, Arg, ArgMatches, Command};

use super::common::{run_sketch_distances, SketchOptions};
use pgr::libs::hash::mash_sketch_distances;

/// Build the clap subcommand for mash.
pub fn make_subcommand() -> Command {
    Command::new("mash")
        .about("Estimates sequence distances with Mash-compatible MinHash sketches")
        .after_help(
            r###"
This command estimates pairwise distances using bottom-k MinHash sketches,
byte-for-byte compatible with Mash (Ondov et al. 2016): canonical k-mers
(min of forward/reverse complement) hashed with MurmurHash3_x64_128
(seed 42), keeping the `--size` smallest unique hashes. The Jaccard index is
Mash's definition (matching hashes in a `--size`-step merge walk of the two
sorted sketches / `--size`, not the standard set Jaccard), and containment is
the full sketch intersection / first-set size (Mash `within` semantics), so
distances match `mash dist` for identical k / sketch size.

* The outputs are printed to stdout as:
    <sequence1> <sequence2> <mash_distance> <jaccard_index> <containment_index>
* With --merge:
    <file1> <file2> <total1> <total2> <inter> <union> <mash> <jaccard> <containment>

* -k k-mer size (default 21); --size sketch size (default 1000, Mash default);
  --seed hash seed (default 42, Mash default).
* MinHash Jaccard is unbiased for similarly-sized sets (Broder 1997);
  containment uses the first input as the denominator and is biased for very
  different-sized sets.
* -p parallelizes sketch loading across input files and pair comparison;
  speedup saturates around the number of input files/sequences (e.g. 4 query
  files -> ~4 threads), so -p beyond that adds little.

Examples:
1. Mash-compatible distances (defaults -k 21 --size 1000 match `mash dist`):
   pgr dist mash input.fa

2. Compare two FA files:
   pgr dist mash file1.fa file2.fa

3. Larger sketch for tighter estimates:
   pgr dist mash file1.fa file2.fa --size 10000

4. Merge all sequences in each file before comparing:
   pgr dist mash file1.fa file2.fa --merge

5. Use 4 threads for parallel processing:
   pgr dist mash input.fa --parallel 4
"###,
        )
        .arg(crate::cmd_pgr::args::pair_infiles_arg())
        .arg(crate::cmd_pgr::args::kmer_arg_with_default("21"))
        .arg(
            Arg::new("size")
                .long("size")
                .num_args(1)
                .default_value("1000")
                .value_parser(value_parser!(usize))
                .help("Sketch size (number of min-hashes)"),
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .num_args(1)
                .default_value("42")
                .value_parser(value_parser!(u32))
                .help("Hash seed (Mash default 42)"),
        )
        .arg(crate::cmd_pgr::args::sim_arg())
        .arg(crate::cmd_pgr::args::zero_arg())
        .arg(crate::cmd_pgr::args::merge_arg())
        .arg(crate::cmd_pgr::args::list_arg())
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the mash command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_size = *args.get_one::<usize>("size").unwrap();
    let opt_seed = *args.get_one::<u32>("seed").unwrap();
    anyhow::ensure!(opt_kmer > 0, "--kmer must be positive: {}", opt_kmer);
    anyhow::ensure!(opt_size > 0, "--size must be positive: {}", opt_size);

    let mut opt = SketchOptions::from_args(args, opt_kmer);
    opt.is_ci = false; // Mash Jaccard semantics differ; CI not applicable
    run_sketch_distances(
        args,
        &opt,
        |infile, is_merge| {
            pgr::libs::hash::load_mash_minhashes(infile, opt_kmer, opt_size, opt_seed, is_merge)
        },
        |s1, s2| mash_sketch_distances(s1, s2, opt_kmer, opt_size),
    )
}
