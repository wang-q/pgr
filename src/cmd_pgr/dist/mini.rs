use clap::{ArgMatches, Command};

use super::common::{run_sketch_distances, SketchOptions};
use pgr::libs::hash::set_distances;

/// Build the clap subcommand for mini.
pub fn make_subcommand() -> Command {
    Command::new("mini")
        .about("Estimates sequence distances with minimizer sketches")
        .after_help(
            r###"
This command estimates pairwise distances between sequences using minimizer
sketches (the minimal-hash k-mer in each window). Output format matches the
other sketch-distance commands (`pgr dist mash` / `pgr dist frac`).

* The outputs are printed to stdout as:
    <sequence1> <sequence2> <mash_distance> <jaccard_index> <containment_index>
* With --merge:
    <file1> <file2> <total1> <total2> <inter> <union> <mash> <jaccard> <containment>

* DNA minimizer defaults: -k 21 -w 5 (Mash convention; applied automatically).
* Protein: -k 7 -w 1 (short sequences).
* --hasher: rapid (default) / fx / murmur / mod. `mod` is a special mode that
  emits canonical k-mers (sequence and reverse complement share the same
  k-mer set) and is DNA-only; combining it with --protein is rejected.
* Increasing the window size speeds up processing.

Note: minimizer Jaccard estimates are biased and inconsistent (Belbasi et al.
2022). Use this command for fast ranking/screening; for unbiased numeric ANI
use `pgr dist frac` (FracMinHash). See notes/benchmarks/dist-cohort-validation.md.

* To get accurate pairwise sequence identities, use clustalo
  https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity

Examples:
1. Calculate distances with default DNA minimizer parameters (-k 21 -w 5):
   pgr dist mini input.fa

2. Use Mod-Minimizer for DNA (canonical k-mers):
   pgr dist mini input.fa --hasher mod -k 21 -w 5

3. Protein sequences (default -k 7 -w 1):
   pgr dist mini proteins.fa --protein

4. Compare two FA files:
   pgr dist mini file1.fa file2.fa

5. Merge all sequences in each file and compare the two sets:
   pgr dist mini file1.fa file2.fa --merge

6. Treat input as a list file (one sequence path per line):
   pgr dist mini list.txt --list-files

7. Use 4 threads for parallel processing:
   pgr dist mini input.fa --parallel 4
"###,
        )
        .arg(crate::cmd_pgr::args::pair_infiles_arg())
        .arg(crate::cmd_pgr::args::kmer_arg_mode_dependent())
        .arg(crate::cmd_pgr::args::window_arg_with_default(
            "5",
            "Window size for minimizers (DNA default 5, protein 1)",
        ))
        .arg(crate::cmd_pgr::args::hasher_arg())
        .arg(crate::cmd_pgr::args::protein_arg())
        .arg(crate::cmd_pgr::args::sim_arg())
        .arg(crate::cmd_pgr::args::zero_arg())
        .arg(crate::cmd_pgr::args::merge_arg())
        .arg(crate::cmd_pgr::args::list_arg())
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the mini command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opt_sampler = "minimizer";
    let opt_hasher = args.get_one::<String>("hasher").unwrap();
    let is_protein = args.get_flag("protein");
    let (opt_kmer, opt_window) =
        crate::cmd_pgr::args::resolve_kmer_window(args, opt_sampler, is_protein);
    anyhow::ensure!(opt_kmer > 0, "--kmer must be positive: {}", opt_kmer);
    anyhow::ensure!(opt_window > 0, "--window must be positive: {}", opt_window);

    if is_protein && opt_hasher == "mod" {
        anyhow::bail!("--hasher mod is DNA-only (canonical reverse complement)");
    }

    let opt = SketchOptions::from_args(args, opt_kmer);
    run_sketch_distances(
        args,
        &opt,
        |infile, is_merge| {
            pgr::libs::hash::load_minimizers(infile, opt_hasher, opt_kmer, opt_window, is_merge)
        },
        |s1, s2| set_distances(s1, s2, opt_kmer),
    )
}
