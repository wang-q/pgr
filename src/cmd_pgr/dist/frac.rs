use clap::{ArgMatches, Command};

use super::common::{run_sketch_distances, SketchOptions};
use pgr::libs::hash::set_distances;

/// Build the clap subcommand for frac.
pub fn make_subcommand() -> Command {
    Command::new("frac")
        .about("Estimates unbiased sequence distances with FracMinHash sketches")
        .after_help(
            r###"
This command estimates pairwise distances using FracMinHash sketches (Irber
et al. 2022): every canonical k-mer is kept with independent probability
1/scale (hash < u64::MAX/scale). Unlike minimizers/syncmers, the Jaccard and
containment estimates are unbiased and comparable across differently-sized
sets, and support ANI bias correction with confidence intervals (Hera et al.
2023). Output format matches the other sketch-distance commands.

* The outputs are printed to stdout as:
    <sequence1> <sequence2> <mash_distance> <jaccard_index> <containment_index>
* With --merge:
    <file1> <file2> <total1> <total2> <inter> <union> <mash> <jaccard> <containment>
* With --ci (unbiased sampler): append the 95% ANI confidence interval.

* DNA default: -k 21; protein default: -k 7 (applied automatically).
* --scale controls density (default 1000; smaller = denser = lower variance).

This is the recommended command for numeric ANI estimation. See
notes/benchmarks/dist-cohort-validation.md for the unbiasedness validation.

Examples:
1. Unbiased numeric ANI with 95% confidence interval:
   pgr dist frac input.fa --ci

2. Compare two FA files:
   pgr dist frac file1.fa file2.fa

3. Denser sampling (lower variance, slower):
   pgr dist frac file1.fa file2.fa --scale 100

4. Protein sequences:
   pgr dist frac proteins.fa --protein --ci

5. Merge all sequences in each file and compare the two sets:
   pgr dist frac file1.fa file2.fa --merge
"###,
        )
        .arg(crate::cmd_pgr::args::pair_infiles_arg())
        .arg(crate::cmd_pgr::args::kmer_arg_mode_dependent())
        .arg(crate::cmd_pgr::args::scale_arg())
        .arg(crate::cmd_pgr::args::protein_arg())
        .arg(crate::cmd_pgr::args::sim_arg())
        .arg(
            clap::Arg::new("ci")
                .long("ci")
                .action(clap::ArgAction::SetTrue)
                .help("Append 95% ANI confidence interval (unbiased sampler)"),
        )
        .arg(crate::cmd_pgr::args::zero_arg())
        .arg(crate::cmd_pgr::args::merge_arg())
        .arg(crate::cmd_pgr::args::list_arg())
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the frac command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opt_sampler = "frachash";
    let is_protein = args.get_flag("protein");
    let (opt_kmer, _) = crate::cmd_pgr::args::resolve_kmer_window(args, opt_sampler, is_protein);
    let opt_scale = *args.get_one::<usize>("scale").unwrap();
    anyhow::ensure!(opt_kmer > 0, "--kmer must be positive: {}", opt_kmer);
    anyhow::ensure!(opt_scale > 0, "--scale must be positive: {}", opt_scale);

    let opt = SketchOptions::from_args(args, opt_kmer);
    run_sketch_distances(
        args,
        &opt,
        |infile, is_merge| {
            pgr::libs::hash::load_fracminhash(infile, opt_kmer, opt_scale, is_protein, is_merge)
        },
        |s1, s2| set_distances(s1, s2, opt_kmer),
    )
}
