use clap::{ArgMatches, Command};

/// Build the clap subcommand for hv.
pub fn make_subcommand() -> Command {
    Command::new("hv")
        .about("Estimates distances between DNA/protein files using hypervectors")
        .after_help(
            r###"
This command calculates pairwise distances between files in FA file(s) using
minimizers or closed syncmers, projected onto hypervectors.

* The outputs are printed to stdout in the following format:
    <file1> <file2> <total1> <total2> <inter> <union> <mash_distance> <jaccard_index> <containment_index>

* Samplers, hash algorithms, --protein, -k/-w semantics are the same as the
  sketch-distance family (`pgr dist mini` / `pgr dist frac`). Syncmer defaults
  (DNA smer=8/window=55, protein smer=7/window=5) are applied automatically
  when --sampler syncmer is used without explicit -k/-w.

* Input Modes:
    * For a single sequence file: Merge all sequences within the file into a single hypervector.
      Note that comparing this set to itself (self-comparison) is not meaningful,
      as the distance will always be 0 and the similarity will always be 1.
    * For two sequence files: Merge all sequences within each file into a single hypervector,
      and calculate distances between the two hypervectors.
    * `.hv` inputs: When the inputs are `.hv` files (produced by `pgr pgi to-hv`),
      they are compared directly; the stored sampling parameters and dimension
      must match between the files.
    * When --list-files is set:
      - For each file listed in the list file, merge all sequences within that file
        into a single hypervector, and calculate distances between these hypervectors.
      - The merging does not span across multiple files listed in the list file.

Examples:
1. Merge all sequences in a file and compare to another:
   pgr dist hv file1.fa file2.fa

2. Use Mod-Minimizer for DNA sequences (canonical k-mers):
   pgr dist hv file1.fa file2.fa --hasher mod -k 21 -w 5

3. Use closed syncmers for DNA:
   pgr dist hv file1.fa file2.fa --sampler syncmer

4. Treat input as a list file and calculate distances:
   pgr dist hv list.txt --list-files

5. Compare two hypervectors from .pgi indexes:
   pgr pgi to-hv a.pgi -o a.hv
   pgr pgi to-hv b.pgi -o b.hv
   pgr dist hv a.hv b.hv

6. Use 4 threads for parallel processing:
   pgr dist hv input.fa --parallel 4

7. Perform six-frame translation on a FA file and match to another
    pgr fa six-frame input.fa |
        pgr dist hv stdin match.fa

"###,
        )
        .arg(crate::cmd_pgr::args::pair_infiles_arg())
        .arg(crate::cmd_pgr::args::sampler_arg())
        .arg(crate::cmd_pgr::args::hasher_arg())
        .arg(crate::cmd_pgr::args::kmer_arg())
        .arg(crate::cmd_pgr::args::window_arg())
        .arg(crate::cmd_pgr::args::protein_arg())
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

/// Execute the hv command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opt_sampler = args.get_one::<String>("sampler").unwrap();
    let opt_hasher = args.get_one::<String>("hasher").unwrap();
    let is_protein = args.get_flag("protein");
    let (opt_kmer, opt_window) =
        crate::cmd_pgr::args::resolve_kmer_window(args, opt_sampler, is_protein);
    let opt_dim = *args.get_one::<usize>("dim").unwrap();
    anyhow::ensure!(opt_kmer > 0, "--kmer must be positive: {}", opt_kmer);
    anyhow::ensure!(opt_window > 0, "--window must be positive: {}", opt_window);
    anyhow::ensure!(opt_dim > 0, "--dim must be positive: {}", opt_dim);

    // mod-minimizer relies on DNA reverse complement and is meaningless on protein.
    if is_protein && opt_sampler == "minimizer" && opt_hasher == "mod" {
        anyhow::bail!(
            "--hasher mod is DNA-only (canonical reverse complement) and cannot be used with --protein"
        );
    }

    let is_sim = args.get_flag("sim");
    let is_list = args.get_flag("list_files");
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    let infiles = crate::cmd_pgr::args::collect_infiles(args);
    if infiles.iter().any(|f| f.ends_with(".hv")) {
        return run_hv_files(
            &infiles,
            is_list,
            is_sim,
            opt_parallel,
            crate::cmd_pgr::args::get_outfile(args),
        );
    }

    let (sender, writer_thread) = pgr::libs::par::spawn_writer_and_pool(
        crate::cmd_pgr::args::get_outfile(args),
        opt_parallel,
    )?;

    let (entries1, entries2) = match opt_sampler.as_str() {
        "syncmer" => pgr::libs::par::load_two_sets(&infiles, is_list, |paths| {
            pgr::libs::par::load_entries(paths, |p| {
                let entry = pgr::libs::hv::load_hv_from_fasta_syncmer(
                    p, opt_kmer, opt_window, is_protein, opt_dim,
                )?;
                Ok(vec![entry])
            })
        })?,
        _ => pgr::libs::par::load_two_sets(&infiles, is_list, |paths| {
            pgr::libs::par::load_entries(paths, |p| {
                let entry = pgr::libs::hv::load_hv_from_fasta(
                    p, opt_hasher, opt_kmer, opt_window, opt_dim,
                )?;
                Ok(vec![entry])
            })
        })?,
    };

    pgr::libs::par::par_run_pairs(&entries1, &entries2, &sender, |e1, e2| {
        let d = pgr::libs::hv::calc_distances(&e1.set, &e2.set, opt_kmer);

        let dist = if is_sim {
            pgr::libs::hash::mash_to_sim(d.mash as f64) as f32
        } else {
            d.mash
        };

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
            e1.name, e2.name, d.card1, d.card2, d.inter, d.union, dist, d.jaccard, d.containment
        );
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

/// Compare `.hv` files directly (produced by `pgr pgi to-hv`).
///
/// Sampling parameters (k) and dimension must match between all compared
/// files; the k-mer size stored in the files drives the Mash distance.
fn run_hv_files(
    infiles: &[&str],
    is_list: bool,
    is_sim: bool,
    opt_parallel: usize,
    outfile: &str,
) -> anyhow::Result<()> {
    let (sender, writer_thread) = pgr::libs::par::spawn_writer_and_pool(outfile, opt_parallel)?;
    let load = |paths: &[String]| -> anyhow::Result<Vec<pgr::libs::pgi::to_hv::HvFile>> {
        pgr::libs::par::load_entries(paths, |p| {
            let mut r = pgr::reader(p)?;
            Ok(vec![pgr::libs::pgi::to_hv::read_hv(&mut r)?])
        })
    };
    let paths1 = pgr::libs::par::resolve_paths(infiles[0], is_list)?;
    let paths2 = if infiles.len() > 1 {
        pgr::libs::par::resolve_paths(infiles[1], is_list)?
    } else {
        paths1.clone()
    };
    let files1 = load(&paths1)?;
    let files2 = load(&paths2)?;
    let k = files1.first().map(|f| f.k).unwrap_or(0);
    let dim = files1.first().map(|f| f.dim).unwrap_or(0);
    let sparse = files1.first().map(|f| f.sparse).unwrap_or(0);
    for f in files1.iter().chain(files2.iter()) {
        anyhow::ensure!(f.k == k, "hv k-mer size mismatch: {} vs {}", f.k, k);
        anyhow::ensure!(f.dim == dim, "hv dimension mismatch: {} vs {}", f.dim, dim);
        anyhow::ensure!(
            f.sparse == sparse,
            "hv sparse-update mismatch: {} vs {}",
            f.sparse,
            sparse
        );
    }
    pgr::libs::par::par_run_pairs(&files1, &files2, &sender, |e1, e2| {
        // Sparse HDC: cosine similarity approximates the k-mer set overlap
        // (shared = cos * sqrt(n1 * n2)); stored k-mer counts give exact
        // set cardinalities.
        let dot: f64 = e1
            .hv
            .iter()
            .zip(&e2.hv)
            .map(|(x, y)| (*x as f64) * (*y as f64))
            .sum();
        let n1 = e1.n_kmer as f64;
        let n2 = e2.n_kmer as f64;
        let na: f64 = e1.hv.iter().map(|x| (*x as f64) * (*x as f64)).sum();
        let nb: f64 = e2.hv.iter().map(|x| (*x as f64) * (*x as f64)).sum();
        let sim = dot / (na.sqrt() * nb.sqrt());
        let inter = (sim * (n1 * n2).sqrt()).round() as usize;
        let inter = inter.min(e1.n_kmer).min(e2.n_kmer);
        let union = e1.n_kmer + e2.n_kmer - inter;
        let jaccard = inter as f32 / union as f32;
        let containment = inter as f32 / e1.n_kmer as f32;
        let mash = pgr::libs::hash::mash_distance(jaccard as f64, e1.k) as f32;
        let dist = if is_sim {
            pgr::libs::hash::mash_to_sim(mash as f64) as f32
        } else {
            mash
        };
        Some(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
            e1.name, e2.name, e1.n_kmer, e2.n_kmer, inter, union, dist, jaccard, containment
        ))
    });
    drop(sender);
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
