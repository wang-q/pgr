//! Shared execution for the sketch-distance command family (`mini` /
//! `mash` / `frac`): input modes, parallelism, and output formatting are
//! identical; each command only supplies its sampler loader and distance
//! function.

use clap::ArgMatches;
use pgr::libs::hash::{ani_ci_from_jaccard, MinimizerEntry, SetDistances};
use pgr::libs::par;
use rapidhash::RapidHashSet;

/// Options shared by `dist mini` / `dist mash` / `dist frac`.
pub struct SketchOptions {
    pub is_merge: bool,
    pub is_list: bool,
    pub is_sim: bool,
    pub is_zero: bool,
    pub is_ci: bool,
    pub parallel: usize,
    pub outfile: String,
    pub kmer: usize,
}

impl SketchOptions {
    pub fn from_args(args: &ArgMatches, kmer: usize) -> Self {
        SketchOptions {
            is_merge: args.get_flag("merge"),
            is_list: args.get_flag("list_files"),
            is_sim: args.get_flag("sim"),
            is_zero: args.get_flag("zero"),
            is_ci: args
                .try_get_one::<bool>("ci")
                .ok()
                .flatten()
                .copied()
                .unwrap_or(false),
            parallel: *args.get_one::<usize>("parallel").unwrap(),
            outfile: crate::cmd_pgr::args::get_outfile(args).to_string(),
            kmer,
        }
    }
}

/// Run pairwise sketch distances. `load(infile, is_merge)` builds the
/// per-file entries; `distance(s1, s2)` computes the metrics.
pub fn run_sketch_distances<F, G>(
    args: &ArgMatches,
    opt: &SketchOptions,
    load: F,
    distance: G,
) -> anyhow::Result<()>
where
    F: Fn(&str, bool) -> anyhow::Result<Vec<MinimizerEntry>> + Send + Sync,
    G: Fn(&RapidHashSet<u64>, &RapidHashSet<u64>) -> SetDistances + Send + Sync,
{
    let infiles = crate::cmd_pgr::args::collect_infiles(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(&opt.outfile, infiles.iter().copied())?;
    let (sender, writer_thread) = par::spawn_writer_and_pool(&opt.outfile, opt.parallel)?;

    let (entries1, entries2) = par::load_two_sets(&infiles, opt.is_list, |paths| {
        par::load_entries(paths, |p| load(p, opt.is_merge))
    })?;

    par::par_run_pairs(&entries1, &entries2, &sender, |e1, e2| {
        let d = distance(&e1.set, &e2.set);

        if !opt.is_zero && d.jaccard == 0. {
            return None;
        }

        let dist = if opt.is_sim {
            pgr::libs::hash::mash_to_sim(d.mash)
        } else {
            d.mash
        };

        let (ci_lo, ci_hi) = if opt.is_ci {
            ani_ci_from_jaccard(d.jaccard, d.union, opt.kmer)
        } else {
            (0.0, 0.0)
        };

        let line = if opt.is_merge {
            if opt.is_ci {
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\n",
                    e1.name,
                    e2.name,
                    d.total1,
                    d.total2,
                    d.inter,
                    d.union,
                    dist,
                    d.jaccard,
                    d.containment,
                    ci_lo,
                    ci_hi
                )
            } else {
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
            }
        } else if opt.is_ci {
            format!(
                "{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\n",
                e1.name, e2.name, dist, d.jaccard, d.containment, ci_lo, ci_hi
            )
        } else {
            format!(
                "{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
                e1.name, e2.name, dist, d.jaccard, d.containment
            )
        };
        Some(line)
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
