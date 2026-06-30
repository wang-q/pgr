use super::common;
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
    pgr six-frame input.fa |
        pgr dist hv stdin match.fa

"###,
        )
        .arg(common::infiles_arg())
        .arg(common::hasher_arg())
        .arg(common::kmer_arg())
        .arg(common::window_arg())
        .arg(
            clap::Arg::new("dim")
                .long("dim")
                .short('d')
                .num_args(1)
                .default_value("4096")
                .value_parser(clap::value_parser!(usize))
                .help("The dimension size should be a multiple of 32."),
        )
        .arg(common::sim_arg())
        .arg(common::list_arg())
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

#[derive(Debug, Default, Clone)]
struct HvEntry {
    name: String,
    set: Vec<i32>,
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
    let is_list = args.get_flag("list");
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    let infiles = common::collect_infiles(args);

    let (sender, writer_thread) =
        common::spawn_writer_and_pool(crate::cmd_pgr::args::get_outfile(args), opt_parallel)?;

    //----------------------------
    // Ops
    //----------------------------
    let (entries1, entries2) = common::load_two_sets(&infiles, is_list, |paths| {
        common::load_entries(paths, |p| {
            load_file(p, opt_hasher, opt_kmer, opt_window, opt_dim)
        })
    })?;

    common::par_run_pairs(&entries1, &entries2, &sender, |e1, e2| {
        let (total1, total2, inter, union, mash, jaccard, containment) =
            calc_distances(&e1.set, &e2.set, opt_kmer);

        let dist = if is_sim { 1.0 - mash } else { mash };

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
            e1.name, e2.name, total1, total2, inter, union, dist, jaccard, containment
        );
        Some(line)
    });

    // Drop the sender to signal the writer thread to exit
    drop(sender);
    // Wait for the writer thread to finish
    writer_thread.join().unwrap();

    Ok(())
}

fn load_file(
    infile: &str,
    opt_hasher: &str,
    opt_kmer: usize,
    opt_window: usize,
    opt_dim: usize,
) -> anyhow::Result<Vec<HvEntry>> {
    let mut fa_in = pgr::libs::fmt::fa::reader(infile)?;

    let mut file_set = rapidhash::RapidHashSet::default();

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;
        let seq = record.sequence();

        let set: rapidhash::RapidHashSet<u64> =
            pgr::libs::hash::seq_mins(&seq[..], opt_hasher, opt_kmer, opt_window)?;

        file_set.extend(set);
    }

    let seed_vec: Vec<u64> = file_set.into_iter().collect();
    let hv: Vec<i32> = pgr::libs::hv::hash_hv_i8(&seed_vec, opt_dim);
    let entry = HvEntry {
        name: infile.to_string(),
        set: hv,
    };

    Ok(vec![entry])
}

// Calculate Jaccard, Containment, and Mash distance between two sets
fn calc_distances(
    s1: &[i32],
    s2: &[i32],
    opt_kmer: usize,
) -> (usize, usize, usize, usize, f32, f32, f32) {
    let card1 = pgr::libs::hv::hv_cardinality(s1);
    let card2 = pgr::libs::hv::hv_cardinality(s2);

    let inter = pgr::libs::hv::hv_dot(s1, s2)
        .min(card1 as f32)
        .min(card2 as f32);
    let union = card1 as f32 + card2 as f32 - inter;

    let jaccard = inter / union;
    let containment = inter / card1 as f32;
    // https://mash.readthedocs.io/en/latest/distances.html#mash-distance-formulation
    let mash = if jaccard == 0.0 {
        1.0
    } else {
        ((-1.0 / opt_kmer as f32) * ((2.0f32 * jaccard) / (1.0f32 + jaccard)).ln()).abs()
    };

    (
        card1,
        card2,
        inter as usize,
        union as usize,
        mash,
        jaccard,
        containment,
    )
}
