use clap::{Arg, ArgMatches, Command};
use std::cell::{Ref, RefCell};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader};
use std::rc::Rc;

use pgr::libs::net::{read_nets, Fill, Gap};

pub fn make_subcommand() -> Command {
    Command::new("class").about("Show stats of net").arg(
        Arg::new("input")
            .required(true)
            .help("Input net file (or stdin if '-')"),
    )
}

struct Stats {
    count: u64,
    bases: u64,
}

impl Default for Stats {
    fn default() -> Self {
        Self { count: 0, bases: 0 }
    }
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("input").unwrap();

    let reader: Box<dyn io::BufRead> = if input_path == "-" {
        Box::new(BufReader::new(io::stdin()))
    } else {
        Box::new(BufReader::new(File::open(input_path)?))
    };

    let chroms = read_nets(reader)?;

    let mut stats_map: HashMap<String, Stats> = HashMap::new();
    let mut total_bases = 0;

    for chrom in chroms {
        total_bases += chrom.size;
        collect_stats_gap(&chrom.root, &mut stats_map);
    }

    // Print results
    // We want to sort by bases desc? Or just alphabetical?
    // UCSC usually sorts by bases or hierarchy.
    // Let's sort by bases desc.

    let mut results: Vec<(String, u64, u64)> = stats_map
        .into_iter()
        .map(|(k, v)| (k, v.count, v.bases))
        .collect();

    results.sort_by(|a, b| b.2.cmp(&a.2)); // Sort by bases descending

    println!(
        "{:<20} {:>10} {:>15} {:>10}",
        "Class", "Count", "Bases", "%"
    );
    println!(
        "{:<20} {:>10} {:>15} {:>10}",
        "-----", "-----", "-----", "-"
    );

    for (class, count, bases) in results {
        let pct = if total_bases > 0 {
            (bases as f64 / total_bases as f64) * 100.0
        } else {
            0.0
        };
        println!("{:<20} {:>10} {:>15} {:>10.2}%", class, count, bases, pct);
    }

    println!("\nTotal bases covered by nets/gaps: {}", total_bases);

    Ok(())
}

fn collect_stats_gap(gap: &Rc<RefCell<Gap>>, stats: &mut HashMap<String, Stats>) {
    let gap_ref: Ref<Gap> = gap.borrow();
    let size = gap_ref.end - gap_ref.start;

    // Gap itself is a "gap" class if we want to count it?
    // Or do we only count fills?
    // UCSC netClass counts "gap" as well.
    // But Gaps contain Fills.
    // The "gap" bases are (size - sum(fills.size)).

    // Actually, usually we count the explicit objects.
    // A Gap object represents a gap in the alignment.
    // But in the net structure, Gap is a container.
    // The "unfilled" part of the Gap is the actual gap.

    let mut fill_bases = 0;
    for fill in &gap_ref.fills {
        let fill_ref: Ref<Fill> = fill.borrow();
        fill_bases += fill_ref.end - fill_ref.start;

        // Count the fill
        let class = if fill_ref.class.is_empty() {
            "unknown".to_string()
        } else {
            fill_ref.class.clone()
        };

        let entry = stats.entry(class).or_default();
        entry.count += 1;
        entry.bases += fill_ref.end - fill_ref.start;

        // Recurse
        collect_stats_fill(fill, stats);
    }

    // The remaining part is gap
    let gap_bases = size - fill_bases;
    if gap_bases > 0 {
        let entry = stats.entry("gap".to_string()).or_default();
        entry.count += 1; // This is tricky. Is it 1 gap? Or multiple fragments?
                          // In this recursive structure, the "gap" is implicitly the background.
                          // We can just add the bases.
                          // Count is hard to define for the implicit background gap.
                          // Let's just track bases for gap.
        entry.bases += gap_bases;
    }
}

fn collect_stats_fill(fill: &Rc<RefCell<Fill>>, stats: &mut HashMap<String, Stats>) {
    let fill_ref: Ref<Fill> = fill.borrow();

    // Fill contains Gaps.
    // The Fill itself covers fill_ref.end - fill_ref.start.
    // This was already added to the stats in the parent.
    // But we need to recurse into its gaps.

    for gap in &fill_ref.gaps {
        collect_stats_gap(gap, stats);
    }
}
