use clap::*;
use pgr::libs::phylo::tree::Tree;
use pgr::libs::phylo::TreeComparison;
use std::io::{self, Write};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("cmp")
        .about("Compare trees (RF distance)")
        .after_help(
            r###"
Compare trees using Robinson-Foulds (RF) distance.

Notes:
* Input:
    * One file: Compares all trees in the file against each other (Pairwise).
    * Two files: Compares each tree in file1 against each tree in file2.

* Output:
    * TSV format: Tree1_ID \t Tree2_ID \t Distance

* IDs are 1-based indices of the trees in the input files.

Examples:
1. Compare all trees in a file:
   pgr nwk cmp trees.nwk

2. Compare trees between two files:
   pgr nwk cmp set1.nwk set2.nwk
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("First input filename (or stdin)"),
        )
        .arg(
            Arg::new("compare_file")
                .num_args(1)
                .index(2)
                .help("Second input filename (optional)"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    // 1. Load first file
    let infile = args.get_one::<String>("infile").unwrap();
    let trees1 = Tree::from_file(infile)?;

    // 2. Load second file (if provided) or point to first
    let trees2 = if let Some(f2) = args.get_one::<String>("compare_file") {
        Tree::from_file(f2)?
    } else {
        // Clone is not ideal for large datasets but simpler for ownership now.
        // Actually, if we just want to compare self-to-self, we can optimize loops,
        // but to keep logic uniform, we'll just treat it as a second list.
        // However, Tree::from_file returns Vec<Tree>.
        // For efficiency, let's just use references/indices logic below.
        Vec::new() 
    };

    let self_compare = args.get_one::<String>("compare_file").is_none();

    // 3. Output writer
    let outfile = args.get_one::<String>("outfile").unwrap();
    let mut writer: Box<dyn Write> = if outfile == "stdout" {
        Box::new(io::stdout())
    } else {
        Box::new(std::fs::File::create(outfile)?)
    };

    // 4. Compare
    // Header
    writeln!(writer, "Tree1\tTree2\tRF_Dist")?;

    if self_compare {
        // Compare trees1 vs trees1 (All pairs? Or just upper triangle?)
        // Standard "cmp" usually implies full list or matrix. 
        // Let's do full pairwise for now (N^2), consistent with "compare A to B".
        for (i, t1) in trees1.iter().enumerate() {
            for (j, t2) in trees1.iter().enumerate() {
                // Skip redundant calculations if desired?
                // But for a full matrix output (even if in list format), users might expect A-B and B-A.
                // RF is symmetric. 
                // Let's optimize: if j < i, use result from i,j? 
                // For simplicity: compute all.
                
                let dist = match t1.robinson_foulds(t2) {
                    Ok(d) => d.to_string(),
                    Err(e) => format!("Error: {}", e), // Handle leaf mismatch gracefully?
                };
                writeln!(writer, "{}\t{}\t{}", i + 1, j + 1, dist)?;
            }
        }
    } else {
        // Compare trees1 vs trees2
        for (i, t1) in trees1.iter().enumerate() {
            for (j, t2) in trees2.iter().enumerate() {
                let dist = match t1.robinson_foulds(t2) {
                    Ok(d) => d.to_string(),
                    Err(e) => format!("Error: {}", e),
                };
                // Use 1-based index, maybe prefix with file? 
                // Just indices 1..N and 1..M
                writeln!(writer, "{}\t{}\t{}", i + 1, j + 1, dist)?;
            }
        }
    }

    Ok(())
}
