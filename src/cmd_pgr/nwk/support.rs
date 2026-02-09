use clap::*;
use pgr::libs::phylo::tree::{self, support};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("support")
        .about("Attribute support values (bootstrap) to a tree")
        .after_help(
            r###"
Attributes bootstrap support values to a target tree based on a set of replicate trees.

Notes:
* The first argument is the target tree (to which support values are attributed).
* The second argument is the replicate trees (e.g., from bootstrap).
* Support values are written as internal node labels.
* Assumes that all trees have the same set of leaves.

Examples:
1. Attribute support values:
   pgr nwk support target.nwk replicates.nwk

2. Output support as percentages:
   pgr nwk support target.nwk replicates.nwk --percent
"###,
        )
        .arg(
            Arg::new("target")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Target tree file"),
        )
        .arg(
            Arg::new("replicates")
                .required(true)
                .num_args(1)
                .index(2)
                .help("Replicate trees file"),
        )
        .arg(
            Arg::new("percent")
                .short('p')
                .long("percent")
                .action(ArgAction::SetTrue)
                .help("Print values as percentages"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let target_file = args.get_one::<String>("target").unwrap();
    let replicates_file = args.get_one::<String>("replicates").unwrap();
    let percent = args.get_flag("percent");

    // 1. Read Replicate Trees
    // We read replicates first to build the leaf map and counts, similar to nw_support logic
    let replicates = tree::io::from_file(replicates_file)?;
    if replicates.is_empty() {
         return Err(anyhow::anyhow!("No replicate trees found"));
    }
    let total_reps = replicates.len();

    // 2. Read Target Trees
    let mut targets = tree::io::from_file(target_file)?;
    if targets.is_empty() {
        return Err(anyhow::anyhow!("No target trees found"));
    }
    
    // 3. Build Leaf Map (from first replicate)
    let leaf_map = support::build_leaf_map(&replicates[0]).map_err(|e| anyhow::anyhow!(e))?;
    
    // 4. Count Clades in Replicates
    let counts = support::count_clades(&replicates, &leaf_map).map_err(|e| anyhow::anyhow!(e))?;
    
    // 5. Annotate Target Trees
    for target in &mut targets {
        let target_bitsets = support::compute_all_bitsets(target, &leaf_map).map_err(|e| anyhow::anyhow!(e))?;
        
        for (id, bs) in target_bitsets {
            let node = target.get_node_mut(id).unwrap();
            
            // Only annotate internal nodes
            if !node.is_leaf() {
                let count = counts.get(&bs).copied().unwrap_or(0);
                
                let label = if percent {
                     if total_reps > 0 {
                        format!("{}", (count * 100) / total_reps)
                     } else {
                        "0".to_string()
                     }
                } else {
                    format!("{}", count)
                };
                
                // Overwrite existing label
                node.name = Some(label);
            }
        }
        
        println!("{}", target.to_newick());
    }
    
    Ok(())
}
