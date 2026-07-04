use clap::{ArgMatches, Command};
use pgr::libs::chain::{chain_psl, GapCalc, ScoreContext, SubMatrix};
use pgr::libs::fmt::twobit::TwoBitFile;
/// Build the clap subcommand for chain.
pub fn make_subcommand() -> Command {
    Command::new("chain")
        .about("Chains PSL alignments")
        .after_help(
            r###"
Processing:
  1. Group PSL blocks by target/query sequence and strand.
  2. Build a KD-tree (k-dimensional tree) for efficient predecessor search.
     - In this context, it's a 2D tree indexing blocks by (query_start, target_start).
     - It allows fast range queries to find candidate predecessor blocks that are "before" the current block in both query and target coordinates.
  3. Connect blocks into chains using dynamic programming:
     - Score = BlockScore + Max(PredecessorScore - GapCost).
     - Block Scoring:
       * Default: Identity matrix (Match: +100, Mismatch: -100).
       * Custom: Use --score-scheme to load a LASTZ format file or preset (hoxd55).
     - Gap Cost (Penalty):
       * Linear (Default): --gap-model loose (suitable for distant species).
                           --gap-model medium (suitable for mouse/human).
       * Affine: Use --align-gap-open and --align-gap-extend to override linear costs.
         (Cost = open + extend * length).
     - Overlaps are trimmed by finding the optimal cut point based on exact sequence scores.
  4. Filter chains by minimum score (controlled by --min-score).
     - Default is 1000 to match UCSC axtChain behavior.

Examples:
1. Chain PSL file with default settings:
   pgr psl chain t.2bit q.2bit in.psl -o out.chain

2. Chain with affine gap costs:
   pgr psl chain t.2bit q.2bit in.psl -o out.chain --align-gap-open 400 --align-gap-extend 30

3. Chain with HoxD55 scoring scheme:
   pgr psl chain t.2bit q.2bit in.psl -o out.chain --score-scheme hoxd55
"###,
        )
        .arg(crate::cmd_pgr::args::target_genome_arg(
            "Path to the target genome 2bit file",
        ))
        .arg(crate::cmd_pgr::args::query_genome_arg(
            "Path to the query genome 2bit file",
        ))
        .arg(crate::cmd_pgr::args::psl_positional_arg("Path to the PSL file"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(crate::cmd_pgr::args::gap_model_arg(
            "loose",
            &["loose", "medium"],
            "Gap model: loose or medium",
        ))
        .arg(crate::cmd_pgr::args::min_score_arg("1000"))
        .arg(crate::cmd_pgr::args::align_gap_open_arg())
        .arg(crate::cmd_pgr::args::align_gap_extend_arg())
        .arg(crate::cmd_pgr::args::score_scheme_arg())
}
/// Execute the chain command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = args.get_one::<String>("psl").unwrap();
    let output = crate::cmd_pgr::args::get_outfile(args);
    let gap_model = args.get_one::<String>("gap_model").unwrap();
    let min_score = *args.get_one::<f64>("min_score").unwrap();
    let target_2bit_path = args.get_one::<String>("target");
    let query_2bit_path = args.get_one::<String>("query");
    let score_scheme_path = args.get_one::<String>("score_scheme");

    let reader = pgr::reader(input)?;
    let mut writer = pgr::writer(output)?;

    let mut t_2bit = if let Some(path) = target_2bit_path {
        Some(TwoBitFile::open(path)?)
    } else {
        None
    };

    let mut q_2bit = if let Some(path) = query_2bit_path {
        Some(TwoBitFile::open(path)?)
    } else {
        None
    };

    let score_matrix = if let Some(path) = score_scheme_path {
        SubMatrix::from_name(path)?
    } else {
        SubMatrix::default()
    };

    let mut score_context = match (t_2bit.as_mut(), q_2bit.as_mut()) {
        (Some(t), Some(q)) => Some(ScoreContext {
            t_2bit: t,
            q_2bit: q,
            matrix: &score_matrix,
        }),
        _ => None,
    };

    let gap_open = args.get_one::<i32>("align_gap_open");
    let gap_extend = args.get_one::<i32>("align_gap_extend");

    if gap_open.is_some() != gap_extend.is_some() {
        anyhow::bail!("--align-gap-open and --align-gap-extend must be provided together");
    }

    let gap_calc = if let (Some(&open), Some(&extend)) = (gap_open, gap_extend) {
        GapCalc::affine(open, extend)
    } else {
        match gap_model.as_str() {
            "loose" => GapCalc::loose(),
            "medium" => GapCalc::medium(),
            _ => anyhow::bail!("invalid gap_model: {}", gap_model),
        }
    };

    chain_psl(reader, &mut writer, &gap_calc, min_score, &mut score_context)?;

    Ok(())
}
