use clap::{ArgMatches, Command};

use pgr::libs::fmt::axt::AxtReader;
use pgr::libs::fmt::psl::Psl;
/// Build the clap subcommand for to-psl.
pub fn make_subcommand() -> Command {
    Command::new("to-psl")
        .about("Converts from axt to psl format")
        .after_help(
            r###"
Where tSizes and qSizes are tab-delimited files with <seqName> <size> columns.

Examples:
  # Convert axt to psl
  pgr axt to-psl in.axt -t t.sizes -q q.sizes -o out.psl
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input AXT file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::t_sizes_arg().required(true))
        .arg(crate::cmd_pgr::args::q_sizes_arg().required(true))
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the to-psl command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let t_sizes_path = args.get_one::<String>("t_sizes").unwrap();
    let q_sizes_path = args.get_one::<String>("q_sizes").unwrap();
    let output = crate::cmd_pgr::args::get_outfile(args);

    // Load sizes
    let t_sizes = pgr::read_sizes::<usize>(t_sizes_path)?;
    let q_sizes = pgr::read_sizes::<usize>(q_sizes_path)?;

    // Open input
    let reader = pgr::reader(input)?;
    let reader = AxtReader::new(reader);

    // Open output
    let mut writer = pgr::writer(output)?;

    for result in reader {
        let axt = result?;

        // Get sizes
        let q_size = *q_sizes
            .get(&axt.q_name)
            .ok_or_else(|| anyhow::anyhow!("Query size not found for {}", axt.q_name))?;
        let t_size = *t_sizes
            .get(&axt.t_name)
            .ok_or_else(|| anyhow::anyhow!("Target size not found for {}", axt.t_name))?;

        // Prepare coordinates
        // libs/axt.rs returns 0-based half-open coordinates
        let mut q_start = i32::try_from(axt.q_start)
            .map_err(|_| anyhow::anyhow!("q_start {} exceeds i32 range", axt.q_start))?;
        let mut q_end = i32::try_from(axt.q_end)
            .map_err(|_| anyhow::anyhow!("q_end {} exceeds i32 range", axt.q_end))?;
        let q_size_i32 = i32::try_from(q_size)
            .map_err(|_| anyhow::anyhow!("q_size {} exceeds i32 range", q_size))?;
        let t_start_i32 = i32::try_from(axt.t_start)
            .map_err(|_| anyhow::anyhow!("t_start {} exceeds i32 range", axt.t_start))?;
        let t_end_i32 = i32::try_from(axt.t_end)
            .map_err(|_| anyhow::anyhow!("t_end {} exceeds i32 range", axt.t_end))?;
        let q_size_u32 = u32::try_from(q_size)
            .map_err(|_| anyhow::anyhow!("q_size {} exceeds u32 range", q_size))?;
        let t_size_u32 = u32::try_from(t_size)
            .map_err(|_| anyhow::anyhow!("t_size {} exceeds u32 range", t_size))?;

        // axtToPsl.c logic: "if (axt->qStrand == '-') reverseIntRange(&qStart, &qEnd, qSize);"
        // This converts strand-relative coordinates (as in AXT) to positive strand coordinates
        // which pslFromAlign expects (so it can reverse them back internally).
        if axt.q_strand == '-' {
            pgr::reverse_range(&mut q_start, &mut q_end, q_size_i32);
        }

        // Construct strand string for PSL (e.g. "-")
        // Note: PSL usually tracks target strand implicitly as +, so strand field is just query strand?
        // axtToPsl.c: strand[0] = axt->qStrand; strand[1] = '\0';
        // So it's just "+" or "-"
        let strand = axt.q_strand.to_string();

        if let Some(psl) = Psl::from_align(
            &axt.q_name,
            q_size_u32,
            q_start,
            q_end,
            &axt.q_sym,
            &axt.t_name,
            t_size_u32,
            t_start_i32,
            t_end_i32,
            &axt.t_sym,
            &strand,
        ) {
            psl.write_to(&mut writer)?;
        } else {
            log::warn!(
                "skipping alignment (invalid coordinates): {} vs {}",
                axt.q_name,
                axt.t_name
            );
        }
    }

    Ok(())
}
