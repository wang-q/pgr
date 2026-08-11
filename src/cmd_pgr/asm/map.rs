use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::map::{map_files, read_fasta, MapOptions};

/// Build the clap subcommand for map.
pub fn make_subcommand() -> Command {
    Command::new("map")
        .about("Maps reads to a reference requiring perfect matches (bbmap perfectmode)")
        .after_help(
            r###"
This command maps reads to a reference (typically an assembly) requiring
every read to match exactly: no mismatches and no gaps, mirroring BBTools
`bbwrap.sh perfectmode maxindel=0 strictmaxindel`. This replaces the
bbwrap call of the anchr `anchors` flow, whose downstream only needs the
mapped/unmapped counts and the per-base coverage.

Mapping is seed-and-verify: the reference's canonical k-mers are indexed
once, each read seeds on its first k-mer, and every candidate position is
verified over the full read length (forward or reverse strand). Reads
matching multiple positions are reported at all of them (`ambiguous=all`).

Notes:
* The reference is FASTA (plain or gzipped); reads are FASTA/FASTQ, one or
  more files (R1/R2 or several single-end files)
* `--outm`/`--outu` write standard SAM (header included); `--basecov`
  writes `RefName Pos Coverage` (0-based) with coverage > 0 only
* Reads shorter than `--kmer` are unmapped
* Processing is parallel (rayon) and the output is deterministic: reads
  are written in input order

Examples:
1. Map reads back to an assembly and report coverage (anchr anchors step):
   pgr asm map UT.fasta R1.fq.gz R2.fq.gz --outm mapped.sam --outu unmapped.sam --basecov basecov.txt

2. Only the coverage profile:
   pgr asm map UT.fasta reads.fq.gz --basecov basecov.txt

3. Use a longer seed k-mer:
   pgr asm map ref.fa reads.fq.gz -k 41 --outm mapped.sam
"###,
        )
        .arg(
            Arg::new("ref")
                .index(1)
                .required(true)
                .help("Reference FASTA file(s) to map against"),
        )
        .arg(
            Arg::new("infiles")
                .index(2)
                .required(true)
                .num_args(1..)
                .help("Read file(s) to map (FASTA/FASTQ)"),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("31")
                .value_parser(value_parser!(usize))
                .help("Seed k-mer length (1..=64)"),
        )
        .arg(
            Arg::new("outm")
                .long("outm")
                .num_args(1)
                .help("SAM output of perfectly matched reads"),
        )
        .arg(
            Arg::new("outu")
                .long("outu")
                .num_args(1)
                .help("SAM output of unmapped reads"),
        )
        .arg(
            Arg::new("basecov")
                .long("basecov")
                .num_args(1)
                .help("Per-base coverage output (RefName Pos Coverage)"),
        )
        .arg(crate::cmd_pgr::args::parallel_arg_with_default("8"))
}

/// Execute the map command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let ref_files: Vec<String> = args.get_many::<String>("ref").unwrap().cloned().collect();
    let read_files: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let opts = MapOptions {
        k: *args.get_one::<usize>("kmer").unwrap(),
        outm: args.get_one::<String>("outm").cloned(),
        outu: args.get_one::<String>("outu").cloned(),
        basecov: args.get_one::<String>("basecov").cloned(),
    };
    let refs = read_fasta(&ref_files).context("failed to read reference")?;
    let stats = map_files(&refs, &read_files, &opts)?;
    eprintln!(
        "Reads in: {}  Mapped: {}  Unmapped: {}  Hits: {}",
        stats.reads_in, stats.mapped, stats.unmapped, stats.hits
    );
    Ok(())
}
