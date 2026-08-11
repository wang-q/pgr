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
mapped/unmapped counts and the per-base coverage (derived from the mapped
SAM with `pgr sam to-rg` and `pgr rg coverage`).

Mapping is seed-and-verify: the reference's canonical k-mers are indexed
once, each read seeds on its first k-mer, and every candidate position is
verified over the full read length (forward or reverse strand). Reads
matching multiple positions are reported at all of them (`ambiguous=all`).

Notes:
* The reference is FASTA (plain or gzipped); reads are FASTA/FASTQ, one or
  more files (R1/R2 or several single-end files)
* `--outm`/`--outu` write standard SAM (header included)
* `--paired` interleaves two read files as R1/R2 pairs; a pair is mapped
  only when both ends match perfectly, and the SAM carries pair flags
  (0x1/0x2/0x40/0x80), mate coordinates and TLEN (for insert-size
  estimation with `pgr sam ihist`)
* Per-base coverage is derived from the mapped SAM, not accumulated here
* Reads shorter than `--kmer` are unmapped
* `--max-reads` stops after processing N read records (pairs count as two)
* Processing is parallel (rayon) and the output is deterministic: reads
  are written in input order

Examples:
1. Map reads back to an assembly (anchr anchors step):
   pgr asm map UT.fasta R1.fq.gz R2.fq.gz --outm mapped.sam --outu unmapped.sam

2. Derive per-base coverage from the mapped SAM (anchr anchors step):
   pgr sam to-rg mapped.sam | pgr rg coverage stdin -m 2 -o cov.json

3. Use a longer seed k-mer:
   pgr asm map ref.fa reads.fq.gz -k 41 --outm mapped.sam

4. Map paired reads and estimate the insert size (anchr 2_insert_size step):
   pgr asm map UT.fasta R1.fq.gz R2.fq.gz --paired \
       --outm mapped.sam --outu unmapped.sam --max-reads 1000000
   pgr sam ihist mapped.sam -o insert_size.ihist.txt
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
            Arg::new("paired")
                .long("paired")
                .action(clap::ArgAction::SetTrue)
                .help("Map reads as R1/R2 pairs (exactly 2 read files)"),
        )
        .arg(
            Arg::new("max_reads")
                .long("max-reads")
                .num_args(1)
                .value_parser(value_parser!(u64))
                .help("Stop after processing this many read records"),
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
        paired: args.get_flag("paired"),
        max_reads: args.get_one::<u64>("max_reads").copied(),
    };
    let refs = read_fasta(&ref_files).context("failed to read reference")?;
    let stats = map_files(&refs, &read_files, &opts)?;
    eprintln!(
        "Reads in: {}  Mapped: {}  Unmapped: {}  Hits: {}",
        stats.reads_in, stats.mapped, stats.unmapped, stats.hits
    );
    Ok(())
}
