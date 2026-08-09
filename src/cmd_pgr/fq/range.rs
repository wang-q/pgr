use anyhow::{bail, Context};
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use pgr::libs::ds::Range;
use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use std::io::{Cursor, Write};

/// Build the clap subcommand for range.
pub fn make_subcommand() -> Command {
    Command::new("range")
        .about("Extracts FASTQ records by name or region")
        .after_help(
            r###"
This command extracts FASTQ records by read name (or a region within a read)
using a `.loc` index that is created automatically. Paired-end input (two
files) writes each mate to its own output with `--outfile-2`.

Notes:
* Index format matches `fa range`: name, plain offset, record size
* Read names with `/1` `/2` suffixes are matched by their pair name
* Interleaved reads with identical names are both returned, in order
* `name:start-end` returns the subsequence of both sequence and quality
* Supports plain text and BGZF (.gz) files (plain gzip is not seekable)
* The index is rebuilt when the input is newer (or with --update)

Examples:
1. Extract whole records:
   pgr fq range in.fq read1 read2

2. Extract a region of a read:
   pgr fq range in.fq "read1:10-100"

3. Extract from a name list with a larger cache:
   pgr fq range in.fq -r names.txt -c 10

4. Force rebuild the index:
   pgr fq range in.fq read1 --update

5. Paired-end extraction (mate via --mate, same ranges, separate outputs):
   pgr fq range R1.fq --mate R2.fq read1 -o r1.out.fq --outfile-2 r2.out.fq
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTQ file to process",
        ))
        .arg(
            Arg::new("mate")
                .long("mate")
                .num_args(1)
                .help("Second mate file (paired-end)"),
        )
        .arg(crate::cmd_pgr::args::ranges_arg())
        .arg(crate::cmd_pgr::args::rgfile_arg())
        .arg(
            Arg::new("cache")
                .long("cache")
                .short('c')
                .value_parser(value_parser!(std::num::NonZeroUsize))
                .num_args(1)
                .default_value("1")
                .help("Set the capacity of the LRU cache"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("outfile_2")
                .long("outfile-2")
                .num_args(1)
                .help("Output filename for the second mate (paired-end)"),
        )
        .arg(
            Arg::new("update")
                .long("update")
                .short('u')
                .action(ArgAction::SetTrue)
                .help("Force update the .loc index file"),
        )
}

/// Writes a FASTQ record, keeping the name/comment and the `+` line plain.
fn write_fq_record<W: Write>(
    w: &mut W,
    rec: &SeqRecord,
    seq: &[u8],
    qual: &[u8],
) -> std::io::Result<()> {
    let comment = rec.comment();
    if comment.is_empty() {
        writeln!(w, "@{}", rec.name())?;
    } else {
        writeln!(w, "@{} {}", rec.name(), comment)?;
    }
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")?;
    Ok(())
}

/// Extracts the requested ranges from one indexed FASTQ file into `out`.
fn extract_file(
    infile: &str,
    ranges: &[String],
    cache_capacity: std::num::NonZeroUsize,
    force_update: bool,
    mut out: impl Write,
) -> anyhow::Result<()> {
    let mut cache: lru::LruCache<String, Vec<u8>> = lru::LruCache::new(cache_capacity);
    let (mut reader, loc_of) = pgr::libs::loc::open_fq_indexed(infile, force_update)?;
    for el in ranges.iter() {
        let rg = Range::from_str(el);
        let hits = pgr::libs::loc::query_fq_locs(&loc_of, rg.chr());
        if hits.is_empty() {
            log::warn!("{} for [{}] not found in the .loc index file", rg.chr(), el);
            continue;
        }
        for (key, offset, size) in hits {
            let data = if let Some(d) = cache.get(key) {
                d.clone()
            } else {
                let d = pgr::libs::loc::read_offset(&mut reader, offset, size)?;
                cache.put(key.to_string(), d.clone());
                d
            };
            let mut fq_in = SeqReader::from_reader(Box::new(Cursor::new(data)));
            let mut rec = SeqRecord::new();
            if !fq_in.read_record(&mut rec)? || !rec.is_fastq() {
                bail!("malformed FASTQ record for {}", key);
            }
            let start = *rg.start() as usize;
            let end = *rg.end() as usize;
            if start == 0 {
                write_fq_record(&mut out, &rec, rec.sequence(), rec.quality_scores())?;
            } else {
                if end < start || end > rec.sequence().len() {
                    bail!("slice error for [{}]", el);
                }
                write_fq_record(
                    &mut out,
                    &rec,
                    &rec.sequence()[start - 1..end],
                    &rec.quality_scores()[start - 1..end],
                )?;
            }
        }
    }
    out.flush()?;
    Ok(())
}

/// Execute the range command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let mate = args.get_one::<String>("mate").map(String::as_str);
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let outfile_2 = args.get_one::<String>("outfile_2").map(String::as_str);
    if mate.is_none() && outfile_2.is_some() {
        bail!("--outfile-2 requires two input files (paired-end)");
    }
    if mate.is_some() && outfile_2.is_none() {
        bail!("--mate requires --outfile-2 (paired-end output)");
    }
    if outfile_2 == Some("stdout") {
        bail!("--outfile-2 must be a file path, not stdout");
    }
    let mut protected: Vec<String> = vec![infile.to_string()];
    if let Some(m) = mate {
        protected.push(m.to_string());
    }
    if let Some(rgfile) = args.get_one::<String>("rgfile") {
        protected.push(rgfile.clone());
    }
    protected.push(format!("{infile}.loc"));
    if let Some(m) = mate {
        protected.push(format!("{m}.loc"));
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, protected.iter().map(|s| s.as_str()))?;
    if let Some(o2) = outfile_2 {
        crate::cmd_pgr::args::ensure_outfile_distinct(o2, protected.iter().map(|s| s.as_str()))?;
        if outfile != "stdout" && pgr::libs::io::same_path(outfile, o2) {
            bail!("output files must be distinct: {} and {}", outfile, o2);
        }
    }

    let ranges = crate::cmd_pgr::args::collect_ranges(args)?;
    let cache_capacity = *args.get_one::<std::num::NonZeroUsize>("cache").unwrap();
    let force_update = args.get_flag("update");
    if let Some(mate) = mate {
        let mut out1 = pgr::writer(outfile)
            .with_context(|| format!("Failed to open writer for {}", outfile))?;
        extract_file(infile, &ranges, cache_capacity, force_update, &mut out1)?;
        let mut out2 = pgr::writer(outfile_2.unwrap())
            .with_context(|| format!("Failed to open writer for {}", outfile_2.unwrap()))?;
        extract_file(mate, &ranges, cache_capacity, force_update, &mut out2)?;
    } else {
        let mut out = pgr::writer(outfile)
            .with_context(|| format!("Failed to open writer for {}", outfile))?;
        extract_file(infile, &ranges, cache_capacity, force_update, &mut out)?;
    }
    Ok(())
}
