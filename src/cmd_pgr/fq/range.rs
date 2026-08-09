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
using a `.loc` index that is created automatically.

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
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTQ file to process",
        ))
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

/// Execute the range command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut protected: Vec<String> = vec![infile.to_string()];
    if let Some(rgfile) = args.get_one::<String>("rgfile") {
        protected.push(rgfile.clone());
    }
    protected.push(format!("{}.loc", infile));
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, protected.iter().map(|s| s.as_str()))?;

    let ranges = crate::cmd_pgr::args::collect_ranges(args)?;
    let opt_cache = *args.get_one::<std::num::NonZeroUsize>("cache").unwrap();
    let mut cache: lru::LruCache<String, Vec<u8>> = lru::LruCache::new(opt_cache);

    let force_update = args.get_flag("update");
    let (mut reader, loc_of) = pgr::libs::loc::open_fq_indexed(infile, force_update)?;

    let mut out =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
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
