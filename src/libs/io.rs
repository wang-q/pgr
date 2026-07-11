use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::str::FromStr;

use anyhow::Context;

/// Random-access reader for subsequence extraction by name and 0-based range.
///
/// Implementors provide `read_sequence(name, start, end)` returning the
/// substring `[start, end)` of sequence `name`. Used by chain/net algorithms
/// to decouple from any specific on-disk format (e.g. 2bit, indexed FASTA).
pub trait SequenceReader {
    /// Read `[start, end)` from sequence `name`. `None` means "from start" /
    /// "to end". Returns the sequence as a `String` (DNA bases, possibly
    /// soft-masked).
    fn read_sequence(
        &mut self,
        name: &str,
        start: Option<usize>,
        end: Option<usize>,
    ) -> anyhow::Result<String>;
}

/// Open a buffered reader for `input` (`stdin` or a file path, `.gz` supported).
///
/// ```ignore
/// # use std::io::BufRead;
/// let reader = pgr::reader("tests/mat/IBPA.list").unwrap();
/// let mut lines = vec![];
/// for line in reader.lines() {
///     lines.push(line);
/// }
/// assert_eq!(lines.len(), 3);
///
/// let reader = pgr::reader("tests/mat/IBPA.list").unwrap();
/// assert_eq!(reader.lines().collect::<Vec<_>>().len(), 3);
/// ```
pub fn reader(input: &str) -> anyhow::Result<Box<dyn BufRead>> {
    if input == "stdin" {
        return Ok(Box::new(BufReader::new(std::io::stdin())));
    }

    let path = Path::new(input);
    let file = File::open(path).with_context(|| format!("could not open {}", path.display()))?;

    if path.extension() == Some(std::ffi::OsStr::new("gz")) {
        Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(
            file,
        ))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Read all lines from `input` (`stdin` or a file path, `.gz` supported).
///
/// ```ignore
/// let lines = pgr::read_lines("tests/mat/IBPA.list").unwrap();
/// assert_eq!(lines.len(), 3);
/// ```
pub fn read_lines(input: &str) -> anyhow::Result<Vec<String>> {
    let mut reader = reader(input)?;
    let mut s = String::new();
    reader.read_to_string(&mut s).context("read error")?;
    Ok(s.lines().map(|s| s.to_string()).collect())
}

/// Safely read a runlist JSON file and convert to IntSpan map.
/// Replaces intspan::read_json + intspan::json2set which panic on errors.
pub fn read_runlist(path: &str) -> anyhow::Result<BTreeMap<String, intspan::IntSpan>> {
    let mut reader = reader(path)?;
    let mut s = String::new();
    reader
        .read_to_string(&mut s)
        .with_context(|| format!("failed to read runlist: {}", path))?;
    let json: BTreeMap<String, serde_json::Value> = serde_json::from_str(&s)
        .with_context(|| format!("failed to parse runlist JSON: {}", path))?;
    let mut set: BTreeMap<String, intspan::IntSpan> = BTreeMap::new();
    for (chr, value) in &json {
        let s = value
            .as_str()
            .ok_or_else(|| anyhow::anyhow!("runlist value for {} is not a string", chr))?;
        set.insert(chr.clone(), intspan::IntSpan::from(s));
    }
    Ok(set)
}

/// Buffered writer that flushes on drop and reports flush errors to stderr.
///
/// Wraps a `BufWriter<Box<dyn Write>>` so that `BufWriter`'s silent flush-on-drop
/// behavior is replaced with a best-effort flush that emits a warning to stderr
/// if flushing fails (e.g. broken pipe, disk full). Callers that need to
/// propagate flush errors should call `flush()?` explicitly before the writer
/// goes out of scope.
pub struct PgrWriter {
    inner: BufWriter<Box<dyn Write>>,
}

impl Write for PgrWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.inner.write(buf)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

impl Drop for PgrWriter {
    fn drop(&mut self) {
        if let Err(e) = self.inner.flush() {
            // Best-effort: report to stderr but never panic from Drop.
            let _ = writeln!(
                std::io::stderr(),
                "pgr: warning: failed to flush writer: {}",
                e
            );
        }
    }
}

/// Open a buffered writer for `output` (`stdout` or a file path).
///
/// Returns a [`PgrWriter`] which flushes on drop with a stderr warning on
/// failure. To get a `Box<dyn Write>` (e.g. for storing in a heterogeneous
/// collection), wrap the result with `Box::new(pgr::writer(...)?)`.
pub fn writer(output: &str) -> anyhow::Result<PgrWriter> {
    let boxed: Box<dyn Write> = if output == "stdout" {
        Box::new(std::io::stdout())
    } else {
        let file = File::create(output).with_context(|| format!("could not create {}", output))?;
        Box::new(file)
    };
    Ok(PgrWriter {
        inner: BufWriter::new(boxed),
    })
}

/// Read a `name<TAB>size` sizes file into a map with the requested value type.
///
/// Lines are split on whitespace; lines with fewer than 2 fields are skipped.
///
/// ```ignore
/// let sizes = pgr::read_sizes::<u64>("tests/pgr/pseudopig.sizes").unwrap();
/// assert_eq!(sizes.len(), 2);
/// assert_eq!(*sizes.get("pig2").unwrap(), 22929);
/// ```
pub fn read_sizes<T>(input: &str) -> anyhow::Result<BTreeMap<String, T>>
where
    T: FromStr,
    T::Err: std::error::Error + Send + Sync + 'static,
{
    let mut sizes: BTreeMap<String, T> = BTreeMap::new();

    for line in read_lines(input)? {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 2 {
            let size: T = fields[1]
                .parse()
                .with_context(|| format!("invalid size value: {}", fields[1]))?;
            sizes.insert(fields[0].to_string(), size);
        }
    }

    Ok(sizes)
}

/// Check whether a file is BGZF-compressed by inspecting the header bytes.
///
/// Returns `false` if the file cannot be read or is too short.
pub fn is_bgzf(path: &str) -> bool {
    let mut f = match File::open(path) {
        Ok(f) => f,
        Err(_) => return false,
    };
    let mut hdr = [0u8; 18];
    if f.read_exact(&mut hdr).is_err() {
        return false;
    }
    // BGZF: gzip magic (1f 8b 08 04), XLEN=6 at [10..12], "BC" at [12..14], SLEN=2 at [14..16]
    hdr[0] == 0x1f
        && hdr[1] == 0x8b
        && hdr[2] == 0x08
        && hdr[3] == 0x04
        && hdr[10] == 0x06
        && hdr[11] == 0x00
        && hdr[12] == b'B'
        && hdr[13] == b'C'
        && hdr[14] == 0x02
        && hdr[15] == 0x00
}

/// List files in `dir` with the given `extension` (non-recursive).
pub fn list_files_ext(dir: &str, extension: &str) -> Vec<String> {
    let mut files = Vec::new();
    let dir_path = Path::new(dir);

    if dir_path.is_dir() {
        if let Ok(entries) = std::fs::read_dir(dir_path) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_file() {
                    if let Some(ext) = path.extension() {
                        if ext == extension {
                            files.push(path.to_string_lossy().into_owned());
                        }
                    }
                }
            }
        }
    }

    files
}

/// Get the basename of `file_path` (the part before the first `.`).
pub fn get_basename(file_path: &str) -> Option<String> {
    let path = Path::new(file_path);
    path.file_stem()
        .and_then(std::ffi::OsStr::to_str)
        .map(|s| s.split('.').next().unwrap_or(s).to_string())
}

/// Like [`get_basename`] but errors with a contextual message on failure.
pub fn basename_or_err(file_path: &str) -> anyhow::Result<String> {
    get_basename(file_path)
        .ok_or_else(|| anyhow::anyhow!("failed to get basename of: {}", file_path))
}

/// Convert a path to `&str`, erroring with a contextual message if not UTF-8.
pub fn path_to_str(path: &std::path::Path) -> anyhow::Result<&str> {
    path.to_str()
        .ok_or_else(|| anyhow::anyhow!("path is not utf-8: {}", path.display()))
}

/// Return the prefix of `name` up to the first separator
/// (`' '`, `'.'`, `','`, or `'-'`); returns `name` unchanged if no separator.
pub fn simplify_name(name: &str) -> &str {
    match name.find(&[' ', '.', ',', '-'][..]) {
        Some(i) => &name[..i],
        None => name,
    }
}

/// Sanitize a string for safe use as a filename by replacing path separators
/// and other special characters with `_`, then collapsing consecutive `_`.
pub fn sanitize_filename(name: &str) -> String {
    name.replace(['/', '\\', '(', ')', ':'], "_")
        .replace("__", "_")
}

/// Get the path of the currently executing binary as a String.
pub fn current_exe_string() -> anyhow::Result<String> {
    Ok(std::env::current_exe()?.display().to_string())
}

/// Read the first column of `path` (one name per line) into a collection.
///
/// Lines are split on whitespace; only the first field is kept. Empty lines
/// are skipped. Order is preserved. Use `read_names::<Vec<String>>` for a
/// vector or `read_names::<HashSet<String>>` for a set.
///
/// ```ignore
/// let names: Vec<String> = pgr::libs::io::read_names("tests/mat/IBPA.list").unwrap();
/// assert_eq!(names.len(), 3);
/// ```
pub fn read_names<T: FromIterator<String>>(path: &str) -> anyhow::Result<T> {
    Ok(read_lines(path)?
        .into_iter()
        .filter_map(|line| {
            let line = line.trim();
            if line.is_empty() {
                None
            } else {
                line.split_whitespace().next().map(|s| s.to_string())
            }
        })
        .collect())
}

/// Read a replacement TSV file into a map of `name -> Vec<replacement_names>`.
///
/// Each line is split on tabs: the first field is the key, remaining fields
/// are replacement names. Multiple lines with the same key accumulate
/// replacements (append semantics). A single-field line (key only) inserts
/// an empty Vec, which callers may interpret as a "delete" directive.
pub fn read_replace_tsv(path: &str) -> anyhow::Result<BTreeMap<String, Vec<String>>> {
    let mut replaces: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in read_lines(path)? {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let key = fields[0].to_string();
        let others: Vec<String> = fields
            .get(1..)
            .map(|rest| rest.iter().map(|s| s.to_string()).collect())
            .unwrap_or_default();
        replaces.entry(key).or_default().extend(others);
    }
    Ok(replaces)
}

/// Borrowed line iterator over a `BufRead`, yielding `String` with the
/// trailing `\n` (and `\r`) stripped. Unlike `BufRead::lines`, the reader
/// is borrowed for the lifetime of the iterator (zero-allocation handle).
pub struct LinesRef<'a, B: 'a> {
    pub(crate) buf: &'a mut B,
}

impl<'a, B: BufRead> Iterator for LinesRef<'a, B> {
    type Item = std::io::Result<String>;

    fn next(&mut self) -> Option<std::io::Result<String>> {
        let mut buf = String::new();
        match self.buf.read_line(&mut buf) {
            Ok(0) => None,
            Ok(_n) => {
                if buf.ends_with('\n') {
                    buf.pop();
                    if buf.ends_with('\r') {
                        buf.pop();
                    }
                }
                Some(Ok(buf))
            }
            Err(e) => Some(Err(e)),
        }
    }
}
