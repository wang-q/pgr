use std::io::{BufRead, BufReader, BufWriter, Write};

/// ```
/// use std::io::BufRead;
/// let reader = pgr::reader("tests/mat/IBPA.list");
/// let mut lines = vec![];
/// for line in reader.lines() {
///     lines.push(line);
/// }
/// assert_eq!(lines.len(), 3);
///
/// let reader = pgr::reader("tests/mat/IBPA.list");
/// assert_eq!(reader.lines().collect::<Vec<_>>().len(), 3);
/// ```
pub fn reader(input: &str) -> Box<dyn BufRead> {
    let reader: Box<dyn BufRead> = if input == "stdin" {
        Box::new(BufReader::new(std::io::stdin()))
    } else {
        let path = std::path::Path::new(input);
        let file = match std::fs::File::open(path) {
            Err(why) => panic!("could not open {}: {}", path.display(), why),
            Ok(file) => file,
        };

        if path.extension() == Some(std::ffi::OsStr::new("gz")) {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        }
    };

    reader
}


/// ```
/// let lines = pgr::read_lines("tests/mat/IBPA.list");
/// assert_eq!(lines.len(), 3);
/// ```
pub fn read_lines(input: &str) -> Vec<String> {
    let mut reader = reader(input);
    let mut s = String::new();
    reader.read_to_string(&mut s).expect("Read error");
    s.lines().map(|s| s.to_string()).collect::<Vec<String>>()
}

pub fn writer(output: &str) -> Box<dyn Write> {
    let writer: Box<dyn Write> = if output == "stdout" {
        Box::new(BufWriter::new(std::io::stdout()))
    } else {
        Box::new(BufWriter::new(std::fs::File::create(output).unwrap()))
    };

    writer
}

/// ```
/// let sizes = intspan::read_sizes("tests/pgr/pseudopig.sizes");
/// assert_eq!(sizes.len(), 2);
/// assert_eq!(*sizes.get("pig2").unwrap(), 22929);
/// ```
pub fn read_sizes(input: &str) -> BTreeMap<String, i32> {
    let mut sizes: BTreeMap<String, i32> = BTreeMap::new();

    for line in read_lines(input) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 2 {
            sizes.insert(fields[0].to_string(), fields[1].parse::<i32>().unwrap());
        }
    }

    sizes
}
