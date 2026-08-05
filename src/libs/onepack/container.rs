//! ONEcode binary container (header, line I/O, object index, footer).
//!
//! This module implements the read and write paths of FastGA's ONEcode
//! container as used by `.1aln` files. It is the P3 building block: readers
//! and writers both operate against a [`Schema`] (built from the embedded
//! header text or from [`aln_schema_text`]).
//!
//! Layout of a binary ONEcode file:
//!
//! ```text
//! <ASCII prolog>       1-line, optional 2-line, '.', '!', '<', '>' lines, '~' schema
//! <$-line>             "$ <isBig>" marks the file as binary
//! <binary data> ...    one binary line per record; ends with a blank line '\n'
//! <footer>             '#','@','+','%' counts, '&' object index, ';' codecs
//! ^                    end-of-footer marker
//! <footer offset>      8-byte off_t pointing at the start of the footer
//! ```

use anyhow::{anyhow, bail, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

use super::ltf::{int_get, int_put};
use super::schema::{
    add_definition, binary_pack, header_schema_text, parse_schema_text, FieldType, LineInfo,
    Schema, MAJOR, MINOR,
};
use super::vc::VcCodec;

/// The `.1aln` schema text (FastGA `alncode.c` `alnSchemaText`, the `aln`
/// primary type portion). The GDB/skeleton section is included so that
/// `.1aln` files carrying a skeleton parse correctly.
pub fn aln_schema_text() -> &'static str {
    "P 3 aln\n\
     D t 1 3 INT\n\
     O g 0\n\
     G S\n\
     O S 1 6 STRING\n\
     D G 1 3 INT\n\
     D C 1 3 INT\n\
     O a 0\n\
     G A\n\
     D p 2 3 INT 3 INT\n\
     O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT\n\
     D L 2 3 INT 3 INT\n\
     D R 0\n\
     D D 1 3 INT\n\
     D T 1 8 INT_LIST\n\
     D X 1 8 INT_LIST\n\
     D Q 1 3 INT\n\
     D E 1 3 INT\n\
     D Z 1 6 STRING\n\
     D U 1 3 INT\n"
}

/// A single field value on a line.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Field {
    Int(i64),
    Real(f64),
    Char(u8),
}

/// The list payload of a line.
#[derive(Clone, Debug, PartialEq)]
pub enum List {
    /// `INT_LIST`: decoded integers.
    Ints(Vec<i64>),
    /// `REAL_LIST`.
    Reals(Vec<f64>),
    /// `STRING`/`DNA`: raw bytes.
    Bytes(Vec<u8>),
    /// `STRING_LIST`: NUL-terminated strings packed together.
    Strings(Vec<Vec<u8>>),
}

/// Provenance entry (`!` line).
#[derive(Clone, Debug)]
pub struct Provenance {
    pub program: String,
    pub version: String,
    pub command: String,
    pub date: String,
}

/// Reference entry (`<` line).
#[derive(Clone, Debug)]
pub struct Reference {
    pub filename: String,
    pub count: i64,
}

/// A decoded line read from a ONEcode file.
#[derive(Debug)]
pub struct Line {
    /// The ASCII line type character.
    pub line_type: u8,
    /// Field values in schema order.
    pub fields: Vec<Field>,
    /// The list payload, if the line type has a list field.
    pub list: Option<List>,
}

impl Line {
    /// The i-th field as an integer (panics if it is not an integer field).
    pub fn int(&self, i: usize) -> i64 {
        match self.fields[i] {
            Field::Int(x) => x,
            _ => panic!("field {i} is not an INT"),
        }
    }

    /// The list length (0 if no list).
    pub fn list_len(&self) -> usize {
        match &self.list {
            Some(List::Ints(v)) => v.len(),
            Some(List::Reals(v)) => v.len(),
            Some(List::Bytes(v)) => v.len(),
            Some(List::Strings(v)) => v.len(),
            None => 0,
        }
    }
}

/// A ONEcode binary file reader.
///
/// Mirrors FastGA's `oneFileOpenRead`/`oneReadLine`: the header (ASCII prolog,
/// schema, provenance, references) and the footer (counts, object index,
/// codecs) are read in one pass, then the file is positioned at the start of
/// the (binary or ASCII) data section.
pub struct Reader {
    file: BufReader<File>,
    schema: Schema,
    binary_unpack: [u8; 256],
    /// Per-line-type list codec (from footer `;` lines).
    codecs: Vec<Option<VcCodec>>,
    /// Per-line-type object index (from footer `&` lines).
    index: Vec<Vec<i64>>,
    /// Per-line-type counts.
    counts: Vec<Counts>,
    /// Byte offset where the data section begins (`startOff` in FastGA).
    data_start: u64,
    /// Byte offset where the footer begins (`footOff` in FastGA); the data
    /// section spans `[data_start, foot_off)`.
    foot_off: u64,
    /// Whether the file is binary.
    pub is_binary: bool,
    /// Provenance entries.
    pub provenance: Vec<Provenance>,
    /// Reference entries.
    pub references: Vec<Reference>,
    /// Current line number.
    pub line: i64,
    /// One-line pushback buffer (used by the record layer to re-examine a
    /// line that starts the next record).
    pending: Option<Line>,
}

#[derive(Clone, Copy, Default, Debug)]
pub struct Counts {
    pub count: i64,
    pub max: i64,
    pub total: i64,
}

impl Reader {
    /// Open a `.1aln`/ONEcode file for reading and parse its header+footer.
    pub fn open(path: &str) -> Result<Reader> {
        let file = BufReader::new(File::open(path)?);
        // Seed the schema with the universal header/footer line types so the
        // footer's `&`/`;`/`#`/`@`/`+`/`%` lines can be decoded.
        let schema = parse_schema_text(header_schema_text())?;
        let mut binary_unpack = [0u8; 256];
        for &code in &schema.defn_order {
            if code & 0x80 == 0 {
                let t = code as u8;
                binary_unpack[binary_pack(t) as usize] = t;
                binary_unpack[(binary_pack(t) + 1) as usize] = t;
            }
        }
        let mut r = Reader {
            file,
            schema,
            binary_unpack,
            codecs: (0..128).map(|_| None).collect(),
            index: vec![Vec::new(); 128],
            counts: vec![Counts::default(); 128],
            data_start: 0,
            foot_off: 0,
            is_binary: false,
            provenance: Vec::new(),
            references: Vec::new(),
            line: 0,
            pending: None,
        };
        r.read_header()?;
        Ok(r)
    }

    fn read_byte(&mut self) -> Result<u8> {
        let mut b = [0u8; 1];
        let n = self.file.read(&mut b)?;
        if n == 0 {
            bail!("unexpected end of file");
        }
        Ok(b[0])
    }

    /// Read the ASCII header (and footer for binary files).
    fn read_header(&mut self) -> Result<()> {
        // First line must be "1 <len> <type> <major> <minor>".
        let mut first = String::new();
        self.file.read_line(&mut first)?;
        let mut it = first.split_whitespace();
        let n = it.next().ok_or_else(|| anyhow!("missing '1' header"))?;
        if n != "1" {
            bail!("not a ONEcode file (first token '{n}')");
        }
        let _dlen: usize = it
            .next()
            .ok_or_else(|| anyhow!("missing type length"))?
            .parse()?;
        let ftype = it
            .next()
            .ok_or_else(|| anyhow!("missing file type"))?
            .to_string();
        let major: i64 = it.next().ok_or_else(|| anyhow!("missing major"))?.parse()?;
        let minor: i64 = it.next().ok_or_else(|| anyhow!("missing minor"))?.parse()?;
        if major != MAJOR {
            bail!("major version file {major} != code {MAJOR}");
        }
        if minor > MINOR {
            bail!("minor version file {minor} > code {MINOR}");
        }
        self.schema.file_type = ftype;

        // Read header lines until an alphabetic (data) line or a blank line.
        loop {
            let peek = self.peek_byte()?;
            if peek == b'\n' {
                // Blank line: end of header/data boundary.
                return Ok(());
            }
            if peek.is_ascii_alphabetic() {
                // Start of data.
                return Ok(());
            }
            let lt = self.read_header_line()?;
            if lt == b'$' {
                // Binary: record where the data section begins, then read the
                // footer (whose offset is stored at EOF) and seek back.
                self.is_binary = true;
                self.data_start = self.file.stream_position()?;
                self.read_footer()?;
                self.file.seek(SeekFrom::Start(self.data_start))?;
                return Ok(());
            }
        }
    }

    fn peek_byte(&mut self) -> Result<u8> {
        let b = self.read_byte()?;
        self.file.seek(SeekFrom::Current(-1))?;
        Ok(b)
    }

    /// Read one ASCII header line and process it. Returns the line type.
    fn read_header_line(&mut self) -> Result<u8> {
        let lt = self.read_byte()?;
        let mut body = String::new();
        self.file.read_line(&mut body)?;
        let mut c = AsciiCursor::new(&body);
        match lt {
            b'2' => {
                let _len = c.read_int()?;
                self.schema.sub_type = Some(c.read_str()?.to_string());
            }
            b'.' => {}
            b'!' => {
                // provenance: ! 4 <plen> <prog> <vlen> <vers> <clen> <cmd> <dlen> <date>
                // Strings may contain spaces, so parse with length prefixes.
                let parts = parse_len_strings(&body)?;
                if parts.len() >= 4 {
                    self.provenance.push(Provenance {
                        program: parts[0].clone(),
                        version: parts[1].clone(),
                        command: parts[2].clone(),
                        date: parts[3].clone(),
                    });
                }
            }
            b'<' => {
                // < <len> filename <count>
                let len = c.read_int()? as usize;
                let filename = String::from_utf8_lossy(c.read_chars(len)?).to_string();
                let count = c.read_int()?;
                self.references.push(Reference { filename, count });
            }
            b'>' => {
                // > <len> filename (deferred)
                let len = c.read_int()? as usize;
                let _ = c.read_chars(len)?;
            }
            b'~' => {
                // ~ O|D|G schema definition line. The field count `n` limits the
                // number of `<len> <TYPE>` pairs; anything after is a comment.
                let kind = c.read_char()?;
                let t = c.read_char()?;
                if kind == b'G' {
                    let _n = c.read_int()?;
                    self.schema.defn_order.push((t as i32) | 0x80);
                } else {
                    let n = c.read_int()? as usize;
                    let mut fts = Vec::with_capacity(n);
                    for _k in 0..n {
                        let _len = c.read_int()?;
                        let name = c.read_str()?;
                        let ft = FieldType::from_name(name)
                            .ok_or_else(|| anyhow!("unknown field type {name}"))?;
                        fts.push(ft);
                    }
                    add_definition(&mut self.schema, kind, t, fts, None)?;
                    self.binary_unpack[binary_pack(t) as usize] = t;
                    self.binary_unpack[(binary_pack(t) + 1) as usize] = t;
                }
            }
            b'&' | b';' | b'#' | b'@' | b'+' | b'%' => {
                // Only valid in the footer; handled by read_footer_line.
                bail!("header line type '{}' outside footer", lt as char);
            }
            _ => {}
        }
        Ok(lt)
    }

    /// Seek to and read the footer (binary files).
    fn read_footer(&mut self) -> Result<()> {
        // The footer offset is stored as 8 little-endian bytes at EOF.
        let end = self.file.seek(SeekFrom::End(0))?;
        let mut off = [0u8; 8];
        self.file.seek(SeekFrom::Start(end - 8))?;
        self.file.read_exact(&mut off)?;
        let foot_off = i64::from_le_bytes(off);
        self.foot_off = foot_off as u64;
        self.file.seek(SeekFrom::Start(foot_off as u64))?;
        // Read footer lines until '^'. ASCII count lines (`#`/`@`/`+`/`%`)
        // are interleaved with binary `&`/`;` lines (high bit set).
        loop {
            let peek = self.peek_byte()?;
            if peek == b'^' {
                self.read_byte()?;
                self.file.read_line(&mut String::new())?;
                return Ok(());
            }
            if peek & 0x80 != 0 {
                let line = self
                    .read_binary_line_raw()?
                    .ok_or_else(|| anyhow!("empty footer line"))?;
                self.process_footer_binary(&line)?;
            } else {
                self.read_footer_line()?;
            }
        }
    }

    /// Apply a binary `&` (object index) or `;` (list codec) footer line.
    fn process_footer_binary(&mut self, line: &Line) -> Result<()> {
        match line.line_type {
            b'&' => {
                let t = match line.fields[0] {
                    Field::Char(c) => c,
                    Field::Int(i) => i as u8,
                    _ => return Ok(()),
                };
                if let Some(List::Ints(idx)) = &line.list {
                    self.index[t as usize] = idx.clone();
                }
            }
            b';' => {
                let t = match line.fields[0] {
                    Field::Char(c) => c,
                    Field::Int(i) => i as u8,
                    _ => return Ok(()),
                };
                if let Some(List::Bytes(bytes)) = &line.list {
                    self.codecs[t as usize] = Some(VcCodec::deserialize(bytes)?);
                }
            }
            _ => {}
        }
        Ok(())
    }

    /// Read one ASCII footer line (counts, index, codec).
    fn read_footer_line(&mut self) -> Result<()> {
        let lt = self.read_byte()?;
        let mut body = String::new();
        self.file.read_line(&mut body)?;
        let fields: Vec<&str> = body.split_whitespace().collect();
        match lt {
            b'#' => {
                let t = fields[0].as_bytes()[0];
                self.counts[t as usize].count = fields[1].parse()?;
            }
            b'@' => {
                let t = fields[0].as_bytes()[0];
                self.counts[t as usize].max = fields[1].parse()?;
            }
            b'+' => {
                let t = fields[0].as_bytes()[0];
                self.counts[t as usize].total = fields[1].parse()?;
            }
            _ => {}
        }
        Ok(())
    }

    /// Read the next data line. Returns `None` at end of data (blank line).
    pub fn read_line(&mut self) -> Result<Option<Line>> {
        if let Some(l) = self.pending.take() {
            return Ok(Some(l));
        }
        if self.is_binary {
            self.read_binary_line()
        } else {
            self.read_ascii_line()
        }
    }

    /// Push a line back so the next `read_line` returns it first.
    ///
    /// Used by the record layer to retain a line ('A') that starts the next
    /// alignment while the current record is being finalized.
    pub fn unread_line(&mut self, line: Line) {
        self.pending = Some(line);
    }

    /// Per-line-type counts (from the footer `#`/`@`/`+` lines).
    pub fn counts(&self) -> &[Counts] {
        self.counts.as_slice()
    }

    fn read_ascii_line(&mut self) -> Result<Option<Line>> {
        let t = self.read_byte()?;
        if t == b'\n' {
            return Ok(None);
        }
        self.line += 1;
        let li = match self.schema.line_info(t) {
            Some(li) => li.clone(),
            None => bail!("unknown line type {}", t as char),
        };
        let mut body = String::new();
        self.file.read_line(&mut body)?;
        let mut toks = body.split_whitespace();
        let mut fields = Vec::with_capacity(li.fields.len());
        let mut list = None;
        for (i, ft) in li.fields.iter().copied().enumerate() {
            if i == li.list_field {
                let len: usize = toks
                    .next()
                    .ok_or_else(|| anyhow!("missing list length"))?
                    .parse()?;
                list = Some(self.read_ascii_list(ft, len, &mut toks)?);
                fields.push(Field::Int(len as i64));
            } else {
                match ft {
                    FieldType::Int => {
                        fields.push(Field::Int(
                            toks.next().ok_or_else(|| anyhow!("missing int"))?.parse()?,
                        ));
                    }
                    FieldType::Real => {
                        fields.push(Field::Real(
                            toks.next()
                                .ok_or_else(|| anyhow!("missing real"))?
                                .parse()?,
                        ));
                    }
                    FieldType::Char => {
                        fields.push(Field::Char(
                            toks.next()
                                .ok_or_else(|| anyhow!("missing char"))?
                                .as_bytes()[0],
                        ));
                    }
                    _ => unreachable!(),
                }
            }
        }
        Ok(Some(Line {
            line_type: t,
            fields,
            list,
        }))
    }

    fn read_ascii_list(
        &mut self,
        ft: super::schema::FieldType,
        len: usize,
        toks: &mut std::str::SplitWhitespace,
    ) -> Result<List> {
        match ft {
            FieldType::IntList => {
                let mut v = Vec::with_capacity(len);
                for _ in 0..len {
                    v.push(
                        toks.next()
                            .ok_or_else(|| anyhow!("missing int list elem"))?
                            .parse()?,
                    );
                }
                Ok(List::Ints(v))
            }
            FieldType::RealList => {
                let mut v = Vec::with_capacity(len);
                for _ in 0..len {
                    v.push(
                        toks.next()
                            .ok_or_else(|| anyhow!("missing real list elem"))?
                            .parse()?,
                    );
                }
                Ok(List::Reals(v))
            }
            FieldType::StrList => {
                let mut v = Vec::with_capacity(len);
                for _ in 0..len {
                    let slen: usize = toks
                        .next()
                        .ok_or_else(|| anyhow!("missing string len"))?
                        .parse()?;
                    let s = toks.next().ok_or_else(|| anyhow!("missing string"))?;
                    v.push(s.chars().take(slen).collect::<String>().into_bytes());
                }
                Ok(List::Strings(v))
            }
            _ => {
                // STRING / DNA: read the raw token.
                let s = toks.next().ok_or_else(|| anyhow!("missing string"))?;
                Ok(List::Bytes(s.as_bytes().to_vec()))
            }
        }
    }

    fn read_binary_line(&mut self) -> Result<Option<Line>> {
        if self.file.stream_position()? >= self.foot_off {
            return Ok(None);
        }
        self.read_binary_line_raw()
    }

    /// Read one binary line without the data/footer boundary check. Used both
    /// for data lines (via [`read_binary_line`]) and for footer `&`/`;` lines.
    fn read_binary_line_raw(&mut self) -> Result<Option<Line>> {
        let x = self.read_byte()?;
        if x == b'\n' {
            return Ok(None);
        }
        self.line += 1;
        let t = self.binary_unpack[x as usize];
        if t == 0 {
            bail!("unknown binary line code {x}");
        }
        let li = match self.schema.line_info(t) {
            Some(li) => li.clone(),
            None => bail!("unknown line type {}", t as char),
        };
        let mut fields = Vec::with_capacity(li.fields.len());
        for ft in &li.fields {
            match ft {
                FieldType::Real => {
                    let mut b = [0u8; 8];
                    self.file.read_exact(&mut b)?;
                    fields.push(Field::Real(f64::from_le_bytes(b)));
                }
                FieldType::Char => {
                    fields.push(Field::Char(self.read_byte()?));
                }
                // INT and all list types (including the list length) are
                // written as an LTF integer in the binary stream.
                _ => {
                    fields.push(Field::Int(self.read_ltf()?));
                }
            }
        }
        let list = if li.list_elt_size > 0 {
            let list_len = match fields[li.list_field] {
                Field::Int(len) => len,
                _ => 0,
            };
            if list_len > 0 {
                Some(self.read_binary_list(&li, t, x, list_len as usize)?)
            } else {
                None
            }
        } else {
            None
        };
        Ok(Some(Line {
            line_type: t,
            fields,
            list,
        }))
    }

    fn read_binary_list(&mut self, li: &LineInfo, t: u8, x: u8, list_len: usize) -> Result<List> {
        let compressed = x & 0x1 != 0;
        if li.list_is_int {
            let first = self.read_ltf()?;
            if list_len == 1 {
                return Ok(List::Ints(vec![first]));
            }
            let used = self.read_byte()? as usize;
            let packed = if compressed {
                let nbits = self.read_ltf()?;
                let nbytes = ((nbits + 7) >> 3) as usize;
                let mut buf = vec![0u8; nbytes];
                self.file.read_exact(&mut buf)?;
                let codec = self
                    .codecs
                    .get(t as usize)
                    .and_then(|c| c.as_ref())
                    .ok_or_else(|| anyhow!("no codec for line type {}", t as char))?;
                codec.decode(&buf, nbits)
            } else {
                let mut buf = vec![0u8; (list_len - 1) * used];
                self.file.read_exact(&mut buf)?;
                buf
            };
            let ints = unpack_ints(&packed, used, first);
            Ok(List::Ints(ints))
        } else if li.is_dna {
            // DNA: read raw bytes (2-bit compression not used for `.1aln`).
            let mut buf = vec![0u8; list_len];
            self.file.read_exact(&mut buf)?;
            Ok(List::Bytes(buf))
        } else {
            // STRING or REAL_LIST.
            let bytes = if compressed {
                let nbits = self.read_ltf()?;
                let nbytes = ((nbits + 7) >> 3) as usize;
                let mut buf = vec![0u8; nbytes];
                self.file.read_exact(&mut buf)?;
                let codec = self
                    .codecs
                    .get(t as usize)
                    .and_then(|c| c.as_ref())
                    .ok_or_else(|| anyhow!("no codec for line type {}", t as char))?;
                codec.decode(&buf, nbits)
            } else {
                let mut buf = vec![0u8; list_len * li.list_elt_size];
                self.file.read_exact(&mut buf)?;
                buf
            };
            if li.fields[li.list_field] == FieldType::RealList {
                let mut v = Vec::with_capacity(list_len);
                for chunk in bytes.chunks_exact(8) {
                    v.push(f64::from_le_bytes(chunk.try_into().unwrap()));
                }
                Ok(List::Reals(v))
            } else {
                Ok(List::Bytes(bytes))
            }
        }
    }

    fn read_ltf(&mut self) -> Result<i64> {
        let mut buf = [0u8; 9];
        buf[0] = self.read_byte()?;
        let n = size_from_first(buf[0]);
        for b in buf.iter_mut().take(n).skip(1) {
            *b = self.read_byte()?;
        }
        let (v, _) = int_get(&buf, 0)?;
        Ok(v)
    }
}

/// Infer the total encoded byte count from the first LTF byte.
///
/// For the `0x00..0x07` and `0x80..0x87` ranges the low 3 bits hold the number
/// of *data* bytes (minus one), so the total is `1 + (h + 1) = h + 2`.
fn size_from_first(b0: u8) -> usize {
    match b0 >> 5 {
        2 | 3 | 6 | 7 => 1,
        1 => 2,
        0 | 4 => {
            let h = (b0 & 0x07) as usize;
            if h == 0 {
                1
            } else {
                h + 2
            }
        }
        _ => 1,
    }
}

/// Unpack `used`-byte little-endian signed integers into differences, then
/// accumulate them onto `first` to recover the original integer list.
fn unpack_ints(packed: &[u8], used: usize, first: i64) -> Vec<i64> {
    let n = packed.len() / used;
    let mut out = Vec::with_capacity(n + 1);
    out.push(first);
    let mut acc = first;
    for i in 0..n {
        let mut v: i64 = 0;
        let base = i * used;
        for k in 0..used {
            v |= (packed[base + k] as i64) << (8 * k);
        }
        if used < 8 && packed[base + used - 1] & 0x80 != 0 {
            v |= !0i64 << (8 * used);
        }
        acc += v;
        out.push(acc);
    }
    out
}

/// A cursor over the ASCII body of a header/footer line.
///
/// Mirrors FastGA's ASCII field readers: whitespace-skipping integer tokens
/// (like `readInt`) and fixed-length byte reads (like `readString`) so that
/// length-prefixed strings may contain spaces.
struct AsciiCursor<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> AsciiCursor<'a> {
    fn new(body: &'a str) -> Self {
        AsciiCursor {
            b: body.as_bytes(),
            i: 0,
        }
    }

    fn skip_ws(&mut self) {
        while self.i < self.b.len() && self.b[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
    }

    /// Read the next whitespace-delimited token.
    fn read_str(&mut self) -> Result<&'a str> {
        self.skip_ws();
        let start = self.i;
        while self.i < self.b.len() && !self.b[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
        if start == self.i {
            bail!("missing token");
        }
        std::str::from_utf8(&self.b[start..self.i]).map_err(|_| anyhow!("bad token"))
    }

    /// Read the next whitespace-delimited integer.
    fn read_int(&mut self) -> Result<i64> {
        self.read_str()?.parse().map_err(|_| anyhow!("bad integer"))
    }

    /// Read exactly `n` bytes, skipping leading whitespace.
    fn read_chars(&mut self, n: usize) -> Result<&'a [u8]> {
        self.skip_ws();
        if self.i + n > self.b.len() {
            bail!("string overruns header line");
        }
        let s = &self.b[self.i..self.i + n];
        self.i += n;
        Ok(s)
    }

    /// Read a single non-whitespace byte.
    fn read_char(&mut self) -> Result<u8> {
        self.skip_ws();
        let c = *self
            .b
            .get(self.i)
            .ok_or_else(|| anyhow!("missing character"))?;
        self.i += 1;
        Ok(c)
    }
}

/// Parse a sequence of length-prefixed strings from a header line body.
///
/// The body starts with the count of strings, then for each string a length
/// followed by that many characters. Strings may contain spaces (e.g. the
/// provenance command), so parsing is byte-cursor based.
fn parse_len_strings(body: &str) -> Result<Vec<String>> {
    let b = body.as_bytes();
    let mut i = 0;
    // Skip leading whitespace and read the count.
    while i < b.len() && b[i].is_ascii_whitespace() {
        i += 1;
    }
    let start = i;
    while i < b.len() && !b[i].is_ascii_whitespace() {
        i += 1;
    }
    let count: usize = std::str::from_utf8(&b[start..i])
        .map_err(|_| anyhow!("bad string count"))?
        .parse()?;
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        while i < b.len() && b[i].is_ascii_whitespace() {
            i += 1;
        }
        let len_start = i;
        while i < b.len() && !b[i].is_ascii_whitespace() {
            i += 1;
        }
        let len: usize = std::str::from_utf8(&b[len_start..i])
            .map_err(|_| anyhow!("bad string length"))?
            .parse()?;
        while i < b.len() && b[i].is_ascii_whitespace() {
            i += 1;
        }
        let content_end = i + len;
        if content_end > b.len() {
            bail!("string overruns header line");
        }
        out.push(String::from_utf8_lossy(&b[i..content_end]).to_string());
        i = content_end;
    }
    Ok(out)
}

/// A ONEcode binary file writer (single-threaded).
pub struct Writer {
    file: BufWriter<File>,
    schema: Schema,
    /// Per-line-type list codecs.
    codecs: Vec<Option<VcCodec>>,
    /// Whether the codec is in use (compressed) for a line type.
    use_codec: Vec<bool>,
    /// Training bytes accumulated per line type.
    list_tack: Vec<i64>,
    /// Object index (byte offsets) per line type.
    index: Vec<Vec<i64>>,
    /// Accumulated counts.
    counts: Vec<Counts>,
    /// Current byte offset.
    byte: i64,
    is_header_out: bool,
    is_last_line_binary: bool,
    provenance: Vec<Provenance>,
    references: Vec<Reference>,
    pub is_binary: bool,
}

impl Writer {
    /// Create a new binary ONEcode file.
    pub fn open(path: &str, schema: Schema, is_binary: bool) -> Result<Writer> {
        Ok(Writer {
            file: BufWriter::new(File::create(path)?),
            schema,
            codecs: (0..128).map(|_| None).collect(),
            use_codec: vec![false; 128],
            list_tack: vec![0; 128],
            index: vec![Vec::new(); 128],
            counts: vec![Counts::default(); 128],
            byte: 0,
            is_header_out: false,
            is_last_line_binary: true,
            provenance: Vec::new(),
            references: Vec::new(),
            is_binary,
        })
    }

    /// Add a reference entry to the header.
    pub fn add_reference(&mut self, filename: &str, count: i64) {
        self.references.push(Reference {
            filename: filename.to_string(),
            count,
        });
    }

    /// Add a provenance entry to the header.
    pub fn add_provenance(&mut self, program: &str, version: &str, command: &str, date: &str) {
        self.provenance.push(Provenance {
            program: program.to_string(),
            version: version.to_string(),
            command: command.to_string(),
            date: date.to_string(),
        });
    }

    /// Write a data line.
    pub fn write_line(&mut self, t: u8, fields: &[Field], list: Option<&List>) -> Result<()> {
        let li = self
            .schema
            .line_info(t)
            .ok_or_else(|| anyhow!("unknown line type {}", t as char))?
            .clone();
        if !self.is_header_out {
            self.write_header()?;
            self.is_header_out = true;
        }
        self.counts[t as usize].count += 1;
        if li.is_object {
            self.index[t as usize].push(self.byte);
        }
        if self.is_binary {
            self.write_binary_line(t, &li, fields, list)?;
        } else {
            self.write_ascii_line(t, &li, fields, list)?;
        }
        Ok(())
    }

    /// Write the ASCII prolog.
    fn write_header(&mut self) -> Result<()> {
        let mut s = String::new();
        s.push_str(&format!(
            "1 {} {} {} {}",
            self.schema.file_type.len(),
            self.schema.file_type,
            MAJOR,
            MINOR
        ));
        if let Some(sub) = &self.schema.sub_type {
            s.push_str(&format!("\n2 {} {}", sub.len(), sub));
        }
        for p in &self.provenance {
            s.push_str(&format!(
                "\n! 4 {} {} {} {} {} {} {} {}",
                p.program.len(),
                p.program,
                p.version.len(),
                p.version,
                p.command.len(),
                p.command,
                p.date.len(),
                p.date
            ));
        }
        s.push_str("\n.");
        if !self.references.is_empty() {
            for r in &self.references {
                s.push_str(&format!(
                    "\n< {} {} {}",
                    r.filename.len(),
                    r.filename,
                    r.count
                ));
            }
            s.push_str("\n.");
        }
        // Schema (each write_info_spec starts with a leading newline).
        for &code in &self.schema.defn_order {
            s.push_str(&write_info_spec(&self.schema, code));
        }
        if self.is_binary {
            // No trailing newline: the first binary line completes the `$` line.
            s.push_str(&format!("\n$ {}", 0)); // little-endian isBig=0
        } else {
            s.push_str("\n.\n");
        }
        self.file.write_all(s.as_bytes())?;
        self.file.flush()?;
        self.is_last_line_binary = false;
        self.byte = self.byte_position()?;
        Ok(())
    }

    fn byte_position(&mut self) -> Result<i64> {
        self.file.flush()?;
        Ok(self.file.stream_position()? as i64)
    }

    fn write_binary_line(
        &mut self,
        t: u8,
        li: &LineInfo,
        fields: &[Field],
        list: Option<&List>,
    ) -> Result<()> {
        // Newline between binary lines.
        if !self.is_last_line_binary {
            self.file.write_all(b"\n")?;
            self.byte += 1;
        }
        let mut x = binary_pack(t);
        let use_codec = li.list_elt_size > 0 && self.use_codec[t as usize];
        if use_codec {
            x |= 0x01;
        }
        self.file.write_all(&[x])?;
        self.byte += 1;
        // Fields.
        for &f in fields {
            match f {
                Field::Int(v) => {
                    let mut buf = Vec::new();
                    int_put(&mut buf, v);
                    self.byte += buf.len() as i64;
                    self.file.write_all(&buf)?;
                }
                Field::Real(v) => {
                    self.byte += 8;
                    self.file.write_all(&v.to_le_bytes())?;
                }
                Field::Char(c) => {
                    self.byte += 1;
                    self.file.write_all(&[c])?;
                }
            }
        }
        // List.
        if li.list_elt_size > 0 {
            if let Some(list) = list {
                let list_len = list_len_of(list);
                if list_len > 0 {
                    self.counts[t as usize].total += list_len as i64;
                    if (list_len as i64) > self.counts[t as usize].max {
                        self.counts[t as usize].max = list_len as i64;
                    }
                    if li.list_is_int {
                        let ints = match list {
                            List::Ints(v) => v,
                            _ => bail!("line {} list is not INT_LIST", t as char),
                        };
                        let first = ints[0];
                        let mut buf = Vec::new();
                        int_put(&mut buf, first);
                        self.byte += buf.len() as i64;
                        self.file.write_all(&buf)?;
                        if ints.len() == 1 {
                            self.is_last_line_binary = true;
                            return Ok(());
                        }
                        let (packed, used) = compact_int_list(ints);
                        self.byte += 1;
                        self.file.write_all(&[used as u8])?;
                        if use_codec {
                            let codec = self.codecs[t as usize].as_ref().unwrap();
                            let (cbuf, nbits) = codec.encode(&packed);
                            let nb = ((nbits + 7) >> 3) as usize;
                            let mut lb = Vec::new();
                            int_put(&mut lb, nbits);
                            self.byte += lb.len() as i64;
                            self.file.write_all(&lb)?;
                            self.byte += nb as i64;
                            self.file.write_all(&cbuf[..nb])?;
                        } else {
                            self.byte += packed.len() as i64;
                            self.file.write_all(&packed)?;
                            self.train_codec(t, &packed)?;
                        }
                    } else {
                        // STRING / REAL_LIST / DNA.
                        let bytes = match list {
                            List::Bytes(b) => b.clone(),
                            List::Reals(r) => {
                                let mut b = Vec::new();
                                for v in r {
                                    b.extend_from_slice(&v.to_le_bytes());
                                }
                                b
                            }
                            _ => bail!("unsupported list for line {}", t as char),
                        };
                        if use_codec {
                            let codec = self.codecs[t as usize].as_ref().unwrap();
                            let (cbuf, nbits) = codec.encode(&bytes);
                            let nb = ((nbits + 7) >> 3) as usize;
                            let mut lb = Vec::new();
                            int_put(&mut lb, nbits);
                            self.byte += lb.len() as i64;
                            self.file.write_all(&lb)?;
                            self.byte += nb as i64;
                            self.file.write_all(&cbuf[..nb])?;
                        } else {
                            self.byte += bytes.len() as i64;
                            self.file.write_all(&bytes)?;
                            self.train_codec(t, &bytes)?;
                        }
                    }
                }
            }
        }
        self.is_last_line_binary = true;
        Ok(())
    }

    fn write_ascii_line(
        &mut self,
        t: u8,
        li: &LineInfo,
        fields: &[Field],
        list: Option<&List>,
    ) -> Result<()> {
        if !self.is_last_line_binary {
            self.file.write_all(b"\n")?;
        }
        self.file.write_all(&[t])?;
        let mut s = String::new();
        for (i, &f) in fields.iter().enumerate() {
            if i == li.list_field {
                let len = list.map(list_len_of).unwrap_or(0);
                s.push_str(&format!(" {}", len));
                if let Some(list) = list {
                    s.push_str(&ascii_list_str(list, len));
                }
            } else {
                match f {
                    Field::Int(v) => s.push_str(&format!(" {v}")),
                    Field::Real(v) => s.push_str(&format!(" {v}")),
                    Field::Char(c) => s.push_str(&format!(" {}", c as char)),
                }
            }
        }
        self.file.write_all(s.as_bytes())?;
        self.is_last_line_binary = false;
        Ok(())
    }

    /// Accumulate histogram bytes to train a list codec.
    fn train_codec(&mut self, t: u8, bytes: &[u8]) -> Result<()> {
        self.list_tack[t as usize] += bytes.len() as i64;
        if self.codecs[t as usize].is_none() {
            self.codecs[t as usize] = Some(VcCodec::new());
        }
        if let Some(codec) = &mut self.codecs[t as usize] {
            codec.add_to_table(bytes);
        }
        if self.list_tack[t as usize] > 100_000 {
            if let Some(codec) = &mut self.codecs[t as usize] {
                codec.create_codec(true)?;
            }
            self.use_codec[t as usize] = true;
            self.list_tack[t as usize] = 0;
        }
        Ok(())
    }

    /// Close the file: write footer and (for binary) offsets.
    pub fn close(mut self) -> Result<()> {
        if self.is_binary {
            let foot_off = self.byte_position()?;
            // Footer.
            let mut s = String::new();
            for &code in &self.schema.defn_order {
                if code & 0x80 != 0 {
                    continue;
                }
                let t = code as u8;
                let c = self.counts[t as usize];
                if c.count > 0 {
                    s.push_str(&format!("# {} {}\n", t as char, c.count));
                    if c.max > 0 {
                        s.push_str(&format!("@ {} {}\n", t as char, c.max));
                    }
                    if c.total > 0 {
                        s.push_str(&format!("+ {} {}\n", t as char, c.total));
                    }
                    if self
                        .schema
                        .line_info(t)
                        .map(|li| li.is_object)
                        .unwrap_or(false)
                        && !self.index[t as usize].is_empty()
                    {
                        let idx = &self.index[t as usize];
                        let mut line = format!("& {} {}", t as char, idx.len() + 1);
                        line.push_str(" 0");
                        for v in idx {
                            line.push_str(&format!(" {}", v));
                        }
                        s.push_str(&line);
                        s.push('\n');
                    }
                    if self.use_codec[t as usize] {
                        if let Some(codec) = &self.codecs[t as usize] {
                            let ser = codec.serialize();
                            s.push_str(&format!("; {} {}", t as char, ser.len()));
                            self.file.write_all(s.as_bytes())?;
                            s.clear();
                            self.file.write_all(&ser)?;
                            self.file.write_all(b"\n")?;
                        }
                    }
                }
            }
            self.file.write_all(s.as_bytes())?;
            self.file.write_all(b"^\n")?;
            self.file.write_all(&foot_off.to_le_bytes())?;
        }
        self.file.flush()?;
        Ok(())
    }
}

fn list_len_of(list: &List) -> usize {
    match list {
        List::Ints(v) => v.len(),
        List::Reals(v) => v.len(),
        List::Bytes(v) => v.len(),
        List::Strings(v) => v.len(),
    }
}

fn ascii_list_str(list: &List, len: usize) -> String {
    match list {
        List::Ints(v) => v.iter().take(len).map(|x| format!(" {x}")).collect(),
        List::Reals(v) => v.iter().take(len).map(|x| format!(" {x}")).collect(),
        List::Bytes(b) => format!(" {}", String::from_utf8_lossy(b)),
        List::Strings(v) => v
            .iter()
            .take(len)
            .map(|s| format!(" {} {}", s.len(), String::from_utf8_lossy(s)))
            .collect(),
    }
}

/// Convert the first bytes of an integer list to differences and pack each
/// into the minimal number of little-endian bytes.
fn compact_int_list(list: &[i64]) -> (Vec<u8>, usize) {
    let len = list.len();
    if len <= 1 {
        return (Vec::new(), 0);
    }
    let mut diffs = vec![0i64; len - 1];
    let mut mask: i64 = 0;
    for i in 1..len {
        let d = list[i] - list[i - 1];
        diffs[i - 1] = d;
        mask |= if d >= 0 { d } else { -(d + 1) };
    }
    let mut used = 8usize;
    let mut m = mask >> 7;
    for d in 1..8usize {
        if m == 0 {
            used = d;
            break;
        }
        m >>= 8;
    }
    let mut packed = Vec::with_capacity((len - 1) * used);
    for d in &diffs {
        let le = d.to_le_bytes();
        packed.extend_from_slice(&le[..used]);
    }
    (packed, used)
}

/// Render one schema definition line as ASCII text.
fn write_info_spec(schema: &Schema, code: i32) -> String {
    let mut s = String::from("\n~ ");
    if code & 0x80 != 0 {
        s.push_str(&format!("G {} 0", (code & 0x7f) as u8 as char));
    } else {
        let t = code as u8;
        if let Some(li) = schema.line_info(t) {
            s.push_str(if li.is_object { "O" } else { "D" });
            s.push_str(&format!(" {} {}", t as char, li.fields.len()));
            for ft in &li.fields {
                s.push_str(&format!(" {} {}", ft.name().len(), ft.name()));
            }
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_path(tag: &str) -> String {
        let mut p = std::env::temp_dir();
        p.push(format!("onecode_{tag}_{}.bin", std::process::id()));
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn compact_unpack_round_trip() {
        let list = [0i64, 100, 200, 150, 1000, 999, -5];
        let (packed, used) = compact_int_list(&list);
        let out = unpack_ints(&packed, used, list[0]);
        assert_eq!(out, list);
    }

    #[test]
    fn int_list_round_trip() {
        let path = temp_path("intlist");
        let schema = parse_schema_text(aln_schema_text()).unwrap();
        {
            let mut w = Writer::open(&path, schema.clone(), true).unwrap();
            w.add_reference("chr1.fa", 1);
            let list = vec![0i64, 100, 200, 150, 1000, 999, -5];
            w.write_line(b'T', &[Field::Int(7)], Some(&List::Ints(list)))
                .unwrap();
            w.close().unwrap();
        }
        {
            let mut r = Reader::open(&path).unwrap();
            let line = r.read_line().unwrap().unwrap();
            assert_eq!(line.line_type, b'T');
            assert_eq!(line.int(0), 7);
            match line.list.unwrap() {
                List::Ints(v) => assert_eq!(v, vec![0, 100, 200, 150, 1000, 999, -5]),
                _ => panic!("expected int list"),
            }
        }
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn multi_line_round_trip() {
        let path = temp_path("multiline");
        let schema = parse_schema_text(aln_schema_text()).unwrap();
        {
            let mut w = Writer::open(&path, schema.clone(), true).unwrap();
            w.add_reference("a.fa", 1);
            w.add_reference("b.fa", 2);
            w.write_line(
                b'A',
                &[
                    Field::Int(0),
                    Field::Int(10),
                    Field::Int(20),
                    Field::Int(1),
                    Field::Int(30),
                    Field::Int(50),
                ],
                None,
            )
            .unwrap();
            w.write_line(b'D', &[Field::Int(3)], None).unwrap();
            w.write_line(b'T', &[Field::Int(2)], Some(&List::Ints(vec![10, 20])))
                .unwrap();
            w.write_line(b'X', &[Field::Int(2)], Some(&List::Ints(vec![1, 2])))
                .unwrap();
            w.close().unwrap();
        }
        {
            let mut r = Reader::open(&path).unwrap();
            assert_eq!(r.references.len(), 2);
            let mut lines = Vec::new();
            while let Some(l) = r.read_line().unwrap() {
                lines.push((l.line_type, l.fields, l.list));
            }
            assert_eq!(lines.len(), 4);
            assert_eq!(lines[0].0, b'A');
            assert_eq!(
                lines[0].1,
                vec![
                    Field::Int(0),
                    Field::Int(10),
                    Field::Int(20),
                    Field::Int(1),
                    Field::Int(30),
                    Field::Int(50)
                ]
            );
            assert_eq!(lines[1].0, b'D');
            assert_eq!(lines[2].0, b'T');
            match &lines[2].2 {
                Some(List::Ints(v)) => assert_eq!(v, &vec![10, 20]),
                _ => panic!(),
            }
            assert_eq!(lines[3].0, b'X');
        }
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn read_golden_file() {
        let path = env!("CARGO_MANIFEST_DIR").to_string() + "/tests/genome/mg1655-sakai.1aln";
        let mut r = Reader::open(&path).unwrap();
        eprintln!("is_binary={} provenance={:?}", r.is_binary, r.provenance);
        eprintln!("references={:?}", r.references);
        eprintln!(
            "counts t={:?} g={:?} S={:?} a={:?} A={:?} D={:?} T={:?} X={:?}",
            r.counts[b't' as usize],
            r.counts[b'g' as usize],
            r.counts[b'S' as usize],
            r.counts[b'a' as usize],
            r.counts[b'A' as usize],
            r.counts[b'D' as usize],
            r.counts[b'T' as usize],
            r.counts[b'X' as usize]
        );
        let mut n = 0;
        while let Some(l) = r.read_line().unwrap() {
            n += 1;
            if n <= 3 {
                eprintln!(
                    "line {} type={} fields={:?}",
                    n, l.line_type as char, l.fields
                );
                if let Some(list) = &l.list {
                    match list {
                        List::Ints(v) => eprintln!(
                            "  ints len={} first10={:?}",
                            v.len(),
                            v.iter().take(10).collect::<Vec<_>>()
                        ),
                        List::Bytes(b) => {
                            eprintln!("  bytes len={} {:?}", b.len(), String::from_utf8_lossy(b))
                        }
                        _ => {}
                    }
                }
            }
        }
        eprintln!("TOTAL LINES={n}");
    }
}
