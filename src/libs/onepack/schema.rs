//! ONEcode schema model and text parser.
//!
//! A ONEcode file embeds its schema in the header as ASCII `~` lines. The
//! schema declares the primary file type plus `O` (object), `D` (data), and
//! `G` (group) record definitions. Each record has an ordered list of field
//! types; at most one of them may be a list type. This module provides the
//! in-memory model ([`Schema`]/[`LineInfo`]) used by both the reader and the
//! writer in [`super::container`].

use anyhow::{anyhow, bail, Result};

/// Container version numbers (must match FastGA's `ONElib.c` `MAJOR`/`MINOR`).
pub const MAJOR: i64 = 2;
pub const MINOR: i64 = 1;

/// A ONEcode field type.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FieldType {
    Int,
    Real,
    Char,
    Str,
    IntList,
    RealList,
    StrList,
    DNA,
}

impl FieldType {
    /// Parse a field type from its schema-text name (e.g. `"INT"`).
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "INT" => Some(Self::Int),
            "REAL" => Some(Self::Real),
            "CHAR" => Some(Self::Char),
            "STRING" => Some(Self::Str),
            "INT_LIST" => Some(Self::IntList),
            "REAL_LIST" => Some(Self::RealList),
            "STRING_LIST" => Some(Self::StrList),
            "DNA" => Some(Self::DNA),
            _ => None,
        }
    }

    /// The canonical type name used in schema text.
    pub fn name(self) -> &'static str {
        match self {
            Self::Int => "INT",
            Self::Real => "REAL",
            Self::Char => "CHAR",
            Self::Str => "STRING",
            Self::IntList => "INT_LIST",
            Self::RealList => "REAL_LIST",
            Self::StrList => "STRING_LIST",
            Self::DNA => "DNA",
        }
    }

    /// Size in bytes of a single element of this field (0 for scalars).
    pub fn elt_size(self) -> usize {
        match self {
            Self::Int | Self::Real | Self::Char => 0,
            Self::Str | Self::StrList | Self::DNA => 1,
            Self::IntList | Self::RealList => 8,
        }
    }

    /// Whether this is a list-bearing field.
    pub fn is_list(self) -> bool {
        matches!(
            self,
            Self::Str | Self::IntList | Self::RealList | Self::StrList | Self::DNA
        )
    }

    /// Whether a list field is a single (possibly length-prefixed) byte blob.
    pub fn is_single(self) -> bool {
        matches!(self, Self::Str | Self::DNA)
    }
}

/// Per-line-type schema information.
#[derive(Clone)]
pub struct LineInfo {
    /// Whether this line type starts an object.
    pub is_object: bool,
    /// Ordered field types.
    pub fields: Vec<FieldType>,
    /// Index into `fields` of the single list field.
    pub list_field: usize,
    /// Byte size of a list element (0 if no list).
    pub list_elt_size: usize,
    /// Whether the list field is an `INT_LIST`.
    pub list_is_int: bool,
    /// Whether the list field is `DNA`.
    pub is_dna: bool,
    /// Object types grouped/contained by this object (for objects).
    pub contains: Vec<u8>,
}

impl LineInfo {
    /// Build a `LineInfo` and derive its list-field metadata.
    pub fn new(is_object: bool, fields: Vec<FieldType>) -> Self {
        let mut info = LineInfo {
            is_object,
            list_field: 0,
            list_elt_size: 0,
            list_is_int: false,
            is_dna: false,
            fields,
            contains: Vec::new(),
        };
        for (i, f) in info.fields.iter().copied().enumerate() {
            if f.is_list() {
                info.list_field = i;
                info.list_elt_size = f.elt_size();
                info.list_is_int = f == FieldType::IntList;
                info.is_dna = f == FieldType::DNA;
            }
        }
        info
    }
}

/// A parsed schema for one primary file type.
#[derive(Clone)]
pub struct Schema {
    /// Primary file type name.
    pub file_type: String,
    /// Optional secondary type.
    pub sub_type: Option<String>,
    /// Per-line-type info, indexed by the ASCII line character.
    pub info: Vec<Option<LineInfo>>,
    /// Definition order; the high bit marks `G` group lines.
    pub defn_order: Vec<i32>,
    /// Maximum number of fields across all line types.
    pub n_field_max: usize,
}

impl Schema {
    /// Look up line info by its ASCII character.
    pub fn line_info(&self, t: u8) -> Option<&LineInfo> {
        self.info.get(t as usize).and_then(|x| x.as_ref())
    }
}

/// Compute the binary pack code for a line type character.
///
/// Two consecutive codes (low bit 0/1) map back to the same line type; the low
/// bit marks whether the list is compressed.
pub fn binary_pack(t: u8) -> u8 {
    if t.is_ascii_uppercase() {
        ((t - b'A') << 1) | 0x80
    } else if t.is_ascii_lowercase() {
        ((26 + (t - b'a')) << 1) | 0x80
    } else {
        match t {
            b';' => (52 << 1) | 0x80,
            b'&' => (53 << 1) | 0x80,
            b'/' => (54 << 1) | 0x80,
            b'.' => (55 << 1) | 0x80,
            _ => 0,
        }
    }
}

/// Build a fresh, empty schema for `file_type`.
pub fn empty_schema(file_type: &str) -> Schema {
    Schema {
        file_type: file_type.to_string(),
        sub_type: None,
        info: vec![None; 128],
        defn_order: Vec::new(),
        n_field_max: 4,
    }
}

/// Add an `O`/`D` definition line to `schema` (groups are handled separately).
///
/// `k` is `b'O'` or `b'D'`, `t` the line character, and `fields` the field
/// types. `current_object` is the currently open object, if any, which owns
/// newly defined `D` line types.
pub fn add_definition(
    schema: &mut Schema,
    k: u8,
    t: u8,
    fields: Vec<FieldType>,
    current_object: Option<u8>,
) -> Result<()> {
    if schema.info[t as usize].is_some() {
        bail!("duplicate schema specification for line type {}", t as char);
    }
    let info = LineInfo::new(k == b'O', fields);
    if info.fields.len() > schema.n_field_max {
        schema.n_field_max = info.fields.len();
    }
    if let Some(obj) = current_object {
        if let Some(li) = &mut schema.info[obj as usize] {
            li.contains.push(t);
        }
    }
    schema.defn_order.push(t as i32);
    schema.info[t as usize] = Some(info);
    Ok(())
}

/// The universal header/footer line-type schema (FastGA `ONElib.c`), shared
/// by every ONEcode file regardless of primary type. Loaded by the reader
/// before reading any `~` schema lines from the file.
pub fn header_schema_text() -> &'static str {
    "P 3 def\n\
     D 1 3 6 STRING 3 INT 3 INT\n\
     D 2 1 6 STRING\n\
     D # 2 4 CHAR 3 INT\n\
     D @ 2 4 CHAR 3 INT\n\
     D + 2 4 CHAR 3 INT\n\
     D % 4 4 CHAR 4 CHAR 4 CHAR 3 INT\n\
     D ! 1 11 STRING_LIST\n\
     D < 2 6 STRING 3 INT\n\
     D > 1 6 STRING\n\
     D ~ 3 4 CHAR 4 CHAR 11 STRING_LIST\n\
     D . 0\n\
     D $ 1 3 INT\n\
     D ^ 0\n\
     D - 1 3 INT\n\
     D & 2 4 CHAR 8 INT_LIST\n\
     D ; 2 4 CHAR 6 STRING\n\
     D / 1 6 STRING\n"
}

/// Parse a schema from its textual form (as embedded in a ONEcode header).
///
/// Lines of the form:
/// - `P <len> <primary>`
/// - `S <len> <secondary>`
/// - `O <char> <n> [<len> <TYPE> ...]`
/// - `D <char> <n> [<len> <TYPE> ...]`
/// - `G <char>`
pub fn parse_schema_text(text: &str) -> Result<Schema> {
    let mut schema: Option<Schema> = None;
    let mut current_object: Option<u8> = None;
    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('.') {
            continue;
        }
        let mut it = line.split_whitespace();
        let first = it.next().unwrap_or("");
        match first {
            "P" => {
                let _len = it.next();
                let name = it.next().unwrap_or("").to_string();
                schema = Some(empty_schema(&name));
            }
            "S" => {
                let _len = it.next();
                if let Some(s) = &mut schema {
                    s.sub_type = Some(it.next().unwrap_or("").to_string());
                }
            }
            "O" | "D" => {
                let s = schema
                    .as_mut()
                    .ok_or_else(|| anyhow!("schema line before P"))?;
                let t = it
                    .next()
                    .ok_or_else(|| anyhow!("missing line type"))?
                    .as_bytes()[0];
                let _n: usize = it
                    .next()
                    .ok_or_else(|| anyhow!("missing field count"))?
                    .parse()?;
                let mut fields = Vec::new();
                while let Some(slen) = it.next() {
                    let _: usize = slen.parse().map_err(|_| anyhow!("bad field length"))?;
                    let name = it.next().ok_or_else(|| anyhow!("missing field type"))?;
                    let ft = FieldType::from_name(name)
                        .ok_or_else(|| anyhow!("unknown field type {name}"))?;
                    fields.push(ft);
                }
                add_definition(s, first.as_bytes()[0], t, fields, current_object)?;
                if let Some(li) = &s.info[t as usize] {
                    if li.is_object {
                        current_object = Some(t);
                    }
                }
            }
            "G" => {
                let s = schema
                    .as_mut()
                    .ok_or_else(|| anyhow!("schema line before P"))?;
                let t = it
                    .next()
                    .ok_or_else(|| anyhow!("missing group type"))?
                    .as_bytes()[0];
                if let Some(obj) = current_object {
                    if let Some(li) = &mut s.info[obj as usize] {
                        li.contains.push(t);
                    }
                }
                s.defn_order.push((t as i32) | 0x80);
            }
            _ => {}
        }
    }
    schema.ok_or_else(|| anyhow!("no P line in schema text"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_aln_schema() {
        let s = parse_schema_text(crate::libs::onepack::container::aln_schema_text()).unwrap();
        assert_eq!(s.file_type, "aln");
        let a = s.line_info(b'A').unwrap();
        assert!(a.is_object);
        assert_eq!(a.fields.len(), 6);
        assert!(a.fields.iter().all(|f| *f == FieldType::Int));
        let t = s.line_info(b'T').unwrap();
        assert!(t.list_is_int);
        assert_eq!(t.list_elt_size, 8);
        let g = s.line_info(b'g').unwrap();
        assert!(g.is_object);
        assert_eq!(g.fields.len(), 0);
    }

    #[test]
    fn binary_pack_letters() {
        assert_eq!(binary_pack(b'A'), 0x80);
        assert_eq!(binary_pack(b'B'), 0x82);
        assert_eq!(binary_pack(b'a'), (26 << 1) | 0x80);
    }
}
