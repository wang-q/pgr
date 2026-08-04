//! Genomic range (`chr:start-end` with optional `name.` prefix and `(strand)`)
//! parsing and manipulation.
//!
//! # Parsing semantics and the original regex
//!
//! [`Range::from_str`] now uses a hand-written byte scanner (see `decode`).
//! The regex it replaced is kept here verbatim for reference — anyone
//! comparing with the external `intspan` crate or debugging parser
//! differences should consult it (the decoder also survives in the test
//! module as an equivalence oracle):
//!
//! ```text
//! r"(?xi)
//!     (?:(?P<name>[\w_]+)\.)?
//!     (?P<chr>[\w/-]+)
//!     (?:\((?P<strand>.+)\))?
//!     [:]                    # spacer
//!     (?P<start>\d+)
//!     [_\-]?                 # spacer
//!     (?P<end>\d+)?
//!     "
//! ```
//!
//! The scanner replicates these semantics exactly:
//! * unanchored **leftmost** match (e.g. `foo I:1-100` matches at `I:1-100`,
//!   `a.b.c:1-2` matches `b.c:1-2` with name `b`);
//! * greedy groups: `name` = longest word run before `.`, `chr` = longest
//!   `[\w/-]+` run, `strand` = longest content between `(` and the last `)`
//!   (at least one char);
//! * a bare `start` (no `-end`) gives `end = start`; an explicit `end = 0` is
//!   treated as "missing" and also defaults to `start` (regex-era behavior);
//! * no match falls back to `chr` = first whitespace token;
//! * digit runs too large for `i32` fail to match (invalid) instead of the
//!   regex path's `parse::<i32>().unwrap()` panic. An oversized `end` must
//!   not be confused with a missing one (which defaults to `start`), or a
//!   malformed line such as `chr1:5-99999999999` would be silently treated
//!   as the point range `chr1:5`. Unicode digits are likewise treated as
//!   invalid (the regex path's `\d+` would match them and then panic in
//!   `parse::<i32>()`).

use super::IntSpan;
use std::fmt;

#[derive(Debug, Default, Clone)]
pub struct Range {
    pub name: String,
    pub chr: String,
    pub strand: String,
    pub start: i32,
    pub end: i32,
}

impl Range {
    // Immutable accessors
    pub fn name(&self) -> &String {
        &self.name
    }
    pub fn chr(&self) -> &String {
        &self.chr
    }
    pub fn strand(&self) -> &String {
        &self.strand
    }
    pub fn start(&self) -> &i32 {
        &self.start
    }
    pub fn end(&self) -> &i32 {
        &self.end
    }

    // Mutable accessors
    pub fn name_mut(&mut self) -> &mut String {
        &mut self.name
    }
    pub fn strand_mut(&mut self) -> &mut String {
        &mut self.strand
    }

    pub fn new() -> Self {
        Self {
            name: "".to_string(),
            chr: "".to_string(),
            strand: "".to_string(),
            start: 0,
            end: 0,
        }
    }

    /// Constructed from chr, start and end
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from("I", 1, 100);
    /// # assert_eq!(*range.chr(), "I");
    /// # assert_eq!(*range.start(), 1);
    /// # assert_eq!(*range.end(), 100);
    /// ```
    pub fn from(chr: &str, start: i32, end: i32) -> Self {
        Self {
            name: "".to_string(),
            chr: chr.to_string(),
            strand: "".to_string(),
            start,
            end,
        }
    }

    /// Constructed from chr, start and end
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from_full("S288c", "I", "-", 1, 100);
    /// # assert_eq!(*range.name(), "S288c");
    /// # assert_eq!(*range.chr(), "I");
    /// # assert_eq!(*range.strand(), "-");
    /// # assert_eq!(*range.start(), 1);
    /// # assert_eq!(*range.end(), 100);
    /// ```
    pub fn from_full(name: &str, chr: &str, strand: &str, start: i32, end: i32) -> Self {
        Self {
            name: name.to_string(),
            chr: chr.to_string(),
            strand: strand.to_string(),
            start,
            end,
        }
    }

    /// Constructed from string
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from_str("I:1-100");
    /// # assert_eq!(*range.chr(), "I");
    /// # assert_eq!(*range.start(), 1);
    /// # assert_eq!(*range.end(), 100);
    /// # assert_eq!(range.to_string(), "I:1-100");
    /// let range = Range::from_str("I:100");
    /// # assert_eq!(*range.chr(), "I");
    /// # assert_eq!(*range.start(), 100);
    /// # assert_eq!(*range.end(), 100);
    /// # assert_eq!(range.to_string(), "I:100");
    /// let range = Range::from_str("S288c.I(-):27070-29557");
    /// # assert_eq!(*range.name(), "S288c");
    /// # assert_eq!(*range.strand(), "-");
    /// # assert_eq!(range.to_string(), "S288c.I(-):27070-29557");
    /// ```
    // Inherent `from_str` kept API-identical to the external intspan crate;
    // call sites use `Range::from_str`, not `str::parse`.
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(range: &str) -> Self {
        let mut new = Self::new();
        new.decode(range);

        new
    }

    /// Valid or not
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from("I", 1, 100);
    /// assert!(range.is_valid());
    /// let range = Range::from_str("I:100");
    /// assert!(range.is_valid());
    /// let range = Range::from_str("invalid");
    /// assert!(!range.is_valid());
    /// ```
    pub fn is_valid(&self) -> bool {
        self.start != 0
    }

    /// IntSpan
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from("I", 1, 100);
    /// assert_eq!(range.intspan().to_string(), "1-100");
    /// let range = Range::from_str("I:100");
    /// assert_eq!(range.intspan().to_string(), "100");
    /// ```
    pub fn intspan(&self) -> IntSpan {
        IntSpan::from_pair(self.start, self.end)
    }

    /// Trim both ends
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from_str("I:100-200");
    /// assert_eq!(range.trim(30).to_string(), "I:130-170");
    /// assert_eq!(range.trim(70).is_valid(), false);
    /// assert_eq!(range.trim(-30).to_string(), "I:70-230");
    /// ```
    pub fn trim(&self, n: i32) -> Self {
        // Saturating arithmetic keeps extreme `n` (e.g. near i32::MAX) from
        // overflowing; `check` turns the resulting out-of-order pair into
        // the invalid (0, 0) range.
        let mut start = self.start.saturating_add(n);
        let mut end = self.end.saturating_sub(n);
        Self::check(&mut start, &mut end);

        Self {
            name: self.name.to_string(),
            chr: self.chr.to_string(),
            strand: self.strand.to_string(),
            start,
            end,
        }
    }

    /// Trim 5p end
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from_str("I(+):100-200");
    /// assert_eq!(range.trim_5p(30).to_string(), "I(+):130-200");
    /// let range = Range::from_str("I(-):100-200");
    /// assert_eq!(range.trim_5p(30).to_string(), "I(-):100-170");
    /// assert_eq!(range.trim_5p(-30).to_string(), "I(-):100-230");
    /// assert_eq!(range.trim_5p(120).is_valid(), false);
    /// ```
    pub fn trim_5p(&self, n: i32) -> Self {
        let mut start = if self.strand == "-" {
            self.start
        } else {
            self.start.saturating_add(n)
        };
        let mut end = if self.strand == "-" {
            self.end.saturating_sub(n)
        } else {
            self.end
        };
        Self::check(&mut start, &mut end);

        Self {
            name: self.name.to_string(),
            chr: self.chr.to_string(),
            strand: self.strand.to_string(),
            start,
            end,
        }
    }

    /// Trim 3p end
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from_str("I(+):100-200");
    /// assert_eq!(range.trim_3p(30).to_string(), "I(+):100-170");
    /// let range = Range::from_str("I(-):100-200");
    /// assert_eq!(range.trim_3p(30).to_string(), "I(-):130-200");
    /// assert_eq!(range.trim_3p(120).is_valid(), false);
    /// ```
    pub fn trim_3p(&self, n: i32) -> Self {
        let mut start = if self.strand == "-" {
            self.start.saturating_add(n)
        } else {
            self.start
        };
        let mut end = if self.strand == "-" {
            self.end
        } else {
            self.end.saturating_sub(n)
        };
        Self::check(&mut start, &mut end);

        Self {
            name: self.name.to_string(),
            chr: self.chr.to_string(),
            strand: self.strand.to_string(),
            start,
            end,
        }
    }

    /// Shift to 5p end
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from_str("I(+):100-200");
    /// assert_eq!(range.shift_5p(30).to_string(), "I(+):70-170");
    /// assert_eq!(range.shift_5p(-30).to_string(), "I(+):130-230");
    /// let range = Range::from_str("I(-):100-200");
    /// assert_eq!(range.shift_5p(30).to_string(), "I(-):130-230");
    /// ```
    pub fn shift_5p(&self, n: i32) -> Self {
        let mut start = if self.strand == "-" {
            self.start.saturating_add(n)
        } else {
            self.start.saturating_sub(n)
        };
        let mut end = if self.strand == "-" {
            self.end.saturating_add(n)
        } else {
            self.end.saturating_sub(n)
        };
        Self::check(&mut start, &mut end);

        Self {
            name: self.name.to_string(),
            chr: self.chr.to_string(),
            strand: self.strand.to_string(),
            start,
            end,
        }
    }

    /// Shift to 3p end
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from_str("I(+):100-200");
    /// assert_eq!(range.shift_3p(30).to_string(), "I(+):130-230");
    /// assert_eq!(range.shift_3p(-30).to_string(), "I(+):70-170");
    /// let range = Range::from_str("I(-):100-200");
    /// assert_eq!(range.shift_3p(30).to_string(), "I(-):70-170");
    /// ```
    pub fn shift_3p(&self, n: i32) -> Self {
        let mut start = if self.strand == "-" {
            self.start.saturating_sub(n)
        } else {
            self.start.saturating_add(n)
        };
        let mut end = if self.strand == "-" {
            self.end.saturating_sub(n)
        } else {
            self.end.saturating_add(n)
        };
        Self::check(&mut start, &mut end);

        Self {
            name: self.name.to_string(),
            chr: self.chr.to_string(),
            strand: self.strand.to_string(),
            start,
            end,
        }
    }

    /// Flanking region of the 5p end.
    /// A negative value for 'n' indicates positions within the range.
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from_str("I(+):100-200");
    /// assert_eq!(range.flank_5p(30).to_string(), "I(+):70-99");
    /// assert_eq!(range.flank_5p(-30).to_string(), "I(+):100-129");
    /// assert_eq!(range.flank_5p(0).is_valid(), false);
    /// let range = Range::from_str("I(-):100-200");
    /// assert_eq!(range.flank_5p(30).to_string(), "I(-):201-230");
    /// assert_eq!(range.flank_5p(-30).to_string(), "I(-):171-200");
    /// assert_eq!(range.flank_5p(0).is_valid(), false);
    /// ```
    pub fn flank_5p(&self, n: i32) -> Self {
        let mut start = if n > 0 {
            if self.strand == "-" {
                self.end.saturating_add(1)
            } else {
                self.start.saturating_sub(n)
            }
        } else if self.strand == "-" {
            self.end.saturating_add(n).saturating_add(1)
        } else {
            self.start
        };
        let mut end = if n > 0 {
            if self.strand == "-" {
                self.end.saturating_add(n)
            } else {
                self.start.saturating_sub(1)
            }
        } else if self.strand == "-" {
            self.end
        } else {
            self.start.saturating_sub(n).saturating_sub(1)
        };
        Self::check(&mut start, &mut end);

        Self {
            name: self.name.to_string(),
            chr: self.chr.to_string(),
            strand: self.strand.to_string(),
            start,
            end,
        }
    }

    /// Flanking region of the 3p end
    ///
    /// ```
    /// # use pgr::libs::ds::Range;
    /// let range = Range::from_str("I(+):100-200");
    /// assert_eq!(range.flank_3p(30).to_string(), "I(+):201-230");
    /// assert_eq!(range.flank_3p(-30).to_string(), "I(+):171-200");
    /// assert_eq!(range.flank_3p(0).is_valid(), false);
    /// let range = Range::from_str("I(-):100-200");
    /// assert_eq!(range.flank_3p(30).to_string(), "I(-):70-99");
    /// assert_eq!(range.flank_3p(-30).to_string(), "I(-):100-129");
    /// assert_eq!(range.flank_3p(0).is_valid(), false);
    /// ```
    pub fn flank_3p(&self, n: i32) -> Self {
        let mut start = if n > 0 {
            if self.strand == "-" {
                self.start.saturating_sub(n)
            } else {
                self.end.saturating_add(1)
            }
        } else if self.strand == "-" {
            self.start
        } else {
            self.end.saturating_add(n).saturating_add(1)
        };
        let mut end = if n > 0 {
            if self.strand == "-" {
                self.start.saturating_sub(1)
            } else {
                self.end.saturating_add(n)
            }
        } else if self.strand == "-" {
            self.start.saturating_sub(n).saturating_sub(1)
        } else {
            self.end
        };
        Self::check(&mut start, &mut end);

        Self {
            name: self.name.to_string(),
            chr: self.chr.to_string(),
            strand: self.strand.to_string(),
            start,
            end,
        }
    }

    /// Parse `header` with the hand-written scanner (see the module docs for
    /// the regex it replaces and the replicated semantics).
    fn decode(&mut self, header: &str) {
        let s = header.as_bytes();
        let n = s.len();
        // Unanchored leftmost match, mirroring `regex::Regex::captures`.
        for p in 0..n {
            if word_or_slash_char_len(s, p) == 0 {
                continue;
            }
            if let Some((name, chr, strand, start, end)) = match_at(s, p, n) {
                self.name = name;
                self.chr = chr;
                self.strand = strand;
                self.start = start;
                self.end = end;
                // Mirror `decode`: an `end` of 0 is treated as "missing"
                // (e.g. `c:911_0` parses end=0, then defaults to start).
                if self.start != 0 && self.end == 0 {
                    self.end = self.start;
                }
                return;
            }
        }
        // Regex-miss fallback: `chr` is the first whitespace token.
        self.chr = header.split(' ').next().unwrap().to_string();
    }

    fn encode(&self) -> String {
        let mut header = String::new();

        if !self.name.is_empty() {
            header += self.name.as_str();
            if !self.chr.is_empty() {
                header += ".";
                header += self.chr.as_str();
            }
        } else if !self.chr.is_empty() {
            header += self.chr.as_str();
        }

        if !self.strand.is_empty() {
            header += "(";
            header += self.strand.as_str();
            header += ")";
        }

        if self.start != 0 {
            header += ":";
            header += self.start.to_string().as_str();
            if self.end != self.start {
                header += "-";
                header += self.end.to_string().as_str();
            }
        }

        header
    }

    fn check(start: &mut i32, end: &mut i32) {
        if *start < 0 {
            *start = 0;
        }
        if *end < 0 {
            *end = 0;
        }
        if *start > *end {
            *start = 0;
            *end = 0;
        }
    }
}

/// Byte length of the UTF-8 character whose first byte is at `i` (the input
/// is valid UTF-8; continuation bytes return 1 so iteration always advances).
fn char_len_at(s: &[u8], i: usize) -> usize {
    let b = s[i];
    if b < 0x80 {
        1
    } else if (0xC2..=0xDF).contains(&b) {
        2
    } else if (0xE0..=0xEF).contains(&b) {
        3
    } else if (0xF0..=0xF4).contains(&b) {
        4
    } else {
        1
    }
}

/// Byte length of the character at byte position `i` when it is a word
/// character (ASCII `[0-9A-Za-z_]` or a Unicode alphanumeric, mirroring the
/// reference regex's Unicode-mode `\w`). Combining marks (e.g. a decomposed
/// accent) are deliberately treated as non-word — std has no mark
/// predicate, and such characters are not expected in names.
fn word_char_len(s: &[u8], i: usize) -> usize {
    let b = s[i];
    if b < 0x80 {
        return if matches!(b, b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z' | b'_') {
            1
        } else {
            0
        };
    }
    let n = char_len_at(s, i);
    if i + n > s.len() {
        return 0;
    }
    match std::str::from_utf8(&s[i..i + n]) {
        Ok(t) if t.chars().next().is_some_and(char::is_alphanumeric) => n,
        _ => 0,
    }
}

/// Word characters plus the ASCII `/` and `-` (`[\w/-]` from the regex).
fn word_or_slash_char_len(s: &[u8], i: usize) -> usize {
    if matches!(s[i], b'/' | b'-') {
        1
    } else {
        word_char_len(s, i)
    }
}

fn is_digit(c: u8) -> bool {
    c.is_ascii_digit()
}

/// Try the whole range pattern at position `p` (leftmost-first, greedy).
fn match_at(s: &[u8], p: usize, n: usize) -> Option<(String, String, String, i32, i32)> {
    // Optional `name.` prefix: the regex's `[\w_]+` is greedy, and a shorter
    // name cannot be followed by '.' (the char before is a word char), so
    // only the full word run plus '.' needs trying; otherwise no name.
    let mut k = p;
    while k < n {
        let len = word_char_len(s, k);
        if len == 0 {
            break;
        }
        k += len;
    }
    if k < n && s[k] == b'.' {
        if let Some((chr, strand, start, end)) = rest_match(s, k + 1, n) {
            let name = String::from_utf8_lossy(&s[p..k]).into_owned();
            return Some((name, chr, strand, start, end));
        }
    }
    rest_match(s, p, n).map(|(chr, strand, start, end)| (String::new(), chr, strand, start, end))
}

/// Match `chr` (greedy `[\w/-]+`), optional `(strand)` (greedy `.+` ending at
/// the last `)`), then the `:start[-end]` tail.
fn rest_match(s: &[u8], q: usize, n: usize) -> Option<(String, String, i32, i32)> {
    let mut r = q;
    while r < n {
        let len = word_or_slash_char_len(s, r);
        if len == 0 {
            break;
        }
        r += len;
    }
    let mut chr_len = r - q;
    while chr_len >= 1 {
        let after = q + chr_len;
        if after < n && s[after] == b'(' {
            // `\(.+\)`: greedy `.+`, so try the last `)` first.
            let mut close = n;
            loop {
                close = close.saturating_sub(1);
                while close > after && s[close] != b')' {
                    close -= 1;
                }
                // `.+` inside `\(.+\)` needs at least one character, so the
                // closing `)` cannot immediately follow the opening `(`.
                if close <= after + 1 {
                    break;
                }
                if let Some((start, end)) = tail_match(s, close + 1, n) {
                    let chr = String::from_utf8_lossy(&s[q..after]).into_owned();
                    let strand = String::from_utf8_lossy(&s[after + 1..close]).into_owned();
                    return Some((chr, strand, start, end));
                }
            }
        } else if let Some((start, end)) = tail_match(s, after, n) {
            let chr = String::from_utf8_lossy(&s[q..after]).into_owned();
            return Some((chr, String::new(), start, end));
        }
        chr_len -= 1;
    }
    None
}

/// Match `:digits` with an optional single `_`/`-` spacer and optional
/// trailing digits (`end` defaults to `start`). Overflowing digit runs do
/// not match (invalid) instead of panicking.
fn tail_match(s: &[u8], t: usize, n: usize) -> Option<(i32, i32)> {
    if t >= n || s[t] != b':' {
        return None;
    }
    let mut u = t + 1;
    while u < n && is_digit(s[u]) {
        u += 1;
    }
    if u == t + 1 {
        return None; // `\d+` needs at least one digit
    }
    let start = parse_i32(&s[t + 1..u])?;
    let mut v = u;
    if v < n && (s[v] == b'_' || s[v] == b'-') {
        v += 1;
    }
    let mut w = v;
    while w < n && is_digit(s[w]) {
        w += 1;
    }
    let end = if w > v { parse_i32(&s[v..w])? } else { start };
    Some((start, end))
}

fn parse_i32(digits: &[u8]) -> Option<i32> {
    let mut v: i64 = 0;
    for &d in digits {
        v = v * 10 + i64::from(d - b'0');
        if v > i64::from(i32::MAX) {
            return None;
        }
    }
    Some(v as i32)
}

/// To string
///
/// ```
/// # use pgr::libs::ds::Range;
/// let range = Range::from("I", 1, 100);
/// assert_eq!(range.to_string(), "I:1-100");
/// let range = Range::from("I", 100, 100);
/// assert_eq!(range.to_string(), "I:100");
/// ```
impl fmt::Display for Range {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.encode())?;
        Ok(())
    }
}

#[test]
fn fa_headers() {
    let tests = vec![
        ("S288c", "S288c"),
        ("S288c The baker's yeast", "S288c"),
        ("1:-100", "1:-100"),
        ("infile_0/1/0_514:19-25", "infile_0/1/0_514:19-25"),
    ];
    for (header, expected) in tests {
        let range = Range::from_str(header);
        assert_eq!(range.to_string(), expected);
    }
}

#[test]
fn extreme_ops_do_not_overflow() {
    // Extreme `n` used to overflow `start + n` / `end - n` (debug panic,
    // release wrap); all ops must stay deterministic and non-panicking.
    // Out-of-order results collapse to the invalid (0, 0) range, while
    // in-order results saturate at the i32 bounds.
    let max = Range::from_str("chr1:2147483645-2147483645");
    assert!(!max.trim(i32::MAX).is_valid());
    assert!(!max.trim_5p(i32::MAX).is_valid());
    assert!(!max.trim_3p(i32::MAX).is_valid());
    assert!(!max.shift_5p(i32::MAX).is_valid());
    assert!(!max.flank_5p(i32::MAX).is_valid());
    assert_eq!(max.shift_3p(i32::MAX).to_string(), "chr1:2147483647");
    assert_eq!(
        max.flank_3p(i32::MAX).to_string(),
        "chr1:2147483646-2147483647"
    );
    // i32::MIN negation paths used to panic on `-n`.
    assert!(!max.trim(i32::MIN).is_valid());
    assert!(!max.trim_5p(i32::MIN).is_valid());
    assert!(!max.shift_3p(i32::MIN).is_valid());
}

#[test]
fn overflow_end_is_invalid_not_start() {
    // An oversized `end` used to be treated as "missing" and defaulted to
    // `start`, silently turning `chr1:5-99999999999` into the point range
    // `chr1:5`. It must be invalid instead.
    let overflow = Range::from_str("chr1:5-99999999999");
    assert!(!overflow.is_valid());
    let overflow = Range::from_str("chr1:5-99999999999999999999");
    assert!(!overflow.is_valid());
    let overflow = Range::from_str("S288c.I(-):5-99999999999");
    assert!(!overflow.is_valid());
    // Literal `0` and missing ends still default to `start` (regex-era
    // behavior), and a plain start without a dash is unaffected.
    assert_eq!(Range::from_str("chr1:5-0").to_string(), "chr1:5");
    assert_eq!(Range::from_str("chr1:5-").to_string(), "chr1:5");
    assert_eq!(Range::from_str("chr1:5").to_string(), "chr1:5");
    // A normal end still parses.
    assert_eq!(Range::from_str("chr1:5-100").to_string(), "chr1:5-100");
}

#[test]
fn regex_and_manual_decoders_agree() {
    // Differential test: the hand-written scanner must produce exactly the
    // same Range as the regex decoder, including the unanchored leftmost
    // match quirks (`foo I:1-100`, `a.b.c:1-2`) and the fallback behavior.
    let corpus = [
        "1:-100",
        "foo I:1-100",
        "S288c.I(-):27070-29557",
        "infile_0/1/0_514:19-25",
        "chr1(+):1-100",
        "I:100",
        "invalid",
        "S288c The baker's yeast",
        "1-100",
        "I:1-100x",
        "I:1x-100",
        "I(+-):1-100",
        "a.b.c:1-2",
        "chr1:1-100",
        "chrM(+):1-16571",
        "NC_000913:100-200",
        "S288c.II(+):1-813184",
        "1:1-23",
        "I:100",
        "I(+):100-200",
        "I(-):100-200",
        "infile_0/1/0_514:19-25",
        "123:456",
        "a-b:1-2",
        "x/y:1-2",
        "a(b)c:1-2",
        "a(b)(c):1-2",
        "I(a:b):1-2",
        ":1-2",
        "chr1:",
        "chr1:1-",
        "chr1:-2",
        "chr1:1_2",
        "chr1:1__2",
        // Unicode word characters must be accepted in name/chr/strand like
        // the reference regex's Unicode-mode `\w` (an ASCII-only scanner
        // used to treat these as non-matches and silently drop the line).
        "chr\u{3b1}:1-5",
        "\u{4e2d}(+):10-20",
        "S288c.\u{30c8}(-):100-200",
        "a.\u{0436}b:1-2",
        "\u{d55c}\u{0627}:3-4",
        "chr\u{00e9}:1-2",
        "",
        " ",
    ];
    for &s in &corpus {
        let regex = from_str_regex(s);
        let manual = Range::from_str(s);
        assert_eq!(
            (regex.name, regex.chr, regex.strand, regex.start, regex.end),
            (
                manual.name,
                manual.chr,
                manual.strand,
                manual.start,
                manual.end
            ),
            "regex vs manual divergence on {s:?}"
        );
    }

    // Randomized fuzz over the grammar's alphabet (coordinates kept small so
    // the regex path never overflows i32). Unicode letters are included;
    // combining marks are deliberately excluded (the scanner treats them as
    // non-word, see `word_char_len`).
    let mut x = 0x9E3779B97F4A7C15u64;
    let alphabet = [
        'a', 'b', 'c', 'Z', '0', '1', '9', '_', '/', '-', ':', '(', ')', '.', ' ', '\u{3b1}',
        '\u{4e2d}', '\u{30c8}', '\u{00e9}', '\u{4e00}', '\u{0436}', '\u{0627}', '\u{05d0}',
        '\u{0928}', '\u{d55c}',
    ];
    for _ in 0..40_000 {
        let len = (next_rand(&mut x) % 24) as usize;
        let s: String = (0..len)
            .map(|_| alphabet[(next_rand(&mut x) % alphabet.len() as u64) as usize])
            .collect();
        let regex = from_str_regex(&s);
        let manual = Range::from_str(&s);
        assert_eq!(
            (regex.name, regex.chr, regex.strand, regex.start, regex.end),
            (
                manual.name,
                manual.chr,
                manual.strand,
                manual.start,
                manual.end
            ),
            "regex vs manual divergence on {s:?}"
        );
    }
}

/// The original regex decoder, kept as the test oracle for the hand-written
/// scanner (see the module docs for the regex pattern itself).
#[cfg(test)]
fn from_str_regex(range: &str) -> Range {
    static RE: std::sync::LazyLock<regex::Regex> = std::sync::LazyLock::new(|| {
        regex::Regex::new(
            r"(?xi)
            (?:(?P<name>[\w_]+)\.)?
            (?P<chr>[\w/-]+)
            (?:\((?P<strand>.+)\))?
            [:]                    # spacer
            (?P<start>\d+)
            [_\-]?                 # spacer
            (?P<end>\d+)?
            ",
        )
        .expect("valid range regex")
    });
    let mut new = Range::new();
    let caps = match RE.captures(range) {
        Some(x) => x,
        None => {
            new.chr = range.split(' ').next().unwrap().to_string();
            return new;
        }
    };
    for name in RE.capture_names().flatten() {
        if let Some(m) = caps.name(name) {
            match name {
                "name" => new.name = m.as_str().to_string(),
                "chr" => new.chr = m.as_str().to_string(),
                "strand" => new.strand = m.as_str().to_string(),
                "start" => new.start = m.as_str().parse::<i32>().unwrap(),
                "end" => new.end = m.as_str().parse::<i32>().unwrap(),
                _ => {}
            }
        }
    }
    if new.start != 0 && new.end == 0 {
        new.end = new.start;
    }
    new
}

#[cfg(test)]
fn next_rand(x: &mut u64) -> u64 {
    *x = x
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *x >> 33
}
