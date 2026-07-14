//! Net data structures: NetNode, Space, Gap, Fill, Chrom.
//!
//! These types model the UCSC Net hierarchy (a tree of Gap/Fill nodes rooted
//! at a chromosome's root Gap). Inherent I/O methods (`write`) emit the
//! unfiltered UCSC Net text format.

use crate::libs::chain::record::Chain;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::io::{self, Write};
use std::rc::Rc;

/// A node in the net tree, either a gap or a fill.
#[derive(Clone, Debug)]
pub enum NetNode {
    /// An unaligned gap node.
    Gap(Rc<RefCell<Gap>>),
    /// An aligned fill node.
    Fill(Rc<RefCell<Fill>>),
}

/// A searchable space on a chromosome (points to its owning gap).
#[derive(Clone, Debug)]
pub struct Space {
    /// Start coordinate on the chromosome (0-based, inclusive).
    pub start: u64,
    /// End coordinate on the chromosome (0-based, exclusive).
    pub end: u64,
    /// Gap that owns this space.
    pub gap: Rc<RefCell<Gap>>,
}

/// A gap region (unaligned stretch) that may contain nested fills.
#[derive(Debug)]
pub struct Gap {
    /// Start coordinate on the target/query chromosome.
    pub start: u64,
    /// End coordinate on the target/query chromosome.
    pub end: u64,
    /// Start coordinate on the other genome's chromosome.
    pub o_start: u64,
    /// End coordinate on the other genome's chromosome.
    pub o_end: u64,
    /// Nested fill regions within this gap.
    pub fills: Vec<Rc<RefCell<Fill>>>,
    /// Count of N bases in the target portion.
    pub t_n: Option<u64>,
    /// Count of N bases in the query portion.
    pub q_n: Option<u64>,
    /// Count of repeat bases in the target portion.
    pub t_r: Option<u64>,
    /// Count of repeat bases in the query portion.
    pub q_r: Option<u64>,
    /// Count of TRF bases in the target portion.
    pub t_trf: Option<u64>,
    /// Count of TRF bases in the query portion.
    pub q_trf: Option<u64>,
}

/// A fill region (aligned stretch) referencing its source chain.
#[derive(Debug)]
pub struct Fill {
    /// Start coordinate on the target/query chromosome.
    pub start: u64,
    /// End coordinate on the target/query chromosome.
    pub end: u64,
    /// Start coordinate on the other genome's chromosome.
    pub o_start: u64,
    /// End coordinate on the other genome's chromosome.
    pub o_end: u64,
    /// Name of the other genome's chromosome.
    pub o_chrom: String,
    /// Strand of the other genome's chromosome.
    pub o_strand: char,
    /// ID of the source chain.
    pub chain_id: u64,
    /// Score of this fill (or source chain score).
    pub score: f64,
    /// Aligned bases in this fill.
    pub ali: u64,
    /// Synteny class (`top`, `syn`, `inv`, `nonSyn`, etc.).
    pub class: String,
    /// Query-side duplication overlap.
    pub q_dup: Option<u64>,
    /// Query overlap beyond the fill's own span.
    pub q_over: Option<u64>,
    /// Distance to the syntenic parent on the query.
    pub q_far: Option<i64>,
    /// Source chain reference (used for subchain scoring).
    pub chain: Option<Rc<Chain>>,
    /// Gaps nested inside this fill.
    pub gaps: Vec<Rc<RefCell<Gap>>>,
    /// Count of N bases in the target portion.
    pub t_n: Option<u64>,
    /// Count of N bases in the query portion.
    pub q_n: Option<u64>,
    /// Count of repeat bases in the target portion.
    pub t_r: Option<u64>,
    /// Count of repeat bases in the query portion.
    pub q_r: Option<u64>,
    /// Count of TRF bases in the target portion.
    pub t_trf: Option<u64>,
    /// Count of TRF bases in the query portion.
    pub q_trf: Option<u64>,
}

/// A chromosome's net tree root and searchable space index.
pub struct Chrom {
    /// Chromosome name.
    pub name: String,
    /// Chromosome size in bases.
    pub size: u64,
    /// Root gap covering the whole chromosome.
    pub root: Rc<RefCell<Gap>>,
    /// Map from space start coordinate to its owning `Space`.
    pub spaces: BTreeMap<u64, Space>,
    /// Header comment lines for this chromosome.
    pub comments: Vec<String>,
}

impl Chrom {
    /// Creates a new chromosome net with a root gap covering `[0, size)`.
    pub fn new(name: &str, size: u64) -> Self {
        let root = Rc::new(RefCell::new(Gap {
            start: 0,
            end: size,
            o_start: 0,
            o_end: 0, // Root gap o_range is 0? UCSC sets it to 0,0
            fills: Vec::new(),
            t_n: None,
            q_n: None,
            t_r: None,
            q_r: None,
            t_trf: None,
            q_trf: None,
        }));

        let space = Space {
            start: 0,
            end: size,
            gap: root.clone(),
        };

        let mut spaces = BTreeMap::new();
        spaces.insert(0, space);

        Chrom {
            name: name.to_string(),
            size,
            root,
            spaces,
            comments: Vec::new(),
        }
    }

    /// Returns all spaces overlapping `[start, end)`.
    pub fn find_spaces(&self, start: u64, end: u64) -> Vec<Space> {
        let mut result = Vec::new();
        for (_, space) in self.spaces.range(..end) {
            if space.end > start {
                result.push(space.clone());
            }
        }
        result
    }

    /// Writes this chromosome net in UCSC Net text format.
    pub fn write<W: Write>(&self, mut writer: W) -> io::Result<()> {
        for comment in &self.comments {
            writeln!(writer, "{}", comment)?;
        }
        writeln!(writer, "net {} {}", self.name, self.size)?;
        for fill in &self.root.borrow().fills {
            fill.borrow().write(&mut writer, 1)?;
        }
        Ok(())
    }
}

impl Fill {
    /// Writes this fill in UCSC Net text format at the given indentation level.
    pub fn write<W: Write>(&self, writer: &mut W, indent: usize) -> io::Result<()> {
        let indent_str = " ".repeat(indent);
        write!(
            writer,
            "{}fill {} {} {} {} {} {} id {} score {} ali {}",
            indent_str,
            self.start,
            self.end - self.start,
            self.o_chrom,
            self.o_strand,
            self.o_start,
            self.o_end - self.o_start,
            self.chain_id,
            self.score,
            self.ali
        )?;

        if let Some(val) = self.q_over {
            write!(writer, " qOver {}", val)?;
        }
        if let Some(val) = self.q_far {
            write!(writer, " qFar {}", val)?;
        }
        if let Some(val) = self.q_dup {
            write!(writer, " qDup {}", val)?;
        }
        if !self.class.is_empty() {
            write!(writer, " type {}", self.class)?;
        }
        if let Some(val) = self.t_n {
            write!(writer, " tN {}", val)?;
        }
        if let Some(val) = self.q_n {
            write!(writer, " qN {}", val)?;
        }
        if let Some(val) = self.t_r {
            write!(writer, " tR {}", val)?;
        }
        if let Some(val) = self.q_r {
            write!(writer, " qR {}", val)?;
        }
        if let Some(val) = self.t_trf {
            write!(writer, " tTrf {}", val)?;
        }
        if let Some(val) = self.q_trf {
            write!(writer, " qTrf {}", val)?;
        }
        writeln!(writer)?;

        for gap in &self.gaps {
            gap.borrow()
                .write(writer, indent + 1, &self.o_chrom, self.o_strand)?;
        }
        Ok(())
    }
}

impl Gap {
    /// Writes this gap in UCSC Net text format at the given indentation level.
    pub fn write<W: Write>(
        &self,
        writer: &mut W,
        indent: usize,
        o_chrom: &str,
        o_strand: char,
    ) -> io::Result<()> {
        let indent_str = " ".repeat(indent);
        write!(
            writer,
            "{}gap {} {} {} {} {} {}",
            indent_str,
            self.start,
            self.end - self.start,
            o_chrom,
            o_strand,
            self.o_start,
            self.o_end - self.o_start
        )?;

        if let Some(val) = self.t_n {
            write!(writer, " tN {}", val)?;
        }
        if let Some(val) = self.q_n {
            write!(writer, " qN {}", val)?;
        }
        if let Some(val) = self.t_r {
            write!(writer, " tR {}", val)?;
        }
        if let Some(val) = self.q_r {
            write!(writer, " qR {}", val)?;
        }
        if let Some(val) = self.t_trf {
            write!(writer, " tTrf {}", val)?;
        }
        if let Some(val) = self.q_trf {
            write!(writer, " qTrf {}", val)?;
        }
        writeln!(writer)?;

        for fill in &self.fills {
            fill.borrow().write(writer, indent + 1)?;
        }
        Ok(())
    }
}
