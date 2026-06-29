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

/// A node in the net tree (gap or fill).
#[derive(Clone, Debug)]
pub enum NetNode {
    Gap(Rc<RefCell<Gap>>),
    Fill(Rc<RefCell<Fill>>),
}

/// A searchable space on a chromosome (points to its owning gap).
#[derive(Clone, Debug)]
pub struct Space {
    pub start: u64,
    pub end: u64,
    pub gap: Rc<RefCell<Gap>>,
}

/// A gap region (unaligned stretch) that may contain nested fills.
#[derive(Debug)]
pub struct Gap {
    pub start: u64,
    pub end: u64,
    pub o_start: u64,
    pub o_end: u64,
    pub fills: Vec<Rc<RefCell<Fill>>>,
    pub t_n: Option<u64>,
    pub q_n: Option<u64>,
    pub t_r: Option<u64>,
    pub q_r: Option<u64>,
    pub t_trf: Option<u64>,
    pub q_trf: Option<u64>,
}

/// A fill region (aligned stretch) referencing its source chain.
#[derive(Debug)]
pub struct Fill {
    pub start: u64,
    pub end: u64,
    pub o_start: u64,
    pub o_end: u64,
    pub o_chrom: String,
    pub o_strand: char,
    pub chain_id: u64,
    pub score: f64,
    pub ali: u64,
    pub class: String,
    pub q_dup: Option<u64>,
    pub q_over: Option<u64>,
    pub q_far: Option<i64>,
    pub chain: Option<Rc<Chain>>,
    pub gaps: Vec<Rc<RefCell<Gap>>>,
    pub t_n: Option<u64>,
    pub q_n: Option<u64>,
    pub t_r: Option<u64>,
    pub q_r: Option<u64>,
    pub t_trf: Option<u64>,
    pub q_trf: Option<u64>,
}

/// A chromosome's net tree root + searchable space index.
pub struct Chrom {
    pub name: String,
    pub size: u64,
    pub root: Rc<RefCell<Gap>>,
    pub spaces: BTreeMap<u64, Space>, // start -> Space
    pub comments: Vec<String>,
}

impl Chrom {
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

    pub fn find_spaces(&self, start: u64, end: u64) -> Vec<Space> {
        let mut result = Vec::new();
        for (_, space) in self.spaces.range(..end) {
            if space.end > start {
                result.push(space.clone());
            }
        }
        result
    }

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
