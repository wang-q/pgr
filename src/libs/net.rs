use crate::libs::chaining::record::{Block, Chain};
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::io::{self, BufRead, Write};
use std::rc::Rc;

#[derive(Clone, Debug)]
pub enum NetNode {
    Gap(Rc<RefCell<Gap>>),
    Fill(Rc<RefCell<Fill>>),
}

#[derive(Clone, Debug)]
pub struct Space {
    pub start: u64,
    pub end: u64,
    pub gap: Rc<RefCell<Gap>>,
}

#[derive(Debug)]
pub struct Gap {
    pub start: u64,
    pub end: u64,
    pub o_start: u64,
    pub o_end: u64,
    pub fills: Vec<Rc<RefCell<Fill>>>,
}

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
    pub q_dup: u64,
    pub q_over: u64,
    pub q_far: i64,
    pub chain: Option<Rc<Chain>>,
    pub gaps: Vec<Rc<RefCell<Gap>>>,
}

pub struct Chrom {
    pub name: String,
    pub size: u64,
    pub root: Rc<RefCell<Gap>>,
    pub spaces: BTreeMap<u64, Space>, // start -> Space
}

impl Chrom {
    pub fn new(name: &str, size: u64) -> Self {
        let root = Rc::new(RefCell::new(Gap {
            start: 0,
            end: size,
            o_start: 0,
            o_end: 0, // Root gap o_range is 0? UCSC sets it to 0,0
            fills: Vec::new(),
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
        }
    }

    pub fn find_spaces(&self, start: u64, end: u64) -> Vec<Space> {
        let mut result = Vec::new();
        // Iterate over spaces that might overlap
        // We can start from the last key <= start, but BTreeMap doesn't support that easily in stable Rust without range
        // range(..end) gives keys < end.
        for (_, space) in self.spaces.range(..end) {
            if space.end > start {
                result.push(space.clone());
            }
        }
        result
    }

    pub fn write<W: Write>(&self, mut writer: W) -> io::Result<()> {
        writeln!(writer, "net {} {}", self.name, self.size)?;
        // The root gap contains the top-level fills.
        // We don't write the root gap itself as a "gap" line, but we write its children.
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

        if self.q_dup > 0 {
            write!(writer, " qDup {}", self.q_dup)?;
        }
        if !self.class.is_empty() {
            write!(writer, " type {}", self.class)?;
        }
        if self.q_over > 0 {
            write!(writer, " qOver {}", self.q_over)?;
        }
        if self.q_far != 0 {
            write!(writer, " qFar {}", self.q_far)?;
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
        writeln!(
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

        for fill in &self.fills {
            fill.borrow().write(writer, indent + 1)?;
        }
        Ok(())
    }
}

pub fn read_nets<R: BufRead>(mut reader: R) -> io::Result<Vec<Chrom>> {
    let mut chroms = Vec::new();
    let mut current_chrom: Option<Chrom> = None;
    let mut stack: Vec<(usize, NetNode)> = Vec::new();

    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        if line.trim().is_empty() {
            line.clear();
            continue;
        }

        let mut indent = 0;
        for c in line.chars() {
            if c == ' ' {
                indent += 1;
            } else {
                break;
            }
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            line.clear();
            continue;
        }

        match parts[0] {
            "net" => {
                if let Some(c) = current_chrom {
                    chroms.push(c);
                }
                let name = parts[1];
                let size = parts[2].parse::<u64>().unwrap();
                let chrom = Chrom::new(name, size);
                stack.clear();
                stack.push((0, NetNode::Gap(chrom.root.clone())));
                current_chrom = Some(chrom);
            }
            "fill" => {
                // fill tStart tLength qName qStrand qStart qLength id chainId score ali [type class]
                let start = parts[1].parse::<u64>().unwrap();
                let len = parts[2].parse::<u64>().unwrap();
                let q_name = parts[3].to_string();
                let q_strand = parts[4].chars().next().unwrap();
                let q_start = parts[5].parse::<u64>().unwrap();
                let q_len = parts[6].parse::<u64>().unwrap();
                // parts[7] is "id"
                let chain_id = parts[8].parse::<u64>().unwrap();
                // parts[9] is "score"
                let score = parts[10].parse::<f64>().unwrap();
                // parts[11] is "ali"
                let ali = parts[12].parse::<u64>().unwrap();

                let mut class = String::new();
                let mut q_dup = 0;
                let mut q_over = 0;
                let mut q_far = 0;

                let mut i = 13;
                while i < parts.len() {
                    match parts[i] {
                        "type" => {
                            if i + 1 < parts.len() {
                                class = parts[i + 1].to_string();
                                i += 2;
                            } else {
                                i += 1;
                            }
                        }
                        "qDup" => {
                            if i + 1 < parts.len() {
                                q_dup = parts[i + 1].parse::<u64>().unwrap_or(0);
                                i += 2;
                            } else {
                                i += 1;
                            }
                        }
                        "qOver" => {
                            if i + 1 < parts.len() {
                                q_over = parts[i + 1].parse::<u64>().unwrap_or(0);
                                i += 2;
                            } else {
                                i += 1;
                            }
                        }
                        "qFar" => {
                            if i + 1 < parts.len() {
                                q_far = parts[i + 1].parse::<i64>().unwrap_or(0);
                                i += 2;
                            } else {
                                i += 1;
                            }
                        }
                        _ => {
                            i += 1;
                        }
                    }
                }

                let fill = Rc::new(RefCell::new(Fill {
                    start,
                    end: start + len,
                    o_start: q_start,
                    o_end: q_start + q_len,
                    o_chrom: q_name,
                    o_strand: q_strand,
                    chain_id,
                    score,
                    ali,
                    class,
                    q_dup,
                    q_over,
                    q_far,
                    chain: None,
                    gaps: Vec::new(),
                }));

                // Find parent gap
                while let Some((parent_indent, parent_node)) = stack.last() {
                    if indent > *parent_indent {
                        if let NetNode::Gap(gap) = parent_node {
                            gap.borrow_mut().fills.push(fill.clone());
                            stack.push((indent, NetNode::Fill(fill)));
                            break;
                        } else {
                            // Parent is Fill, but Fill cannot have Fill children directly (must be via Gap)
                            // But maybe the file format allows implicit gaps? No, UCSC net format is strict.
                            // However, if we see indent > parent_indent and parent is Fill, something is wrong or I misunderstood.
                            // In Net format: fill -> gap -> fill
                            // So parent of fill must be gap.
                            stack.pop();
                        }
                    } else {
                        stack.pop();
                    }
                }
            }
            "gap" => {
                // gap tStart tLength qName qStrand qStart qLength
                let start = parts[1].parse::<u64>().unwrap();
                let len = parts[2].parse::<u64>().unwrap();
                let _q_name = parts[3].to_string();
                let _q_strand = parts[4].chars().next().unwrap();
                let q_start = parts[5].parse::<u64>().unwrap();
                let q_len = parts[6].parse::<u64>().unwrap();

                let gap = Rc::new(RefCell::new(Gap {
                    start,
                    end: start + len,
                    o_start: q_start,
                    o_end: q_start + q_len,
                    fills: Vec::new(),
                }));

                // Find parent fill
                while let Some((parent_indent, parent_node)) = stack.last() {
                    if indent > *parent_indent {
                        if let NetNode::Fill(fill) = parent_node {
                            fill.borrow_mut().gaps.push(gap.clone());
                            stack.push((indent, NetNode::Gap(gap)));
                            break;
                        } else {
                            stack.pop();
                        }
                    } else {
                        stack.pop();
                    }
                }
            }
            _ => {}
        }
        line.clear();
    }
    if let Some(c) = current_chrom {
        chroms.push(c);
    }
    Ok(chroms)
}

pub struct ChainNet {
    pub chroms: HashMap<String, RefCell<Chrom>>, // Use RefCell to allow mutation of Chroms
    pub chains: Vec<Rc<Chain>>,
}

impl ChainNet {
    pub fn new(target_sizes: &HashMap<String, u64>) -> Self {
        let mut chroms = HashMap::new();
        for (name, size) in target_sizes {
            chroms.insert(name.clone(), RefCell::new(Chrom::new(name, *size)));
        }
        Self {
            chroms,
            chains: Vec::new(),
        }
    }

    pub fn add_chain(&mut self, chain: Chain, min_space: u64, min_fill: u64) {
        let chain_rc = Rc::new(chain);
        self.chains.push(chain_rc.clone());

        // Add to target net
        if let Some(chrom) = self.chroms.get(&chain_rc.header.t_name) {
            let mut chrom = chrom.borrow_mut();
            let blocks = chain_rc.to_blocks();
            add_chain_core(
                &mut chrom,
                chain_rc.clone(),
                blocks,
                false,
                min_space,
                min_fill,
            );
        }
    }

    pub fn add_chain_as_q(&mut self, chain: Chain, min_space: u64, min_fill: u64) {
        let chain_rc = Rc::new(chain);
        self.chains.push(chain_rc.clone());

        if let Some(chrom) = self.chroms.get(&chain_rc.header.q_name) {
            let mut chrom = chrom.borrow_mut();
            let mut blocks = chain_rc.to_blocks();

            if chain_rc.header.q_strand == '-' {
                reverse_blocks_q(&mut blocks, chain_rc.header.q_size);
            }

            add_chain_core(
                &mut chrom,
                chain_rc.clone(),
                blocks,
                true,
                min_space,
                min_fill,
            );
        }
    }
}

// Helper to calculate intersection
pub fn range_intersection(start1: u64, end1: u64, start2: u64, end2: u64) -> u64 {
    let s = start1.max(start2);
    let e = end1.min(end2);
    e.saturating_sub(s)
}

fn reverse_range(start: &mut u64, end: &mut u64, size: u64) {
    let tmp = *start;
    *start = size - *end;
    *end = size - tmp;
}

fn chain_base_count_sub_t(chain: &Chain, t_min: u64, t_max: u64) -> u64 {
    let mut total = 0;
    // We need block list. chain.to_blocks() returns blocks with absolute coords.
    // However, recreating blocks every time is expensive.
    // Ideally Chain should store blocks or we iterate data.
    // ChainData: size, dt, dq.
    let mut t_curr = chain.header.t_start;
    for d in &chain.data {
        let t_start = t_curr;
        let t_end = t_curr + d.size;
        total += range_intersection(t_start, t_end, t_min, t_max);
        t_curr += d.size + d.dt;
    }
    total
}

fn chain_base_count_sub_q(chain: &Chain, q_min: u64, q_max: u64) -> u64 {
    let mut total = 0;
    let mut q_curr = chain.header.q_start;
    for d in &chain.data {
        let q_start = q_curr;
        let q_end = q_curr + d.size;
        total += range_intersection(q_start, q_end, q_min, q_max);
        q_curr += d.size + d.dq;
    }
    total
}

fn subchain_info(chain: &Chain, start: u64, end: u64, is_q: bool) -> (u64, f64) {
    // let full_size = chain.header.score as u64; // Approx? No, chain.header.score is score.
    // UCSC chainBaseCount calculates bases in gap-free alignments.
    // We need to calculate full_size first.
    let mut full_ali_size = 0;
    for d in &chain.data {
        full_ali_size += d.size;
    }

    if full_ali_size == 0 {
        return (0, 0.0);
    }

    let sub_size = if is_q {
        let mut s = start;
        let mut e = end;
        if chain.header.q_strand == '-' {
            reverse_range(&mut s, &mut e, chain.header.q_size);
        }
        if s <= chain.header.q_start && e >= chain.header.q_end {
            full_ali_size
        } else {
            chain_base_count_sub_q(chain, s, e)
        }
    } else if start <= chain.header.t_start && end >= chain.header.t_end {
        full_ali_size
    } else {
        chain_base_count_sub_t(chain, start, end)
    };

    let sub_score = chain.header.score * (sub_size as f64) / (full_ali_size as f64);
    (sub_size, sub_score)
}

fn reverse_blocks_q(blocks: &mut [Block], size: u64) {
    blocks.reverse();
    for b in blocks {
        reverse_range(&mut b.q_start, &mut b.q_end, size);
    }
}

fn add_chain_core(
    chrom: &mut Chrom,
    chain: Rc<Chain>,
    blocks: Vec<Block>,
    is_q: bool,
    min_space: u64,
    min_fill: u64,
) {
    let (start, end) = if is_q {
        let mut s = chain.header.q_start;
        let mut e = chain.header.q_end;
        if chain.header.q_strand == '-' {
            reverse_range(&mut s, &mut e, chain.header.q_size);
        }
        (s, e)
    } else {
        (chain.header.t_start, chain.header.t_end)
    };

    let spaces = chrom.find_spaces(start, end);
    let mut start_block_idx = 0;

    for space in spaces {
        let mut first_idx = None;
        let mut last_idx = None;
        let mut s = u64::MAX;
        let mut e = 0;

        for i in start_block_idx..blocks.len() {
            let b = &blocks[i];
            let (b_start, b_end) = if is_q {
                (b.q_start, b.q_end)
            } else {
                (b.t_start, b.t_end)
            };

            if b_end <= space.start {
                continue;
            }
            if b_start >= space.end {
                break;
            }

            if first_idx.is_none() {
                first_idx = Some(i);
            }
            last_idx = Some(i);

            let curr_s: u64 = b_start.max(space.start);
            let curr_e: u64 = b_end.min(space.end);

            if curr_s < s {
                s = curr_s;
            }
            if curr_e > e {
                e = curr_e;
            }
        }

        if let Some(idx) = first_idx {
            start_block_idx = idx;
        } else {
            continue;
        }

        if s >= e || (e - s) < min_fill {
            continue;
        }

        fill_space(
            chrom,
            space,
            chain.clone(),
            &blocks,
            first_idx.unwrap(),
            last_idx.unwrap(),
            s,
            e,
            min_space,
            is_q,
        );
    }
}

// Logic to add chain to Q side would be similar but swapping coords and handling strand
// For now, let's assume we use add_chain_t for both, by constructing a "proxy" chain with swapped coords if needed.
// Or we implement add_chain_q separately.

fn fill_space(
    chrom: &mut Chrom,
    space: Space,
    chain: Rc<Chain>,
    blocks: &[Block],
    first_idx: usize,
    last_idx: usize,
    fill_start: u64,
    fill_end: u64,
    min_space: u64,
    is_q: bool,
) {
    // Remove old space
    chrom.spaces.remove(&space.start);

    // Calculate other side coords for the fill
    let (o_start, o_end) = if !is_q {
        let b1 = &blocks[first_idx];
        let offset1 = fill_start - b1.t_start;
        let q1 = b1.q_start + offset1;

        let b2 = &blocks[last_idx];
        let offset2 = fill_end - b2.t_start;
        let q2 = b2.q_start + offset2;

        (q1, q2)
    } else {
        let b1 = &blocks[first_idx];
        let offset1 = fill_start - b1.q_start;

        let t1 = if chain.header.q_strand == '-' {
            b1.t_end - offset1
        } else {
            b1.t_start + offset1
        };

        let b2 = &blocks[last_idx];
        let offset2 = fill_end - b2.q_start;
        let t2 = if chain.header.q_strand == '-' {
            b2.t_end - offset2
        } else {
            b2.t_start + offset2
        };

        if t1 > t2 {
            (t2, t1)
        } else {
            (t1, t2)
        }
    };

    let o_chrom = if is_q {
        &chain.header.t_name
    } else {
        &chain.header.q_name
    };
    let o_strand = if is_q { '+' } else { chain.header.q_strand };

    // Create Fill
    let fill = Rc::new(RefCell::new(Fill {
        start: fill_start,
        end: fill_end,
        o_start,
        o_end,
        o_chrom: o_chrom.clone(),
        o_strand,
        chain_id: chain.header.id,
        score: 0.0,
        ali: 0,
        class: String::new(),
        q_dup: 0,
        q_over: 0,
        q_far: 0,
        chain: Some(chain.clone()),
        gaps: Vec::new(),
    }));

    // Add Left Space
    if fill_start > space.start && (fill_start - space.start) >= min_space {
        chrom.spaces.insert(
            space.start,
            Space {
                start: space.start,
                end: fill_start,
                gap: space.gap.clone(),
            },
        );
    }

    // Add Right Space
    if fill_end < space.end && (space.end - fill_end) >= min_space {
        chrom.spaces.insert(
            fill_end,
            Space {
                start: fill_end,
                end: space.end,
                gap: space.gap.clone(),
            },
        );
    }

    // Internal gaps
    for i in first_idx..last_idx {
        let b1 = &blocks[i];
        let b2 = &blocks[i + 1];

        let (gap_start, gap_end) = if is_q {
            (b1.q_end, b2.q_start)
        } else {
            (b1.t_end, b2.t_start)
        };

        if gap_start > fill.borrow().start
            && gap_end < fill.borrow().end
            && (gap_end - gap_start) >= min_space
        {
            let (os, oe) = if !is_q {
                (b1.q_end, b2.q_start)
            } else if chain.header.q_strand == '-' {
                (b2.t_end, b1.t_start)
            } else {
                (b1.t_end, b2.t_start)
            };

            let new_gap = Rc::new(RefCell::new(Gap {
                start: gap_start,
                end: gap_end,
                o_start: os,
                o_end: oe,
                fills: Vec::new(),
            }));

            chrom.spaces.insert(
                gap_start,
                Space {
                    start: gap_start,
                    end: gap_end,
                    gap: new_gap.clone(),
                },
            );

            fill.borrow_mut().gaps.push(new_gap);
        }
    }

    // Add fill to parent gap
    space.gap.borrow_mut().fills.push(fill);
}

// Calculate o_start/o_end for fills
pub fn finalize_net(chrom: &mut Chrom, is_q: bool) {
    // Sort fills/gaps and calculate other ranges
    sort_net(&chrom.root);
    calc_other_fill(&chrom.root, is_q);
}

fn sort_net(gap: &Rc<RefCell<Gap>>) {
    let mut gap_borrow = gap.borrow_mut();
    gap_borrow.fills.sort_by_key(|f| f.borrow().start);

    for fill in &gap_borrow.fills {
        let mut fill_borrow = fill.borrow_mut();
        fill_borrow.gaps.sort_by_key(|g| g.borrow().start);
        for g in &fill_borrow.gaps {
            sort_net(g);
        }
    }
}

fn calc_other_fill(gap: &Rc<RefCell<Gap>>, is_q: bool) {
    let gap_borrow = gap.borrow();
    for fill in &gap_borrow.fills {
        let mut fill_borrow = fill.borrow_mut();

        if let Some(chain) = fill_borrow.chain.clone() {
            let clip_start = fill_borrow.start;
            let clip_end = fill_borrow.end;

            if !is_q {
                let mut q_min = u64::MAX;
                let mut q_max = 0;

                let mut t_curr = chain.header.t_start;
                let mut q_curr = chain.header.q_start;

                for d in &chain.data {
                    let t_s = t_curr;
                    let t_e = t_curr + d.size;
                    let q_s = q_curr;

                    let start = t_s.max(clip_start);
                    let end = t_e.min(clip_end);

                    if start < end {
                        let offset = start - t_s;
                        let len = end - start;
                        let qs = q_s + offset;
                        let qe = qs + len;

                        if qs < q_min {
                            q_min = qs;
                        }
                        if qe > q_max {
                            q_max = qe;
                        }
                    }

                    t_curr += d.size + d.dt;
                    q_curr += d.size + d.dq;
                }

                if q_min < q_max {
                    if chain.header.q_strand == '-' {
                        let mut s = q_min;
                        let mut e = q_max;
                        reverse_range(&mut s, &mut e, chain.header.q_size);
                        fill_borrow.o_start = s;
                        fill_borrow.o_end = e;
                    } else {
                        fill_borrow.o_start = q_min;
                        fill_borrow.o_end = q_max;
                    }
                }
            } else {
                let mut t_min = u64::MAX;
                let mut t_max = 0;

                let mut t_curr = chain.header.t_start;
                let mut q_curr = chain.header.q_start;

                for d in &chain.data {
                    let t_s = t_curr;
                    let q_s = q_curr;
                    let q_e = q_curr + d.size;

                    let (c_start, c_end) = if chain.header.q_strand == '-' {
                        let mut s = clip_start;
                        let mut e = clip_end;
                        reverse_range(&mut s, &mut e, chain.header.q_size);
                        (s, e)
                    } else {
                        (clip_start, clip_end)
                    };

                    let start = q_s.max(c_start);
                    let end = q_e.min(c_end);

                    if start < end {
                        let offset = start - q_s;
                        let len = end - start;
                        let ts = t_s + offset;
                        let te = ts + len;

                        if ts < t_min {
                            t_min = ts;
                        }
                        if te > t_max {
                            t_max = te;
                        }
                    }

                    t_curr += d.size + d.dt;
                    q_curr += d.size + d.dq;
                }

                if t_min < t_max {
                    fill_borrow.o_start = t_min;
                    fill_borrow.o_end = t_max;
                }
            }
        }

        drop(fill_borrow);
        for g in &fill.borrow().gaps {
            calc_other_fill(g, is_q);
        }
    }
}

// Writing
pub fn write_net<W: Write>(
    chrom: &Chrom,
    writer: &mut W,
    is_q: bool,
    min_score: f64,
    min_fill: u64,
) -> io::Result<()> {
    if chrom.root.borrow().fills.is_empty() {
        return Ok(());
    }
    writeln!(writer, "net {} {}", chrom.name, chrom.size)?;

    for fill in &chrom.root.borrow().fills {
        write_fill(fill, writer, 1, is_q, min_score, min_fill)?;
    }
    Ok(())
}

fn write_fill<W: Write>(
    fill: &Rc<RefCell<Fill>>,
    writer: &mut W,
    depth: usize,
    is_q: bool,
    min_score: f64,
    min_fill: u64,
) -> io::Result<()> {
    let f = fill.borrow();

    // Calculate subscore/subsize if chain is available, otherwise use stored
    let (sub_size, sub_score) = if let Some(chain) = &f.chain {
        subchain_info(chain, f.start, f.end, is_q)
    } else {
        (f.ali, f.score)
    };

    if sub_score >= min_score && sub_size >= min_fill {
        write_indent(writer, depth)?;

        let mut line = format!(
            "fill {} {} {} {} {} {} id {} score {:.0} ali {}",
            f.start,
            f.end - f.start,
            f.o_chrom,
            f.o_strand,
            f.o_start,
            f.o_end - f.o_start,
            f.chain_id,
            sub_score,
            sub_size
        );

        if !f.class.is_empty() {
            line.push_str(" type ");
            line.push_str(&f.class);
        }

        if f.q_dup > 0 {
            line.push_str(" qDup ");
            line.push_str(&f.q_dup.to_string());
        }

        if f.q_over > 0 {
            line.push_str(" qOver ");
            line.push_str(&f.q_over.to_string());
        }

        if f.q_far != 0 {
            line.push_str(" qFar ");
            line.push_str(&f.q_far.to_string());
        }

        writeln!(writer, "{}", line)?;

        for gap in &f.gaps {
            write_gap(gap, fill, writer, depth + 1, is_q, min_score, min_fill)?;
        }
    }
    Ok(())
}

fn write_gap<W: Write>(
    gap: &Rc<RefCell<Gap>>,
    parent: &Rc<RefCell<Fill>>,
    writer: &mut W,
    depth: usize,
    is_q: bool,
    min_score: f64,
    min_fill: u64,
) -> io::Result<()> {
    let g = gap.borrow();
    let p = parent.borrow();
    let o_chrom = &p.o_chrom;
    let o_strand = p.o_strand;

    write_indent(writer, depth)?;
    writeln!(
        writer,
        "gap {} {} {} {} {} {}",
        g.start,
        g.end - g.start,
        o_chrom,
        o_strand,
        g.o_start,
        g.o_end - g.o_start
    )?;

    for fill in &g.fills {
        write_fill(fill, writer, depth + 1, is_q, min_score, min_fill)?;
    }
    Ok(())
}

fn write_indent<W: Write>(writer: &mut W, depth: usize) -> io::Result<()> {
    for _ in 0..depth {
        write!(writer, " ")?;
    }
    Ok(())
}
