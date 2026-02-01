use clap::{Arg, ArgMatches, Command};
use pgr::libs::net::{range_intersection, read_nets, write_net, Chrom, Fill, Gap};
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::rc::Rc;

pub fn make_subcommand() -> Command {
    Command::new("syntenic")
        .about("Add synteny info to net")
        .arg(Arg::new("in_net").required(true).help("Input net file"))
        .arg(Arg::new("out_net").required(true).help("Output net file"))
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let in_file = matches.get_one::<String>("in_net").unwrap();
    let out_file = matches.get_one::<String>("out_net").unwrap();

    let reader = BufReader::new(File::open(in_file)?);
    let nets = read_nets(reader)?;

    // Build DupeTrees for all query chromosomes
    let mut q_chrom_map: HashMap<String, DupeTree> = HashMap::new();

    for net in &nets {
        r_calc_dupes(net, &mut q_chrom_map);
    }

    // Process DupeTrees
    for dt in q_chrom_map.values_mut() {
        dt.build();
    }

    // Classify
    for net in &nets {
        r_net_syn(net, &q_chrom_map);
    }

    // Write output
    let mut writer = BufWriter::new(File::create(out_file)?);
    for net in &nets {
        write_net(net, &mut writer, false, 0.0, 0)?;
    }

    Ok(())
}

fn r_calc_dupes(chrom: &Chrom, map: &mut HashMap<String, DupeTree>) {
    r_calc_dupes_gap(&chrom.root, map);
}

fn r_calc_dupes_gap(gap: &Rc<RefCell<Gap>>, map: &mut HashMap<String, DupeTree>) {
    let g = gap.borrow();
    for fill in &g.fills {
        r_calc_dupes_fill(fill, map);
    }
}

fn r_calc_dupes_fill(fill: &Rc<RefCell<Fill>>, map: &mut HashMap<String, DupeTree>) {
    let f = fill.borrow();
    let q_name = &f.o_chrom;
    let start = f.o_start;
    let end = f.o_end;

    if !q_name.is_empty() {
        let dt = map.entry(q_name.clone()).or_insert_with(DupeTree::new);
        dt.add(start, end);
    }

    // Recursively process gaps inside fill
    for gap in &f.gaps {
        let g = gap.borrow();
        // Gap inside Fill shares query chrom with Fill
        // But Gap subtracts coverage
        let q_name = &f.o_chrom;
        let start = g.o_start;
        let end = g.o_end;

        if !q_name.is_empty() {
            let dt = map.entry(q_name.clone()).or_insert_with(DupeTree::new);
            dt.subtract(start, end);
        }

        // Recurse into fills inside gap
        r_calc_dupes_gap(gap, map);
    }
}

fn r_net_syn(chrom: &Chrom, map: &HashMap<String, DupeTree>) {
    r_net_syn_gap(&chrom.root, map, None);
}

fn r_net_syn_gap(
    gap: &Rc<RefCell<Gap>>,
    map: &HashMap<String, DupeTree>,
    parent_fill: Option<&Rc<RefCell<Fill>>>,
) {
    let g = gap.borrow();
    for fill in &g.fills {
        r_net_syn_fill(fill, map, parent_fill);
    }
}

fn r_net_syn_fill(
    fill: &Rc<RefCell<Fill>>,
    map: &HashMap<String, DupeTree>,
    parent: Option<&Rc<RefCell<Fill>>>,
) {
    // Need to borrow_mut to update fields
    // But we also need to pass `fill` (Rc) to children.
    // So we borrow mut, update, drop borrow, then recurse.

    let (q_name, start, end, strand) = {
        let f = fill.borrow();
        (f.o_chrom.clone(), f.o_start, f.o_end, f.o_strand)
    };

    let mut q_dup = 0;
    if let Some(dt) = map.get(&q_name) {
        q_dup = dt.count_over(start, end, 2);
    }

    let type_str;
    let mut q_over = 0;
    let mut q_far = 0;

    if parent.is_none() {
        type_str = "top".to_string();
    } else {
        let p = parent.unwrap().borrow();
        if q_name != p.o_chrom {
            type_str = "nonSyn".to_string();
        } else {
            let p_start = p.o_start;
            let p_end = p.o_end;
            let intersection = range_intersection(start, end, p_start, p_end);

            if intersection > 0 {
                q_over = intersection;
                q_far = 0;
            } else {
                q_over = 0;
                q_far = -(intersection as i64); // Always 0 if intersection returns 0
            }

            if p.o_strand == strand {
                type_str = "syn".to_string();
            } else {
                type_str = "inv".to_string();
            }
        }
    }

    {
        let mut f = fill.borrow_mut();
        f.class = type_str;
        f.q_dup = q_dup;
        f.q_over = q_over;
        f.q_far = q_far;
    }

    // Recurse
    // Children of fill are in `f.gaps`
    // We need to access `f.gaps` without holding mutable borrow on `f`
    let gaps = fill.borrow().gaps.clone();
    for gap in gaps {
        r_net_syn_gap(&gap, map, Some(fill));
    }
}

// DupeTree implementation
struct DupeTree {
    intervals: Vec<(u64, u64, i32)>,
    segments: Vec<Segment>,
}

struct Segment {
    start: u64,
    end: u64,
    depth: i32,
}

impl DupeTree {
    fn new() -> Self {
        Self {
            intervals: Vec::new(),
            segments: Vec::new(),
        }
    }

    fn add(&mut self, start: u64, end: u64) {
        if start < end {
            self.intervals.push((start, end, 1));
        }
    }

    fn subtract(&mut self, start: u64, end: u64) {
        if start < end {
            self.intervals.push((start, end, -1));
        }
    }

    fn build(&mut self) {
        if self.intervals.is_empty() {
            return;
        }

        let mut events = Vec::new();
        for (s, e, d) in &self.intervals {
            events.push((*s, *d));
            events.push((*e, -*d));
        }
        // Sort by position, then by delta (to process all updates at same pos)
        // Actually, order of processing at same position doesn't matter for final segments between positions.
        events.sort_by(|a, b| a.0.cmp(&b.0));

        let mut current_depth = 0;
        let mut segments = Vec::new();

        for i in 0..events.len() - 1 {
            let (pos, delta) = events[i];
            current_depth += delta;

            let next_pos = events[i + 1].0;
            if next_pos > pos {
                segments.push(Segment {
                    start: pos,
                    end: next_pos,
                    depth: current_depth,
                });
            }
        }

        self.segments = segments;
    }

    fn count_over(&self, start: u64, end: u64, threshold: i32) -> u64 {
        if self.segments.is_empty() {
            return 0;
        }

        // Binary search for first segment ending after start
        let idx = self.segments.partition_point(|seg| seg.end <= start);

        let mut count = 0;
        for seg in &self.segments[idx..] {
            if seg.start >= end {
                break;
            }

            if seg.depth >= threshold {
                let s = seg.start.max(start);
                let e = seg.end.min(end);
                if s < e {
                    count += e - s;
                }
            }
        }
        count
    }
}
