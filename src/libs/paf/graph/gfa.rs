//! GFA v1.0 emission for [`super::PafGraph`].

use super::PafGraph;
use std::io::Write;

impl PafGraph {
    /// Write GFA v1.0 (rGFA with SN/SO/SR tags on S lines) to a writer.
    pub fn write_gfa<W: Write>(&self, mut w: W) -> std::io::Result<()> {
        // S lines (1-based node ids in GFA convention) with rGFA SN/SO/SR tags.
        // SN: source sequence name; SO: 0-based start offset; SR: rank (0 = primary).
        for (i, seq) in self.node_seqs.iter().enumerate() {
            let id = (i + 1) as u32;
            let len = self.node_lens[i];
            let (sn, so) = &self.node_origins[i];
            if seq.is_empty() {
                // Topology-only mode: emit '*' with LN:i: tag.
                if sn.is_empty() {
                    writeln!(w, "S\t{id}\t*\tLN:i:{len}")?;
                } else {
                    writeln!(w, "S\t{id}\t*\tLN:i:{len}\tSN:Z:{sn}\tSO:i:{so}\tSR:i:0")?;
                }
            } else {
                let s = String::from_utf8_lossy(seq);
                if sn.is_empty() {
                    writeln!(w, "S\t{id}\t{s}")?;
                } else {
                    writeln!(w, "S\t{id}\t{s}\tSN:Z:{sn}\tSO:i:{so}\tSR:i:0")?;
                }
            }
        }
        // L lines.
        for e in &self.edges {
            writeln!(
                w,
                "L\t{}\t{}\t{}\t{}\t0M",
                e.from + 1,
                e.from_orient,
                e.to + 1,
                e.to_orient
            )?;
        }
        // P lines.
        for (name, steps) in &self.paths {
            let path_str: Vec<String> = steps
                .iter()
                .map(|s| format!("{}{}", s.node + 1, s.orient))
                .collect();
            let overlaps = vec!["0M"; steps.len().saturating_sub(1)];
            writeln!(
                w,
                "P\t{name}\t{}\t{}",
                path_str.join(","),
                overlaps.join(",")
            )?;
        }
        w.flush()?;
        Ok(())
    }
}
