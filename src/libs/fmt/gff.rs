//! GFF reading via noodles, shared by the `gff` subcommands.

use anyhow::Context;
use std::collections::BTreeMap;
use std::io::BufRead;

/// Core fields of one GFF record, with attributes flattened to strings
/// (array values are joined with `,`).
#[derive(Debug, Clone)]
pub struct GffRecord {
    pub seqid: String,
    pub ty: String,
    pub strand: char,
    pub start: u64,
    pub end: u64,
    pub attributes: BTreeMap<String, String>,
}

/// Read all GFF records from `reader` with the noodles parser.
pub fn read_records<R: BufRead>(reader: R) -> anyhow::Result<Vec<GffRecord>> {
    let mut gff = noodles_gff::io::Reader::new(reader);
    let mut out = Vec::new();
    for result in gff.record_bufs() {
        let record = result.context("parsing GFF record")?;
        let strand = match record.strand() {
            noodles_gff::feature::record::Strand::Reverse => '-',
            _ => '+',
        };
        let attributes = record
            .attributes()
            .as_ref()
            .iter()
            .map(|(tag, value)| {
                let key = String::from_utf8_lossy(tag.as_ref()).into_owned();
                let value = match value {
                    noodles_gff::feature::record_buf::attributes::field::Value::String(s) => {
                        String::from_utf8_lossy(s).into_owned()
                    }
                    noodles_gff::feature::record_buf::attributes::field::Value::Array(vs) => vs
                        .iter()
                        .map(|s| String::from_utf8_lossy(s).into_owned())
                        .collect::<Vec<_>>()
                        .join(","),
                };
                (key, value)
            })
            .collect();
        out.push(GffRecord {
            seqid: record.reference_sequence_name().to_string(),
            ty: record.ty().to_string(),
            strand,
            start: usize::from(record.start()) as u64,
            end: usize::from(record.end()) as u64,
            attributes,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_basic_gff() {
        let gff = "##gff-version 3\n\
                   chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;Name=gene one\n\
                   chr1\tsrc\tmRNA\t50\t80\t.\t-\t.\tID=m1\n";
        let recs = read_records(std::io::Cursor::new(gff)).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].seqid, "chr1");
        assert_eq!(recs[0].ty, "gene");
        assert_eq!((recs[0].start, recs[0].end), (1, 100));
        assert_eq!(recs[0].strand, '+');
        assert_eq!(recs[0].attributes["ID"], "g1");
        assert_eq!(recs[0].attributes["Name"], "gene one");
        assert_eq!(recs[1].strand, '-');
    }
}
