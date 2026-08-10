//! Streaming pair reader for FASTQ inputs.

use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use anyhow::Result;

/// Yields `(R1, R2)` pairs from one interleaved file or two files (R1, R2).
///
/// Owned and `Send`, so it can feed parallel pipelines; reads are produced
/// on demand rather than loaded wholesale.
pub struct PairReader {
    reader1: SeqReader<'static>,
    reader2: Option<SeqReader<'static>>,
    rec1: SeqRecord,
    rec2: SeqRecord,
}

impl PairReader {
    /// Opens the input stream(s).
    pub fn new(infiles: &[String]) -> Result<Self> {
        let reader1 = SeqReader::new(&infiles[0])?;
        let reader2 = if infiles.len() == 2 {
            Some(SeqReader::new(&infiles[1])?)
        } else {
            None
        };
        Ok(Self {
            reader1,
            reader2,
            rec1: SeqRecord::new(),
            rec2: SeqRecord::new(),
        })
    }
}

impl Iterator for PairReader {
    type Item = Result<(SeqRecord, Option<SeqRecord>)>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader1.read_record(&mut self.rec1) {
            Err(e) => return Some(Err(e)),
            Ok(false) => return None,
            Ok(true) => {}
        }
        let has2 = match self.reader2.as_mut() {
            Some(r2) => match r2.read_record(&mut self.rec2) {
                Err(e) => return Some(Err(e)),
                Ok(b) => b,
            },
            None => match self.reader1.read_record(&mut self.rec2) {
                Err(e) => return Some(Err(e)),
                Ok(b) => b,
            },
        };
        Some(Ok((self.rec1.clone(), has2.then(|| self.rec2.clone()))))
    }
}
