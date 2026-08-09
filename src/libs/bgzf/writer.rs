//! BGZF writers: single-threaded and parallel block compression.

use std::collections::BTreeMap;
use std::io::{self, Write};
use std::thread::{self, JoinHandle};

use crc32fast::Hasher;
use crossbeam::channel::{bounded, unbounded, Receiver, Sender};
use libdeflater::{CompressionLvl, Compressor};

use crate::libs::bgzf::{BGZF_HEADER_SIZE, MAX_ISIZE, TRAILER_SIZE};

/// Maximum uncompressed input per block: leaves room for a stored (uncompressed)
/// fallback within the 64 KiB block limit (htslib uses the same 0xff00 bound).
const MAX_BLOCK_INPUT: usize = 0xff00;

/// Compression level used by default.
const DEFAULT_LEVEL: i32 = 6;

fn compress_block(
    compressor: &mut Compressor,
    stored: &mut Compressor,
    data: &[u8],
    out: &mut Vec<u8>,
) -> io::Result<()> {
    out.resize(compressor.deflate_compress_bound(data.len()), 0);
    let n = compressor
        .deflate_compress(data, out)
        .map_err(io::Error::other)?;
    out.truncate(n);

    // Fall back to a stored block when the compressed block would exceed the
    // 64 KiB limit (defeats compression, but keeps the file valid BGZF).
    if BGZF_HEADER_SIZE + out.len() + TRAILER_SIZE > MAX_ISIZE {
        out.resize(stored.deflate_compress_bound(data.len()), 0);
        let n = stored
            .deflate_compress(data, out)
            .map_err(io::Error::other)?;
        out.truncate(n);
    }
    Ok(())
}

fn write_bgzf_block<W: Write>(w: &mut W, data: &[u8], cdata: &[u8]) -> io::Result<()> {
    let block_size = BGZF_HEADER_SIZE + cdata.len() + TRAILER_SIZE;
    debug_assert!(block_size <= MAX_ISIZE);
    let bsize = (block_size - 1) as u16;

    let mut header = [0u8; BGZF_HEADER_SIZE];
    header[..2].copy_from_slice(&[0x1f, 0x8b]);
    header[2] = 8; // CM = deflate
    header[3] = 4; // FLG = FEXTRA
    header[9] = 255; // OS
    header[10..12].copy_from_slice(&6u16.to_le_bytes()); // XLEN
    header[12..16].copy_from_slice(&[b'B', b'C', 2, 0]); // BC subfield
    header[16..18].copy_from_slice(&bsize.to_le_bytes());

    w.write_all(&header)?;
    w.write_all(cdata)?;

    let mut crc = Hasher::new();
    crc.update(data);
    let mut trailer = [0u8; TRAILER_SIZE];
    trailer[..4].copy_from_slice(&crc.finalize().to_le_bytes());
    trailer[4..8].copy_from_slice(&(data.len() as u32).to_le_bytes());
    w.write_all(&trailer)?;
    Ok(())
}

/// Writes the 28-byte BGZF EOF marker.
fn write_eof<W: Write>(w: &mut W) -> io::Result<()> {
    const EOF_BLOCK: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
    w.write_all(&EOF_BLOCK)
}

/// A single-threaded BGZF writer.
pub struct BgzfWriter<W: Write> {
    inner: W,
    buf: Vec<u8>,
    compressor: Compressor,
    stored: Compressor,
}

impl<W: Write> BgzfWriter<W> {
    /// Creates a writer over `inner` at the default compression level.
    pub fn new(inner: W) -> io::Result<Self> {
        Self::with_level(inner, DEFAULT_LEVEL)
    }

    /// Creates a writer with a specific compression level.
    pub fn with_level(inner: W, level: i32) -> io::Result<Self> {
        let level = CompressionLvl::new(level)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("{e:?}")))?;
        let stored_level = CompressionLvl::new(0).expect("level 0 valid");
        Ok(Self {
            inner,
            buf: Vec::with_capacity(MAX_BLOCK_INPUT),
            compressor: Compressor::new(level),
            stored: Compressor::new(stored_level),
        })
    }

    /// Flushes pending data as a partial block and returns the inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        if !self.buf.is_empty() {
            let data = std::mem::take(&mut self.buf);
            let mut cdata = Vec::new();
            compress_block(&mut self.compressor, &mut self.stored, &data, &mut cdata)?;
            write_bgzf_block(&mut self.inner, &data, &cdata)?;
        }
        write_eof(&mut self.inner)?;
        Ok(self.inner)
    }

    /// Access to the underlying writer.
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut rest = buf;
        while !rest.is_empty() {
            let room = MAX_BLOCK_INPUT - self.buf.len();
            let take = rest.len().min(room);
            self.buf.extend_from_slice(&rest[..take]);
            rest = &rest[take..];
            if self.buf.len() == MAX_BLOCK_INPUT {
                let data = std::mem::take(&mut self.buf);
                let mut cdata = Vec::new();
                compress_block(&mut self.compressor, &mut self.stored, &data, &mut cdata)?;
                write_bgzf_block(&mut self.inner, &data, &cdata)?;
            }
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        // BGZF blocks are independent; nothing to flush between blocks.
        Ok(())
    }
}

struct BlockJob {
    seq: u64,
    data: Vec<u8>,
}

/// A BGZF writer that compresses blocks on a worker pool while preserving
/// output order.
pub struct ParallelBgzfWriter<W: Write + Send + 'static> {
    inner: Option<W>,
    buf: Vec<u8>,
    tx: Sender<BlockJob>,
    rx: Receiver<(u64, Vec<u8>)>,
    workers: Vec<JoinHandle<()>>,
    next_seq: u64,
    next_out: u64,
    results: BTreeMap<u64, Vec<u8>>,
}

impl<W: Write + Send + 'static> ParallelBgzfWriter<W> {
    /// Creates a writer with `worker_count` compression threads.
    pub fn new(inner: W, worker_count: usize) -> io::Result<Self> {
        Self::with_level(inner, worker_count, DEFAULT_LEVEL)
    }

    /// Creates a writer with a specific compression level.
    pub fn with_level(inner: W, worker_count: usize, level: i32) -> io::Result<Self> {
        let level = CompressionLvl::new(level)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("{e:?}")))?;
        let stored_level = CompressionLvl::new(0).expect("level 0 valid");
        let (job_tx, job_rx) = bounded::<BlockJob>(worker_count * 2);
        let (res_tx, res_rx) = unbounded::<(u64, Vec<u8>)>();

        let mut workers = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let job_rx = job_rx.clone();
            let res_tx = res_tx.clone();
            workers.push(thread::spawn(move || {
                let mut compressor = Compressor::new(level);
                let mut stored = Compressor::new(stored_level);
                for job in job_rx {
                    let mut cdata = Vec::new();
                    if compress_block(&mut compressor, &mut stored, &job.data, &mut cdata).is_ok() {
                        let mut block = Vec::new();
                        if write_bgzf_block(&mut block, &job.data, &cdata).is_ok() {
                            let _ = res_tx.send((job.seq, block));
                            continue;
                        }
                    }
                    let _ = res_tx.send((job.seq, Vec::new()));
                }
            }));
        }

        Ok(Self {
            inner: Some(inner),
            buf: Vec::with_capacity(MAX_BLOCK_INPUT),
            tx: job_tx,
            rx: res_rx,
            workers,
            next_seq: 0,
            next_out: 0,
            results: BTreeMap::new(),
        })
    }

    fn submit_pending(&mut self) -> io::Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let data = std::mem::take(&mut self.buf);
        self.tx
            .send(BlockJob {
                seq: self.next_seq,
                data,
            })
            .map_err(|_| io::Error::other("compression worker exited"))?;
        self.next_seq += 1;
        Ok(())
    }

    /// Flushes pending data and returns the inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        self.submit_pending()?;
        drop(self.tx);
        let mut writer = self.inner.take().expect("inner present");
        for handle in self.workers.drain(..) {
            let _ = handle.join();
        }
        while let Ok((seq, cdata)) = self.rx.recv() {
            self.results.insert(seq, cdata);
        }
        while let Some(cdata) = self.results.remove(&self.next_out) {
            writer.write_all(&cdata)?;
            self.next_out += 1;
        }
        if self.next_out != self.next_seq {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "missing compressed blocks",
            ));
        }
        write_eof(&mut writer)?;
        Ok(writer)
    }
}

impl<W: Write + Send + 'static> Write for ParallelBgzfWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut rest = buf;
        while !rest.is_empty() {
            let room = MAX_BLOCK_INPUT - self.buf.len();
            let take = rest.len().min(room);
            self.buf.extend_from_slice(&rest[..take]);
            rest = &rest[take..];
            if self.buf.len() == MAX_BLOCK_INPUT {
                self.submit_pending()?;
            }
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::bgzf::reader::GzReader;
    use std::io::Read;

    #[test]
    fn single_threaded_roundtrip() {
        let data: Vec<u8> = (0..300_000u32).map(|i| (i % 251) as u8).collect();
        let mut out = Vec::new();
        {
            let mut w = BgzfWriter::new(&mut out).unwrap();
            w.write_all(&data).unwrap();
            w.finish().unwrap();
        }
        let mut decoded = Vec::new();
        GzReader::new(out.as_slice())
            .unwrap()
            .read_to_end(&mut decoded)
            .unwrap();
        assert_eq!(decoded, data);
    }
}
