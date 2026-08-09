//! Streaming multi-member gzip reader (replaces flate2's MultiGzDecoder).

use std::io::{self, Read};

use std::collections::BTreeMap;
use std::fs::File;
use std::thread::{self, JoinHandle};

use crc32fast::Hasher;
use crossbeam::channel::{bounded, unbounded, Receiver};
use libz_ng_sys::{Z_BUF_ERROR, Z_DATA_ERROR, Z_NO_FLUSH, Z_OK, Z_STREAM_END, Z_STREAM_ERROR};

use crate::libs::bgzf::{BlockInflater, LibdeflaterInflater, MAX_ISIZE};

const OUT_BUF_SIZE: usize = 1 << 16;
const TRAILER_SIZE: usize = 8;

// zlib-ng's z_stream must be zero-initialized (NULL zalloc/zfree select the
// default allocator); the struct contains `extern "C" fn` fields, which trip
// the invalid_value lint even though null fn pointers are valid.
#[allow(invalid_value)]
fn zeroed_z_stream() -> libz_ng_sys::z_stream {
    unsafe { std::mem::MaybeUninit::<libz_ng_sys::z_stream>::zeroed().assume_init() }
}

/// A streaming reader over gzip data with any number of concatenated members
/// (ordinary gzip or BGZF), inflating with zlib-ng.
pub struct GzReader<R: Read> {
    inner: R,
    // Box keeps the z_stream address stable: zlib's inflate_state records the
    // address of the stream at init time and rejects a moved one.
    strm: Box<libz_ng_sys::z_stream>,
    in_buf: Vec<u8>,
    in_pos: usize,
    out_buf: Vec<u8>,
    out_pos: usize,
    out_len: usize,
    eof: bool,
    need_header: bool,
    need_trailer: bool,
    crc: Hasher,
    member_len: u64,
}

// z_stream contains raw pointers, but the inner state is heap-allocated and
// next_in/next_out are re-pointed at self buffers on every pump, so moving the
// reader across threads is safe.
unsafe impl<R: Read> Send for GzReader<R> {}

impl<R: Read> GzReader<R> {
    /// Creates a gzip stream reader over `inner`.
    pub fn new(inner: R) -> io::Result<Self> {
        let mut strm = Box::new(zeroed_z_stream());
        let rc = unsafe {
            libz_ng_sys::inflateInit2_(
                strm.as_mut(),
                -15, // raw deflate; gzip header/trailer handled here
                libz_ng_sys::zlibVersion(),
                std::mem::size_of::<libz_ng_sys::z_stream>() as i32,
            )
        };
        if rc != Z_OK {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "inflateInit2 failed",
            ));
        }
        Ok(Self {
            inner,
            strm,
            in_buf: Vec::new(),
            in_pos: 0,
            out_buf: vec![0; OUT_BUF_SIZE],
            out_pos: 0,
            out_len: 0,
            eof: false,
            need_header: true,
            need_trailer: false,
            crc: Hasher::new(),
            member_len: 0,
        })
    }

    fn fill_more(&mut self) -> io::Result<bool> {
        let mut tmp = [0u8; 1 << 16];
        let n = self.inner.read(&mut tmp)?;
        if n == 0 {
            return Ok(false);
        }
        self.in_buf.extend_from_slice(&tmp[..n]);
        Ok(true)
    }

    fn pending(&self) -> &[u8] {
        &self.in_buf[self.in_pos..]
    }

    fn consume(&mut self, n: usize) {
        self.in_pos += n;
        if self.in_pos == self.in_buf.len() {
            self.in_buf.clear();
            self.in_pos = 0;
        }
    }

    fn parse_header(&mut self) -> io::Result<bool> {
        loop {
            if let Some(header) = crate::libs::bgzf::parse_gz_header(self.pending())? {
                self.consume(header.header_len);
                self.crc = Hasher::new();
                self.member_len = 0;
                self.need_header = false;
                return Ok(true);
            }
            // Need more input for the header.
            if self.in_pos > 0 {
                self.in_buf.drain(..self.in_pos);
                self.in_pos = 0;
            }
            if !self.fill_more()? {
                if self.pending().is_empty() {
                    return Ok(false);
                }
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated gzip header",
                ));
            }
        }
    }

    fn read_trailer(&mut self) -> io::Result<bool> {
        while self.pending().len() < TRAILER_SIZE {
            if self.in_pos > 0 {
                self.in_buf.drain(..self.in_pos);
                self.in_pos = 0;
            }
            if !self.fill_more()? {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated gzip trailer",
                ));
            }
        }
        let trailer = &self.pending()[..TRAILER_SIZE];
        let expected_crc = u32::from_le_bytes(trailer[..4].try_into().unwrap());
        let expected_len = u32::from_le_bytes(trailer[4..].try_into().unwrap());
        if self.crc.clone().finalize() != expected_crc || self.member_len != u64::from(expected_len)
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "gzip CRC/ISIZE mismatch",
            ));
        }
        self.consume(TRAILER_SIZE);
        self.need_trailer = false;
        self.need_header = true;
        Ok(true)
    }

    fn pump(&mut self) -> io::Result<bool> {
        loop {
            if self.need_header {
                if !self.parse_header()? {
                    self.eof = true;
                    return Ok(false);
                }
                continue;
            }
            if self.need_trailer {
                self.read_trailer()?;
                continue;
            }
            if self.pending().is_empty() && !self.fill_more()? {
                self.eof = true;
                return Ok(false);
            }
            self.out_buf.resize(OUT_BUF_SIZE, 0);
            let out_cap = self.out_buf.len();
            self.strm.next_in = self.pending().as_ptr().cast_mut();
            self.strm.avail_in = self.pending().len() as u32;
            self.strm.next_out = self.out_buf.as_mut_ptr();
            self.strm.avail_out = out_cap as u32;
            let rc = unsafe { libz_ng_sys::inflate(self.strm.as_mut(), Z_NO_FLUSH) };
            let consumed = self.pending().len() - self.strm.avail_in as usize;
            self.consume(consumed);
            self.out_len = out_cap - self.strm.avail_out as usize;
            self.out_pos = 0;
            self.crc.update(&self.out_buf[..self.out_len]);
            self.member_len += self.out_len as u64;

            match rc {
                Z_STREAM_END => {
                    let reset = unsafe { libz_ng_sys::inflateReset(self.strm.as_mut()) };
                    if reset != Z_OK {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "inflateReset failed",
                        ));
                    }
                    self.need_trailer = true;
                    return Ok(self.out_len > 0);
                }
                Z_OK | Z_BUF_ERROR => {
                    if self.out_len > 0 {
                        return Ok(true);
                    }
                    if self.pending().is_empty() {
                        // No progress possible without more input.
                        continue;
                    }
                    continue;
                }
                Z_DATA_ERROR | Z_STREAM_ERROR => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "inflate failed: rc={rc} avail_in={} avail_out={}",
                            self.strm.avail_in, self.strm.avail_out
                        ),
                    ));
                }
                other => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("inflate error {other} avail_out={}", self.strm.avail_out),
                    ));
                }
            }
        }
    }
}

impl<R: Read> Read for GzReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            if self.out_pos < self.out_len {
                let n = (self.out_len - self.out_pos).min(buf.len());
                buf[..n].copy_from_slice(&self.out_buf[self.out_pos..self.out_pos + n]);
                self.out_pos += n;
                return Ok(n);
            }
            if self.eof {
                return Ok(0);
            }
            if !self.pump()? {
                return Ok(0);
            }
        }
    }
}

impl<R: Read> Drop for GzReader<R> {
    fn drop(&mut self) {
        unsafe {
            libz_ng_sys::inflateEnd(self.strm.as_mut());
        }
    }
}

struct BlockJob {
    seq: u64,
    cdata: Vec<u8>,
    isize: usize,
    crc: u32,
}

/// Parses one BGZF block at `cpos`, returning its payload and metadata.
///
/// Parsed block payload: (compressed data, uncompressed size, CRC32, block size).
type ParsedBlock = (Vec<u8>, usize, u32, u64);

/// Returns `Ok(None)` at the EOF marker or end of file.
fn parse_block(file: &mut File, cpos: u64) -> io::Result<Option<ParsedBlock>> {
    use std::io::{Seek, SeekFrom};
    file.seek(SeekFrom::Start(cpos))?;
    let mut header = [0u8; 12];
    match file.read_exact(&mut header) {
        Ok(_) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }
    if header[0] != 0x1f || header[1] != 0x8b || header[3] & 4 == 0 {
        return Ok(None);
    }
    let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;
    let mut extra = vec![0u8; xlen];
    file.read_exact(&mut extra)?;
    let mut bsize = 0u16;
    let mut found = false;
    let mut cursor = 0;
    while cursor + 4 <= extra.len() {
        let si1 = extra[cursor];
        let si2 = extra[cursor + 1];
        let slen = u16::from_le_bytes([extra[cursor + 2], extra[cursor + 3]]) as usize;
        if si1 == b'B' && si2 == b'C' && slen == 2 && cursor + 6 <= extra.len() {
            bsize = u16::from_le_bytes([extra[cursor + 4], extra[cursor + 5]]);
            found = true;
            break;
        }
        cursor += 4 + slen;
    }
    if !found {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing BC subfield",
        ));
    }
    let block_size = u64::from(bsize) + 1;
    let data_start = (12 + xlen) as u64;
    if data_start + 8 > block_size {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "malformed block size",
        ));
    }
    let cdata_len = (block_size - data_start - 8) as usize;
    let mut cdata = vec![0u8; cdata_len];
    file.read_exact(&mut cdata)?;
    let mut trailer = [0u8; 8];
    file.read_exact(&mut trailer)?;
    let isize = u32::from_le_bytes(trailer[4..].try_into().unwrap()) as usize;
    if isize == 0 || isize > MAX_ISIZE {
        return Ok(None); // empty block / EOF marker
    }
    let crc = u32::from_le_bytes(trailer[..4].try_into().unwrap());
    Ok(Some((cdata, isize, crc, block_size)))
}

/// A BGZF reader that inflates blocks on a worker pool, delivering output in
/// file order.
pub struct ParallelBgzfReader {
    rx: Receiver<(u64, Result<Vec<u8>, String>)>,
    handles: Vec<JoinHandle<()>>,
    results: BTreeMap<u64, Vec<u8>>,
    next_seq: u64,
    cur: Option<(usize, Vec<u8>)>,
    eof: bool,
}

impl ParallelBgzfReader {
    /// Opens a BGZF file with `worker_count` inflate threads.
    pub fn open(path: impl AsRef<std::path::Path>, worker_count: usize) -> anyhow::Result<Self> {
        let file = File::open(path.as_ref())?;
        let (job_tx, job_rx) = bounded::<BlockJob>(worker_count * 2);
        let (res_tx, res_rx) = unbounded::<(u64, Result<Vec<u8>, String>)>();

        let mut handles = Vec::with_capacity(worker_count + 1);
        {
            let res_tx = res_tx.clone();
            handles.push(thread::spawn(move || {
                let mut file = file;
                let mut cpos = 0u64;
                let mut seq = 0u64;
                loop {
                    match parse_block(&mut file, cpos) {
                        Ok(Some((cdata, isize, crc, block_size))) => {
                            if job_tx
                                .send(BlockJob {
                                    seq,
                                    cdata,
                                    isize,
                                    crc,
                                })
                                .is_err()
                            {
                                break;
                            }
                            seq += 1;
                            cpos += block_size;
                        }
                        Ok(None) => break,
                        Err(e) => {
                            let _ = res_tx.send((seq, Err(e.to_string())));
                            break;
                        }
                    }
                }
            }));
        }
        for _ in 0..worker_count {
            let job_rx = job_rx.clone();
            let res_tx = res_tx.clone();
            handles.push(thread::spawn(move || {
                let mut inflater = LibdeflaterInflater::default();
                for job in job_rx {
                    let mut out = vec![0u8; job.isize];
                    let result = inflater.inflate(&job.cdata, &mut out).and_then(|()| {
                        let mut crc = Hasher::new();
                        crc.update(&out);
                        if crc.finalize() != job.crc {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "CRC32 mismatch",
                            ));
                        }
                        Ok(out)
                    });
                    let _ = res_tx.send((job.seq, result.map_err(|e| e.to_string())));
                }
            }));
        }
        Ok(Self {
            rx: res_rx,
            handles,
            results: BTreeMap::new(),
            next_seq: 0,
            cur: None,
            eof: false,
        })
    }
}

impl Read for ParallelBgzfReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            if let Some((pos, data)) = &mut self.cur {
                if *pos < data.len() {
                    let n = (data.len() - *pos).min(buf.len());
                    buf[..n].copy_from_slice(&data[*pos..*pos + n]);
                    *pos += n;
                    return Ok(n);
                }
                self.cur = None;
                self.next_seq += 1;
                continue;
            }
            if self.eof {
                return Ok(0);
            }
            if let Some(data) = self.results.remove(&self.next_seq) {
                self.cur = Some((0, data));
                continue;
            }
            match self.rx.recv() {
                Ok((seq, Ok(data))) => {
                    if seq == self.next_seq {
                        self.cur = Some((0, data));
                    } else {
                        self.results.insert(seq, data);
                    }
                }
                Ok((_, Err(msg))) => {
                    return Err(io::Error::new(io::ErrorKind::InvalidData, msg));
                }
                Err(_) => {
                    self.eof = true;
                    return Ok(0);
                }
            }
        }
    }
}

impl Drop for ParallelBgzfReader {
    fn drop(&mut self) {
        for handle in self.handles.drain(..) {
            let _ = handle.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn gz_member(data: &[u8]) -> Vec<u8> {
        crate::libs::bgzf::gzip_compress(data, 6).unwrap()
    }

    #[test]
    fn single_member() {
        let data = b"hello world\n".repeat(1000);
        let compressed = gz_member(&data);
        let mut out = Vec::new();
        GzReader::new(compressed.as_slice())
            .unwrap()
            .read_to_end(&mut out)
            .unwrap();
        assert_eq!(out, data);
    }

    #[test]
    fn multi_member() {
        let a = b"first member ".repeat(500);
        let b = b"second member ".repeat(700);
        let mut compressed = gz_member(&a);
        compressed.extend_from_slice(&gz_member(&b));
        let mut out = Vec::new();
        GzReader::new(compressed.as_slice())
            .unwrap()
            .read_to_end(&mut out)
            .unwrap();
        let mut expected = a;
        expected.extend_from_slice(&b);
        assert_eq!(out, expected);
    }

    #[test]
    fn partial_reads() {
        let data = b"partial reads ".repeat(3000);
        let compressed = gz_member(&data);
        let mut reader = GzReader::new(compressed.as_slice()).unwrap();
        let mut out = Vec::new();
        let mut chunk = vec![0u8; 4093];
        loop {
            let n = reader.read(&mut chunk).unwrap();
            if n == 0 {
                break;
            }
            out.extend_from_slice(&chunk[..n]);
        }
        assert_eq!(out, data);
    }

    #[test]
    fn parallel_reader_matches_plain() {
        let data: Vec<u8> = (0..500_000u32).map(|i| (i % 251) as u8).collect();
        let path =
            std::env::temp_dir().join(format!("pgr_bgzf_parallel_{}.gz", std::process::id()));
        {
            let mut w =
                crate::libs::bgzf::BgzfWriter::new(std::fs::File::create(&path).unwrap()).unwrap();
            w.write_all(&data).unwrap();
            w.finish().unwrap();
        }
        for workers in [1usize, 2, 4] {
            let mut reader = ParallelBgzfReader::open(&path, workers).unwrap();
            let mut out = Vec::new();
            reader.read_to_end(&mut out).unwrap();
            assert_eq!(out, data, "workers={workers}");
        }
        std::fs::remove_file(&path).ok();
    }
}
