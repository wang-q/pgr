//! Block-cached BGZF reader for random access.

use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::Context;
use flate2::{Crc, Decompress, FlushDecompress};
use lru::LruCache;
use noodles_bgzf as bgzf;

const HEADER_SIZE: usize = 12;
const TRAILER_SIZE: usize = 8;
const MAX_ISIZE: usize = 1 << 16;

/// Decompresses a raw DEFLATE stream (one BGZF block) into `dst`.
pub trait BlockInflater {
    /// Inflates `cdata` into `dst`, which must be large enough for the output.
    fn inflate(&mut self, cdata: &[u8], dst: &mut [u8]) -> io::Result<()>;
}

/// Default inflater backed by `flate2` (zlib-rs), reused across blocks.
pub struct Flate2Inflater {
    decompress: Decompress,
}

impl Default for Flate2Inflater {
    fn default() -> Self {
        Self {
            decompress: Decompress::new(false),
        }
    }
}

impl BlockInflater for Flate2Inflater {
    fn inflate(&mut self, cdata: &[u8], dst: &mut [u8]) -> io::Result<()> {
        self.decompress.reset(false);
        let status = self
            .decompress
            .decompress(cdata, dst, FlushDecompress::Finish)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if status != flate2::Status::StreamEnd {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "incomplete inflate",
            ));
        }
        Ok(())
    }
}

/// libdeflate-backed inflater (C), the default block inflater.
pub struct LibdeflaterInflater {
    decompressor: libdeflater::Decompressor,
}

impl Default for LibdeflaterInflater {
    fn default() -> Self {
        Self {
            decompressor: libdeflater::Decompressor::new(),
        }
    }
}

impl BlockInflater for LibdeflaterInflater {
    fn inflate(&mut self, cdata: &[u8], dst: &mut [u8]) -> io::Result<()> {
        let n = self
            .decompressor
            .deflate_decompress(cdata, dst)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if n != dst.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unexpected output size",
            ));
        }
        Ok(())
    }
}

/// A decompressed BGZF block plus its compressed size (for advancing to the
/// next block during sequential reads).
struct CachedBlock {
    data: Box<[u8]>,
    size: u64,
}

/// A BGZF reader with an LRU cache of decompressed blocks.
///
/// Repeated seeks into the same block hit the cache instead of re-reading and
/// re-inflating the block. The inflater instance is reused across blocks.
pub struct CachedBgzfReader {
    file: File,
    index: bgzf::gzi::Index,
    cache: LruCache<u64, CachedBlock>,
    inflater: Box<dyn BlockInflater>,
    /// Current block: (compressed offset, compressed size, uncompressed cursor, uncompressed len).
    current: Option<(u64, u64, usize, usize)>,
}

impl CachedBgzfReader {
    /// Opens a BGZF file (with its sibling `.gzi` index) and a block cache.
    pub fn open(path: impl AsRef<Path>, capacity: NonZeroUsize) -> anyhow::Result<Self> {
        Self::open_with_inflater(path, capacity, Box::<LibdeflaterInflater>::default())
    }

    /// Opens a BGZF file with a custom block inflater (for benchmarking).
    pub fn open_with_inflater(
        path: impl AsRef<Path>,
        capacity: NonZeroUsize,
        inflater: Box<dyn BlockInflater>,
    ) -> anyhow::Result<Self> {
        let path = path.as_ref();
        let index_path = format!("{}.gzi", path.display());
        let index = bgzf::gzi::fs::read(&index_path)
            .with_context(|| format!("failed to read index {}", index_path))?;
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        Ok(Self {
            file,
            index,
            cache: LruCache::new(capacity),
            inflater,
            current: None,
        })
    }

    /// Ensures the block at `cpos` is in the cache, returning its size and
    /// uncompressed length.
    fn load_block(&mut self, cpos: u64) -> io::Result<(u64, usize)> {
        if let Some(block) = self.cache.get(&cpos) {
            return Ok((block.size, block.data.len()));
        }
        let block = read_block(&mut self.file, cpos, self.inflater.as_mut())?;
        let (size, len) = (block.size, block.data.len());
        if !block.data.is_empty() {
            self.cache.put(cpos, block);
        }
        Ok((size, len))
    }
}

impl Read for CachedBgzfReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            let Some((cpos, bsize, upos, blen)) = self.current else {
                return Ok(0);
            };
            if upos >= blen {
                let next_cpos = cpos + bsize;
                match self.load_block(next_cpos) {
                    Ok((next_size, next_len)) => {
                        self.current = Some((next_cpos, next_size, 0, next_len));
                        continue;
                    }
                    Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => {
                        self.current = None;
                        return Ok(0);
                    }
                    Err(e) => return Err(e),
                }
            }
            let n = (blen - upos).min(buf.len());
            {
                let block = self
                    .cache
                    .get(&cpos)
                    .ok_or_else(|| io::Error::other("block not cached"))?;
                buf[..n].copy_from_slice(&block.data[upos..upos + n]);
            }
            let (_, _, u, _) = self
                .current
                .as_mut()
                .ok_or_else(|| io::Error::other("no current block"))?;
            *u += n;
            return Ok(n);
        }
    }
}

impl Seek for CachedBgzfReader {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let SeekFrom::Start(pos) = pos else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "only SeekFrom::Start is supported",
            ));
        };
        let vpos = self.index.query(pos)?;
        let (cpos, upos): (u64, u16) = vpos.into();
        let (bsize, blen) = self.load_block(cpos)?;
        self.current = Some((cpos, bsize, usize::from(upos), blen));
        Ok(pos)
    }
}

fn read_block(
    file: &mut File,
    cpos: u64,
    inflater: &mut dyn BlockInflater,
) -> io::Result<CachedBlock> {
    file.seek(SeekFrom::Start(cpos))?;

    let mut header = [0u8; HEADER_SIZE];
    file.read_exact(&mut header)?;
    if header[0] != 0x1f || header[1] != 0x8b || header[3] & 4 == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a BGZF block",
        ));
    }

    let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;
    let mut extra = vec![0u8; xlen];
    file.read_exact(&mut extra)?;

    let mut bsize = 0u16;
    let mut found_bc = false;
    let mut cursor = 0;
    while cursor + 4 <= extra.len() {
        let si1 = extra[cursor];
        let si2 = extra[cursor + 1];
        let slen = u16::from_le_bytes([extra[cursor + 2], extra[cursor + 3]]) as usize;
        if si1 == b'B' && si2 == b'C' && slen == 2 {
            if cursor + 6 <= extra.len() {
                bsize = u16::from_le_bytes([extra[cursor + 4], extra[cursor + 5]]);
                found_bc = true;
            }
            break;
        }
        cursor += 4 + slen;
    }
    if !found_bc {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing BC subfield",
        ));
    }

    let block_size = u64::from(bsize) + 1;
    let data_start = (HEADER_SIZE + xlen) as u64;
    if data_start + TRAILER_SIZE as u64 > block_size {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "malformed block size",
        ));
    }

    let cdata_len = (block_size - data_start - TRAILER_SIZE as u64) as usize;
    let mut cdata = vec![0u8; cdata_len];
    file.read_exact(&mut cdata)?;
    let mut trailer = [0u8; TRAILER_SIZE];
    file.read_exact(&mut trailer)?;

    let isize = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]) as usize;
    if isize > MAX_ISIZE {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid ISIZE"));
    }

    if isize == 0 {
        // Empty block (e.g. the BGZF EOF marker): skip it, nothing to decode.
        return Ok(CachedBlock {
            data: Box::new([]),
            size: block_size,
        });
    }

    let mut out = vec![0u8; isize];
    inflater.inflate(&cdata, &mut out)?;

    let expected_crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
    let mut crc = Crc::new();
    crc.update(&out);
    if crc.sum() != expected_crc {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "CRC32 mismatch"));
    }

    Ok(CachedBlock {
        data: out.into_boxed_slice(),
        size: block_size,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn make_bgzf(data: &[u8], tag: &str) -> std::path::PathBuf {
        let path =
            std::env::temp_dir().join(format!("pgr_bgzf_test_{}_{}.gz", std::process::id(), tag));
        let mut writer = bgzf::io::Writer::new(File::create(&path).unwrap());
        writer.write_all(data).unwrap();
        writer.finish().unwrap();
        crate::libs::fmt::fa::build_gzi_index(path.to_str().unwrap()).unwrap();
        path
    }

    #[test]
    fn random_access_matches_plain() {
        let data: Vec<u8> = (0..300_000u32).map(|i| (i % 251) as u8).collect();
        let path = make_bgzf(&data, "random_access");
        let capacity = NonZeroUsize::new(4).expect("non-zero");
        let mut reader = CachedBgzfReader::open(&path, capacity).unwrap();

        let probes = [
            (0u64, 100usize),
            (65_000, 200),
            (65_536, 300),
            (131_000, 150),
            (299_800, 200),
            (50_000, 5000),
        ];
        for &(off, len) in &probes {
            let mut got = vec![0u8; len];
            reader.seek(SeekFrom::Start(off)).unwrap();
            reader.read_exact(&mut got).unwrap();
            assert_eq!(got, data[off as usize..off as usize + len], "offset {off}");
        }

        // Sequential read crossing block boundaries after a seek.
        let mut got = Vec::new();
        reader.seek(SeekFrom::Start(0)).unwrap();
        reader.read_to_end(&mut got).unwrap();
        assert_eq!(got, data);

        std::fs::remove_file(&path).ok();
        std::fs::remove_file(format!("{}.gzi", path.display())).ok();
    }

    #[test]
    fn cache_reuses_blocks() {
        let data: Vec<u8> = (0..200_000u32).map(|i| (i % 251) as u8).collect();
        let path = make_bgzf(&data, "cache_reuse");
        let capacity = NonZeroUsize::new(1).expect("non-zero");
        let mut reader = CachedBgzfReader::open(&path, capacity).unwrap();
        let mut got = vec![0u8; 50];
        reader.seek(SeekFrom::Start(10_000)).unwrap();
        reader.read_exact(&mut got).unwrap();
        assert_eq!(got, data[10_000..10_050]);
        // Same block again: still correct after a fresh seek.
        reader.seek(SeekFrom::Start(10_100)).unwrap();
        reader.read_exact(&mut got).unwrap();
        assert_eq!(got, data[10_100..10_150]);

        std::fs::remove_file(&path).ok();
        std::fs::remove_file(format!("{}.gzi", path.display())).ok();
    }
}
