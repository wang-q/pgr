//! Self-contained BGZF/gzip support: block parsing, inflate/deflate,
//! random-access reads, and sequential/parallel readers and writers.

pub mod index;
pub mod reader;
pub mod writer;

pub use self::index::GziIndex;
pub use self::reader::{GzReader, ParallelBgzfReader};
pub use self::writer::{BgzfWriter, ParallelBgzfWriter};

use std::fs::File;
use std::io::{self, BufRead, Read, Seek, SeekFrom};
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::Context;
use crc32fast::Hasher;
use libdeflater::{CompressionLvl, Compressor, Decompressor};
use lru::LruCache;

/// Fixed gzip header size (before any FLG optional fields).
const GZIP_FIXED_HEADER: usize = 10;
/// BGZF extra field header: gzip fixed (10) + XLEN (2) + BC subfield (6).
pub(crate) const BGZF_HEADER_SIZE: usize = 18;
const TRAILER_SIZE: usize = 8;
/// Maximum uncompressed size of a BGZF block (2^16).
pub const MAX_ISIZE: usize = 1 << 16;

/// A BGZF virtual position: compressed offset (high 48 bits) + offset within
/// the uncompressed block (low 16 bits).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct VirtualPos(u64);

impl VirtualPos {
    /// The minimum virtual position (start of the file).
    pub const MIN: Self = Self(0);

    /// Creates a virtual position.
    pub fn new(cpos: u64, upos: u16) -> Self {
        Self((cpos << 16) | u64::from(upos))
    }

    /// Compressed offset of the block.
    pub fn compressed(self) -> u64 {
        self.0 >> 16
    }

    /// Offset within the uncompressed block.
    pub fn uncompressed(self) -> u16 {
        self.0 as u16
    }
}

impl From<u64> for VirtualPos {
    fn from(value: u64) -> Self {
        Self(value)
    }
}

impl From<VirtualPos> for u64 {
    fn from(pos: VirtualPos) -> Self {
        pos.0
    }
}

impl From<VirtualPos> for (u64, u16) {
    fn from(pos: VirtualPos) -> Self {
        (pos.compressed(), pos.uncompressed())
    }
}

impl Default for VirtualPos {
    fn default() -> Self {
        Self::MIN
    }
}

/// Parsed gzip member header.
#[derive(Debug)]
pub struct GzHeader {
    /// BGZF `BSIZE` (total block size - 1) when the `BC` subfield is present.
    pub bsize: Option<u16>,
    /// Total header length (fixed + optional fields).
    pub header_len: usize,
}

/// Parses a gzip member header from `buf`, which must contain the full header.
///
/// Handles FEXTRA (including the BGZF `BC` subfield), FNAME, FCOMMENT and
/// FHCRC. Returns `None` when `buf` is too short to contain the whole header.
pub fn parse_gz_header(buf: &[u8]) -> io::Result<Option<GzHeader>> {
    if buf.len() < GZIP_FIXED_HEADER {
        return Ok(None);
    }
    if buf[0] != 0x1f || buf[1] != 0x8b || buf[2] != 8 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a gzip stream",
        ));
    }

    let flg = buf[3];
    let mut cursor = GZIP_FIXED_HEADER;
    let mut bsize = None;

    if flg & 4 != 0 {
        // FEXTRA
        if buf.len() < cursor + 2 {
            return Ok(None);
        }
        let xlen = u16::from_le_bytes([buf[cursor], buf[cursor + 1]]) as usize;
        cursor += 2;
        if buf.len() < cursor + xlen {
            return Ok(None);
        }
        let extra = &buf[cursor..cursor + xlen];
        let mut e = 0;
        while e + 4 <= extra.len() {
            let si1 = extra[e];
            let si2 = extra[e + 1];
            let slen = u16::from_le_bytes([extra[e + 2], extra[e + 3]]) as usize;
            if si1 == b'B' && si2 == b'C' && slen == 2 && e + 6 <= extra.len() {
                bsize = Some(u16::from_le_bytes([extra[e + 4], extra[e + 5]]));
            }
            e += 4 + slen;
        }
        cursor += xlen;
    }
    if flg & 8 != 0 {
        // FNAME: NUL-terminated
        match buf[cursor..].iter().position(|&b| b == 0) {
            Some(p) => cursor += p + 1,
            None => return Ok(None),
        }
    }
    if flg & 16 != 0 {
        // FCOMMENT: NUL-terminated
        match buf[cursor..].iter().position(|&b| b == 0) {
            Some(p) => cursor += p + 1,
            None => return Ok(None),
        }
    }
    if flg & 2 != 0 {
        // FHCRC
        if buf.len() < cursor + 2 {
            return Ok(None);
        }
        cursor += 2;
    }

    Ok(Some(GzHeader {
        bsize,
        header_len: cursor,
    }))
}

/// Compresses `data` as a single gzip member (full-buffer).
pub fn gzip_compress(data: &[u8], level: i32) -> io::Result<Vec<u8>> {
    let level = CompressionLvl::new(level)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("{e:?}")))?;
    let mut compressor = Compressor::new(level);
    let mut out = vec![0u8; compressor.gzip_compress_bound(data.len())];
    let n = compressor
        .gzip_compress(data, &mut out)
        .map_err(io::Error::other)?;
    out.truncate(n);
    Ok(out)
}

/// Decompresses a single-member gzip buffer, rejecting output above `max_size`.
///
/// The output size is taken from the gzip ISIZE trailer; multi-member input
/// must go through [`GzReader`] instead.
pub fn gzip_decompress(data: &[u8], max_size: usize) -> io::Result<Vec<u8>> {
    if data.len() < 18 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "gzip buffer too short",
        ));
    }
    let isize = u32::from_le_bytes(data[data.len() - 4..].try_into().unwrap()) as usize;
    if isize > max_size {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("decompressed size exceeds maximum {max_size} bytes"),
        ));
    }
    let mut decompressor = Decompressor::new();
    let mut out = vec![0u8; isize];
    let n = decompressor
        .gzip_decompress(data, &mut out)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    out.truncate(n);
    Ok(out)
}

/// Decompresses a raw DEFLATE stream (one BGZF block payload) into `dst`.
pub trait BlockInflater {
    /// Inflates `cdata` into `dst`, which must be large enough for the output.
    fn inflate(&mut self, cdata: &[u8], dst: &mut [u8]) -> io::Result<()>;
}

/// libdeflate-backed inflater (C), reused across blocks.
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
pub(crate) struct CachedBlock {
    data: Box<[u8]>,
    size: u64,
}

/// A BGZF reader with an LRU cache of decompressed blocks for random access.
///
/// Repeated seeks into the same block hit the cache instead of re-reading and
/// re-inflating the block.
pub struct CachedBgzfReader {
    file: File,
    index: Option<GziIndex>,
    cache: LruCache<u64, CachedBlock>,
    inflater: Box<dyn BlockInflater>,
    /// Current block: (compressed offset, compressed size, uncompressed cursor, uncompressed len).
    current: Option<(u64, u64, usize, usize)>,
}

impl CachedBgzfReader {
    /// Opens a BGZF file (with its sibling `.gzi` index) and a block cache.
    pub fn open(path: impl AsRef<Path>, capacity: NonZeroUsize) -> anyhow::Result<Self> {
        let path = path.as_ref();
        let index_path = format!("{}.gzi", path.display());
        let index = GziIndex::read(&index_path)
            .with_context(|| format!("failed to read index {}", index_path))?;
        Self::open_inner(
            path,
            Some(index),
            capacity,
            Box::<LibdeflaterInflater>::default(),
        )
    }

    /// Opens a BGZF file in virtual-position mode (no `.gzi` index required).
    ///
    /// Seeks use pre-computed virtual positions, e.g. PAF lazy CIGAR loading.
    pub fn open_virtual(path: impl AsRef<Path>, capacity: NonZeroUsize) -> anyhow::Result<Self> {
        Self::open_inner(
            path.as_ref(),
            None,
            capacity,
            Box::<LibdeflaterInflater>::default(),
        )
    }

    /// Opens a BGZF file with a custom block inflater (for benchmarking).
    pub fn open_with_inflater(
        path: impl AsRef<Path>,
        capacity: NonZeroUsize,
        inflater: Box<dyn BlockInflater>,
    ) -> anyhow::Result<Self> {
        let path = path.as_ref();
        let index_path = format!("{}.gzi", path.display());
        let index = GziIndex::read(&index_path)
            .with_context(|| format!("failed to read index {}", index_path))?;
        Self::open_inner(path, Some(index), capacity, inflater)
    }

    fn open_inner(
        path: &Path,
        index: Option<GziIndex>,
        capacity: NonZeroUsize,
        inflater: Box<dyn BlockInflater>,
    ) -> anyhow::Result<Self> {
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

    /// Returns the current virtual position (compressed offset + block offset).
    pub fn virtual_position(&self) -> Option<VirtualPos> {
        self.current
            .map(|(cpos, _, upos, _)| VirtualPos::new(cpos, upos as u16))
    }

    /// Seeks to a pre-computed virtual position.
    pub fn seek_virtual(&mut self, vpos: VirtualPos) -> io::Result<()> {
        let (cpos, upos) = vpos.into();
        let (bsize, blen) = self.load_block(cpos)?;
        self.current = Some((cpos, bsize, usize::from(upos), blen));
        Ok(())
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
            if self.current.is_none() {
                match self.load_block(0) {
                    Ok((bsize, blen)) => self.current = Some((0, bsize, 0, blen)),
                    Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
                    Err(e) => return Err(e),
                }
            }
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

impl BufRead for CachedBgzfReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        loop {
            if self.current.is_none() {
                match self.load_block(0) {
                    Ok((bsize, blen)) => self.current = Some((0, bsize, 0, blen)),
                    Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(&[]),
                    Err(e) => return Err(e),
                }
            }
            let Some((cpos, bsize, upos, blen)) = self.current else {
                return Ok(&[]);
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
                        return Ok(&[]);
                    }
                    Err(e) => return Err(e),
                }
            }
            let block = self
                .cache
                .get(&cpos)
                .ok_or_else(|| io::Error::other("block not cached"))?;
            return Ok(&block.data[upos..blen]);
        }
    }

    fn consume(&mut self, amt: usize) {
        if let Some((_, _, upos, _)) = &mut self.current {
            *upos += amt;
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
        let (cpos, upos) = self
            .index
            .as_ref()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing index"))?
            .query(pos)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid index position"))?;
        let (bsize, blen) = self.load_block(cpos)?;
        self.current = Some((cpos, bsize, usize::from(upos), blen));
        Ok(pos)
    }
}

/// Reads and inflates the BGZF block at `cpos`, skipping empty blocks.
pub(crate) fn read_block(
    file: &mut File,
    cpos: u64,
    inflater: &mut dyn BlockInflater,
) -> io::Result<CachedBlock> {
    file.seek(SeekFrom::Start(cpos))?;

    // Fixed gzip header (10) + XLEN (2); the FEXTRA payload follows.
    let mut header = [0u8; 12];
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
    let data_start = (12 + xlen) as u64;
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
    let mut crc = Hasher::new();
    crc.update(&out);
    if crc.finalize() != expected_crc {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "CRC32 mismatch"));
    }

    Ok(CachedBlock {
        data: out.into_boxed_slice(),
        size: block_size,
    })
}

/// Builds a `.gzi` index for a BGZF file.
///
/// The first BGZF block (offset 0, 0) is implicitly skipped and empty blocks
/// (like the EOF marker, ISIZE = 0) are excluded, matching `bgzip -i`.
pub fn build_gzi_index(path: &str) -> anyhow::Result<()> {
    let mut file = File::open(path)?;
    let mut entries = Vec::new();
    let mut uncompressed_offset = 0u64;
    let mut compressed_offset = 0u64;

    loop {
        file.seek(SeekFrom::Start(compressed_offset))?;
        let mut header_fixed = [0u8; 12];
        match file.read_exact(&mut header_fixed) {
            Ok(_) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e.into()),
        }
        if header_fixed[0] != 0x1f || header_fixed[1] != 0x8b || header_fixed[3] & 4 == 0 {
            break;
        }
        let xlen = u16::from_le_bytes([header_fixed[10], header_fixed[11]]) as usize;
        if xlen == 0 {
            break;
        }
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
            anyhow::bail!("missing BC subfield at offset {}", compressed_offset);
        }
        let block_size = u64::from(bsize) + 1;
        if block_size < 4 {
            anyhow::bail!("malformed BGZF block at offset {}", compressed_offset);
        }
        file.seek(SeekFrom::Start(compressed_offset + block_size - 4))?;
        let mut isize_buf = [0u8; 4];
        file.read_exact(&mut isize_buf)?;
        let isize = u64::from(u32::from_le_bytes(isize_buf));

        if compressed_offset > 0 && isize > 0 {
            entries.push((compressed_offset, uncompressed_offset));
        }
        compressed_offset += block_size;
        uncompressed_offset += isize;
    }

    GziIndex::from_entries(entries).write(format!("{}.gzi", path))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::bgzf::writer::BgzfWriter;
    use std::io::Write;

    fn make_bgzf(data: &[u8], tag: &str) -> std::path::PathBuf {
        let path =
            std::env::temp_dir().join(format!("pgr_bgzf_test_{}_{}.gz", std::process::id(), tag));
        let mut writer = BgzfWriter::new(File::create(&path).unwrap()).unwrap();
        writer.write_all(data).unwrap();
        writer.finish().unwrap();
        build_gzi_index(path.to_str().unwrap()).unwrap();
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
        reader.seek(SeekFrom::Start(10_100)).unwrap();
        reader.read_exact(&mut got).unwrap();
        assert_eq!(got, data[10_100..10_150]);

        std::fs::remove_file(&path).ok();
        std::fs::remove_file(format!("{}.gzi", path.display())).ok();
    }

    #[test]
    fn gzip_roundtrip_large() {
        let data: Vec<u8> = (0..1_000_000u32).map(|i| (i % 251) as u8).collect();
        let packed = gzip_compress(&data, 6).unwrap();
        let isize = u32::from_le_bytes(packed[packed.len() - 4..].try_into().unwrap()) as usize;
        assert_eq!(isize, data.len());
        let out = gzip_decompress(&packed, data.len()).unwrap();
        assert_eq!(out, data);
    }

    #[test]
    fn gzip_roundtrip_small() {
        for data in [
            b"hello".to_vec(),
            b"hello world hello world".to_vec(),
            vec![0u8; 100],
        ] {
            let packed = gzip_compress(&data, 6).unwrap();
            let isize = u32::from_le_bytes(packed[packed.len() - 4..].try_into().unwrap()) as usize;
            eprintln!(
                "data_len={} packed_len={} isize={}",
                data.len(),
                packed.len(),
                isize
            );
            assert_eq!(isize, data.len());
            let out = gzip_decompress(&packed, 1000).unwrap();
            assert_eq!(out, data);
        }
    }
}
