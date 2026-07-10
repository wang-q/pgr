# pbit 设计笔记

> 本文档覆盖 `pgr pbit` 压缩格式，包括 LZ-diff 与 PAF 驱动的 CIGAR delta 编码。
> `pgr paf` 泛基因组命令族（query / graph / to-gfa 等）见 [`paf-pangenome.md`](../paf-pangenome.md)。

## 目标

借鉴 AGC (Assembled Genomes Compressor) v3.2.3 的 C++ 算法（LZ-diff、段级参考压缩、k-mer minimizer
参考选择，源码分析见 [agc-cpp.md](../references/agc-cpp.md)），设计 pgr 原生的群体基因组压缩格式
`pbit`（p = population 或 plus；2bit 记录 + delta），集成到 pgr 作为 `pgr pbit` 命令族。**不兼容**
C++ AGC 的 `.agc` 文件格式。

**关键约束**：

- **不兼容 `.agc` 文件格式**：设计原生 `pbit` 格式（扩展名 `.pbit`），仅借鉴 AGC 的算法（LZ-diff、
  段级参考压缩、k-mer minimizer 参考选择），不移植 CArchive 容器 / varint / 前缀编码
- 不引入 zstd 依赖，使用已有的 `flate2`（gzip/DEFLATE）替代
- **参考层直接复用 2bit 记录格式**（`dna_size + n_blocks + mask_blocks + packed_dna`），保留 N
  blocks 和 mask blocks（AGC 丢失了 mask），提取 `TwoBitFile::read_sequence` 核心逻辑为共享函数
  （见 §复用-2）
- 深度复用 pgr 现有基础设施（`libs/fmt/twobit`、`intspan::Range`、`lru::LruCache`、`libs/fmt/fa`）
- 遵循 pgr 分层原则：复杂逻辑放 `libs/pbit/`，`cmd_pgr/pbit/` 保持薄壳

## pgr 现有基础设施复用

### 1. 随机访问模式 — `libs/loc.rs`

pbit 的随机访问基于**文件偏移**（不再有 CArchive 的 stream/part 抽象）： `Decompressor` 解析
Header + Footer + Index 后，得到各参考段和 delta 的文件偏移，按需 seek + read。

| pbit 组件                                  | pgr Rust                                            | 说明             |
|--------------------------------------------|-----------------------------------------------------|------------------|
| `ref_groups[i].segment_offset`（u64 偏移） | `loc::read_offset(reader, offset, size)`            | 按偏移读取数据块 |
| `delta_data_offset` + delta `packed_size`  | seek + `read_exact`                                 | delta 随机读取   |
| `IndexMap<String, ...>`（参考组/样本索引） | `IndexMap<String, (u64, usize)>` (name→offset,size) | 名称→位置索引    |

**复用方式**：`Decompressor::read_sequence` 中 seek 到参考段偏移后调用 `read_2bit_record`， seek 到
delta 偏移后按 `packed_size` 读取——与 `loc::read_offset` 的"按偏移读取"模式一致。

参考：[fa/range.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/fa/range.rs) 的
`loc::open_indexed` + `loc::fetch_record` 模式。

### 2. 直接复用 2bit 记录格式 — `libs/fmt/twobit`（参考层底层）

pbit 的**参考层直接复用标准 2bit 记录**（`dna_size + n_blocks + mask_blocks + packed_dna`），
不是"借鉴模式"而是"调用同一函数"。这样 2bit 和 pbit 共享参考段的读写代码，pbit 仅在 2bit 记录
之上叠加 delta 层。

#### 2.1 提取共享函数 `read_2bit_record` / `write_2bit_record`

[twobit.rs:384-451](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L384) 的
`TwoBitFile::read_sequence` 内嵌了"从当前 reader 位置读取一个 2bit 记录"的逻辑（dna_size →n_blocks
→ mask_blocks → reserved → packed_dna → unpack）。将其核心抽取为模块级共享函数：

```rust
// libs/fmt/twobit.rs (新增 pub 函数)
/// Read a single 2bit record from the current reader position and return the
/// decoded DNA string with masks applied. Reused by TwoBitFile and pbit::Decompressor.
pub fn read_2bit_record<R: Read + Seek>(
    reader: &mut R,
    is_swapped: bool,
    start: Option<usize>,
    end: Option<usize>,
    no_mask: bool,
) -> Result<String>;

/// Write a single 2bit record (dna_size + n_blocks + mask_blocks + reserved +
/// packed_dna) to the writer. Reused by TwoBitWriter and pbit::Compressor.
pub fn write_2bit_record<W: Write>(
    writer: &mut W,
    dna: &str,
    do_mask: bool,
) -> Result<()>;
```

`TwoBitFile::read_sequence` 改为 `seek(offset) → read_2bit_record(...)` 的薄壳，
`TwoBitWriter::write` 改为循环调用 `write_2bit_record`。pbit 的 `Decompressor` 读参考段时也 调用
`read_2bit_record`，`Compressor` 写参考段时调用 `write_2bit_record`——**参考段字节级兼容 2bit 记录**。

#### 2.2 `SequenceReader` trait — 统一的随机访问接口

[libs/io.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/io.rs) 定义了 pgr 的序列随机访问抽象：

```rust
pub trait SequenceReader {
    fn read_sequence(
        &mut self,
        name: &str,
        start: Option<usize>,
        end: Option<usize>,
    ) -> anyhow::Result<String>;
}
```

`TwoBitFile<R>` 在 [twobit.rs:500](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L500-L510)
实现了此 trait（委托给固有的 5 参数 `read_sequence`）。**pbit 的 `Decompressor<R>` 也实现此 trait**，
但语义限定为**读参考层序列**（参考 contig 名唯一），供未来 chain/net 按参考坐标查询
pbit 归档——与查询 2bit 完全等价。样本层提取（`getctg` 语义，遍历多样本输出多 FASTA）不走
`SequenceReader`，改调 `Decompressor::get_contig`（见 §2.6）。

#### 2.3 `TwoBitFile<R: Read + Seek>` 结构模式 → `Decompressor<R>`

twobit.rs 的 `TwoBitFile<R>` 结构
（[twobit.rs:277-284](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L277-L284)）展示了
pgr 处理二进制基因组格式的标准模式，pbit 读取器逐项镜像：

- **泛型边界** `R: Read + Seek`
    - `TwoBitFile<R>` — 2bit 读取器
    - `Decompressor<R>` — pbit 读取器（直接持有 `R`，无中间容器层）
- **三构造器模式**
  （[twobit.rs:286-360](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L286-L360)）：
    - `open<P: AsRef<Path>>(path)` → 文件支持的 `BufReader<File>`
    - `open_and_read<P: AsRef<Path>>(path)` → 内存 `Cursor<Vec<u8>>`（小文件/测试用）
    - `new(reader)` → 解析 header + 索引（2bit）/ 解析 header + ref_index + sample_index（pbit）
- **名称→位置索引**：
    - 2bit 单层：`sequence_offsets: HashMap<String, u64>`（name → 文件偏移，
      [twobit.rs:280](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L280)）
    - pbit 双层：参考段层 `ref_groups: Vec<RefGroupEntry>`（每个 group = 一个参考段，含
      `contig_name` + `segment_offset`）+ `contig_groups: IndexMap<String, Vec<u32>>`（contig
      名 →ref_group_id 列表，供 `SequenceReader` 按 contig 名拼接多段）+ 样本层
      `collection: Collection`（sample → 各 contig 的段描述，供 `get_contig`/`get_sample` 遍历样本）
- **坐标提取 API**：
    - `TwoBitFile::read_sequence(name, start, end, no_mask)` — seek 到偏移，`read_2bit_record` 解包，
      应用 mask（[twobit.rs:384](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L384)）
    - `Decompressor::read_sequence(name, start, end)`（`SequenceReader`，读**参考层**）— 查
      `contig_groups[name]` 得 `Vec<ref_group_id>` → 遍历各 id: seek `ref_groups[id].segment_offset`
      → `read_2bit_record` → 拼接 → 按 `[start, end)` 切片
    - `Decompressor::get_contig(contig, start, end, strand, out)`（读**样本层**）—
      遍历 `collection.samples` → 对含此 contig 的样本：查 `SegmentDescs` → 每段按
      `ref_groups[ref_group_id].segment_offset` seek → `read_2bit_record` 读参考 → seek delta →
      flate2 解压 → LZ-diff 解码 →（if `is_rev_comp`: 反向互补）→ 拼接 → 切片 →（if `strand == "-"`:
      `rev_comp`）→ 写一条 FASTA

#### 2.4 `TwoBitWriter<W: Write>` 模式 → `Compressor<W>`

twobit.rs 的 `TwoBitWriter<W>`（
[twobit.rs:168-275](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L168-L275)）
展示了写入器模式，pbit 写入器对应：

- `TwoBitWriter<W: Write>` / `new(writer)` / `write(sequences, do_mask)` → 写 header + 索引 + 循环
  `write_2bit_record`
- `Compressor<W: Write + Seek>` / `new(writer, ...)` / `finish()` → 写 header → 写参考段
  （`write_2bit_record`）→ 写 delta 数据 → 写 sample_index → 写 footer（footer 在文件末尾，需
  `Seek`）。**参考段字节级与 2bit 记录一致**

#### 2.5 二进制解析工具函数

twobit.rs 的模块级 helper（`read_u32`、`read_u64`、`read_u32_vec`，
[twobit.rs:512-540](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L512-L540)）展示了
pgr 读取固定大小小端整数的模式。Phase 0 已将这些函数改为 `pub`（供未来复用），但 pbit 实际实现中
为避免冗余的 `is_swapped` 参数（pbit 统一小端，无需字节序切换），在 `format.rs` 内定义了独立的
`read_u32_le` / `read_u64_le` / `write_u32_le` / `write_u64_le`。记录级的 `read_2bit_record` /
`write_2bit_record` 则被 pbit 直接复用（见 §复用-2.4）。

#### 2.6 命令层复用 — `pbit/range.rs` 与 `twobit/range.rs` 共享区间工具

[twobit/range.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/twobit/range.rs) 和
`pbit/range.rs` 共享**区间解析与输出工具**（`ranges_arg()` / `collect_ranges()` /`intspan::Range`/
`nt::rev_comp` / `outfile_arg()` / `get_outfile()`），但提取逻辑不同——twobit 每个 range 读 1 条序列
（`SequenceReader::read_sequence`），pbit 每个 range 读所有样本 （`get_contig`，`getctg`语义，
输出多 FASTA）：

```rust
// twobit/range.rs (现有)                     // pbit/range.rs (待实现)
let mut tb = TwoBitFile::open(infile)?;        let mut pbit = Decompressor::open(infile)?;
for el in ranges.iter() {                      for el in ranges.iter() {
    let rg = intspan::Range::from_str(el);         let rg = intspan::Range::from_str(el);
    let seq_id = rg.chr();                         let seq_id = rg.chr();
    if !tb.sequence_offsets.contains_key(seq_id)   if !pbit.contains_contig(seq_id)
        { warn; continue; }                            { warn; continue; }
    let (start, end) = ...;  // 1-based→0-based     let (start, end) = ...;  // 同
    let mut seq =                                  // 遍历所有样本，每个含此 contig 的样本写一条 FASTA；
        tb.read_sequence(seq_id, start, end, false)?;  // strand 由 get_contig 内部逐条 rev_comp
    if rg.strand() == "-" {                        pbit.get_contig(seq_id, start, end,
        seq = rev_comp(seq);                                              rg.strand(), &mut writer)?;
    }                                              }
    writeln!(writer, ">{}", rg)?;              }
    writeln!(writer, "{}", seq)?;
}                                              }
```

### 3. 坐标/区间基础设施

| 组件                                    | 位置                                                               | pbit 用途                               |
|-----------------------------------------|--------------------------------------------------------------------|-----------------------------------------|
| `intspan::Range`                        | crate                                                              | `range` 命令的区间解析（`chr1:1-1000`） |
| `args::ranges_arg()`                    | [args.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/args.rs) | CLI 区间参数                            |
| `args::collect_ranges()`                | args.rs                                                            | 收集区间列表                            |
| `args::infile_arg_required_with_help()` | args.rs                                                            | 输入文件参数                            |
| `args::outfile_arg()` / `get_outfile()` | args.rs                                                            | 输出文件参数                            |

**复用方式**：`pgr pbit range` 可支持 `chr1:1-1000` 形式的区间提取，复用 `ranges_arg` +
`collect_ranges`，与 `fa/range.rs`、`twobit/range.rs` 保持一致。

### 4. FASTA I/O — `libs/fmt/fa`

- `pgr::libs::fmt::fa::reader(infile)` — noodles-based FASTA reader（支持 gzip）
- `pgr::libs::fmt::fa::writer(outfile)` — FASTA writer
- `libs::fasta` — chunk/dedup/filter/stat 等操作

**复用方式**：AGC 的 `genome_io.cpp`（FASTA 读写）完全不需要移植。压缩端用 `fa::reader` 读输入，
解压端用 `fa::writer` 写输出。

### 5. LRU 缓存 — `lru::LruCache`

[fa/range.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/fa/range.rs#L77) 已用
`lru::LruCache<String, noodles_fasta::Record>` 缓存序列。

**复用方式**：AGC 解压器缓存解压的参考段和 delta pack，用 `LruCache<u32, String>` （ref_group_id →
解码参考段；group = segment，单 u32 键即可）。

### 6. 压缩后端 — `flate2`（已有依赖）

AGC 用 zstd 的两处，均改为 flate2：

| AGC 用途                       | zstd API                                    | flate2 替代                                            |
|--------------------------------|---------------------------------------------|--------------------------------------------------------|
| 段级 delta 压缩（`CSegment`）  | `ZSTD_compressCCtx` / `ZSTD_decompressDCtx` | `flate2::write::GzEncoder` / `flate2::read::GzDecoder` |
| 元数据压缩（`CCollection_V3`） | 同上                                        | 同上                                                   |

**权衡**：gzip 压缩率略低于 zstd（通常差 5-15%），但无需新依赖。gzip 的压缩速度与 zstd level 3-6
相当，解压速度 gzip 略快。

**实现**：

```rust
use flate2::{write::GzEncoder, read::GzDecoder, Compression};

fn compress(data: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

fn decompress(data: &[u8]) -> Vec<u8> {
    let mut decoder = GzDecoder::new(data);
    let mut decoded = Vec::new();
    decoder.read_to_end(&mut decoded).unwrap();
    decoded
}
```

### 7. 其他可复用

| 组件                 | 位置                                                          | pbit 用途                                     |
|----------------------|---------------------------------------------------------------|-----------------------------------------------|
| `libs::nt::rev_comp` | [libs/nt](file:///Volumes/ExtHome/Scripts/pgr/src/libs/nt.rs) | LZ-diff 的反向互补参考                        |
| `indexmap::IndexMap` | crate                                                         | 保序 HashMap（pgr 统一模式）                  |
| `minimizer-iter`     | crate                                                         | 参考选择的 minimizer 计算                     |
| `murmurhash3`        | crate                                                         | LZ-diff 哈希表（与 C++ AGC 的 MurMur64 对齐） |
| `rayon`              | crate                                                         | 并行压缩（多 contig 同时编码）                |

## Rust 架构建议

### 模块结构

```
src/libs/pbit/
├── mod.rs              # 模块导出
├── format.rs           # Header/Footer/Index 序列化、版本常量、二进制 helper
├── lz_diff.rs          # CLZDiffBase/V1/V2 → LzDiff (LZ-diff 算法)
├── segment.rs          # CSegment → Segment (段级 delta 压缩/解压)
├── collection.rs       # CCollection_V3 → Collection (sample/contig/segment 元数据)
├── compressor.rs       # CAGCCompressor → Compressor (写 header + 参考段 + delta + index + footer)
└── decompressor.rs     # CAGCDecompressor → Decompressor (读 header + index + 按需解压段)

src/cmd_pgr/pbit/
├── mod.rs              # 子命令注册
├── create.rs           # pgr pbit create
├── append.rs           # pgr pbit append
├── to_fa.rs            # pgr pbit to-fa
├── some.rs             # pgr pbit some
├── range.rs            # pgr pbit range
└── stat.rs             # pgr pbit stat
```

> **无 `archive.rs`**：pbit 不用 CArchive 多流容器，无需单独模块。Compressor/Decompressor 直接持有
> `W: Write + Seek` / `R: Read + Seek`，按文件偏移读写参考段（2bit 记录）和 delta 数据块。

### 核心数据结构

```rust
// format.rs — Header/Footer + binary helpers (independent read_u32_le / read_u64_le etc.;
//             reuses twobit's read_2bit_record / write_2bit_record at record level)
pub const PBIT_MAGIC: u32 = 0x54494250;  // 'PBIT' — native "2bit + delta" format
pub const PBIT_VERSION_MAJOR: u32 = 1;
pub const PBIT_VERSION_MINOR: u32 = 0;

/// File header (fixed 36 bytes, at file start).
pub struct PbitHeader {
    pub magic: u32,           // PBIT_MAGIC
    pub version: u32,         // major*1000 + minor
    pub segment_size: u32,
    pub kmer_len: u32,
    pub min_match_len: u32,   // LZ-diff min match length (CLI -l, default 18)
    pub ref_group_count: u32,
    pub sample_count: u32,
    pub ref_records_offset: u64,  // offset to Reference Records (usually 36)
}

/// Footer (fixed 24 bytes, at end of file). Located by seeking to file_size - 24.
pub struct PbitFooter {
    pub ref_index_offset: u64,    // offset to reference group index
    pub delta_data_offset: u64,   // offset to delta data section
    pub sample_index_offset: u64, // offset to sample index
}

/// Reference group index entry (in ref_index section). Each group = one reference
/// segment (one 2bit record), mirroring C++ AGC where group = segment.
pub struct RefGroupEntry {
    pub contig_name: String,  // contig this segment belongs to
    pub segment_offset: u64,  // file offset to one 2bit record
}

/// Delta entry (in delta data section).
pub struct DeltaEntry {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_data: Vec<u8>,       // flate2(LZ-diff encoded)
    // packed_size is not stored as a field; it's derived via `meta()` → packed_data.len()
}

// lz_diff.rs — LZ-diff encoder/decoder (V2), mirrors C++ CLZDiff_V2.
pub struct LzDiff {
    reference: Vec<u8>,       // 2-bit encoded (A=0,C=1,G=2,T=3,N=4, other=31)
    ht: HashTable,            // u16 or u32 positions (built lazily by prepare_index)
    min_match_len: u32,       // default 18
    key_len: u32,             // = min_match_len - hashing_step + 1 (= 15)
    key_mask: u64,            // 2 * key_len bits
    hashing_step: u32,        // 4 (USE_SPARSE_HT)
    ht_size: u64,             // hash table size
    ht_mask: u64,             // hash table mask (= ht_size - 1)
    index_ready: bool,        // lazy index build flag
    min_nrun_len: u32,        // 4
}

enum HashTable {
    Short(Vec<u16>),   // ref < 65535 * hashing_step
    Long(Vec<u32>),
}

impl LzDiff {
    pub fn new(min_match_len: u32) -> Self;
    /// 2-bit encode + store reference. Required for both encode and decode.
    pub fn prepare(&mut self, reference: &[u8]);
    /// Build hash table over the stored reference. Required only for encode
    /// (decode does not need the hash table). Split from `prepare` so that
    /// `Decompressor` can skip the expensive index build.
    pub fn prepare_index(&mut self);
    pub fn encode(&mut self, text: &[u8], encoded: &mut Vec<u8>);
    /// Decode using `self.reference` (stored by `prepare`). No external
    /// reference parameter — avoids redundancy with the internal field.
    pub fn decode(&self, encoded: &[u8], decoded: &mut Vec<u8>);
    pub fn estimate(&mut self, text: &[u8], bound: u32) -> usize;
}

// segment.rs — delta encoding/decoding per reference group (no archive dependency)
/// Segment-level delta compression/decompression.
/// Holds a single `LzDiff` (which stores the 2-bit encoded reference internally);
/// does NOT duplicate the reference as ASCII.
pub struct Segment {
    lz_diff: LzDiff,
}

impl Segment {
    /// Prepare reference (call before add/get). Delegates to `lz_diff.prepare`.
    pub fn prepare(&mut self, ref_dna: &[u8]);
    /// Build LZ-diff hash table (call before add, not needed for get).
    pub fn prepare_index(&mut self);
    /// Encode `seq` against the prepared reference, return uncompressed delta bytes.
    pub fn add(&mut self, seq: &[u8]) -> Vec<u8>;
    /// Decode `delta` (uncompressed) against the prepared reference, return raw sequence.
    pub fn get(&self, delta: &[u8]) -> anyhow::Result<Vec<u8>>;
}

// collection.rs — sample/contig/segment metadata (fixed u32 fields, flate2 compressed)
pub struct SegmentDesc {
    pub ref_group_id: u32,    // index into ref_groups
    pub delta_id: u32,        // index into ref_group's delta list
    // is_rev_comp / raw_length live in DeltaEntry (properties of the encoding,
    // shared by all segments pointing to the same delta).
}

pub struct ContigSegs {
    pub contig_name: String,
    pub segments: Vec<SegmentDesc>,
}

pub struct Collection {
    pub samples: IndexMap<String, Vec<ContigSegs>>,  // sample → contigs
    pub cmd_line: String,
}

impl Collection {
    pub fn register_sample_contig(&mut self, sample: &str, contig: &str);
    pub fn add_segment(&mut self, sample: &str, contig: &str,
                       ref_group_id: u32, delta_id: u32);
    pub fn get_contig_segments(&self, sample: &str, contig: &str) -> Option<&[SegmentDesc]>;
    pub fn list_samples(&self) -> Vec<&str>;
    pub fn list_contigs(&self, sample: &str) -> Vec<&str>;
    /// Serialize with fixed-size u32 LE fields + flate2 (no prefix coding).
    pub fn serialize(&self) -> Vec<u8>;
    pub fn deserialize(data: &[u8]) -> anyhow::Result<Self>;
}

// compressor.rs — holds W: Write + Seek directly (no archive wrapper)
pub struct Compressor<W: Write + Seek> {
    writer: W,
    header: PbitHeader,
    ref_groups: Vec<RefGroupEntry>,         // one per reference segment (group = segment)
    deltas: Vec<Vec<DeltaEntry>>,           // deltas[ref_group_id][delta_id]
    collection: Collection,
    segments: Vec<Segment>,                 // one per ref group (each Segment holds its own min_match_len)
    segment_size: usize,
    kmer_len: usize,
    // min_match_len is not stored as a field; it's passed into Segment::new() and
    // held there. The header still records it (PbitHeader.min_match_len).
}

impl Compressor<BufWriter<File>> {
    /// Create from output path + reference FASTA.
    /// Writes header (placeholder) + reference records (via write_2bit_record).
    /// The CLI `create` command calls this, then `append_sample` for each `-i`
    /// input FASTA, then `finish`.
    pub fn create<P: AsRef<Path>>(out_path: P, ref_fasta: &str, segment_size: usize,
                                  kmer_len: usize, min_match_len: u32) -> anyhow::Result<Self>;
    /// Open an existing `.pbit` for appending samples (powers `pgr pbit append`).
    /// Uses Decompressor to read existing header / ref_groups / deltas (with
    /// packed_data) / collection, reads each reference segment via
    /// `read_2bit_record` and rebuilds `Segment`s (prepare + prepare_index),
    /// positions the writer at `footer.ref_index_offset` (start of Reference
    /// Index). `finish` then rewrites Reference Index + Delta Data + Sample Index
    /// + Footer as a contiguous block. No reference FASTA needed — reuses
    /// reference records already embedded in the archive.
    pub fn open_for_append<P: AsRef<Path>>(in_path: P) -> anyhow::Result<Self>;
}

impl<W: Write + Seek> Compressor<W> {
    /// Append a sample: read FASTA → segment → k-mer minimizer reference selection
    /// → LZ-diff encode → flate2 compress → delta dedup → store. Sample name is
    /// derived from the FASTA basename by the caller (CLI layer).
    pub fn append_sample(&mut self, sample_name: &str, fasta_path: &str) -> anyhow::Result<()>;
    /// Finalize: write Reference Index → Delta Data → Sample Index → Footer →
    /// patch Header offsets. Called once after all samples are appended.
    pub fn finish(self) -> anyhow::Result<()>;
}

// decompressor.rs — holds R: Read + Seek directly (no archive wrapper)
pub struct Decompressor<R: Read + Seek> {
    reader: R,
    header: PbitHeader,
    footer: PbitFooter,
    ref_groups: Vec<RefGroupEntry>,
    contig_groups: IndexMap<String, Vec<u32>>,  // contig name → ref_group_id list (for SequenceReader)
    contig_set: HashSet<String>,                // all contig names in collection (for contains_contig)
    collection: Collection,
    ref_cache: LruCache<u32, String>,           // ref_group_id → decoded ref segment
    // Delta metadata loaded during `new` by scanning Delta Data headers (9 bytes each).
    // Stores is_rev_comp / raw_length / packed_size per delta so that `get_contig`
    // can compute segment coordinates without re-reading headers from disk.
    delta_meta: Vec<Vec<DeltaMeta>>,            // [ref_group_id][delta_id] → (is_rc, raw_len, packed_size)
    delta_offsets: Vec<Vec<u64>>,              // delta_offsets[ref_group_id][delta_id] → file offset of header
    delta_cache: LruCache<(u32, u32), Vec<u8>>,  // (ref_group_id, delta_id) → decoded raw seq
}

/// In-memory delta header (loaded during `new`, no packed_data).
pub struct DeltaMeta {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_size: u32,
}

impl Decompressor<BufReader<File>> {
    /// Open from file path (mirrors TwoBitFile::open).
    pub fn open<P: AsRef<Path>>(path: P) -> anyhow::Result<Self>;
}

impl Decompressor<std::io::Cursor<Vec<u8>>> {
    /// Open and read entire file into memory (mirrors TwoBitFile::open_and_read).
    pub fn open_and_read<P: AsRef<Path>>(path: P) -> anyhow::Result<Self>;
}

impl<R: Read + Seek> Decompressor<R> {
    /// Construct from an already-opened reader: parse header + footer + indexes
    /// (mirrors TwoBitFile::new parsing header + sequence_offsets). Builds
    /// `contig_groups`, `contig_set`, and scans Delta Data headers (9 bytes each)
    /// to populate `delta_meta` + `delta_offsets` without decompressing data.
    pub fn new(reader: R) -> anyhow::Result<Self>;

    /// Check if a contig name exists in any sample's collection (getctg/getset
    /// semantics). Mirrors tb.sequence_offsets.contains_key but queries the
    /// sample layer, not the reference layer. Built once in `new` as a
    /// `HashSet<String>` from `collection.samples[*][*].contig_name`.
    pub fn contains_contig(&self, name: &str) -> bool;
    pub fn list_samples(&self) -> Vec<&str>;
    pub fn list_contigs(&self, sample: Option<&str>) -> Vec<&str>;
    pub fn get_sample(&mut self, sample: &str, out: &mut impl Write) -> anyhow::Result<()>;
    /// Extract a contig from ALL samples (getctg semantics), optionally sliced to
    /// [start, end). Writes one FASTA entry per sample that has this contig.
    /// Uses `delta_meta` to compute segment coordinates and only decodes segments
    /// overlapping [start, end) (smart selection, like the reference layer).
    /// Output header: `>{sample_name} {contig}:{start}-{end}({strand})`.
    /// If strand is "-", each sequence is reverse-complemented before writing.
    pub fn get_contig(&mut self, contig: &str, start: Option<usize>, end: Option<usize>,
                      strand: &str, out: &mut impl Write) -> anyhow::Result<()>;
}

// SequenceReader reads the REFERENCE layer (reference contig names are unique),
// NOT sample sequences. This lets future chain/net modules query pbit by reference
// coordinate just like they query 2bit. Sample extraction goes through get_contig /
// get_sample (multi-sample output, doesn't fit SequenceReader's single-return contract).
//
// **Assumption**: reference contig names and sample contig names share the same
// namespace (e.g. both use "chr1"). This holds when reference and samples are from
// the same species assembly. `range` (getctg) queries the sample layer via
// `contains_contig`; `SequenceReader` queries the reference layer via
// `contig_groups` — both keyed by contig name, so the assumption must hold for
// a contig name to be resolvable in both layers.
impl<R: Read + Seek> crate::libs::io::SequenceReader for Decompressor<R> {
    fn read_sequence(
        &mut self,
        name: &str,
        start: Option<usize>,
        end: Option<usize>,
    ) -> anyhow::Result<String> {
        // name = reference contig name → contig_groups[name] = Vec<ref_group_id> →
        //   walk segments, accumulate lengths, skip those before `start`, read only
        //   segments overlapping [start, end) via read_2bit_record, concat the slice.
        // (Mirrors TwoBitFile::read_sequence: seek → read_2bit_record → slice,
        //  but spans multiple segments per contig via contig_groups without
        //  materializing the full contig.)
        ...
    }
}
```

### CLI 设计

遵循 twobit 的单字子命令 + 分组惯例（`2bit` 有 info/subset/transform 三组，见
[twobit/mod.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/twobit/mod.rs)），pbit 命令设计如下
（括注对应的 C++ AGC 命令，源码分析见 [agc-cpp.md](../references/agc-cpp.md)）：

- `create` → `create`（build 组，无 twobit 对应）
- `append` → `append`（build 组，无 twobit 对应）
- `getcol` → `to-fa`（transform 组，对应 `2bit to-fa`）
- `getset` → `some`（subset 组，对应 `2bit some`）
- `getctg` → `range`（subset 组，对应 `2bit range`）
- `info` + `listref` + `listset` + `listctg` → `stat`（info 组，对应 `2bit size`，用 flag 区分输出）

```
pgr pbit create    -i input.fa [-i input2.fa ...] -o out.pbit -r ref.fa [-s segment_size] [-k kmer_len] [-l min_match_len]
pgr pbit append    in.pbit -i input.fa [-i input2.fa ...] [-o out.pbit]  # 归档为位置参数输入，-o 省略时原地修改
pgr pbit to-fa     in.pbit -o out_dir/          # 提取所有样本为 FASTA，每样本一个文件 out_dir/{sample}.fa
pgr pbit some      in.pbit sample_list.txt [-o out.fa] [--invert]  # 按样本名列表提取，输出多 FASTA
pgr pbit range     in.pbit "chr1" "chr2:1-1000" [-o out.fa]  # 按 contig/区间提取（遍历所有样本，输出多 FASTA）
pgr pbit stat      in.pbit [--samples | --refs | --contigs [-s sample]]  # 统计/列表
```

> **参数说明**：

> - `-s segment_size`：分段大小（bp，默认 4096）
> - `-k kmer_len`：k-mer minimizer 参考选择的 k-mer 长度（默认 15）
> - `-l min_match_len`：LZ-diff 最小匹配长度（默认 18，对应 key_len=15）
> - 样本名从 `-i` 指定的 FASTA 文件 basename 派生（`libs::io::get_basename`），如 `path/to/sample1.fa`
  → 样本名 `sample1`
> - `append` 的归档文件为位置参数（输入），`-o` 可选（省略时原地修改 `in.pbit`， 指定时先复制再追加）；
  参考已嵌入归档，无需 `-r`
> - `some` / `range` 输出 FASTA header 格式：`>{sample_name} {contig}:{start}-{end}({strand})`
  （`some` 无区间时为 `>{sample_name} {contig}`）

> **stat flag 语义**：

> - 默认（无 flag）：输出归档总览统计（ref_group_count、sample_count、segment_size、kmer_len、
  min_match_len、各参考 contig 段数、各样本 contig 数）。
> - `--samples`：列出所有样本名（对应 C++ `listset`）。
> - `--refs`：列出参考 contig 名（来自 `ref_groups[*].contig_name`，去重保序；对应 C++ `listref`， 但
  C++ 列"参考样本名"——pbit 无参考样本实体，参考由 contig 名标识）。
> - `--contigs [-s sample]`：无 `-s` 列出所有样本的 `sample\tcontig` 对（对应 C++ `listctg`）； 有
  `-s` 仅列指定样本的 contig 名。

子命令分组（镜像 `2bit` mod.rs 结构）：

- build: `create` / `append`
- info: `stat`
- subset: `range` / `some`
- transform: `to-fa`

`range` 复用 `ranges_arg` 支持区间提取：`pgr pbit range in.pbit "chr1:1-1000" "chr2(-):500-1000"`。
区间解析与 `pgr 2bit range` 一致（共享 `ranges_arg` / `intspan::Range`），但提取语义不同——pbit
遍历所有样本输出多 FASTA（`getctg` 语义），2bit 读单序列（见 §复用-2.6）。

## 文件格式规范

pbit 格式（原生"2bit + delta"，扩展名 `.pbit`，区别于 C++ AGC 的 `.agc`）。所有整数使用
**固定大小 小端序**（u32 = 4 字节，u64 = 8 字节），不使用 varint / 前缀编码。字符串采用**长度前缀**
（u32 len + UTF-8 bytes），不使用 null 终止。

### 文件结构总览

```
┌─────────────────────────────────────┐  ← offset 0
│ Header (固定 36 字节)               │
├─────────────────────────────────────┤
│ Reference Records                   │  ← 参考段，每段为标准 2bit 记录
│   ref_group 0: seg 0, seg 1, ...    │     (dna_size + n_blocks + mask_blocks
│   ref_group 1: seg 0, seg 1, ...    │      + reserved + packed_dna)
│   ...                               │
├─────────────────────────────────────┤  ← footer.ref_index_offset
│ Reference Index                     │  ← 每个参考组的名称 + 段文件偏移列表
├─────────────────────────────────────┤  ← footer.delta_data_offset
│ Delta Data                          │  ← 每个参考组的 delta 列表
│   ref_group 0: delta 0, delta 1,... │     (is_rc + raw_len + packed_size +
│   ref_group 1: delta 0, delta 1,... │      flate2(LZ-diff encoded))
├─────────────────────────────────────┤  ← footer.sample_index_offset
│ Sample Index (collection)           │  ← flate2(序列化的 samples/contigs/segments)
├─────────────────────────────────────┤
│ Footer (固定 24 字节)               │  ← 三个 section 的偏移
└─────────────────────────────────────┘  ← EOF
```

### Header（36 字节，文件起始）

```
offset  size  field              说明
0       4     magic              0x54494250 ('PBIT', 小端)
4       4     version            major*1000 + minor (当前 1000)
8       4     segment_size       分段大小（bp，如 4096）
12      4     kmer_len           k-mer 长度（如 15）
16      4     min_match_len      LZ-diff 最小匹配长度（如 18，对应 -l 参数）
20      4     ref_group_count    参考段总数（每段 = 一个 group，非 contig 数）
24      4     sample_count       样本数（不含参考）
28      8     ref_records_offset Reference Records 起始偏移（通常 = 36）
```

### Reference Records（标准 2bit 记录，连续存储）

每个参考段为一个**标准 2bit 记录**（与
[twobit.rs:384](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L384)`read_sequence`
读取的格式字节级一致），由共享函数 `write_2bit_record` / `read_2bit_record`读写：

```
offset  size           field
0       4              dna_size                 序列长度
4       4              n_block_count            N-block 数
8       4*n_block_count n_starts                N-block 起始位置
        4*n_block_count n_sizes                 N-block 长度
        4              mask_block_count         Mask-block 数
        4*mask_block_count mask_starts          Mask-block 起始位置
        4*mask_block_count mask_sizes           Mask-block 长度
        4              reserved                 保留字段（0）
        (dna_size+3)/4 packed_dna               2-bit 打包 DNA
```

> **与 C++ AGC 的关键差异**：保留 mask blocks（C++ AGC 丢失 mask）；参考段不经过 zstd/flate2
> 二次压缩（2bit 打包已压缩 4 倍，delta 层负责进一步压缩）。

### Reference Index（参考组索引）

每个 ref_group = 一个参考段（一个 2bit 记录）：

```
offset  size           field
0       4              ref_group_count          (= header.ref_group_count)
for each ref_group:
  4      name_len       contig 名字节长度
  name_len name_bytes   contig 名（UTF-8，该段所属的参考 contig）
  8      segment_offset 该段的文件偏移（u64 LE，指向 Reference Records 中的一个 2bit 记录）
```

### Delta Data（delta 数据区）

每个参考组有一组唯一 delta（去重后的 LZ-diff 编码）。同一参考组的 delta 连续存储：

```
offset  size           field
0       4              ref_group_count
for each ref_group:
  4      delta_count    唯一 delta 数
  for each delta:
    1      is_rev_comp   0 = 正向, 1 = 反向互补
    4      raw_length    原始序列长度
    4      packed_size   flate2 压缩后字节数
    packed_size packed_data  flate2(LZ-diff encoded bytes)
```

> **随机访问**：delta 为变长记录（`packed_size` 各异），不在文件中存储 per-delta offset 表。
> `Decompressor::new` 顺序扫描一次 Delta Data 区，按每条 delta 头部的 `packed_size` 累加偏移，
> 构建内存中的 `delta_offsets: Vec<Vec<u64>>`（`[ref_group_id][delta_id]` → 文件偏移）。
> 扫描仅需读取每条 delta 的 9 字节头部（`is_rev_comp + raw_length + packed_size`），不解压数据。

### Sample Index（collection，flate2 压缩）

序列化为固定大小字段后用 flate2 压缩：

```
[未压缩的逻辑结构]
4      sample_count
for each sample:
  4      name_len
  name_len name_bytes     样本名
  4      contig_count
  for each contig:
    4      contig_name_len
    contig_name_len contig_name_bytes
    4      segment_count
    for each segment:
      4      ref_group_id    参考组索引
      4      delta_id        该参考组内的 delta 索引
      (is_rev_comp / raw_length 由 DeltaEntry 提供，此处不重复)
4      cmd_line_len
cmd_line_len cmd_line_bytes  命令行（记录创建参数）
```

### Footer（24 字节，文件末尾）

```
offset  size  field
0       8     ref_index_offset     Reference Index 起始偏移
8       8     delta_data_offset    Delta Data 起始偏移
16      8     sample_index_offset  Sample Index 起始偏移
```

> **设计要点**：

> - **无 magic number 在末尾**：Header 已含 magic，Footer 仅 24 字节偏移
> - **无变长整数**：所有字段固定大小，简化解析
> - **参考段 = 2bit 记录**：可直接用 `read_2bit_record` 读取，与 twobit.rs 共享代码
> - **delta 层独立压缩**：每个 delta 单独 flate2 压缩，支持随机访问单个样本/contig
> - **collection 整块压缩**：元数据通常较小，一次 flate2 压缩即可
> - **冗余计数字段**：`ref_group_count`/`sample_count` 在 Header、Reference Index、Delta Data、
  Sample Index 中重复出现，这是有意为之——读取各 section 时可就地校验一致性（不匹配则报错），
  避免损坏文件导致越界访问

## PAF 驱动的 CIGAR delta 编码

> **状态**：Phase 8a–8e 已全部实现。PAF 驱动模式是 pbit 格式的正式组成部分，用 PAF 比对结果驱动
> 参考选择并直接存储 CIGAR，与 LZ-diff 共存于同一归档。
>
> **实现索引**：CIGAR 编解码 [cigar_delta.rs](../../src/libs/pbit/cigar_delta.rs)；
> 格式扩展 [format.rs](../../src/libs/pbit/format.rs)（`DeltaEncoding`/10 字节 `DeltaMeta`/16 字节 `SegmentDesc`，版本 1001）；
> query-side 索引 [paf_index.rs](../../src/libs/pbit/paf_index.rs)（`PafQueryIndex`）；
> 压缩端 [compressor.rs](../../src/libs/pbit/compressor.rs)（`append_sample_with_paf`）；
> 解压端 [decompressor.rs](../../src/libs/pbit/decompressor.rs)（`decode_delta(&SegmentDesc)`）；
> CLI [create.rs](../../src/cmd_pgr/pbit/create.rs) / [append.rs](../../src/cmd_pgr/pbit/append.rs)（`--paf`/`-p`、3 列 TSV）；
> 集成测试 [cli_pbit_paf.rs](../../tests/cli_pbit_paf.rs)。

### 1. 动机

#### 当前局限

pbit 当前的参考选择是**按段位置索引匹配**（非 minimizer；[compressor.rs:280-336](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/compressor.rs#L280-L336)）：

- 按 contig 名匹配参考段（`contig_ref_groups[contig_name]`）
- 按段位置索引匹配（`ref_group_ids[seg_idx]`，clamped 到最后一段）
- 首段 k-mer 采样**仅用于方向检测**（`detect_rev_comp`，非参考选择）+ 逐段 delta 大小回退

**无法处理**：重排、转位、跨 contig 比对、大段 indel 导致的位置偏移。对标准群体基因组
（同源染色体对齐）足够，但对含结构变异的多样性基因组压缩率差。

#### PAF 驱动的优势

用户用 minimap2/wfmash 等工具将样本比对到参考，生成 PAF（含精确 `=/X` CIGAR）。pbit 可利用
该比对结果：

1. **精确参考选择**：PAF 直接给出每段样本对应哪个参考段（含跨 contig/重排）
2. **精确方向**：PAF 的 strand 字段比 k-mer 采样检测更可靠
3. **精确差异**：CIGAR 的 `=/X/I/D` 操作已完整描述样本与参考的差异，可直接存储
4. **可复用比对信息**：解压后可还原比对关系，支持变异分析，无需重新跑比对工具

### 2. pgr 已有的 PAF + CIGAR 基础设施

pgr 已有完整的 PAF 处理栈，pbit 可直接复用：

| 组件 | 位置 | 说明 |
|------|------|------|
| `CigarOp` | [libs/paf/cigar.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/cigar.rs) | bit-packed u32（3位op + 29位length），已有 `from_raw`/`op()`/`len()`/`target_delta`/`query_delta` |
| `CigarStore` | [libs/paf/index/mod.rs:22](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/index/mod.rs#L22) | `Owned(Vec<CigarOp>)` / `Lazy(u64)` / `LazyReversed(u64)` |
| `PafIndex` | [libs/paf/index/](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/index) | coitrees 区间树索引，按坐标查询比对 |
| `build_pairwise_block` | [libs/paf/msa_build.rs:191](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/msa_build.rs#L191) | 从 CIGAR + FastaStore 重建比对序列（含正负链处理） |
| `extract_cigar` | libs/paf/cigar.rs | 从 PAF tags 解析 `cg:Z:` CIGAR 字符串为 `Vec<CigarOp>` |
| `reverse_cigar` | libs/paf/cigar.rs | 反转 CIGAR（用于负链） |
| `FastaStore` | libs/paf/fasta.rs | 序列存储，按名+区间提取 |

**关键洞察**：`CigarOp` 的 bit-packed u32 格式天然适合 pbit 的二进制存储。pgr 已有从 CIGAR +
参考序列重建样本序列的完整逻辑（`build_pairwise_block` → `build_maf_block`），pbit 解压时可复用。

### 3. 方案选项

#### 方案 A：PAF 仅指导参考选择（LZ-diff 不变）

PAF 用于选择参考段和方向，delta 仍用 LZ-diff 编码。PAF 不存储在文件中。

```
压缩：PAF → 确定段级对应 → LZ-diff 编码差异 → flate2 压缩
解压：LZ-diff 解码 → 重建样本（与当前一致）
```

- **优点**：改动最小，文件格式不变，解压端无需改动
- **缺点**：不保留比对信息；LZ-diff 对含 indel 的段编码效率不如 CIGAR（indel 导致后续位置错位，
  LZ-diff 的 match 失效，退化为大量 literal）

#### 方案 B：PAF 指导 + 存储比对元数据

PAF 指导参考选择，在 Collection 中存储比对摘要（ref_group_id, strand, query/target 坐标），
delta 仍用 LZ-diff 编码。

```
压缩：PAF → 确定段级对应 → LZ-diff 编码差异 → flate2 压缩；存储比对摘要
解压：LZ-diff 解码 → 重建样本；比对摘要可供查询
```

- **优点**：保留比对关系，解压后可查询；格式变化小
- **缺点**：仍用 LZ-diff 编码（indel 效率问题）；比对摘要与 LZ-diff delta 有冗余

#### 方案 C：CIGAR 替代 LZ-diff 作为 delta 编码（推荐）

用 PAF 的 CIGAR 直接作为 delta 编码，替代 LZ-diff。CIGAR 的 `=/X/I/D` 操作完整描述差异，
解压时按 CIGAR 从参考序列重建样本序列。

```
压缩：解析 PAF → 提取每段 CIGAR + X/I 碱基 → bit-packed → flate2 压缩
解压：flate2 解压 → Vec<CigarOp> + X/I 碱基流 → 按 CIGAR 从参考段重建样本段
```

> `packed_data` 的完整二进制格式（CIGAR ops + X/I 碱基流）见 §8 "X/I 操作的碱基存储问题"。

- **优点**：
  - 复用 `CigarOp` bit-packed 编码（4 字节/op，高效）
  - 不需要重新计算 LZ-diff（省压缩时间）
  - 天然处理 indel（I/D 操作，不受位置错位影响）
  - 解压逻辑与 `build_pairwise_block` 一致
  - 保留比对信息（CIGAR 本身就是比对结果）
- **缺点**：
  - 依赖外部比对工具（用户需先跑 minimap2/wfmash）
  - 需处理未覆盖区域（PAF 没覆盖的样本序列）
  - 文件格式需扩展（delta 编码类型标志）

#### 方案 D：混合编码（LZ-diff + CIGAR 可选）

delta 层支持两种编码，每段独立选择：
- 有 PAF 覆盖的段 → CIGAR 编码
- 无 PAF 覆盖的段 → LZ-diff 编码（回退到当前 minimizer 逻辑）

```
DeltaEntry {
    encoding: u8,  // 0 = LZ-diff, 1 = CIGAR
    packed_data: Vec<u8>,  // flate2(编码数据)
}
```

- **优点**：兼顾两种场景；未覆盖区域不丢数据
- **缺点**：解压端需支持两种解码路径；复杂度增加

### 4. 推荐方案：C + D 混合

**推荐方案 C 为主，D 为回退**：

- 有 PAF 输入时：PAF 覆盖的段用 CIGAR 编码，未覆盖段回退 LZ-diff
- 无 PAF 输入时：全部用 LZ-diff（当前行为，完全向后兼容）

理由：
1. CIGAR 对含 indel/重排的段编码效率优于 LZ-diff
2. `CigarOp` bit-packed 格式与 pbit 二进制风格一致
3. pgr 已有完整的 CIGAR 重建逻辑可复用
4. 混合模式保证未覆盖区域不丢数据，且与现有格式向后兼容

### 5. 格式扩展

#### 5.1 DeltaEntry / DeltaMeta 扩展

pbit 现有**两个**结构体共享同一在盘头部：
- `DeltaMeta`（9 字节，仅头部）—— `Decompressor::new` 扫描所有 delta 头部构建
  `delta_meta: Vec<Vec<DeltaMeta>>`（[decompressor.rs:42](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/decompressor.rs#L42)），**不读 packed_data**
- `DeltaEntry`（完整，含 packed_data）—— 按需读取，通过 `meta()` 派生 `DeltaMeta`

`encoding` 字段**必须加入 `DeltaMeta`**：`get_contig` 决定走 LZ-diff 还是 CIGAR 解码路径时只看
`delta_meta`（不读 packed_data）。若 `encoding` 不在 `DeltaMeta`，则必须先读完整 delta 才能判断编码
——破坏随机访问。

```rust
pub struct DeltaMeta {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_size: u32,
    pub encoding: DeltaEncoding,  // 新增（在盘第 10 字节）
}

pub struct DeltaEntry {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_data: Vec<u8>,
    pub encoding: DeltaEncoding,  // 新增；meta() 自动带上
}

pub enum DeltaEncoding {
    LzDiff = 0,
    Cigar = 1,
}
```

**在盘格式**（Delta Data 区，每条 delta）：

```
offset  size  field
0       1     is_rev_comp
1       4     raw_length
5       4     packed_size
9       1     encoding        ← 新增（0 = LZ-diff, 1 = CIGAR）
10      packed_size packed_data
```

> encoding 字段是新增的。现有文件的 delta 头部是 9 字节（无 encoding），新文件为 10 字节。
> `DeltaMeta::read_header` / `write_header` 改为读写 10 字节。
> 本项目尚未正式发布，无需考虑向后兼容，直接升级 Header 版本号（1000 → 1001）。
>
> **raw_length 在 CIGAR 模式下的语义**：LZ-diff 的 `raw_length` 是样本段长度；
> CIGAR 模式下 `raw_length = Σ op.query_delta()`（query 轴长度，含 `X`/`I`，不含 `D`），
> 即样本段长度。`get_contig` 用 `raw_length` 做切片计算，CIGAR 段必须提供正确的 query 轴长度。

#### 5.2 SegmentDesc 扩展（CIGAR 模式，固定大小存储）

CIGAR 模式下，每段需要记录比对坐标（样本段在参考段中的偏移）。采用**固定大小存储**：
SegmentDesc 新增 `ref_start` / `ref_end` 两个 u32 字段，LZ-diff 段填 0。

```rust
pub struct SegmentDesc {
    pub ref_group_id: u32,
    pub delta_id: u32,
    /// Segment-relative start offset within the reference 2bit record
    /// (= `target_start - seg_idx * segment_size`). 0 for LZ-diff segments.
    pub ref_start: u32,
    /// Segment-relative end offset (exclusive) within the reference 2bit
    /// record (= `target_end - seg_idx * segment_size`). 0 for LZ-diff segments.
    pub ref_end: u32,
    // is_rev_comp / raw_length / encoding 由 DeltaEntry 提供
}
```

**在盘格式**（Sample Index 区，每段，固定 16 字节）：

```
offset  size  field
0       4     ref_group_id
4       4     delta_id
8       4     ref_start         ← 新增（CIGAR 段：参考段内相对偏移；LZ-diff 段：0）
12      4     ref_end           ← 新增（CIGAR 段：参考段内相对结束；LZ-diff 段：0）
```

> **坐标语义（相对坐标）**：`ref_start`/`ref_end` 存储的是**相对于参考段（2bit 记录）起始的
> 偏移**，而非 PAF 的 contig 绝对坐标。计算方式：`seg_idx = target_start / segment_size`，
> `ref_start = target_start - seg_idx * segment_size`，
> `ref_end = target_end - seg_idx * segment_size`。解压时 `read_2bit_record` 读出的参考段
> 从 `seg_idx * segment_size` 开始，直接用 `ref_start..ref_end` 切片即可定位 CIGAR 对应的
> 参考区间，无需反查 `contig_ref_groups`。
>
> **为什么需要 ref_start/ref_end**：CIGAR 描述的是"样本段 vs 参考段某区间"的差异。
> 同一参考段可能被多个样本段引用（不同区间），需要 ref_start/ref_end 定位。LZ-diff 段
> 参考段是整条 2bit 记录，ref_start=0、ref_end=记录长度，但为保持固定大小统一填 0
> （解压时 LZ-diff 路径不读这两个字段）。
>
> **为什么不用变长存储**：变长方案（按 encoding 决定是否后跟 ref_start/ref_end）会破坏
> `SegmentDesc` 的 `Copy` trait，并让 serialize/deserialize 复杂化。固定大小方案虽然 LZ-diff
> 段多占 8 字节，但 Collection 整体经 flate2 压缩，连续零值压缩率极高，实际开销可忽略。
>
> **encoding 位置**：encoding 只放 `DeltaMeta`（§5.1），不放 SegmentDesc。一个 delta 的
> encoding 是其 packed_data 的固有属性，不会因引用它的 segment 不同而变化。`Decompressor::new`
> 扫描 delta 头时即可获知编码类型，`get_contig` 据
> `delta_meta[ref_group_id][delta_id].encoding` 分支解压路径，无需读 SegmentDesc。

#### 5.3 CIGAR 的 bit-packed 存储

`Vec<CigarOp>` 直接序列化为 `Vec<u32>`（每个 CigarOp 是一个 u32），配合 X/I 碱基流，
然后 flate2 压缩。完整 `packed_data` 格式见 §8。

```
packed_data = flate2( u32 op_count + [CigarOp; op_count] + u32 base_count + [u8; base_count] )
```

解压时：`flate2_decompress(packed_data) → (Vec<CigarOp>, Vec<u8> X/I bases)`

> **空间效率**：一个 4096 bp 的纯 match 段 → `[=4096]` → 1 个 CigarOp = 4 字节 → flate2 后更小。
> 含 40 个 SNP 的段 → `=99 X1` × 40 + `=96` ≈ 82 个 CigarOp = 328 字节 → flate2 后约 100-200 字节。
>
> **压缩率 caveat**：上述估算需基准验证。CigarOp 是 `(op << 29) | len` 的 bit-packed u32，op 占
> 高位 3 bit、len 占低 29 bit。相邻 CigarOp 的 op 字段可能有局部重复（如 `=99 X1 =99 X1...`），
> 但 len 字段差异大。DEFLATE 对这种半结构化 u32 数组的压缩率通常不如对文本/LZ-diff 字节流。
> **SoA 优化候选**（Phase 8f）：将 ops 与 lengths 分离存储（两个独立数组），可提升 DEFLATE 对
> op 字段（高度重复）的压缩效率。作为可选优化，默认采用 interleaved（AoS）布局。

### 6. CLI 设计

#### 6.1 新增 `--paf` 参数

```
pgr pbit create -r ref.fa -i sample.fa --paf sample.paf -o out.pbit
pgr pbit append in.pbit -i sample.fa --paf sample.paf
```

- `--paf` 可选，指定样本比对到参考的 PAF 文件
- 一个 `-i` 对应一个 `--paf`（或支持一个 PAF 含多个样本的比对）
- 省略 `--paf` 时使用当前 LZ-diff 模式（按段位置索引匹配参考，见 §1）

#### 6.2 参数语义

```
pgr pbit create -r ref.fa \
    -i sample1.fa --paf sample1.paf \    # sample1 用 CIGAR 模式
    -i sample2.fa                        # sample2 用 LZ-diff 模式（无 PAF）
```

> 同一归档可混用两种模式（每样本独立选择）。delta 层的 encoding 字段区分。

#### 6.3 `--name` TSV 模式与 `--paf` 的关系

现有 `--name` TSV 接受两列 `sample_name<TAB>fasta_path`（[create.rs:84-88](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/pbit/create.rs#L84-L88) /
[append.rs:48-52](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/pbit/append.rs#L48-L52)）。为支持 PAF，
扩展为**可选三列**：

```
# samples.tsv
sample1    sample1.fa    sample1.paf    # CIGAR 模式
sample2    sample2.fa                   # LZ-diff 模式（第三列缺失）
sample3    sample3.fa                   # LZ-diff 模式
```

**互斥规则**：
- `--paf` 仅与 `-i` 配套（`-i` 与 `--paf` 一一对应，按出现顺序配对）
- `--name` 与 `--paf` **互斥**（`--name` 模式下 PAF 由 TSV 第三列指定，不允许同时用 `--paf`）
- `--name` TSV 第三列缺失时该样本走 LZ-diff（向后兼容现有两列 TSV）

> **推荐**：多样本 + PAF 场景优先使用 `--name` TSV 三列模式。`-i`/`--paf` 按顺序配对
> 在样本数较多时容易出错（如 `-i s1.fa -i s2.fa --paf s1.paf --paf s2.paf` 与
> `-i s1.fa --paf s1.paf -i s2.fa --paf s2.paf` 顺序不同但意图相同，易混淆）。TSV 模式
> 将样本名、FASTA、PAF 显式绑定在同一行，无歧义。CLI 实现时对 `-i`/`--paf` 配对应做
> 计数校验（数量必须相等，否则报错）。

### 7. 压缩流程（CIGAR 模式）

```
append_sample_with_paf(sample_name, fasta_path, paf_path):
  1. 解析 PAF → 构建 query-side 区间树（按 query 坐标索引，见下方"PAF query-side
     索引"说明）
  2. 读取样本 FASTA → 按 segment_size 分段
  3. 对每个样本段 [seg_start, seg_end):
     a. 查 query-side 区间树：哪些 PAF alignment 覆盖此样本段区间？
     b. 选最佳 alignment（按 identity / 覆盖度）
     c. 提取该段对应的 CIGAR 片段（从 PAF CIGAR 中按坐标截取）
     d. 若 CIGAR 含 `M` 操作：对比样本碱基与参考碱基，拆分为 `=/X`
        （fallback：`pgr maf to-paf` 已在源头区分 `=/X`，此步骤仅处理
        minimap2/fastga 等直接输出的 `M` CIGAR）。参考碱基从
        `self.segments[ref_group_id].reference_dna()` 按 CIGAR target 坐标切片获取
        （target 坐标减去 `seg_idx * segment_size` 转为参考段内偏移）。
     d'. 提取 X/I 碱基：按 CIGAR 正向遍历顺序，收集 `X`/`I` 操作对应的样本碱基。
        **`-` 链记录**：CIGAR 描述的是 RC(query) vs forward(target) 比对，因此 X/I 碱基
        必须从 **RC(sample)** 提取（与 CIGAR 描述的比对方向一致），而非原始正向样本序列。
        解压时先正向应用 CIGAR（用存储的 RC 样本碱基）→ 得到 RC(sample) → 再 RC 得到 sample
        （见 §8 负链语义）
     e. 确定 ref_group_id（PAF target → 参考段映射）+ strand + ref_start/ref_end
     f. bit-pack CIGAR → flate2 压缩 → DeltaEntry(encoding=CIGAR)
     g. 记录 SegmentDesc(ref_group_id, delta_id, ref_start, ref_end)
     h. Delta 去重按 packed_data 字节比较（与 LZ-diff 一致），相同 CIGAR+XI 的段共享 delta_id
  4. 未覆盖段：回退 LZ-diff（对到同名参考段，当前逻辑）
```

#### 关键问题：PAF 坐标→段映射

PAF alignment 的 query_start/end 是连续坐标，pbit 按 segment_size 切段。一条 PAF alignment
可能横跨多段。处理方式：

1. **保持固定段切分**：样本仍按 segment_size 切分
2. **CIGAR 截取**：对每段，从覆盖它的 PAF alignment 的 CIGAR 中截取对应区间
   - 用 `CigarOp::target_delta` / `query_delta` 做坐标投影
   - 一段可能被多条 alignment 覆盖 → 选最佳（identity 最高 + 覆盖度最大）
   - 一段可能部分被覆盖或跨多条 alignment 衔接点 → **整段回退 LZ-diff**（见 §11 决策 3，
     不拆段、不混用 CIGAR 与 LZ-diff）
3. **CIGAR target 投影跨参考段边界检查**：样本段按 query 坐标切分，但其 CIGAR 投影到
   target 轴时，因 indel 偏移，target 区间 `[target_start, target_end)` 可能跨越参考段边界
   （即 `target_start / segment_size != (target_end - 1) / segment_size`）。此时单个
   `ref_group_id` 无法覆盖完整 target 区间，`ref_end` 会超出参考段长度。处理方式：
   **整段回退 LZ-diff**（与决策 3 一致，不拆段、不跨段引用多 ref_group）。
   此情况在 indel 密集或段尾接近边界时可能发生，但 segment_size（默认 4096）远大于
   典型 indel，实际命中率影响小。

#### 关键问题：PAF target → ref_group_id 映射

PAF 的 `target_name` + `target_start`/`target_end` 是 contig 内坐标，需映射到 pbit 的
`ref_group_id`。ref_groups 按 contig 分组、按 segment_size 切段顺序排列
（[compressor.rs:233-263](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/compressor.rs#L233-L263)），
`contig_ref_groups: IndexMap<String, Vec<u32>>` 提供 contig 名 → ref_group_id 列表。

映射规则：
```
ref_group_id = contig_ref_groups[target_name][target_start / segment_size]
```

参考段 i 覆盖 contig 内区间 `[i * segment_size, (i+1) * segment_size)`（最后一段可能短于
segment_size，需 clamped）。`ref_start`/`ref_end` 存储**相对参考段起始的偏移**（见 §5.2 坐标
语义），即 `ref_start = target_start - seg_idx * segment_size`、
`ref_end = target_end - seg_idx * segment_size`，其中 `seg_idx = target_start / segment_size`。
解压时 `read_2bit_record` 读出的参考段从 `seg_idx * segment_size` 开始，直接用
`ref_start..ref_end` 切片定位 CIGAR 对应的参考区间。

#### 关键问题：未覆盖区域

PAF 没覆盖的样本序列区域：
- **回退 LZ-diff**：对到同名参考段（当前逻辑），encoding=LzDiff
- **存 raw**：直接存储原始序列（无参考），encoding=Raw（需第三种 encoding）
- **推荐**：回退 LZ-diff（复用现有逻辑，格式变化最小）

#### 关键问题：PAF query-side 索引

pbit 压缩时需按**样本（query）段坐标**查询覆盖该段的 PAF alignment。现有
`PafIndex`（[libs/paf/index/mod.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/index/mod.rs)）
的 `reverse_trees` 虽按 query 坐标索引，但**不适合 pbit 直接使用**：

1. **仅覆盖 `+` 链记录**：`insert_record`（[mod.rs:179](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/index/mod.rs#L179)）
   的 mirror entry 只在 `rec.strand == '+'` 时插入。`-` 链记录在 `reverse_trees` 中无条目。
2. **元数据被交换**：mirror entry 的 query_id = 原 target_id，target_start/end = 原 query
   坐标，strand 强制为 `'+'`，CIGAR 被 `reverse_cigar`（I/D 交换）。这是为 BFS 设计的
   "角色互换"视图，与 pbit 需要的原始 PafMetadata 不一致。
3. **无公开查询方法**：`query()`（[query.rs:44](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/index/query.rs#L44)）
   只查 `self.trees`（target 侧）。`reverse_trees` 是 `pub(crate)`，无公开 query 接口。

**方案**：pbit 自建 query-side 区间树。读取 PAF 后，构建独立的
`BasicCOITree<PafMetadata, u32>`（key = query_id，interval = query 坐标），存储**原始未交换**
的 PafMetadata（含原始 strand、原始 target 坐标、原始 CIGAR）。不依赖 `reverse_trees`，
不改现有 PafIndex 的行为。

> **query_name → query_id 映射**：pbit 压缩端按 contig **名**（非 id）处理样本段，PAF 的
> `query_name` 是字符串。构建区间树时先建立 `IndexMap<String, u32>` 映射 query_name →
> query_id（按 PAF 出现顺序分配），区间树以 `query_id` 为 key。压缩端遍历样本段时按
> contig_name 查此映射获取 query_id，再查对应的区间树。

> 实现位置：`libs/pbit/` 内新增 query-side 索引构建逻辑（复用 `libs/paf/` 的 PAF 解析 +
> `coitrees` 区间树，但不依赖 `PafIndex` 结构体）。PAF 记录若无 `cg:Z:` CIGAR 标签
> （`extract_cigar` 返回空 Vec），该记录不插入索引（视为未覆盖，相关段回退 LZ-diff）。

#### append 兼容性

`pgr pbit append` 复用 `Compressor::open_for_append`（[compressor.rs:198-267](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/compressor.rs#L198-L267)），
CIGAR 模式仅影响新追加的样本，已有数据不变。`open_for_append` 重建的 `Segment` 对象
（调用 `prepare` / `prepare_index`）供 LZ-diff 回退段使用，CIGAR 段不依赖它。新追加样本
若带 `--paf` 则走 CIGAR 路径，不带则走 LZ-diff 路径，与已有样本的编码类型互不影响。

### 8. 解压流程（CIGAR 模式）

```
get_contig / get_sample (CIGAR 段):
  1. 读取 DeltaEntry → 检查 encoding（从 DeltaMeta 获取，无需读 packed_data）
  2. 若 encoding == CIGAR:
     a. flate2 解压 packed_data → Vec<u32> → Vec<CigarOp> + X/I 碱基流
        （X/I 碱基按 CIGAR 正向遍历顺序连续存储，解压时按相同顺序消费）
     b. 读取参考段（ref_group_id → seek → read_2bit_record）
     c. 按 SegmentDesc.ref_start/ref_end 截取参考段区间（正向坐标）
     d. 按 CIGAR 重建样本段（正向应用 CIGAR，is_rev_comp 见下方语义说明）:
        - '=' : 从参考取 len 个碱基
        - 'X' : 跳过参考 len 个碱基，从 X/I 碱基流取 len 个样本碱基
        - 'I' : 从 X/I 碱基流取 len 个样本碱基（参考不前进）
        - 'D' : 跳过参考 len 个碱基
     e. 若 is_rev_comp: 对重建结果做反向互补
     f. 拼接所有段
  3. 若 encoding == LZ-diff: 当前逻辑（LZ-diff 解码）
```

> **`decode_delta` 签名扩展**：当前 `decode_delta(ref_group_id, delta_id)`
> （[decompressor.rs:260](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/decompressor.rs#L260)）
> 不接收 ref_start/ref_end。CIGAR 模式下步骤 2c 需要 SegmentDesc 的 ref_start/ref_end 来切片
> 参考段。实现时将 `decode_delta` 签名扩展为接收 `&SegmentDesc`（或额外传入 ref_start/ref_end），
> 由 `get_contig`/`get_sample` 在调用时从 SegmentDesc 提取。LZ-diff 路径忽略这两个字段。

#### 负链 `is_rev_comp` 语义（CIGAR 模式）

`is_rev_comp` 来源于 PAF strand 字段，与 [build_pairwise_block](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/msa_build.rs#L191-L220)
的处理一致：

- **压缩时**：`is_rev_comp = (PAF strand == '-')`。CIGAR 直接取自 PAF（不重写为正向）。
  ref_start/ref_end 始终是**正向参考坐标**（PAF target 坐标，target 在 PAF 中总是正向）。
- **解压时**：先按 CIGAR 在正向参考区间上重建（得到 rev-comp 后的样本段），再对结果做
  `nt::rev_comp` 得到原始样本段。等价于 `build_pairwise_block` 中对 '-' 链先 rev-comp query
  再正向应用 CIGAR 的逆操作。

> **不变量**：ref_start/ref_end 永远是正向参考坐标；is_rev_comp 仅影响样本段方向。
>
> 此操作是 `build_pairwise_block` 的逆——后者对 '-' 链先 RC query 再正向应用 CIGAR，
> 前者先正向应用 CIGAR 再 RC 结果。

#### X/I 操作的碱基存储问题

CIGAR 的 `X`（mismatch）和 `I`（insertion）操作需要存储样本特有的碱基（参考中没有）。
两种处理方式：

**方式 1**：CIGAR 之外额外存储 X/I 碱基
```
packed_data = flate2( Vec<CigarOp> + Vec<u8>(X/I 碱基) )
```
解压时按 CIGAR 遍历，遇到 X/I 从碱基流中取。

**方式 2**：用 CIGAR 的 `M` 操作替代 `=/X`
- `M`（match/mismatch）操作不区分 match 和 mismatch
- 样本碱基全部存储（match 区域也存），仅用 CIGAR 定位 indel
- 压缩率差（match 区域冗余存储），不推荐

**推荐方式 1**：CIGAR ops + X/I 碱基流。存储格式：
```
packed_data = flate2(
    u32 cigar_op_count,
    [CigarOp; cigar_op_count],
    u32 x_i_base_count,
    [u8; x_i_base_count]  // X/I 操作的碱基（编码见下）
)
```

> **N 碱基编码**：pbit 支持 ACGTN（[create.rs:21-22](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/pbit/create.rs#L21-L22)），
> 2-bit 编码（A=0,C=1,G=2,T=3）无法表示 N。采用与 LZ-diff 一致的 5 状态编码
> （LZ-diff 的 5 状态编码（`A=0,C=1,G=2,T=3,N=4`） `A=0,C=1,G=2,T=3,N=4`，
> 需 3 bit）。实现上每碱基存 1 字节（ASCII，简单），或打包为 3-bit 流（紧凑但复杂）。
> 默认采用 1 字节/碱基（ASCII），3-bit 打包作为 Phase 8f 优化候选。

### 9. LZ-diff vs CIGAR 压缩率对比

| 场景 | LZ-diff | CIGAR (bit-packed) | 胜出 |
|------|---------|---------------------|------|
| 纯 match（高相似） | `!` back-ref + match-to-end，极紧凑 | `[=4096]`，1 op = 4B | LZ-diff 略优 |
| 稀疏 SNP（每 100bp 1 SNP） | literal + match，~每 SNP 2B | `=99 X1` ≈ 每 SNP 8B | LZ-diff 优 |
| 密集 SNP | 大量 literal | `Xn`，每段几个 op | CIGAR 优 |
| 含 indel | 位置错位→退化为 literal | `In Dn`，直接描述 | CIGAR 显著优 |
| 重排/转位 | 无法处理（参考选择错误） | PAF 坐标直接描述 | CIGAR 唯一可行 |

> **结论**：CIGAR 在含 indel/重排的场景下显著优于 LZ-diff；在纯 SNP 场景下 LZ-diff 略优。
> 混合模式（方案 D）可让用户按需选择。实际压缩率取决于数据特征，需基准测试验证。
>
> **CIGAR 命中率 caveat**：上述对比假设段被单条 PAF alignment 完整覆盖。实际中若一段跨多条
> alignment 的衔接点（如 alignment1 覆盖 [0,3000)、alignment2 覆盖 [3000,6000)，段为 [1000,5000)），
> 则该段被视为部分覆盖，按 §11 决策 3 整段回退 LZ-diff。高连续性 PAF（minimap2 `--paf-no-hit`、
> wfmash）少见；低连续性 PAF（mashmap 多 mapping、分段输出）会显著降低 CIGAR 命中率。
> Phase 8f 基准测试应包含"PAF 连续性 vs CIGAR 命中率"维度。

### 10. 实施阶段（草案）

#### Phase 8a: CIGAR 编解码基础设施
- `libs/pbit/cigar_delta.rs`：CIGAR ↔ bit-packed 存储 + X/I 碱基流编码
- 复用 `libs/paf/cigar.rs` 的 `CigarOp`
- 公共 API（与 `segment.rs` 的 `Segment::add/get` 对应）：
  ```rust
  /// Pack CIGAR ops + X/I bases into a flate2-compressed byte buffer.
  pub fn pack_cigar(ops: &[CigarOp], xi_bases: &[u8]) -> Vec<u8>;
  /// Unpack a flate2-compressed buffer into (CIGAR ops, X/I bases).
  pub fn unpack_cigar(packed: &[u8]) -> anyhow::Result<(Vec<CigarOp>, Vec<u8>)>;
  /// Apply CIGAR to a reference slice, consuming X/I bases, producing sample seq.
  /// Used by Decompressor. Simplified variant of build_pairwise_block: produces
  /// raw sample sequence (no '-' gap insertion, no coordinate trimming), logic
  /// follows build_maf_block's =/X/M/I/D branches.
  pub fn apply_cigar(ref_seq: &[u8], ops: &[CigarOp], xi_bases: &[u8]) -> anyhow::Result<Vec<u8>>;
  ```
- 单元测试：CIGAR 往返、X/I 碱基流往返、空 CIGAR、纯 match、含 N 的 X/I 碱基

#### Phase 8b: 格式扩展
- `format.rs`：`DeltaMeta`/`DeltaEntry` 新增 `encoding` 字段（10 字节头）；
  `SegmentDesc` 固定大小存储（新增 ref_start/ref_end，LZ-diff 段填 0，见 §5.2）
- Header 版本号升级（1000 → 1001），不保留旧格式读取路径

#### Phase 8c: 压缩端（PAF 驱动）
- `compressor.rs`：新增 `append_sample_with_paf` 方法
  ```rust
  /// Append a sample using PAF-driven CIGAR encoding. Segments covered by PAF
  /// alignments are CIGAR-encoded; uncovered segments fall back to LZ-diff.
  pub fn append_sample_with_paf(
      &mut self,
      sample_name: &str,
      fasta_path: &str,
      paf_path: &str,
  ) -> anyhow::Result<()>;
  ```
- PAF → query-side 区间树构建（见 §7 "PAF query-side 索引"）→ 段级 CIGAR 提取 → CIGAR delta 存储
- 未覆盖段回退 LZ-diff

#### Phase 8d: 解压端
- `decompressor.rs`：`get_contig` / `get_sample` 支持 CIGAR 解码
- `decode_delta` 签名扩展为接收 `&SegmentDesc`（CIGAR 段需 ref_start/ref_end，见 §8 说明）
- 复用 `libs/paf/msa_build.rs` 的 CIGAR 应用逻辑

#### Phase 8e: CLI 集成
- `create.rs` / `append.rs`：新增 `--paf` 参数；`--name` TSV 扩展为可选三列（见 §6.3）
- 测试场景：
  - PAF 驱动往返（`=/X/I/D` CIGAR，`+` 链）
  - `-` 链 PAF 往返（验证 RC 语义：X/I 碱基从 RC(sample) 提取，解压后正向还原）
  - `M` 操作拆分（模拟 minimap2 未加 `--eqx` 的 CIGAR，验证 `M` → `=/X` 拆分正确性）
  - 混合模式（同归档内部分样本 CIGAR、部分 LZ-diff）
  - 未覆盖段回退（PAF 未覆盖的样本段走 LZ-diff）
  - CIGAR target 投影跨参考段边界 → 回退 LZ-diff（§7 决策 3c）
  - 重排场景（样本段比对到不同参考 contig）
  - PAF 无 CIGAR 标签 → 全部回退（决策 7）
  - PAF 文件为空 → 全部回退（决策 7）
  - PAF 单记录解析错误 → 跳过该记录 + warn（决策 8）
  - `--name` 三列 TSV（CIGAR + LZ-diff 混用）
  - `--paf` 与 `--name` 互斥校验
  - `-i`/`--paf` 数量不匹配 → 报错

#### Phase 8f: 基准
- 压缩率对比：LZ-diff vs CIGAR vs 混合（不同数据特征）
- 压缩速度对比（CIGAR 省去 LZ-diff 计算，预期更快）
- 解压速度对比
- **PAF 连续性 vs CIGAR 命中率**：用 minimap2/wfmash/mashmap 分别生成 PAF，统计 CIGAR 模式
  命中段比例（部分覆盖回退 LZ-diff 的段比例）
- **优化候选验证**：SoA 布局（ops/lens 分离）对 flate2 压缩率的提升；3-bit X/I 碱基打包
  对空间的节省

### 11. 决策记录

1. **PAF 粒度**：采用 per-sample PAF，`--paf` 与 `-i` 一一对应。

2. **多 alignment 重叠**：选最佳（identity 最高 + 覆盖度最大），不合并。

3. **部分覆盖**：整段回退 LZ-diff，不拆段。触发回退的三种情况：
   (a) 最佳 alignment 未完整覆盖该段（段有未覆盖的 flank）；
   (b) 该段跨多条 alignment 的衔接点（需合并 CIGAR，不合并则覆盖不全）；
   (c) CIGAR 的 target 投影区间跨越参考段边界（`target_start / segment_size !=
       (target_end - 1) / segment_size`），单个 `ref_group_id` 无法覆盖完整 target 区间。
   判定标准：选最佳 alignment 后，检查其 query 覆盖区间是否完整包含 `[seg_start, seg_end)`，
   且 target 投影区间不跨参考段边界。任一条件不满足 → 整段回退 LZ-diff。

4. **`M` 操作处理**：`pgr maf to-paf` 已改为在源头区分 `=/X`（[cigar.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/cigar.rs)
   的 `cigar_from_alignment` 对比 MAF 对齐碱基，case-insensitive）。因此本项目通过
   UCSC pipeline → `maf to-paf` 路径生成的 PAF 已经是 `=/X`。pbit 压缩时仍保留 `M` → `=/X`
   拆分逻辑作为 fallback，处理 minimap2（未用 `--eqx`）、fastga 等直接输出的 `M` CIGAR。
   **注意**：`M` → `=/X` 拆分需读取样本序列（比对 M 区域的样本碱基与参考碱基）。压缩端有样本
   序列（从 FASTA 读取），可行；解压端无需此逻辑（CIGAR 已是 `=/X`）。`M` → `=/X` 拆分仅在
   CIGAR 含 `M` 操作时触发（如 minimap2 未加 `--eqx`）；纯 `=/X` CIGAR（来自 `maf to-paf`、
   minimap2 `--eqx`）无此开销。

5. **版本兼容**：无需考虑。本项目尚未正式发布，直接升级 Header 版本号（1000 → 1001），
   不保留旧格式读取路径。

6. **query-side 索引方案**：pbit 自建 query-side 区间树，不复用 `PafIndex.reverse_trees`。
   理由：(1) `reverse_trees` 仅覆盖 `+` 链记录，`-` 链完全缺失；(2) mirror entry 的元数据
   被交换（query↔target 角色互换、strand 强制为 `+`、CIGAR 被 reverse），不适合 pbit 需要
   原始 PafMetadata 的场景；(3) `reverse_trees` 无公开查询接口。pbit 在 `libs/pbit/` 内构建
   独立的 `BasicCOITree<PafMetadata, u32>`（key = query_id，interval = query 坐标），存储原始
   未交换的 PafMetadata。详见 §7 "PAF query-side 索引"。

7. **PAF 无 CIGAR 处理**：若 PAF 记录无 `cg:Z:` CIGAR 标签（`extract_cigar` 返回空 Vec），
   该记录不插入 query-side 索引（视为未覆盖，相关样本段回退 LZ-diff）。若整个 PAF 文件无
   任何 CIGAR，所有样本段回退 LZ-diff（等价于无 PAF 输入）。实现时不报错，仅 log 警告。

8. **PAF 解析错误处理**：单条 PAF 记录解析错误（格式错误、坐标非法、字段缺失等）时跳过
   该记录并 `log::warn`，相关样本段视为未覆盖，回退 LZ-diff（与决策 7 一致）。若整个 PAF
   文件无法打开或非 PAF 格式（如文件不存在、首行解析失败且无有效记录），返回错误终止压缩
   （`anyhow::bail!`），避免用户误以为 PAF 已生效但实际全部回退。

## 参考资料

- AGC C++ 源码分析: [notes/references/agc-cpp.md](../references/agc-cpp.md)
- AGC GitHub: https://github.com/refresh-bio/agc
- AGC 论文: Deorowicz et al., "AGC: Assembly Genomes Compressor", Bioinformatics (2024)
- pgr loc 模块: [libs/loc.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/loc.rs)
- pgr 2bit 模块: [libs/fmt/twobit](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs)
- pgr FASTA I/O: [libs/fmt/fa](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/fa.rs)

